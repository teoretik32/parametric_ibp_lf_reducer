"""Offline Method.11 reconstruction audit from the persistent modular record cache.

This script performs NO RREF solves.  It detects support-changing sample points,
excludes them conservatively from fitting, reconstructs the generic rational
coefficient functions, and checks every modular record exactly.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from pathlib import Path

import sympy as sp

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from parametric_ibp_lf_reducer.reconstruction import reconstruct_coefficients  # noqa: E402

DEFAULT_CACHE = REPO_ROOT / "validation" / "external_int2_method11_records.jsonl"
DEFAULT_OUT = REPO_ROOT / "validation" / "external_int2_method11_reconstruction_audit.json"


def _load_runner():
    path = REPO_ROOT / "scripts" / "run_external_int2_method11.py"
    spec = importlib.util.spec_from_file_location("external_int2_method11_audit_load", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _sample_key(sample) -> str:
    return ",".join(f"{k}={Fraction(v)}" for k, v in sorted(sample.items()))


def _eval_mod(expr: sp.Expr, sample, prime: int) -> int:
    syms = {sp.Symbol(str(k)): sp.Rational(Fraction(v).numerator, Fraction(v).denominator) for k, v in sample.items()}
    value = sp.cancel(expr.subs(syms))
    num, den = sp.fraction(value)
    num = sp.Rational(num)
    den = sp.Rational(den)
    nmod = int(num.p) * pow(int(num.q), -1, prime) % prime
    dmod = int(den.p) * pow(int(den.q), -1, prime) % prime
    return nmod * pow(dmod, -1, prime) % prime


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args(argv)

    runner = _load_runner()
    lines = args.cache.read_text(encoding="utf-8").splitlines()
    header = json.loads(lines[0])
    records = [runner._record_from_json(json.loads(line)) for line in lines[1:] if line.strip()]

    grouped = defaultdict(list)
    for rec in records:
        grouped[_sample_key(rec.sample)].append(rec)

    rank_hist = Counter(rec.rank for rec in records if rec.formal_success)
    generic_rank = max(rank_hist)
    generic = [rec for rec in records if rec.formal_success and rec.rank == generic_rank]
    support_counts = Counter(tuple(sorted(rec.coeffs)) for rec in generic)
    modal_support, modal_count = support_counts.most_common(1)[0]
    modal_support = tuple(modal_support)

    changed_samples = []
    stable_records = []
    for key, recs in sorted(grouped.items()):
        supports = {tuple(sorted(rec.coeffs)) for rec in recs}
        ranks = {rec.rank for rec in recs}
        if ranks == {generic_rank} and supports == {modal_support}:
            stable_records.extend(recs)
        else:
            changed_samples.append({
                "sample": key,
                "ranks": sorted(ranks),
                "supports": [[list(lab) for lab in support] for support in sorted(supports)],
                "n_records": len(recs),
            })

    report = {}
    coeffs = reconstruct_coefficients(stable_records, ["ep", "r"], report=report)

    stable_mismatches = []
    changed_mismatches = []
    for rec in records:
        key = _sample_key(rec.sample)
        target = changed_mismatches if any(x["sample"] == key for x in changed_samples) else stable_mismatches
        for lab, expr in coeffs.items():
            pred = _eval_mod(expr, rec.sample, rec.prime)
            actual = rec.coeffs.get(lab, 0)
            if pred != actual:
                target.append({
                    "sample": key,
                    "prime": rec.prime,
                    "label": list(lab),
                    "predicted": pred,
                    "recorded": actual,
                })

    payload = {
        "schema": "external-int2-method11-reconstruction-audit/v1",
        "cache_schema": header.get("schema"),
        "n_records": len(records),
        "n_samples": len(grouped),
        "n_primes": len({rec.prime for rec in records}),
        "generic_rank": generic_rank,
        "rank_histogram": {str(k): v for k, v in sorted(rank_hist.items())},
        "modal_support": [list(lab) for lab in modal_support],
        "modal_support_record_count": modal_count,
        "n_stable_samples": len({ _sample_key(r.sample) for r in stable_records }),
        "n_changed_samples": len(changed_samples),
        "changed_samples": changed_samples,
        "classification": "basis_specialization_not_special_zero" if changed_samples and changed_mismatches else "no_basis_specialization_detected",
        "coefficients": {
            str(lab): {
                "expression": str(sp.factor(expr)),
                "num_deg": report[str(lab)].get("num_deg"),
                "den_deg": report[str(lab)].get("den_deg"),
                "n_points": report[str(lab)].get("n_points"),
                "n_fit": report[str(lab)].get("n_fit"),
                "n_hold": report[str(lab)].get("n_hold"),
                "status": report[str(lab)].get("status"),
            }
            for lab, expr in coeffs.items()
        },
        "stable_record_mismatches": stable_mismatches,
        "changed_record_mismatches": changed_mismatches,
        "stable_records_exactly_matched": not stable_mismatches,
        "changed_sample_is_not_generic_zero": bool(changed_mismatches),
        "recommendation": (
            "Exclude support-changing samples from fitting. Reconstruct the generic coefficient functions from full-support, max-rank records; "
            "use support-changing points only after reconstruction to distinguish a genuine coefficient zero from a basis/pivot specialization."
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    return 0 if not stable_mismatches else 2


if __name__ == "__main__":
    raise SystemExit(main())
