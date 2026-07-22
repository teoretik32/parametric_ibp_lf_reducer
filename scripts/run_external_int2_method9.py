#!/usr/bin/env python
"""External Int2 Method.9: surface-policy audit (limit vs convergence-chamber signs).

The production surface filters (``surface.coordinate_primitive_surface_free`` /
``surface.vector_field_surface_free``) decide boundary-term vanishing via ``regulated_sign``,
i.e. in the one-sided limit ``ep -> 0^-`` (value at 0, marginal cases by first-order slope).
The source integral representation, however, converges in a chamber around ``ep = -3/5``
(Part A of the audit verified the integrand conventions exactly at that point; ``r`` never
enters any exponent, so the chamber condition depends on ``ep`` alone).

This script re-runs the SAME library filter code under both sign policies on a small label
box inside the recorded Method.7 witness box and classifies every candidate
(coordinate-IBP multiplier and tangent-field seed) as:

    agree_keep / agree_reject / chamber_only (kept at ep=-3/5, rejected in the limit)
    / limit_only (kept in the limit, rejected at ep=-3/5)

``limit_only`` is the soundness-critical direction: such rows are admitted into production
systems although their boundary terms are NOT justified at the point where the integral
actually converges. Disagreement rows are built as primitive rows and paired against all
recorded Method.7 obstruction witnesses (NO RREF; pairing only, exactly as Method.7).

Diagnostic only: production modules and recorded artifacts are never modified.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from contextlib import contextmanager
from itertools import product
from pathlib import Path

import sympy as sp

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

T2_SCRIPT = Path(__file__).resolve().parent / "run_external_int2_t2_rankrepair.py"


def _load_t2_module():
    """Load the Method.6 runner as a module (scripts/ is not a package)."""
    name = "run_external_int2_t2_rankrepair"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, T2_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


T2 = _load_t2_module()

from parametric_ibp_lf_reducer import surface  # noqa: E402
from parametric_ibp_lf_reducer.lf_obstruction_witness import (  # noqa: E402
    STATUS_WITNESS,
    test_rows_against_obstruction_witness,
    witness_from_payload,
)
from parametric_ibp_lf_reducer.row_generation import (  # noqa: E402
    coordinate_ibp_primitive_row,
    monomials_up_to,
    tangent_ibp_primitive_row,
)
from parametric_ibp_lf_reducer.tangent_fields import generate_tangent_fields  # noqa: E402

INPUT_FILE = REPO_ROOT / "examples" / "external_int2_dimensionless_input.wl.txt"
DEFAULT_OUT = REPO_ROOT / "validation" / "external_int2_method9.json"
DEFAULT_WITNESS = REPO_ROOT / "validation" / "external_int2_method7_witness.json"
DEFAULT_EP = sp.Rational(-3, 5)

# Small audit box: n inside the recorded witness n-box, m in {-1,0}^4 (witness m-box [-3,0]).
N_RANGES = ((-1, 1), (0, 1), (0, 1))
M_RANGES = ((-1, 0),) * 4
MAX_IBP_DEGREE = int(T2.BASELINE["max_ibp_degree"])
TANGENT_BLOCKS = tuple(tuple(b) for b in T2.BASELINE["tangent_degree_blocks"])


def chamber_sign_factory(ep_value: sp.Rational):
    """A ``regulated_sign``-compatible policy: exact rational sign at ``ep = ep_value``."""

    def chamber_sign(expr, regulators, direction: str = "minus") -> str:
        del direction  # the chamber is a point; no one-sided limit is involved
        val = sp.sympify(expr)
        for reg in regulators:
            val = val.subs(sp.Symbol(str(reg)), ep_value)
        val = sp.nsimplify(sp.expand(val))
        if val.free_symbols:
            return "Unknown"
        if val > 0:
            return "pos"
        if val < 0:
            return "neg"
        return "zero"

    return chamber_sign


@contextmanager
def sign_policy(fn):
    """Temporarily swap ``surface.regulated_sign`` (the ONLY policy difference)."""
    orig = surface.regulated_sign
    surface.regulated_sign = fn
    try:
        yield
    finally:
        surface.regulated_sign = orig


def audit_labels():
    """Seed labels of the small audit box (inside the recorded witness label box)."""
    n_axes = [range(lo, hi + 1) for lo, hi in N_RANGES]
    m_axes = [range(lo, hi + 1) for lo, hi in M_RANGES]
    return [tuple(n) + tuple(m) for n in product(*n_axes) for m in product(*m_axes)]


def classify(family, labels, fields, ep_value):
    """Run both sign policies over every candidate; return verdict table + disagreements."""
    chamber = chamber_sign_factory(ep_value)
    monos = monomials_up_to(family.nvars, MAX_IBP_DEGREE)
    table = []  # (kind, provenance, verdict_limit, verdict_chamber)

    def _verdicts(check):
        v_lim = check()
        with sign_policy(chamber):
            v_ch = check()
        return v_lim, v_ch

    for label in labels:
        for var_index in range(family.nvars):
            for p in monos:
                v_lim, v_ch = _verdicts(
                    lambda: surface.coordinate_primitive_surface_free(
                        family, label, var_index, multiplier_exps=p, eps_direction="minus"
                    )
                )
                table.append(("coordinate_ibp", {"seed": label, "var": var_index, "P": p}, v_lim, v_ch))
        for fi, tf in enumerate(fields):
            v_lim, v_ch = _verdicts(
                lambda: surface.vector_field_surface_free(
                    family, label, list(tf.components), eps_direction="minus"
                )
            )
            table.append(("tangent_ibp", {"seed": label, "field_index": fi}, v_lim, v_ch))
    return table


def build_rows(family, fields, entries):
    """Primitive rows for disagreement entries (dedup; trivial rows dropped)."""
    rows, seen, dropped = [], set(), 0
    for kind, prov, _vl, _vc in entries:
        if kind == "coordinate_ibp":
            row = coordinate_ibp_primitive_row(family, prov["seed"], prov["var"], prov["P"])
        else:
            row = tangent_ibp_primitive_row(family, prov["seed"], fields[prov["field_index"]])
        if row.is_trivial():
            dropped += 1
            continue
        key = row.dedup_key()
        if key in seen:
            continue
        seen.add(key)
        rows.append((row, prov))
    return rows, dropped


def screen(rows_with_prov, witness_payloads):
    """Pair rows against every recorded witness at its own (sample, prime) — Method.7 style."""
    rows = [rp[0] for rp in rows_with_prov]
    out = []
    for wp in witness_payloads:
        witness = witness_from_payload(wp)
        sample = T2._sample_from_repr([list(pair) for pair in witness.sample])
        pairings = test_rows_against_obstruction_witness(rows, witness, sample, witness.prime)
        breaking = [p for p in pairings if p.breaks]
        out.append(
            {
                "witness_prime": witness.prime,
                "witness_sample": [list(pair) for pair in witness.sample],
                "n_rows": len(rows),
                "n_breaks": len(breaking),
                "breaking_examples": [
                    {"kind": rows[p.index].kind,
                     "provenance": {k: str(v) for k, v in rows_with_prov[p.index][1].items()}}
                    for p in breaking[:8]
                ],
            }
        )
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--witness-file", type=Path, default=DEFAULT_WITNESS)
    parser.add_argument("--ep", type=str, default=str(DEFAULT_EP),
                        help="chamber point for the exact sign policy (rational, default -3/5)")
    args = parser.parse_args(argv)

    t0 = time.time()
    ep_value = sp.Rational(args.ep)
    family = T2.parse_family_text(INPUT_FILE.read_text(encoding="utf-8"))
    labels = audit_labels()
    fields = generate_tangent_fields(family, TANGENT_BLOCKS)
    print(f"[m9] box: {len(labels)} seeds, {len(fields)} tangent fields, "
          f"ibp deg<={MAX_IBP_DEGREE}, chamber ep={ep_value}", flush=True)

    table = classify(family, labels, fields, ep_value)
    counts = {"agree_keep": 0, "agree_reject": 0, "chamber_only": 0, "limit_only": 0, "unknown": 0}
    chamber_only, limit_only = [], []
    for entry in table:
        _kind, _prov, v_lim, v_ch = entry
        if v_lim == "Unknown" or v_ch == "Unknown":
            counts["unknown"] += 1
        elif v_lim is True and v_ch is True:
            counts["agree_keep"] += 1
        elif v_lim is not True and v_ch is not True:
            counts["agree_reject"] += 1
        elif v_ch is True:
            counts["chamber_only"] += 1
            chamber_only.append(entry)
        else:
            counts["limit_only"] += 1
            limit_only.append(entry)
    print(f"[m9] verdicts: {counts}", flush=True)

    witness_payloads = [
        wp for wp in json.loads(args.witness_file.read_text(encoding="utf-8"))["points"]
        if wp["status"] == STATUS_WITNESS
    ]

    report = {
        "script": "scripts/run_external_int2_method9.py",
        "method": "External Int2 Method.9: surface-policy audit (limit vs chamber)",
        "scope_note": (
            "Two sign policies through the IDENTICAL library filter code: production "
            "regulated_sign at ep->0^- vs exact rational sign at the convergence-chamber point "
            f"ep={ep_value}. r never enters exponents, so the chamber test depends on ep only. "
            "Pairing-only against recorded Method.7 witnesses; NO RREF; diagnostic only."
        ),
        "chamber_ep": str(ep_value),
        "box": {"n_ranges": [list(r) for r in N_RANGES], "m_ranges": [list(r) for r in M_RANGES],
                "n_seeds": len(labels), "max_ibp_degree": MAX_IBP_DEGREE,
                "tangent_degree_blocks": [list(b) for b in TANGENT_BLOCKS],
                "n_tangent_fields": len(fields)},
        "n_checks": len(table),
        "verdict_counts": counts,
        "families": {},
    }

    for name, entries in (("chamber_only", chamber_only), ("limit_only", limit_only)):
        rows_with_prov, dropped = build_rows(family, fields, entries)
        by_kind: dict[str, int] = {}
        for kind, _p, _vl, _vc in entries:
            by_kind[kind] = by_kind.get(kind, 0) + 1
        fam_report = {
            "n_entries": len(entries),
            "by_kind": by_kind,
            "n_rows": len(rows_with_prov),
            "n_trivial_dropped": dropped,
            "per_witness": screen(rows_with_prov, witness_payloads) if rows_with_prov else [],
        }
        fam_report["n_breaks_total"] = sum(pw["n_breaks"] for pw in fam_report["per_witness"])
        report["families"][name] = fam_report
        print(f"[m9] {name}: {len(entries)} entries -> {len(rows_with_prov)} rows, "
              f"breaks_total={fam_report['n_breaks_total']}", flush=True)

    report["elapsed_seconds"] = round(time.time() - t0, 3)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=1), encoding="utf-8")
    print(f"[m9] report -> {args.out} ({report['elapsed_seconds']}s)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
