"""External Int2 Method.11 Phase C: independent certificate for the reconstructed LF relation.

Rebuilds the Method.11 Level-0 chamber row system ONCE (same construction as
``scripts/run_external_int2_method11.py``: 3072 labels, 49439 rows, chamber ``ep = -3/5``) and
then uses the first-class certificate API to check that

    J[0,0,0,0,0,0,0] - sum_i C_i(ep, r) * J[label_i]

lies in the row span at NEW exact rational ``(ep, r)`` points, for >= 2 large primes per point.

Independence from Phase B: no cached ``NormalFormRecord`` is read; the claimed coefficients
``C_i`` are the closed-form rational functions from the reconstruction audit artifact
(``validation/external_int2_method11_reconstruction_audit.json``), and the check points are
required to be disjoint from every ``(ep, r)`` sample recorded in the Phase B record cache.

Success criteria (all required for ``verdict: certified``):
  * every requested ``(point, prime)`` pair is informative (no BadSpecialization poles),
  * every informative pair reports ``InSpan``,
  * the rank at every pair equals the recorded generic rank (26984),
  * at least ``--min-points`` distinct points were certified.
"""

from __future__ import annotations

import argparse
import ast
import importlib.util
import json
import sys
import time
from dataclasses import replace
from fractions import Fraction
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPTS_DIR.parent
sys.path.insert(0, str(REPO_ROOT / "src"))

from parametric_ibp_lf_reducer import parse_family_text  # noqa: E402
from parametric_ibp_lf_reducer.certificate import (  # noqa: E402
    STATUS_BAD_SPECIALIZATION,
    STATUS_IN_SPAN,
    verify_reduction_relations_mod_p,
)


def _load_script(stem: str):
    """Load a sibling runner as a module (``scripts/`` is not a package)."""
    if stem in sys.modules:
        return sys.modules[stem]
    spec = importlib.util.spec_from_file_location(stem, SCRIPTS_DIR / f"{stem}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[stem] = mod
    spec.loader.exec_module(mod)
    return mod


M11 = _load_script("run_external_int2_method11")
T2 = M11.T2

DEFAULT_POINTS = "5/2,7/3;8/3,13/4;11/4,17/5"
DEFAULT_PRIMES = "2147483647,2147483629"
DEFAULT_AUDIT = REPO_ROOT / "validation" / "external_int2_method11_reconstruction_audit.json"
DEFAULT_RECORDS = REPO_ROOT / "validation" / "external_int2_method11_records.jsonl"
DEFAULT_OUT = REPO_ROOT / "validation" / "external_int2_lf_certificate.json"

EXPECTED_LABELS = 3072
EXPECTED_ROWS = 49439
GENERIC_RANK = 26984

SCOPE_NOTE = (
    "Method.11 Phase C: independent row-span certificate for the reconstructed four-term LF "
    "reduction relation at NEW off-sample (ep, r) points. 'Certified' claims only that the "
    "claimed relation vector lies in the chamber-policy row span at the listed exact rational "
    "points modulo the listed primes, with the generic rank reproduced. It does NOT claim the "
    "analytic Laurent value; that cross-check is a separate later phase."
)


def parse_points(text: str) -> list[tuple[Fraction, Fraction]]:
    """Parse ``ep,r;ep,r;...`` into exact rational pairs (floats rejected by Fraction)."""
    points: list[tuple[Fraction, Fraction]] = []
    for chunk in text.split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        parts = chunk.split(",")
        if len(parts) != 2:
            raise ValueError(f"bad point {chunk!r}: expected 'ep,r'")
        points.append((Fraction(parts[0]), Fraction(parts[1])))
    if len(set(points)) != len(points):
        raise ValueError("duplicate points requested")
    return points


def load_claimed_coefficients(audit_path: Path) -> dict[tuple, str]:
    """The four reconstructed coefficients; every entry must have status ``validated``."""
    data = json.loads(audit_path.read_text(encoding="utf-8"))
    coeffs = data.get("coefficients")
    if not isinstance(coeffs, dict) or not coeffs:
        raise SystemExit(f"no coefficients in {audit_path}")
    out: dict[tuple, str] = {}
    for key, entry in coeffs.items():
        if entry.get("status") != "validated":
            raise SystemExit(f"coefficient {key} is not validated: {entry.get('status')!r}")
        label = tuple(ast.literal_eval(key))
        out[label] = entry["expression"]
    return out


def recorded_sample_points(records_path: Path) -> set[tuple[str, str]]:
    """All distinct ``(ep, r)`` sample coordinates in the Phase B record cache (as strings)."""
    seen: set[tuple[str, str]] = set()
    if not records_path.exists():
        return seen
    with records_path.open(encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            data = json.loads(line)
            sample = data.get("sample")
            if isinstance(sample, dict) and "ep" in sample and "r" in sample:
                seen.add((str(sample["ep"]), str(sample["r"])))
    return seen


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--points", default=DEFAULT_POINTS, help="semicolon-separated ep,r pairs")
    parser.add_argument("--primes", default=DEFAULT_PRIMES, help="comma-separated odd primes")
    parser.add_argument("--min-points", type=int, default=3)
    parser.add_argument("--audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    parser.add_argument("--input", type=Path, default=T2.DEFAULT_INPUT)
    parser.add_argument("--ep", default="-3/5", help="chamber point for the surface policy")
    parser.add_argument("--rref-backend", choices=("auto", "python", "numba"), default="auto")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args(argv)

    points = parse_points(args.points)
    primes = tuple(int(p) for p in args.primes.split(",") if p.strip())
    if len(primes) < 2:
        parser.error("need >= 2 primes")
    if len(points) < args.min_points:
        parser.error(f"need >= {args.min_points} points")

    claimed = load_claimed_coefficients(args.audit)
    print(f"[cert] claimed relation: 1 target -> {len(claimed)} RHS labels", flush=True)
    for label, expr in sorted(claimed.items()):
        print(f"[cert]   C{list(label)} = {expr}", flush=True)

    recorded = recorded_sample_points(args.records)
    for ep, r in points:
        if (str(ep), str(r)) in recorded:
            raise SystemExit(f"point (ep={ep}, r={r}) is NOT off-sample: it appears in "
                             f"{args.records}")
    print(f"[cert] {len(points)} points, all disjoint from {len(recorded)} recorded samples",
          flush=True)

    t0 = time.time()
    ep_value = Fraction(args.ep)
    policy = M11.chamber_policy(ep_value)
    family = parse_family_text(args.input.read_text(encoding="utf-8"))
    _target, base_config = T2.build_reducer_config(family)
    labels, target_label, config = T2.build_level_config(family, base_config, T2.LEVELS[0])
    config = replace(config, surface_policy=policy)
    print(f"[cert] Level 0 box: {len(labels)} labels; chamber ep={ep_value}", flush=True)

    t_asm = time.time()
    asm = M11.assemble_chamber_rows(family, labels, config)
    rows = asm["merged"]
    print(f"[cert] chamber-policy rows: {asm['n_rows_total']} {asm['by_kind']} "
          f"({time.time() - t_asm:.1f}s)", flush=True)
    if len(labels) != EXPECTED_LABELS or asm["n_rows_total"] != EXPECTED_ROWS:
        raise SystemExit(f"row system mismatch: {len(labels)} labels / {asm['n_rows_total']} "
                         f"rows, expected {EXPECTED_LABELS} / {EXPECTED_ROWS}")

    unknown_rhs = [lab for lab in claimed if lab not in set(labels)]
    if unknown_rhs:
        raise SystemExit(f"claimed RHS labels not in the Level-0 box: {unknown_rhs}")

    relations = {tuple(target_label): claimed}
    checks: list[dict] = []
    n_in_span = 0
    n_bad = 0
    for ep, r in points:
        sample = {"ep": ep, "r": r}
        for prime in primes:
            t_chk = time.time()
            res = verify_reduction_relations_mod_p(
                family, rows, relations, sample, prime, rref_backend=args.rref_backend,
            )[tuple(target_label)]
            elapsed = round(time.time() - t_chk, 1)
            entry = {
                "ep": str(ep),
                "r": str(r),
                "prime": prime,
                "status": res.status,
                "in_span": res.in_span,
                "rank": res.rank,
                "nrows": res.nrows,
                "rank_matches_generic": res.rank == GENERIC_RANK,
                "n_residual_terms": len(res.residual or {}),
                "elapsed_seconds": elapsed,
            }
            checks.append(entry)
            if res.status == STATUS_IN_SPAN:
                n_in_span += 1
            elif res.status == STATUS_BAD_SPECIALIZATION:
                n_bad += 1
            print(f"[cert] (ep={ep}, r={r}) p={prime}: {res.status} rank={res.rank} "
                  f"({elapsed}s)", flush=True)

    informative = [c for c in checks if c["status"] != STATUS_BAD_SPECIALIZATION]
    certified_points = {
        (c["ep"], c["r"])
        for c in checks
        if c["status"] == STATUS_IN_SPAN and c["rank_matches_generic"]
    }
    all_pairs_informative = len(informative) == len(checks)
    all_informative_in_span = all(c["status"] == STATUS_IN_SPAN for c in informative)
    all_ranks_generic = all(c["rank_matches_generic"] for c in informative)
    fully_certified_points = {
        pt for pt in certified_points
        if all(c["status"] == STATUS_IN_SPAN and c["rank_matches_generic"]
               for c in checks if (c["ep"], c["r"]) == pt)
    }
    certified = (
        all_pairs_informative
        and all_informative_in_span
        and all_ranks_generic
        and len(fully_certified_points) >= args.min_points
    )

    report = {
        "schema": "external-int2-lf-certificate/v1",
        "script": "scripts/certify_external_int2_lf.py",
        "scope_note": SCOPE_NOTE,
        "chamber_ep": str(ep_value),
        "level": 0,
        "n_labels": len(labels),
        "target_label": list(target_label),
        "claimed_coefficients": {str(list(k)): v for k, v in sorted(claimed.items())},
        "coefficient_source": str(args.audit),
        "off_sample_guard": {"recorded_samples": len(recorded), "points_disjoint": True},
        "generic_rank_expected": GENERIC_RANK,
        "rows": {k: v for k, v in asm.items() if k not in ("merged", "row_diagnostics")},
        "primes": list(primes),
        "checks": checks,
        "n_checks": len(checks),
        "n_in_span": n_in_span,
        "n_bad_specialization": n_bad,
        "n_certified_points": len(fully_certified_points),
        "min_points": args.min_points,
        "verdict": "certified" if certified else "not-certified",
        "elapsed_seconds": round(time.time() - t0, 1),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[cert] verdict: {report['verdict']} ({len(fully_certified_points)} certified "
          f"points, {n_in_span}/{len(checks)} InSpan); wrote {args.out}", flush=True)
    return 0 if certified else 1


if __name__ == "__main__":
    raise SystemExit(main())
