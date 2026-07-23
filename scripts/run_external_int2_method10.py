#!/usr/bin/env python
"""External Int2 Method.10: chamber-policy controlled LF-feasibility solve (Phase E gate).

Method.9 found the first policy-dependent lever in the implementation: under the exact
convergence-chamber sign policy (``ep = -3/5``) strictly more rows pass the IDENTICAL library
surface-filter code than under the production one-sided limit ``ep -> 0^-``, and the chamber-only
rows break the recorded Method.7 dual obstruction witnesses. Per the audit protocol a witness
break is NECESSARY but not sufficient; the next authorized step is exactly one small controlled
solve — this runner:

1. assembles the recorded Level-0 system under the production limit policy (semantics unchanged);
2. re-assembles the SAME system inside a diagnostic ``sign_policy(chamber_sign_factory(ep))``;
3. diffs by ``dedup_key`` -> chamber-only rows (and checks ``limit_only == 0``: chamber acceptance
   strictly extends production acceptance, so no production row loses its justification);
4. pairs the chamber-only rows against every recorded Method.7 witness (cheap, NO RREF);
5. ONLY if a chamber-only row breaks a witness: runs ``lf_reduction_feasible_mod_p`` on the merged
   (limit + chamber-only) system at generic (sample, prime) points, gated by ``--allow-heavy``.

Diagnostic only: ``surface.regulated_sign`` is swapped inside a context manager local to this
script; reducer core, production filters, certificates and LF gates are untouched. ``Feasible``
here is modular span evidence through LF-True labels only — NOT an LF basis; rational
reconstruction + row-span certificate + per-term LF export checks remain mandatory before any
Success claim. ``Obstructed`` stays an honest per-(sample, prime), per-box negative.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path

import sympy as sp

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

SCRIPTS_DIR = Path(__file__).resolve().parent


def _load_script(stem: str):
    """Load a sibling runner as a module (``scripts/`` is not a package)."""
    if stem in sys.modules:
        return sys.modules[stem]
    spec = importlib.util.spec_from_file_location(stem, SCRIPTS_DIR / f"{stem}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[stem] = mod
    spec.loader.exec_module(mod)
    return mod


T2 = _load_script("run_external_int2_t2_rankrepair")
M9 = _load_script("run_external_int2_method9")

DEFAULT_OUT = REPO_ROOT / "validation" / "external_int2_method10.json"
DEFAULT_WITNESS = REPO_ROOT / "validation" / "external_int2_method7_witness.json"
DEFAULT_EP = sp.Rational(-3, 5)
DEFAULT_PRIMES = T2.DEFAULT_WITNESS_PRIMES

VERDICT_SCREENING_ONLY = "ScreeningOnly"
VERDICT_FEASIBLE = "Feasible(modular)"
VERDICT_OBSTRUCTED = "Obstructed"
VERDICT_MIXED = "Mixed"

SCOPE_NOTE = (
    "Method.10 controlled solve: Level-0 system re-assembled under the exact convergence-chamber "
    "sign policy (diagnostic swap of surface.regulated_sign inside this script only). The span "
    "test is constrained to LF-True labels; a span through divergent labels is NOT an LF basis "
    "and is never reported as one. 'Feasible(modular)' is per-(sample, prime) span evidence "
    "only -- reconstruction, row-span certificate and per-term LF export checks are still "
    "required before Success. 'Obstructed' remains an honest per-point, per-box negative, never "
    "a global claim. Reducer core, production surface policy, certificates and LF gates are "
    "untouched."
)


def diff_by_dedup_key(base_rows, other_rows):
    """Rows of ``other_rows`` not in ``base_rows`` and vice versa (first occurrence per key)."""
    base_keys = {row.dedup_key() for row in base_rows}
    other_keys = set()
    only_other = []
    for row in other_rows:
        key = row.dedup_key()
        if key in other_keys:
            continue
        other_keys.add(key)
        if key not in base_keys:
            only_other.append(row)
    only_base = [row for row in base_rows if row.dedup_key() not in other_keys]
    return only_other, only_base


def overall_verdict(statuses) -> str:
    """Aggregate per-point feasibility statuses (empty -> screening only)."""
    if not statuses:
        return VERDICT_SCREENING_ONLY
    feas = sum(1 for s in statuses if s == T2.STATUS_FEASIBLE)
    obst = sum(1 for s in statuses if s == T2.STATUS_OBSTRUCTED)
    if feas == len(statuses):
        return VERDICT_FEASIBLE
    if obst and not feas:
        return VERDICT_OBSTRUCTED
    return VERDICT_MIXED


def _row_key_hash(row) -> str:
    return hashlib.sha256(repr(row.dedup_key()).encode()).hexdigest()[:16]


def screen_against_witnesses(rows, witness_payloads) -> list[dict]:
    """Pair ``rows`` against every recorded witness at its own (sample, prime). NO RREF."""
    out = []
    for wp in witness_payloads:
        witness = T2.witness_from_payload(wp)
        sample = T2._sample_from_repr([list(pair) for pair in witness.sample])
        pairings = T2.test_rows_against_obstruction_witness(rows, witness, sample, witness.prime)
        breaking = [p for p in pairings if p.breaks]
        out.append(
            {
                "witness_prime": witness.prime,
                "witness_sample": [list(pair) for pair in witness.sample],
                "n_rows": len(rows),
                "n_breaks": len(breaking),
                "breaking_examples": [
                    {"kind": rows[p.index].kind, "dedup_sha16": _row_key_hash(rows[p.index])}
                    for p in breaking[:8]
                ],
            }
        )
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", type=Path, default=T2.DEFAULT_INPUT)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--witness-file", type=Path, default=DEFAULT_WITNESS)
    parser.add_argument("--ep", type=str, default=str(DEFAULT_EP),
                        help="chamber point for the exact sign policy (rational, default -3/5)")
    parser.add_argument("--allow-heavy", action="store_true",
                        help="required for any RREF solve (--max-solves > 0)")
    parser.add_argument("--max-solves", type=int, default=0,
                        help="max (sample, prime) feasibility solves; 0 = screening only")
    parser.add_argument("--primes", type=str,
                        default=",".join(str(p) for p in DEFAULT_PRIMES))
    parser.add_argument("--no-stop-on-obstructed", action="store_true",
                        help="keep solving remaining points after an Obstructed generic point")
    args = parser.parse_args(argv)

    if args.max_solves > 0 and not args.allow_heavy:
        parser.error("--max-solves > 0 requires --allow-heavy (Level-0 RREF is expensive)")
    if args.max_solves < 0:
        parser.error("--max-solves must be >= 0")
    primes = tuple(int(p) for p in args.primes.split(",") if p.strip())
    if any(p < 3 for p in primes):
        parser.error("primes must be odd primes > 2")

    t0 = time.time()
    ep_value = sp.Rational(args.ep)
    family = T2.parse_family_text(args.input.read_text(encoding="utf-8"))
    _target, base_config = T2.build_reducer_config(family)
    spec = T2.LEVELS[0]
    labels, target_label, config = T2.build_level_config(family, base_config, spec)
    print(f"[m10] Level 0 box: {len(labels)} labels; chamber ep={ep_value}", flush=True)

    t_asm = time.time()
    asm_limit = T2.assemble_level_rows(family, labels, config)
    print(f"[m10] limit-policy rows: {asm_limit['n_rows_total']} "
          f"({time.time() - t_asm:.1f}s)", flush=True)

    t_asm = time.time()
    with M9.sign_policy(M9.chamber_sign_factory(ep_value)):
        asm_chamber = T2.assemble_level_rows(family, labels, config)
    print(f"[m10] chamber-policy rows: {asm_chamber['n_rows_total']} "
          f"({time.time() - t_asm:.1f}s)", flush=True)

    chamber_only, limit_only = diff_by_dedup_key(asm_limit["merged"], asm_chamber["merged"])
    by_kind: dict[str, int] = {}
    for row in chamber_only:
        by_kind[row.kind] = by_kind.get(row.kind, 0) + 1
    print(f"[m10] chamber_only: {len(chamber_only)} {by_kind}; "
          f"limit_only: {len(limit_only)}", flush=True)

    witness_payloads = [
        wp for wp in json.loads(args.witness_file.read_text(encoding="utf-8"))["points"]
        if wp["status"] == T2.STATUS_WITNESS
    ]
    screening = screen_against_witnesses(chamber_only, witness_payloads)
    n_breaks_total = sum(s["n_breaks"] for s in screening)
    rerun_justified = n_breaks_total > 0 and not limit_only
    print(f"[m10] screening: {n_breaks_total} breaks over {len(screening)} witnesses; "
          f"rerun_justified={rerun_justified}", flush=True)

    solve_points: list[dict] = []
    if args.max_solves > 0 and rerun_justified:
        merged = list(asm_limit["merged"]) + chamber_only
        lf_map = {lab: T2.is_locally_finite(family, lab) for lab in labels}
        if target_label not in lf_map:
            lf_map[target_label] = T2.is_locally_finite(family, target_label)
        samples = T2.default_scattered_samples(family.parameters)
        generic = [s for i, s in enumerate(samples) if i != T2.SPECIAL_SAMPLE_INDEX]
        n_done = 0
        stop = False
        for sample in generic:
            for prime in primes:
                if n_done >= args.max_solves or stop:
                    break
                t_solve = time.time()
                res = T2.lf_reduction_feasible_mod_p(
                    merged, labels, target_label, lf_map, sample, prime
                )
                payload = T2.feasibility_to_payload(res)
                payload["solve_seconds"] = round(time.time() - t_solve, 3)
                solve_points.append(payload)
                n_done += 1
                print(f"[m10] solve {n_done}: prime={prime} "
                      f"sample={T2._sample_repr(sample)} -> {res.status} "
                      f"rank={res.rank} ({payload['solve_seconds']}s)", flush=True)
                if res.status == T2.STATUS_OBSTRUCTED and not args.no_stop_on_obstructed:
                    stop = True
            if stop or n_done >= args.max_solves:
                break
    elif args.max_solves > 0:
        print("[m10] solve skipped: screening gate not passed", flush=True)

    verdict = overall_verdict([p["status"] for p in solve_points])
    report = {
        "script": "scripts/run_external_int2_method10.py",
        "method": "External Int2 Method.10: chamber-policy controlled LF-feasibility solve",
        "scope_note": SCOPE_NOTE,
        "chamber_ep": str(ep_value),
        "level": 0,
        "n_labels": len(labels),
        "target_label": list(target_label),
        "rows": {
            "n_limit": asm_limit["n_rows_total"],
            "n_chamber": asm_chamber["n_rows_total"],
            "n_chamber_only": len(chamber_only),
            "n_limit_only": len(limit_only),
            "chamber_only_by_kind": by_kind,
        },
        "screening": {
            "n_witnesses": len(screening),
            "n_breaks_total": n_breaks_total,
            "per_witness": screening,
            "rerun_justified": rerun_justified,
        },
        "solve": {
            "requested_max_solves": args.max_solves,
            "primes": list(primes),
            "points": solve_points,
            "verdict": verdict,
        },
        "elapsed_seconds": round(time.time() - t0, 3),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"[m10] verdict: {verdict}; wrote {args.out} "
          f"({report['elapsed_seconds']}s total)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
