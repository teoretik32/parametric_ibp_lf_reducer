"""External Int2 Method.1 diagnostics: directional LF audit + LF-constrained feasibility.

Phases (``--phases a|ab|abc``):

* ``a`` — directional local-finiteness audit (``explain_local_finiteness``): LF verdict
  counts over the label box, a full report for the target, and capped per-label detail
  (failing rays, strictly improving unit shifts) for the first ``--audit-cap`` non-LF-True
  labels;
* ``b`` — LF-constrained span feasibility mod p (``lf_reduction_feasible_mod_p``) on a
  small deterministic (sample, prime) grid: can the target be reduced through LF-True
  labels only, independent of the reducer's normal-form ranking?
* ``c`` — one explicit LF-supported reduction per Feasible point
  (``lf_reduction_coefficients_mod_p``) plus a support-stability summary; per-point
  statuses must agree with phase b (exit 1 otherwise).

Levels (``--level A|B``) are the first two boxes of ``default_search_levels``: A = the
document's base box (``max_ibp_degree=1``, no tangent rows); B = every m-range deepened
by one (``max_ibp_degree=2``, tangent blocks ``((1, 1),)``).

Caps by design: ``--samples`` (default 3) x ``--primes`` (default 2) modular points and
``--audit-cap`` detailed audits. This script NEVER launches the production multi-sample
multi-prime reduction and NEVER modifies certificate or LF gates. Every verdict is a
statement about THIS row system and label box at the tested points: ``Obstructed`` means
"the projected target unit vector is outside the projected row span there" — never a
global impossibility claim.

Diagnostics deliberately reach one level below the public API (``_enumerate_labels`` /
``_generate_rows``) — read-only reuse, no reducer state is touched.

Exit codes: 0 = diagnostics completed and artifacts written (``Obstructed`` is a finding,
not a failure); 1 = internal inconsistency (span test vs coefficient extraction);
2 = usage/input problem.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import replace
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:  # allow running without installation
    sys.path.insert(0, str(REPO_ROOT / "src"))

from parametric_ibp_lf_reducer import ParserError, parse_family_text  # noqa: E402
from parametric_ibp_lf_reducer.adaptive import default_search_levels  # noqa: E402
from parametric_ibp_lf_reducer.api import (  # noqa: E402
    build_reducer_config,
    default_scattered_samples,
)
from parametric_ibp_lf_reducer.lf_feasibility import (  # noqa: E402
    STATUS_FEASIBLE,
    STATUS_OBSTRUCTED,
    feasibility_to_payload,
    lf_reduction_coefficients_mod_p,
    lf_reduction_feasible_mod_p,
)
from parametric_ibp_lf_reducer.reducer import _enumerate_labels, _generate_rows  # noqa: E402
from parametric_ibp_lf_reducer.valuations import (  # noqa: E402
    explain_local_finiteness,
    is_locally_finite,
)

DEFAULT_INPUT = REPO_ROOT / "examples" / "external_int2_dimensionless_input.wl.txt"
DEFAULT_OUT = REPO_ROOT / "validation" / "external_int2_method1.json"

VERDICT_FEASIBLE = "Feasible"
VERDICT_OBSTRUCTED = "Obstructed"
VERDICT_MIXED = "Mixed"
VERDICT_INCONCLUSIVE = "Inconclusive"

SCOPE_NOTE = (
    "Per-point statements about this row system and label box only; 'Obstructed' means the "
    "projected target unit vector is outside the projected row span at the tested "
    "(sample, prime) points — never a global impossibility claim. Certificate and LF gates "
    "are untouched."
)


def _verdict_str(verdict) -> object:
    """True/False stay booleans (JSON true/false); anything else becomes a string."""
    return verdict if isinstance(verdict, bool) else str(verdict)


def _serialize_ray_verdict(rv) -> dict:
    return {
        "direction": list(rv.ray.direction),
        "kind": rv.ray.kind,
        "score": None if rv.score is None else str(rv.score),
        "classification": rv.classification,
        "detail": rv.detail,
    }


def _serialize_shift(shift) -> dict:
    return {
        "shift": list(shift.shift),
        "deltas_on_failing": list(shift.deltas_on_failing),
        "improves_all": bool(shift.improves_all),
    }


def _serialize_report(report, max_shifts: int = 8) -> dict:
    """JSON-safe LocalFinitenessReport (failing rays in full, shifts capped)."""
    return {
        "label": list(report.label),
        "verdict": _verdict_str(report.verdict),
        "n_rays": len(report.rays),
        "n_failing_rays": len(report.failing_rays),
        "n_unknown_rays": len(report.unknown_rays),
        "bulk_safe": report.bulk_safe,
        "failing_rays": [_serialize_ray_verdict(rv) for rv in report.failing_rays],
        "n_recommended_shifts": len(report.recommended_shifts),
        "recommended_shifts": [
            _serialize_shift(s) for s in report.recommended_shifts[:max_shifts]
        ],
        "notes": list(report.notes),
    }


def _aggregate(statuses) -> str:
    """Fold per-point feasibility statuses into one honest verdict."""
    statuses = list(statuses)
    n_feasible = statuses.count(STATUS_FEASIBLE)
    n_obstructed = statuses.count(STATUS_OBSTRUCTED)
    if n_feasible and not n_obstructed:
        return VERDICT_FEASIBLE
    if n_obstructed and not n_feasible:
        return VERDICT_OBSTRUCTED
    if n_feasible and n_obstructed:
        return VERDICT_MIXED
    return VERDICT_INCONCLUSIVE


def _phase_a(family, target_label, labels, lf_map, *, audit_cap: int, audit_trials: int) -> dict:
    """Directional audit: verdict counts + full target report + capped non-LF details."""
    counts: dict[str, int] = {"True": 0, "False": 0, "Unknown": 0}
    for lab in labels:
        counts[str(lf_map[lab])] = counts.get(str(lf_map[lab]), 0) + 1
    target_report = explain_local_finiteness(family, target_label, random_trials=audit_trials)
    detailed = []
    n_non_lf = 0
    for lab in sorted(labels):
        if lab == target_label or lf_map[lab] is True:
            continue
        n_non_lf += 1
        if len(detailed) < audit_cap:
            detailed.append(
                _serialize_report(explain_local_finiteness(family, lab, random_trials=audit_trials))
            )
    return {
        "audit_trials": audit_trials,
        "lf_counts": counts,
        "target": _serialize_report(target_report),
        "n_non_lf_labels": n_non_lf,
        "detailed_non_lf_labels": detailed,
        "n_detailed_omitted": max(0, n_non_lf - len(detailed)),
    }


def _phase_b(rows, labels, target_label, lf_map, samples, primes) -> dict:
    """LF-constrained span feasibility at every (sample, prime) point."""
    points = []
    statuses = []
    for sample in samples:
        for prime in primes:
            res = lf_reduction_feasible_mod_p(rows, labels, target_label, lf_map, sample, prime)
            points.append(feasibility_to_payload(res))
            statuses.append(res.status)
    return {
        "n_points": len(points),
        "n_feasible": statuses.count(STATUS_FEASIBLE),
        "n_obstructed": statuses.count(STATUS_OBSTRUCTED),
        "n_bad_specialization": len(statuses)
        - statuses.count(STATUS_FEASIBLE)
        - statuses.count(STATUS_OBSTRUCTED),
        "verdict": _aggregate(statuses),
        "points": points,
    }


def _phase_c(rows, labels, target_label, lf_map, samples, primes, *, coeff_cap: int) -> dict:
    """Explicit LF-supported coefficients at each point (values stored when support is small)."""
    points = []
    supports = []
    for sample in samples:
        for prime in primes:
            res, coeffs = lf_reduction_coefficients_mod_p(
                rows, labels, target_label, lf_map, sample, prime
            )
            entry = feasibility_to_payload(res)
            if res.status == STATUS_FEASIBLE:
                support = sorted(coeffs)
                supports.append(tuple(support))
                entry["n_support"] = len(support)
                entry["support"] = [list(lab) for lab in support]
                if len(support) <= coeff_cap:
                    entry["coefficients"] = {
                        ",".join(map(str, lab)): int(coeffs[lab]) for lab in support
                    }
            points.append(entry)
    support_stable = len(set(supports)) <= 1
    return {
        "n_points": len(points),
        "n_feasible": len(supports),
        "support_stable": support_stable,
        "common_support_size": len(supports[0]) if supports and support_stable else None,
        "points": points,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="run_external_int2_method1",
        description="Method.1 diagnostics: directional LF audit + LF-constrained feasibility.",
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--phases", choices=("a", "ab", "abc"), default="abc")
    parser.add_argument("--level", choices=("A", "B"), default="A")
    parser.add_argument("--samples", type=int, default=3, help="scattered samples (cap)")
    parser.add_argument("--primes", type=int, default=2, help="primes from the document list (cap)")
    parser.add_argument("--audit-cap", type=int, default=24, help="detailed non-LF audits (cap)")
    parser.add_argument("--audit-trials", type=int, default=64)
    parser.add_argument("--coeff-cap", type=int, default=64, help="store values iff support <= cap")
    args = parser.parse_args(argv)

    if args.samples < 1 or args.primes < 1:
        print("error: --samples and --primes must be >= 1", file=sys.stderr)
        return 2
    try:
        text = args.input.read_text(encoding="utf-8")
    except OSError as exc:
        print(f"error: cannot read {args.input}: {exc}", file=sys.stderr)
        return 2
    try:
        family = parse_family_text(text)
        target_label, base_config = build_reducer_config(family)
        levels = default_search_levels(family, base_config)
    except (ParserError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    lvl = levels[0] if args.level == "A" else levels[1]
    config = replace(
        base_config,
        labels=None,
        label_box=lvl.label_box,
        max_ibp_degree=lvl.max_ibp_degree,
        tangent_degree_blocks=lvl.tangent_degree_blocks,
    )

    t0 = time.time()
    labels = list(_enumerate_labels(family, config))
    lf_map = {lab: is_locally_finite(family, lab) for lab in labels}
    if target_label not in lf_map:
        lf_map[target_label] = is_locally_finite(family, target_label)

    payload: dict = {
        "script": "run_external_int2_method1.py",
        "method": "Method.1: directional LF audit + LF-constrained feasibility mod p",
        "scope_note": SCOPE_NOTE,
        "input": str(args.input),
        "level": args.level,
        "level_name": lvl.name,
        "phases": args.phases,
        "target": list(target_label),
        "target_lf_verdict": _verdict_str(lf_map[target_label]),
        "label_box": lvl.label_box,
        "n_labels": len(labels),
    }
    payload["phase_a"] = _phase_a(
        family,
        target_label,
        labels,
        lf_map,
        audit_cap=args.audit_cap,
        audit_trials=args.audit_trials,
    )

    inconsistent = False
    if len(args.phases) >= 2:
        rows, row_diag = _generate_rows(family, labels, config, None)
        samples = default_scattered_samples(family.parameters, args.samples)
        primes = list(base_config.primes)[: args.primes]
        payload["n_rows"] = len(rows)
        payload["row_diagnostics"] = row_diag
        payload["samples"] = [sorted((k, str(v)) for k, v in s.items()) for s in samples]
        payload["primes"] = primes
        payload["phase_b"] = _phase_b(rows, labels, target_label, lf_map, samples, primes)
        if len(args.phases) >= 3:
            phase_c = _phase_c(
                rows, labels, target_label, lf_map, samples, primes, coeff_cap=args.coeff_cap
            )
            statuses_b = [p["status"] for p in payload["phase_b"]["points"]]
            statuses_c = [p["status"] for p in phase_c["points"]]
            phase_c["consistent_with_phase_b"] = statuses_b == statuses_c
            inconsistent = statuses_b != statuses_c
            payload["phase_c"] = phase_c
    payload["elapsed_seconds"] = round(time.time() - t0, 3)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")

    print(f"Method.1 diagnostics: level {args.level} ({lvl.name}), phases '{args.phases}'")
    print(f"  labels: {len(labels)}  LF counts: {payload['phase_a']['lf_counts']}")
    print(f"  target LF verdict: {payload['target_lf_verdict']}")
    if "phase_b" in payload:
        pb = payload["phase_b"]
        print(
            f"  feasibility: {pb['verdict']} "
            f"(Feasible {pb['n_feasible']}/{pb['n_points']}, Obstructed {pb['n_obstructed']})"
        )
    if "phase_c" in payload:
        pc = payload["phase_c"]
        print(
            f"  coefficients: support_stable={pc['support_stable']} "
            f"common_support_size={pc['common_support_size']}"
        )
    print(f"  wrote {args.out}")
    if inconsistent:
        print("error: span test and coefficient extraction disagree", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
