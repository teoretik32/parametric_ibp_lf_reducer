"""External Int2 Method.5: label-box geometry audit (per-variable n ranges).

Hypothesis under test (Method.4 follow-up): the generic ``Obstructed`` verdict is
caused by label-box *geometry* — too-narrow numerator (``n``) ranges truncating IBP
chains — rather than by row generation. Method.4 already showed that adding richer
tangent blocks ``(3,3)``/``(4,4)`` at fixed box does not flip any point.

Boxes (labels are SHIFTS against the Int2 base, offset convention):

* Box A: ``n2 in [-1, 1]``, ``n5, n7 in [0, 1]``, ``m0..m3 in [-3, 0]`` (3072 labels)
* Box B: ``n2 in [-1, 2]``, ``n5, n7 in [0, 1]``, ``m0..m3 in [-3, 0]`` (4096 labels)

Rows: baseline (coordinate IBP ``max_ibp_degree=2`` + tangent blocks
``((1,1),(2,2))``) **plus** richer tangent blocks ``(3,3)``/``(4,4)`` — i.e. the
enriched Method.4 row system, rebuilt on each box.

Genericity guard (binding): for each box and prime, a (sample, prime) point whose
rank falls below the maximal rank observed at that prime is classified as
``rank-deficient special`` and is EXCLUDED from the verdict. The verdict requires
at least ``--min-generic`` generic points, otherwise it is ``Inconclusive``.
Feasibility that appears only at rank-deficient special points is rejected.

Scope: read-only modular LF-feasibility diagnostics via
``lf_reduction_feasible_mod_p``; the span test is constrained to LF-True labels
only; no reconstruction is run; reducer core, certificates and LF gates are
untouched. ``Obstructed`` is an honest negative for this row system and label box
at the tested points — never a global impossibility claim.
"""

from __future__ import annotations

import argparse
import itertools
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
    lf_reduction_feasible_mod_p,
)
from parametric_ibp_lf_reducer.reducer import _enumerate_labels, _generate_rows  # noqa: E402
from parametric_ibp_lf_reducer.row_generation import generate_tangent_ibp_rows  # noqa: E402
from parametric_ibp_lf_reducer.tangent_fields import generate_tangent_fields  # noqa: E402
from parametric_ibp_lf_reducer.valuations import is_locally_finite  # noqa: E402

DEFAULT_INPUT = REPO_ROOT / "examples" / "external_int2_dimensionless_input.wl.txt"
DEFAULT_OUT = REPO_ROOT / "validation" / "external_int2_method5.json"

VERDICT_FEASIBLE = "Feasible"
VERDICT_OBSTRUCTED = "Obstructed"
VERDICT_MIXED = "Mixed"
VERDICT_INCONCLUSIVE = "Inconclusive"

RICHER_BLOCKS = ((3, 3), (4, 4))

# Per-variable numerator ranges (vars x2, x5, x7) + uniform poly ranges (G0..G3).
BOXES: dict[str, dict] = {
    "A": {"n_ranges": ((-1, 1), (0, 1), (0, 1)), "m_range": (-3, 0)},
    "B": {"n_ranges": ((-1, 2), (0, 1), (0, 1)), "m_range": (-3, 0)},
}

SCOPE_NOTE = (
    "Label-box geometry audit (Method.4 follow-up). Modular LF-feasibility only: the "
    "span test is constrained to LF-True labels; a span passing through divergent "
    "labels is never reported as an LF basis. Rank-deficient special points are "
    "excluded from verdicts (generic feasibility required). Per-point statements "
    "about this row system and label box only; 'Obstructed' is an honest negative at "
    "the tested (sample, prime) points — never a global impossibility claim. No "
    "reconstruction is run; reducer core, certificates and LF gates are untouched."
)


def build_box_labels(box: dict) -> list[tuple[int, ...]]:
    """Explicit per-variable label enumeration (config.labels wins over label_box)."""
    n_axes = [range(lo, hi + 1) for lo, hi in box["n_ranges"]]
    m_lo, m_hi = box["m_range"]
    m_axes = [range(m_lo, m_hi + 1)] * 4
    return [tuple(parts) for parts in itertools.product(*(n_axes + m_axes))]


def _analyze_points(points: list[dict], min_generic: int) -> dict:
    """Split points into generic vs rank-deficient special; verdict over generic only."""
    max_rank_by_prime: dict[int, int] = {}
    for p in points:
        r = p["result"]["rank"]
        prime = p["prime"]
        max_rank_by_prime[prime] = max(max_rank_by_prime.get(prime, 0), r)
    generic, special = [], []
    for p in points:
        deficit = max_rank_by_prime[p["prime"]] - p["result"]["rank"]
        p["rank_deficit"] = deficit
        p["classification"] = "generic" if deficit == 0 else "rank_deficient_special"
        (generic if deficit == 0 else special).append(p)

    statuses = [p["result"]["status"] for p in generic]
    real = [s for s in statuses if s in (STATUS_FEASIBLE, STATUS_OBSTRUCTED)]
    if len(real) < min_generic:
        verdict = VERDICT_INCONCLUSIVE
    elif all(s == STATUS_FEASIBLE for s in real):
        verdict = VERDICT_FEASIBLE
    elif all(s == STATUS_OBSTRUCTED for s in real):
        verdict = VERDICT_OBSTRUCTED
    else:
        verdict = VERDICT_MIXED

    special_only_feasible = (
        any(p["result"]["status"] == STATUS_FEASIBLE for p in special)
        and not any(s == STATUS_FEASIBLE for s in real)
    )

    supports = [tuple(map(tuple, p["result"]["residual_support"])) for p in generic]
    support_sets = [set(s) for s in supports] or [set()]
    common = set.intersection(*support_sets) if supports else set()
    return {
        "n_points": len(points),
        "n_generic": len(generic),
        "n_special": len(special),
        "max_rank_by_prime": {str(k): v for k, v in sorted(max_rank_by_prime.items())},
        "generic_verdict": verdict,
        "min_generic_required": min_generic,
        "special_only_feasible_rejected": special_only_feasible,
        "residual_support_stable": bool(supports) and len(set(supports)) == 1,
        "common_residual_support_size": len(common),
        "common_residual_support": sorted(list(lab) for lab in common),
    }


def run_box(family, base_config, deep, name: str, box: dict, richer_fields,
            samples, primes) -> dict:
    labels = build_box_labels(box)
    target_label = tuple([0] * (len(box["n_ranges"]) + 4))
    config = replace(
        base_config,
        labels=tuple(labels),
        label_box=None,
        max_ibp_degree=deep.max_ibp_degree,
        tangent_degree_blocks=deep.tangent_degree_blocks,
    )
    assert list(_enumerate_labels(family, config)) == labels

    t0 = time.time()
    lf_map = {lab: is_locally_finite(family, lab) for lab in labels}
    if target_label not in lf_map:
        lf_map[target_label] = is_locally_finite(family, target_label)

    base_rows, row_diag = _generate_rows(family, labels, config, None)
    extra = generate_tangent_ibp_rows(family, labels, richer_fields)
    rejected: dict[str, int] = {}
    for rej in extra.rejected:
        rejected[rej.reason] = rejected.get(rej.reason, 0) + 1

    seen = {row.dedup_key() for row in base_rows}
    merged = list(base_rows)
    n_dup = 0
    for row in extra.rows:
        key = row.dedup_key()
        if key in seen:
            n_dup += 1
            continue
        seen.add(key)
        merged.append(row)

    points = []
    for sample in samples:
        for prime in primes:
            res = lf_reduction_feasible_mod_p(
                merged, labels, target_label, lf_map, sample, prime
            )
            points.append(
                {
                    "prime": prime,
                    "sample": sorted((k, str(v)) for k, v in sample.items()),
                    "result": feasibility_to_payload(res),
                }
            )

    payload = {
        "box_name": name,
        "n_ranges": [list(r) for r in box["n_ranges"]],
        "m_range": list(box["m_range"]),
        "target": list(target_label),
        "target_lf_verdict": str(lf_map[target_label]),
        "n_labels": len(labels),
        "n_lf_true": sum(1 for lab in labels if lf_map[lab] is True),
        "n_base_rows": len(base_rows),
        "n_richer_rows_offered": len(extra.rows),
        "n_richer_rows_new": len(merged) - len(base_rows),
        "n_richer_rows_duplicate": n_dup,
        "n_rows_total": len(merged),
        "row_diagnostics": row_diag,
        "richer_row_rejections": rejected,
        "points": points,
        "elapsed_s": round(time.time() - t0, 3),
    }
    payload["analysis"] = _analyze_points(points, payload_min_generic(payload))
    return payload


_MIN_GENERIC = 3


def payload_min_generic(_payload: dict) -> int:
    return _MIN_GENERIC


def main(argv=None) -> int:
    global _MIN_GENERIC
    parser = argparse.ArgumentParser(
        prog="run_external_int2_method5",
        description="Method.5: label-box geometry audit (per-variable n ranges).",
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--boxes", type=str, default="A,B", help="comma list from {A,B}")
    parser.add_argument("--samples", type=int, default=4, help="scattered samples (cap)")
    parser.add_argument("--primes", type=int, default=2, help="primes from the document list")
    parser.add_argument(
        "--min-generic", type=int, default=3,
        help="minimum generic (sample, prime) points for a verdict",
    )
    args = parser.parse_args(argv)

    if args.samples < 1 or args.primes < 1 or args.min_generic < 1:
        print("error: --samples, --primes, --min-generic must be >= 1", file=sys.stderr)
        return 2
    box_names = [b.strip() for b in args.boxes.split(",") if b.strip()]
    unknown = [b for b in box_names if b not in BOXES]
    if unknown or not box_names:
        print(f"error: unknown boxes {unknown!r} (choose from {sorted(BOXES)})", file=sys.stderr)
        return 2
    _MIN_GENERIC = args.min_generic

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
    deep = levels[2]

    t0 = time.time()
    fields_t0 = time.time()
    richer_fields = generate_tangent_fields(family, list(RICHER_BLOCKS))
    fields_elapsed = round(time.time() - fields_t0, 3)

    samples = default_scattered_samples(family.parameters, args.samples)
    primes = list(base_config.primes)[: args.primes]

    payload: dict = {
        "script": "run_external_int2_method5.py",
        "method": "Method.5: label-box geometry audit (per-variable n ranges)",
        "scope_note": SCOPE_NOTE,
        "input": str(args.input),
        "baseline_tangent_blocks": [list(b) for b in deep.tangent_degree_blocks],
        "richer_blocks": [list(b) for b in RICHER_BLOCKS],
        "n_richer_fields": len(richer_fields),
        "richer_fields_elapsed_s": fields_elapsed,
        "samples": [sorted((k, str(v)) for k, v in s.items()) for s in samples],
        "primes": primes,
        "min_generic": args.min_generic,
        "boxes": {},
    }

    for name in box_names:
        print(f"[method5] box {name} ...", flush=True)
        payload["boxes"][name] = run_box(
            family, base_config, deep, name, BOXES[name], richer_fields, samples, primes
        )
        b = payload["boxes"][name]
        a = b["analysis"]
        print(
            f"[method5] box {name}: labels={b['n_labels']} rows={b['n_rows_total']} "
            f"generic_verdict={a['generic_verdict']} "
            f"(generic {a['n_generic']}/{a['n_points']}, special {a['n_special']}) "
            f"elapsed={b['elapsed_s']}s",
            flush=True,
        )

    payload["elapsed_s"] = round(time.time() - t0, 3)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
