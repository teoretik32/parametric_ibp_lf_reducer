"""External Int2 Method.4: same-dimension LF-basis completeness audit with richer rows.

Phases (``--phases {s, sd}``):

* ``s`` — setup feasibility for richer tangent (syzygy) degree blocks: time
  ``generate_tangent_fields`` at each requested block (default ``(3,3)`` then ``(4,4)``)
  under a soft time budget (``--field-budget``). The symbolic nullspace
  (``sympy.linear_eq_to_matrix`` + ``Matrix.nullspace``) is the suspected bottleneck, so
  it is measured before any audit work; blocks are skipped (recorded, not silently
  dropped) once the budget is exhausted.
* ``d`` — completeness audit on the deep-level label box: baseline rows
  (coordinate IBP ``max_ibp_degree=2`` + tangent blocks ``((1,1),(2,2))``) versus the
  enriched system (baseline + surface-filtered tangent IBP rows from the richer fields),
  compared point-by-point via ``lf_reduction_feasible_mod_p`` at the same scattered
  ``(sample, prime)`` points.

Scope correction (binding): the span test is constrained to LF-True labels ONLY.
A combination that passes through divergent labels — i.e. a sum of separately divergent
integrals — is NEVER reported as an LF basis. ``Obstructed`` outcomes are reported as
honest negatives for this row system and label box; they are not global impossibility
claims and are not reinterpreted.

Reuses Method.1 scaffolding read-only (``_enumerate_labels``, ``_generate_rows``);
no reducer state is touched.
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
    lf_reduction_feasible_mod_p,
)
from parametric_ibp_lf_reducer.reducer import _enumerate_labels, _generate_rows  # noqa: E402
from parametric_ibp_lf_reducer.row_generation import generate_tangent_ibp_rows  # noqa: E402
from parametric_ibp_lf_reducer.tangent_fields import generate_tangent_fields  # noqa: E402
from parametric_ibp_lf_reducer.valuations import is_locally_finite  # noqa: E402

DEFAULT_INPUT = REPO_ROOT / "examples" / "external_int2_dimensionless_input.wl.txt"
DEFAULT_OUT = REPO_ROOT / "validation" / "external_int2_method4.json"

VERDICT_FEASIBLE = "Feasible"
VERDICT_OBSTRUCTED = "Obstructed"
VERDICT_MIXED = "Mixed"
VERDICT_INCONCLUSIVE = "Inconclusive"

SCOPE_NOTE = (
    "Same-dimension LF-basis completeness audit. Reduction is constrained to LF-True "
    "labels only; a span passing through divergent labels (a sum of separately divergent "
    "integrals) is never reported as an LF basis. Per-point statements about this row "
    "system and label box only; 'Obstructed' is an honest negative at the tested "
    "(sample, prime) points — never a global impossibility claim. Certificate and LF "
    "gates are untouched."
)


def _parse_blocks(text: str) -> list[tuple[int, int]]:
    """Parse ``"3,3;4,4"`` into ``[(3, 3), (4, 4)]``."""
    blocks: list[tuple[int, int]] = []
    for part in text.split(";"):
        part = part.strip()
        if not part:
            continue
        pieces = part.split(",")
        if len(pieces) != 2:
            raise ValueError(f"bad block {part!r}: expected 'd_q,d_h'")
        d_q, d_h = int(pieces[0]), int(pieces[1])
        if d_q < 0 or d_h < 0:
            raise ValueError(f"bad block {part!r}: degrees must be >= 0")
        blocks.append((d_q, d_h))
    if not blocks:
        raise ValueError("no degree blocks given")
    return blocks


def _phase_s(family, blocks, budget_s: float) -> tuple[dict, list]:
    """Time tangent-field generation per block; stop starting new blocks once over budget."""
    entries = []
    fields = []
    spent = 0.0
    for block in blocks:
        if spent > budget_s:
            entries.append(
                {
                    "block": list(block),
                    "status": "skipped_budget",
                    "budget_s": budget_s,
                    "spent_s": round(spent, 3),
                }
            )
            continue
        t0 = time.time()
        found = generate_tangent_fields(family, [block])
        dt = time.time() - t0
        spent += dt
        entries.append(
            {
                "block": list(block),
                "status": "ok",
                "n_fields": len(found),
                "elapsed_s": round(dt, 3),
            }
        )
        fields.extend(found)
    payload = {
        "budget_s": budget_s,
        "total_elapsed_s": round(spent, 3),
        "blocks": entries,
        "n_richer_fields": len(fields),
    }
    return payload, fields


def _phase_d(family, labels, target_label, lf_map, base_rows, extra_rows, samples, primes) -> dict:
    """Baseline vs enriched LF-span feasibility at identical (sample, prime) points."""
    seen = {row.dedup_key() for row in base_rows}
    merged = list(base_rows)
    n_dup = 0
    for row in extra_rows:
        key = row.dedup_key()
        if key in seen:
            n_dup += 1
            continue
        seen.add(key)
        merged.append(row)

    points = []
    flipped = 0
    for sample in samples:
        for prime in primes:
            base = lf_reduction_feasible_mod_p(
                base_rows, labels, target_label, lf_map, sample, prime
            )
            rich = lf_reduction_feasible_mod_p(
                merged, labels, target_label, lf_map, sample, prime
            )
            if base.status == STATUS_OBSTRUCTED and rich.status == STATUS_FEASIBLE:
                flipped += 1
            points.append(
                {
                    "prime": prime,
                    "sample": sorted((k, str(v)) for k, v in sample.items()),
                    "baseline": feasibility_to_payload(base),
                    "enriched": feasibility_to_payload(rich),
                }
            )

    def _aggregate(key: str) -> str:
        statuses = [p[key]["status"] for p in points]
        real = [s for s in statuses if s in (STATUS_FEASIBLE, STATUS_OBSTRUCTED)]
        if not real:
            return VERDICT_INCONCLUSIVE
        if all(s == STATUS_FEASIBLE for s in real):
            return VERDICT_FEASIBLE
        if all(s == STATUS_OBSTRUCTED for s in real):
            return VERDICT_OBSTRUCTED
        return VERDICT_MIXED

    return {
        "n_base_rows": len(base_rows),
        "n_extra_rows_offered": len(extra_rows),
        "n_extra_rows_new": len(merged) - len(base_rows),
        "n_extra_rows_duplicate": n_dup,
        "points": points,
        "baseline_verdict": _aggregate("baseline"),
        "enriched_verdict": _aggregate("enriched"),
        "n_points_flipped_to_feasible": flipped,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="run_external_int2_method4",
        description="Method.4: same-dimension LF-basis completeness audit with richer rows.",
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--phases", choices=("s", "sd"), default="sd")
    parser.add_argument("--blocks", type=str, default="3,3;4,4", help="richer blocks 'dq,dh;...'")
    parser.add_argument(
        "--field-budget", type=float, default=600.0, help="soft budget (s) for phase s"
    )
    parser.add_argument("--samples", type=int, default=3, help="scattered samples (cap)")
    parser.add_argument("--primes", type=int, default=2, help="primes from the document list (cap)")
    args = parser.parse_args(argv)

    if args.samples < 1 or args.primes < 1:
        print("error: --samples and --primes must be >= 1", file=sys.stderr)
        return 2
    try:
        blocks = _parse_blocks(args.blocks)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
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

    deep = levels[2]
    config = replace(
        base_config,
        labels=None,
        label_box=deep.label_box,
        max_ibp_degree=deep.max_ibp_degree,
        tangent_degree_blocks=deep.tangent_degree_blocks,
    )

    t0 = time.time()
    payload: dict = {
        "script": "run_external_int2_method4.py",
        "method": "Method.4: same-dimension LF-basis completeness audit (richer rows)",
        "scope_note": SCOPE_NOTE,
        "input": str(args.input),
        "level_name": deep.name,
        "label_box": deep.label_box,
        "baseline_tangent_blocks": [list(b) for b in deep.tangent_degree_blocks],
        "richer_blocks_requested": [list(b) for b in blocks],
        "phases": args.phases,
        "target": list(target_label),
    }

    phase_s, richer_fields = _phase_s(family, blocks, args.field_budget)
    payload["phase_s"] = phase_s

    if "d" in args.phases:
        labels = list(_enumerate_labels(family, config))
        lf_map = {lab: is_locally_finite(family, lab) for lab in labels}
        if target_label not in lf_map:
            lf_map[target_label] = is_locally_finite(family, target_label)
        payload["n_labels"] = len(labels)
        payload["n_lf_true"] = sum(1 for lab in labels if lf_map[lab] is True)
        payload["target_lf_verdict"] = str(lf_map[target_label])

        base_rows, row_diag = _generate_rows(family, labels, config, None)
        payload["row_diagnostics"] = row_diag

        extra = generate_tangent_ibp_rows(family, labels, richer_fields)
        rejected: dict[str, int] = {}
        for rej in extra.rejected:
            rejected[rej.reason] = rejected.get(rej.reason, 0) + 1
        payload["richer_row_rejections"] = rejected

        samples = default_scattered_samples(family.parameters, args.samples)
        primes = list(base_config.primes)[: args.primes]
        payload["samples"] = [sorted((k, str(v)) for k, v in s.items()) for s in samples]
        payload["primes"] = primes
        payload["phase_d"] = _phase_d(
            family, labels, target_label, lf_map, base_rows, extra.rows, samples, primes
        )

    payload["elapsed_s"] = round(time.time() - t0, 3)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    print(f"wrote {args.out}")
    if "d" in args.phases:
        d = payload["phase_d"]
        print(
            f"baseline={d['baseline_verdict']} enriched={d['enriched_verdict']} "
            f"new_rows={d['n_extra_rows_new']} flipped={d['n_points_flipped_to_feasible']}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
