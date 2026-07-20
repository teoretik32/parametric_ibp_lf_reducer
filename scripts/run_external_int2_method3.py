"""External Int2 Method.3 diagnostics: composite locally-finite master feasibility.

Given the non-locally-finite terms of the certified Method.2 corrected normal form,
build a deterministic candidate pool (integer label shifts + finite numerators) and
solve for "composite masters" M = sum_i c_i * J(label_i) whose bad Laurent layers
cancel identically on the primary ray, then verify every kernel vector on all
candidate rays plus a deterministic random safety net
(``composite_master_feasibility``).

Coefficients live in the field of rational functions of the non-regulator
parameters; a fixed-sample rank cross-check guards against special-locus artifacts
(``BadSpecialization``).  ``FeasibleCompositeBasis`` means "these explicit integer
linear combinations are locally finite on every checked ray" — a statement about
THIS pool and ray set, never a completeness claim about all composites.

Caps by design: the pool is bounded by ``--max-pool`` (default 512) and the safety
net by ``--random-trials`` (default 64, cap 256).  This script never modifies
reducer state, certificates, or LF gates.

Exit codes: 0 = diagnostics completed and artifacts written (``NoComposite...`` is
a finding, not a failure); 2 = usage/input problem or cap violation.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:  # allow running without installation
    sys.path.insert(0, str(REPO_ROOT / "src"))

from parametric_ibp_lf_reducer import ParserError, parse_family_text  # noqa: E402
from parametric_ibp_lf_reducer.composite_masters import (  # noqa: E402
    build_candidate_pool,
    composite_master_feasibility,
    feasibility_to_payload,
)

DEFAULT_INPUT = REPO_ROOT / "examples" / "external_int2_dimensionless_input.wl.txt"
DEFAULT_JSON = REPO_ROOT / "validation" / "external_int2_composite_feasibility.json"

# Non-LF terms of the certified Method.2 corrected normal form (see
# notes/EXTERNAL_INT2_AUDIT.md and validation/external_int2_corrected_reduction.json).
NF_LABELS = (
    (-1, 0, 0, -1, -1, 0, 0),
    (-1, 0, 0, 0, -1, 0, -1),
    (0, 0, 0, -1, 0, 0, -1),
    (0, 0, 0, 0, -1, 0, 0),
    (0, 0, 0, 0, 0, -1, 0),
    (0, 0, 1, -1, 0, 0, -1),
)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--primary-ray", default="-1,0,0", help="comma-separated ray")
    parser.add_argument("--numerator-degree", type=int, default=2)
    parser.add_argument("--random-trials", type=int, default=64)
    parser.add_argument("--max-pool", type=int, default=512)
    args = parser.parse_args(argv)

    if not (0 < args.random_trials <= 256):
        print("cap violation: --random-trials must be in 1..256", file=sys.stderr)
        return 2
    try:
        primary_ray = tuple(int(t) for t in args.primary_ray.split(","))
    except ValueError:
        print(f"bad --primary-ray: {args.primary_ray!r}", file=sys.stderr)
        return 2
    try:
        family = parse_family_text(args.input.read_text(encoding="utf-8"))
    except (OSError, ParserError) as exc:
        print(f"cannot load family: {exc}", file=sys.stderr)
        return 2

    t0 = time.time()
    pool = build_candidate_pool(
        family,
        list(NF_LABELS),
        var_shift_axes=("x2",),
        poly_shift_axes=("G0", "G3"),
        shift_depths=(-1, -2),
        numerator_vars=("x5", "x7"),
        numerator_degree=args.numerator_degree,
    )
    if len(pool) > args.max_pool:
        print(f"cap violation: pool size {len(pool)} > --max-pool {args.max_pool}", file=sys.stderr)
        return 2
    result = composite_master_feasibility(
        family, pool, primary_ray, random_trials=args.random_trials
    )
    elapsed = time.time() - t0

    payload = {
        "script": Path(__file__).name,
        "input": str(args.input),
        "nf_labels": [list(label) for label in NF_LABELS],
        "elapsed_seconds": round(elapsed, 2),
        "feasibility": feasibility_to_payload(result),
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8"
    )

    print(
        f"status={result.status} pool={result.pool_size} "
        f"participants={len(result.participants)} kernel={result.kernel_dimension} "
        f"full={result.full_dimension} checked_rays={result.checked_rays} "
        f"elapsed={elapsed:.1f}s -> {args.json}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
