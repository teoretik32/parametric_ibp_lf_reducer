# Perf.8: profile the RREF kernel on a REAL matrix from corrected Example 4*.
#
# Builds ONE representative row system from the explicit family in
# examples/example4_star_corrected_input.wl.txt (real rows, real ranking order,
# first configured prime/sample), assembles it mod p, and runs rref_mod_p with
# collect_stats=True for each requested backend. Prints one JSON object.
#
# Default mode is "medium": every label-box range is shaved by one unit off the
# top, so the matrix is real (same family, same row generators, same ranking)
# but small enough to profile in minutes. Pass --full for the untouched 972-label
# box (the full RREF alone is tens of minutes on this machine). The mode is
# recorded in the JSON so numbers are never mislabeled.
#
# Read-only: no repo artifacts are written unless --out is given. No stdout
# except the final JSON. Deterministic (no RNG anywhere on this path).
from __future__ import annotations

import argparse
import dataclasses
import json
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:  # allow running without installation
    sys.path.insert(0, str(REPO_ROOT / "src"))

from parametric_ibp_lf_reducer import parse_family_text  # noqa: E402
from parametric_ibp_lf_reducer.api import build_reducer_config  # noqa: E402
from parametric_ibp_lf_reducer.modular_normal_form import assemble_rows_mod_p  # noqa: E402
from parametric_ibp_lf_reducer.ranking import rank_labels  # noqa: E402
from parametric_ibp_lf_reducer.reducer import _enumerate_labels, _generate_rows  # noqa: E402
from parametric_ibp_lf_reducer.sparse_rref import RREF_BACKENDS, rref_mod_p  # noqa: E402

INPUT_PATH = REPO_ROOT / "examples" / "example4_star_corrected_input.wl.txt"


def _shrink(rng):
    """Shave one unit off the top of every (lo, hi) range; handles nested range lists."""
    if isinstance(rng, (list, tuple)) and len(rng) == 2 and all(isinstance(x, int) for x in rng):
        lo, hi = rng
        return (lo, max(lo, hi - 1))
    return tuple(_shrink(r) for r in rng)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--full", action="store_true", help="use the untouched label box (slow)")
    ap.add_argument(
        "--backends",
        default=",".join(RREF_BACKENDS),
        help="comma-separated backend list (default: all registered)",
    )
    ap.add_argument("--out", default=None, help="also write the JSON to this path")
    args = ap.parse_args()
    backends = [b.strip() for b in args.backends.split(",") if b.strip()]

    family = parse_family_text(INPUT_PATH.read_text(encoding="utf-8"))
    _, config = build_reducer_config(family)
    if config.label_box is None:
        raise SystemExit("document did not configure a label box")

    mode = "full-real" if args.full else "medium-real (box shaved by 1 per range)"
    n_range, m_range = config.label_box
    if not args.full:
        n_range, m_range = _shrink(n_range), _shrink(m_range)
    cfg = dataclasses.replace(config, labels=None, label_box=(n_range, m_range))

    labels = _enumerate_labels(family, cfg)
    t0 = time.perf_counter()
    rows, row_diag = _generate_rows(family, labels, cfg)
    t_rows = time.perf_counter() - t0

    t0 = time.perf_counter()
    ranked = rank_labels(family, labels)
    t_rank = time.perf_counter() - t0

    prime = config.primes[0]
    sample = dict(config.samples[0])
    t0 = time.perf_counter()
    matrix = assemble_rows_mod_p(family, rows, sample, prime)
    t_asm = time.perf_counter() - t0

    per_backend: dict = {}
    ref = None
    for backend in backends:
        res = rref_mod_p(
            [dict(r) for r in matrix],  # fresh copies: backends must not share row dicts
            prime,
            column_order=ranked.ordered,
            backend=backend,
            collect_stats=True,
        )
        per_backend[backend] = res.stats
        if ref is None:
            ref = res
        elif (res.pivots, res.pivot_order, res.free_cols) != (
            ref.pivots,
            ref.pivot_order,
            ref.free_cols,
        ):
            print(f"BACKEND MISMATCH: {backend}", file=sys.stderr)
            return 1

    payload = {
        "mode": mode,
        "input": str(INPUT_PATH.relative_to(REPO_ROOT)),
        "label_box_used": [n_range, m_range],
        "n_labels": len(labels),
        "rows_by_kind": row_diag,
        "row_generation_time_s": round(t_rows, 3),
        "ranking_time_s": round(t_rank, 3),
        "assemble_time_s": round(t_asm, 3),
        "prime": prime,
        "sample": {str(k): str(v) for k, v in sample.items()},
        "backends": per_backend,
    }
    text = json.dumps(payload, indent=2, default=str)
    print(text)
    if args.out:
        Path(args.out).write_text(text + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
