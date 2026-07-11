# Perf.7: benchmark the RREF backends on synthetic label-tuple-keyed matrices.
#
# Shapes are chosen to mimic the real reducer runs (sparse rows, label-tuple
# columns, prime ~2^31). Single process, wall clock only, no stdout except the
# final table. Usage:  python scripts/bench_rref_backends.py [--fast]
from __future__ import annotations

import argparse
import random
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from parametric_ibp_lf_reducer.sparse_rref import RREF_BACKENDS, rref_mod_p  # noqa: E402

P = 2_147_483_629

SHAPES = [
    # (name, n_rows, n_cols, nnz_lo, nnz_hi)
    ("tiny 200x150, ~6", 200, 150, 4, 8),
    ("medium 1000x800, ~8", 1000, 800, 6, 10),
    ("D4-like 3000x2400, ~8", 3000, 2400, 6, 10),
]


def build(seed: int, n_rows: int, n_cols: int, lo: int, hi: int):
    rng = random.Random(seed)
    cols = [("T1", i, (i * 7 + 3) % 5) for i in range(n_cols)]
    rows = [
        {c: rng.randrange(1, P) for c in rng.sample(cols, rng.randint(lo, hi))}
        for _ in range(n_rows)
    ]
    order = sorted(cols, key=lambda c: (c[2], -c[1]))
    return rows, order


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true", help="skip the largest shape")
    args = ap.parse_args()
    shapes = SHAPES[:-1] if args.fast else SHAPES

    print(f"| Matrix (rows x cols, ~nnz/row) | {' | '.join(RREF_BACKENDS)} | ratio |")
    print("|---|" + "---|" * (len(RREF_BACKENDS) + 1))
    for name, n_rows, n_cols, lo, hi in shapes:
        rows, order = build(2026, n_rows, n_cols, lo, hi)
        times = {}
        ref = None
        for backend in RREF_BACKENDS:
            t0 = time.perf_counter()
            res = rref_mod_p(rows, P, column_order=order, backend=backend)
            times[backend] = time.perf_counter() - t0
            if ref is None:
                ref = res
            elif res != ref:
                print(f"BACKEND MISMATCH on {name}", file=sys.stderr)
                return 1
        base = times[RREF_BACKENDS[0]]
        cells = " | ".join(f"{times[b]:.2f}s" for b in RREF_BACKENDS)
        ratio = times[RREF_BACKENDS[1]] / base if base else float("nan")
        print(f"| {name} (rank {ref.rank}) | {cells} | {ratio:.2f}x |")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
