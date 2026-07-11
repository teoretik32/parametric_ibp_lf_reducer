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

from parametric_ibp_lf_reducer.sparse_rref import (  # noqa: E402
    RREF_BACKENDS,
    rref_backend_available,
    rref_mod_p,
)

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

    backends = [b for b in RREF_BACKENDS if rref_backend_available(b)]
    # Warm up once per backend so one-time costs (numba JIT compile) are not billed
    # to the first timed shape.
    for backend in backends:
        rref_mod_p([{0: 1, 1: 2}, {1: 3}], P, column_order=[0, 1], backend=backend)

    speedups = [f"{b} speedup" for b in backends[1:]]
    print(f"| Matrix (rows x cols, ~nnz/row) | {' | '.join(backends + speedups)} |")
    print("|---|" + "---|" * (len(backends) + len(speedups)))
    for name, n_rows, n_cols, lo, hi in shapes:
        rows, order = build(2026, n_rows, n_cols, lo, hi)
        times = {}
        ref = None
        for backend in backends:
            t0 = time.perf_counter()
            res = rref_mod_p(rows, P, column_order=order, backend=backend)
            times[backend] = time.perf_counter() - t0
            if ref is None:
                ref = res
            elif res != ref:
                print(f"BACKEND MISMATCH on {name}", file=sys.stderr)
                return 1
        base = times[backends[0]]
        cells = " | ".join(f"{times[b]:.2f}s" for b in backends)
        ratios = " | ".join(
            f"{(base / times[b]):.2f}x" if times[b] else "inf" for b in backends[1:]
        )
        print(f"| {name} (rank {ref.rank}) | {cells} | {ratios} |")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
