"""Tests for the minimal sparse modular RREF (Pass 2E)."""

from __future__ import annotations

import random

from parametric_ibp_lf_reducer import rref_mod_p


def _dense_rank_mod_p(rows, cols, p):
    """Reference rank via dense Gaussian elimination over GF(p)."""
    mat = [[r.get(c, 0) % p for c in cols] for r in rows]
    rank = 0
    ncols = len(cols)
    for col in range(ncols):
        piv = next((i for i in range(rank, len(mat)) if mat[i][col] % p != 0), None)
        if piv is None:
            continue
        mat[rank], mat[piv] = mat[piv], mat[rank]
        inv = pow(mat[rank][col], p - 2, p)
        mat[rank] = [(x * inv) % p for x in mat[rank]]
        for i in range(len(mat)):
            if i != rank and mat[i][col] % p != 0:
                f = mat[i][col]
                mat[i] = [(x - f * y) % p for x, y in zip(mat[i], mat[rank])]
        rank += 1
    return rank


def _reduce_row_by_pivots(row, pivots, p):
    """Subtract pivot rows until all pivot columns are cleared; return the residual."""
    r = {c: v % p for c, v in row.items() if v % p}
    for pc, prow in pivots.items():
        f = r.get(pc, 0)
        if f:
            for c, v in prow.items():
                nv = (r.get(c, 0) - f * v) % p
                if nv:
                    r[c] = nv
                elif c in r:
                    del r[c]
    return r


def test_rref_rank_matches_dense_reference():
    p = 10_007
    rng = random.Random(7)
    for _ in range(30):
        cols = list(range(6))
        rows = [{c: rng.randrange(0, p) for c in cols if rng.random() < 0.6} for _ in range(5)]
        res = rref_mod_p(rows, p, column_order=cols)
        assert res.rank == _dense_rank_mod_p(rows, cols, p)


def test_pivot_rows_are_fully_reduced():
    p = 97
    rows = [{0: 1, 1: 2, 2: 3}, {0: 2, 1: 1, 2: 0}, {1: 1, 2: 5}]
    res = rref_mod_p(rows, p, column_order=[0, 1, 2])
    # Each pivot column appears only in its own pivot row (RREF property).
    for pc, prow in res.pivots.items():
        assert prow[pc] == 1
        for other_pc, other in res.pivots.items():
            if other_pc != pc:
                assert other.get(pc, 0) == 0


def test_every_original_row_is_in_the_pivot_span():
    p = 101
    rng = random.Random(11)
    cols = list(range(5))
    rows = [{c: rng.randrange(0, p) for c in cols if rng.random() < 0.7} for _ in range(6)]
    res = rref_mod_p(rows, p, column_order=cols)
    for row in rows:
        assert _reduce_row_by_pivots(row, res.pivots, p) == {}  # reduces to zero


def test_column_order_controls_free_column():
    # r1: a + b = 0, r2: b + c = 0  (rank 2, 3 columns -> exactly one free column).
    p = 13
    rows = [{"a": 1, "b": 1}, {"b": 1, "c": 1}]
    # Prefer eliminating a, b -> c should be the free (master) column.
    res_c_free = rref_mod_p(rows, p, column_order=["a", "b", "c"])
    assert res_c_free.free_cols == ["c"]
    # Prefer eliminating b, c -> a should be the free column instead.
    res_a_free = rref_mod_p(rows, p, column_order=["c", "b", "a"])
    assert res_a_free.free_cols == ["a"]
    assert res_c_free.rank == res_a_free.rank == 2


def test_dependent_rows_do_not_inflate_rank():
    p = 1009
    base = {0: 1, 1: 2, 2: 3}
    rows = [base, {k: (2 * v) % p for k, v in base.items()}, {0: 0, 1: 1, 2: 1}]
    res = rref_mod_p(rows, p, column_order=[0, 1, 2])
    assert res.rank == 2  # the second row is 2x the first
