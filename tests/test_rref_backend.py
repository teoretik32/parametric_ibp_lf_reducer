# Perf.7: backend equivalence + counters for the sparse modular RREF kernel.
from __future__ import annotations

import json
import random
from fractions import Fraction

import pytest

from parametric_ibp_lf_reducer.sparse_rref import (
    DEFAULT_RREF_BACKEND,
    NUMBA_RREF_BACKEND,
    RREF_BACKENDS,
    rref_backend_available,
    rref_mod_p,
)

P = 2_147_483_629

# Backends that can actually run here; the numba one drops out cleanly when not installed.
AVAILABLE_BACKENDS = [b for b in RREF_BACKENDS if rref_backend_available(b)]


def _random_rows(rng, n_rows, cols, lo=2, hi=8):
    return [
        {c: rng.randrange(1, P) for c in rng.sample(cols, rng.randint(lo, min(hi, len(cols))))}
        for _ in range(n_rows)
    ]


def _label_cols(n):
    # Label tuples like the real ranking keys (family, exponent-vector index, shift).
    return [("T1", i, (i * 7 + 3) % 5) for i in range(n)]


def _assert_identical(a, b):
    assert a.pivots == b.pivots
    assert a.pivot_order == b.pivot_order
    assert a.free_cols == b.free_cols
    assert a.all_cols == b.all_cols
    assert a.rank == b.rank
    assert a == b  # dataclass equality (stats excluded via compare=False)


# ---------------------------------------------------------------- equivalence


def test_default_backend_is_dict():
    assert DEFAULT_RREF_BACKEND == "dict"
    assert "int_sparse_experimental" in RREF_BACKENDS
    assert NUMBA_RREF_BACKEND in RREF_BACKENDS
    assert rref_backend_available("dict") and rref_backend_available("int_sparse_experimental")
    assert not rref_backend_available("numpy")


def test_unknown_backend_rejected():
    with pytest.raises(ValueError, match="unknown rref_backend"):
        rref_mod_p([{0: 1}], P, backend="numpy")


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
def test_backends_identical_random_label_matrices(seed):
    rng = random.Random(seed)
    cols = _label_cols(60)
    rows = _random_rows(rng, 90, cols)
    order = sorted(cols, key=lambda c: (c[2], -c[1]))
    a = rref_mod_p(rows, P, column_order=order)
    b = rref_mod_p(rows, P, column_order=order, backend="int_sparse_experimental")
    _assert_identical(a, b)


def test_backends_identical_int_columns_and_partial_order():
    rng = random.Random(99)
    cols = list(range(40))
    rows = _random_rows(rng, 55, cols)
    # Partial column_order: unlisted columns must be pivoted last, in sorted order.
    order = [37, 5, 11, 2]
    a = rref_mod_p(rows, P, column_order=order)
    b = rref_mod_p(rows, P, column_order=order, backend="int_sparse_experimental")
    _assert_identical(a, b)


def test_backends_identical_rank_deficient_with_duplicate_rows():
    rng = random.Random(5)
    cols = _label_cols(25)
    base = _random_rows(rng, 12, cols)
    # Duplicates, scalar multiples, and zero rows: rank must stay <= 12.
    rows = base + [{c: (3 * v) % P for c, v in r.items()} for r in base] + [{}, {cols[0]: P}]
    a = rref_mod_p(rows, P, column_order=cols)
    b = rref_mod_p(rows, P, column_order=cols, backend="int_sparse_experimental")
    _assert_identical(a, b)
    assert a.rank <= 12


def test_backend_none_matches_default():
    rng = random.Random(11)
    cols = _label_cols(20)
    rows = _random_rows(rng, 30, cols)
    _assert_identical(rref_mod_p(rows, P), rref_mod_p(rows, P, backend="dict"))


# --------------------------------------------------------------- correctness


def test_pivot_rows_annihilate_original_system():
    # Every original row must reduce to zero against the RREF pivot rows.
    rng = random.Random(21)
    cols = _label_cols(30)
    rows = _random_rows(rng, 45, cols)
    for backend in AVAILABLE_BACKENDS:
        res = rref_mod_p(rows, P, column_order=cols, backend=backend)
        for row in rows:
            r = {c: v % P for c, v in row.items() if v % P}
            for pc in res.pivot_order:
                f = r.get(pc, 0)
                if f:
                    for c, v in res.pivots[pc].items():
                        nv = (r.get(c, 0) - f * v) % P
                        if nv:
                            r[c] = nv
                        else:
                            r.pop(c, None)
            assert r == {}, f"row not in span under backend {backend}"


def test_matches_fraction_rref_small():
    # 3x4 known system, checked against exact rational elimination.
    rows = [{0: 1, 1: 2, 2: 3}, {1: 1, 2: 1, 3: 4}, {0: 2, 1: 5, 2: 7, 3: 4}]
    exact = [[Fraction(x) for x in (1, 2, 3, 0)], [0, 1, 1, 4], [2, 5, 7, 4]]
    # Rational RREF (columns 0..3, natural order).
    exact = [[Fraction(x) for x in row] for row in exact]
    pr = 0
    for c in range(4):
        piv = next((i for i in range(pr, 3) if exact[i][c]), None)
        if piv is None:
            continue
        exact[pr], exact[piv] = exact[piv], exact[pr]
        exact[pr] = [x / exact[pr][c] for x in exact[pr]]
        for i in range(3):
            if i != pr and exact[i][c]:
                f = exact[i][c]
                exact[i] = [a - f * b for a, b in zip(exact[i], exact[pr])]
        pr += 1
    for backend in AVAILABLE_BACKENDS:
        res = rref_mod_p(rows, P, column_order=[0, 1, 2, 3], backend=backend)
        assert res.rank == pr
        for i, pc in enumerate(sorted(res.pivot_order)):
            got = res.pivots[pc]
            want = exact[i]
            for c in range(4):
                num, den = want[c].numerator, want[c].denominator
                assert got.get(c, 0) == (num * pow(den, P - 2, P)) % P


# ------------------------------------------------------------------ counters


def test_stats_absent_by_default_and_json_safe_when_requested():
    rng = random.Random(31)
    cols = _label_cols(15)
    rows = _random_rows(rng, 20, cols)
    assert rref_mod_p(rows, P).stats is None
    for backend in AVAILABLE_BACKENDS:
        res = rref_mod_p(rows, P, column_order=cols, backend=backend, collect_stats=True)
        s = res.stats
        assert s is not None and s["backend"] == backend
        json.dumps(s)  # JSON-safe
        assert s["n_rows"] == 20 and s["n_cols"] == len(res.all_cols)
        assert s["rank"] == res.rank == s["pivot_count"] == s["inversions"]
        assert s["nnz_initial"] >= s["row_nnz_before_max"] >= 1
        assert s["nnz_final"] == sum(len(r) for r in res.pivots.values())
        assert s["fill_in_ratio"] == pytest.approx(s["nnz_final"] / s["nnz_initial"])
        assert s["elimination_time_s"] >= 0.0


def test_stats_do_not_affect_result():
    rng = random.Random(41)
    cols = _label_cols(18)
    rows = _random_rows(rng, 25, cols)
    _assert_identical(
        rref_mod_p(rows, P, column_order=cols),
        rref_mod_p(rows, P, column_order=cols, collect_stats=True),
    )


def test_stats_empty_system():
    res = rref_mod_p([], P, collect_stats=True)
    assert res.rank == 0 and res.stats["nnz_initial"] == 0
    assert res.stats["fill_in_ratio"] is None
    json.dumps(res.stats)
