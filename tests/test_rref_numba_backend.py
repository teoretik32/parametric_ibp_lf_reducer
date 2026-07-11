# Perf.10: numba int-array RREF backend — equivalence vs the dict backend + guards.
#
# The whole module skips cleanly when numba is not installed; availability of the
# backend itself is covered (backend-list membership, clear error paths) in
# test_rref_backend.py, which does not require numba.
from __future__ import annotations

import random

import pytest

pytest.importorskip("numba")

from parametric_ibp_lf_reducer.sparse_rref import (
    NUMBA_RREF_BACKEND,
    rref_backend_available,
    rref_mod_p,
)

P = 2_147_483_629


def _random_rows(rng, n_rows, cols, lo=2, hi=8):
    return [
        {c: rng.randrange(1, P) for c in rng.sample(cols, rng.randint(lo, min(hi, len(cols))))}
        for _ in range(n_rows)
    ]


def _label_cols(n):
    return [("T1", i, (i * 7 + 3) % 5) for i in range(n)]


def _assert_identical(a, b):
    assert a.pivots == b.pivots
    assert a.pivot_order == b.pivot_order
    assert a.free_cols == b.free_cols
    assert a.all_cols == b.all_cols
    assert a.rank == b.rank
    assert a == b


def test_backend_reports_available():
    assert rref_backend_available(NUMBA_RREF_BACKEND)


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4, 5, 6])
def test_identical_to_dict_on_random_label_matrices(seed):
    rng = random.Random(seed)
    cols = _label_cols(60)
    rows = _random_rows(rng, 90, cols)
    order = sorted(cols, key=lambda c: (c[2], -c[1]))
    a = rref_mod_p(rows, P, column_order=order)
    b = rref_mod_p(rows, P, column_order=order, backend=NUMBA_RREF_BACKEND)
    _assert_identical(a, b)


def test_identical_on_dense_high_fill_in_matrix():
    # Denser rows force heavy fill-in during elimination (the axpy/merge hot path).
    rng = random.Random(17)
    cols = list(range(40))
    rows = _random_rows(rng, 80, cols, lo=10, hi=16)
    a = rref_mod_p(rows, P, column_order=cols)
    b = rref_mod_p(rows, P, column_order=cols, backend=NUMBA_RREF_BACKEND)
    _assert_identical(a, b)
    assert a.rank == 40  # dense random system: full column rank, lots of fill-in


def test_identical_rank_deficient_duplicates_and_zero_rows():
    rng = random.Random(5)
    cols = _label_cols(25)
    base = _random_rows(rng, 12, cols)
    rows = base + [{c: (3 * v) % P for c, v in r.items()} for r in base] + [{}, {cols[0]: P}]
    a = rref_mod_p(rows, P, column_order=cols)
    b = rref_mod_p(rows, P, column_order=cols, backend=NUMBA_RREF_BACKEND)
    _assert_identical(a, b)
    assert a.rank <= 12


def test_identical_partial_column_order():
    rng = random.Random(99)
    cols = list(range(40))
    rows = _random_rows(rng, 55, cols)
    order = [37, 5, 11, 2]  # unlisted columns must be pivoted last, in sorted order
    a = rref_mod_p(rows, P, column_order=order)
    b = rref_mod_p(rows, P, column_order=order, backend=NUMBA_RREF_BACKEND)
    _assert_identical(a, b)


def test_small_primes_and_tiny_systems():
    for prime in (5, 97, 2_147_483_647):  # 2**31 - 1 is the largest accepted prime
        rows = [{0: 1, 1: prime - 1}, {1: 2, 2: 3}, {0: 1, 1: 1, 2: 3}]
        a = rref_mod_p(rows, prime, column_order=[0, 1, 2])
        b = rref_mod_p(rows, prime, column_order=[0, 1, 2], backend=NUMBA_RREF_BACKEND)
        _assert_identical(a, b)
    # Empty system.
    _assert_identical(rref_mod_p([], P), rref_mod_p([], P, backend=NUMBA_RREF_BACKEND))


def test_prime_at_or_above_2_pow_31_is_rejected():
    for prime in (1 << 31, (1 << 31) + 11):
        with pytest.raises(ValueError, match="2\\*\\*31"):
            rref_mod_p([{0: 1}], prime, backend=NUMBA_RREF_BACKEND)


def test_stats_collected_and_json_safe():
    import json

    rng = random.Random(31)
    cols = _label_cols(15)
    rows = _random_rows(rng, 20, cols)
    a = rref_mod_p(rows, P, column_order=cols, collect_stats=True)
    b = rref_mod_p(rows, P, column_order=cols, backend=NUMBA_RREF_BACKEND, collect_stats=True)
    _assert_identical(a, b)
    s = b.stats
    assert s is not None and s["backend"] == NUMBA_RREF_BACKEND
    json.dumps(s)
    # Counters that must not depend on the representation:
    for key in ("n_rows", "n_cols", "nnz_initial", "nnz_final", "rank", "inversions"):
        assert s[key] == a.stats[key], key
    # Result values must be plain Python ints (no numpy leakage into pivots).
    pc = b.pivot_order[0]
    assert all(type(v) is int for v in b.pivots[pc].values())
    assert all(type(c) is type(a_c) for c, a_c in zip(b.pivot_order, a.pivot_order))


# ------------------------------------------------------------------ integration


ONE_VAR = """
IBPInput = <|
  "Variables" -> {x}, "Parameters" -> {ep}, "Regulators" -> {ep},
  "Polynomials" -> <| "G0" -> 1 + x |>,
  "MonomialExponents" -> <| x -> 0 |>,
  "PolynomialExponents" -> <| "G0" -> -2 |>
|>
"""


def test_modular_normal_form_identical_under_numba_backend(monkeypatch):
    # Same tiny reduction as test_modular_normal_form, once per backend, by swapping
    # the module-level default (rref_mod_p reads it at call time).
    import parametric_ibp_lf_reducer.sparse_rref as sr
    from parametric_ibp_lf_reducer import (
        algebraic_row,
        modular_normal_form,
        parse_family_text,
        zero_label,
    )

    fam = parse_family_text(ONE_VAR)
    row = algebraic_row(fam, zero_label(1, 1), 0)
    base = modular_normal_form(fam, [row], (0, 0), {"ep": 5}, 2_147_483_629)
    monkeypatch.setattr(sr, "DEFAULT_RREF_BACKEND", NUMBA_RREF_BACKEND)
    alt = modular_normal_form(fam, [row], (0, 0), {"ep": 5}, 2_147_483_629)
    assert alt.status == base.status == "Reduced"
    assert alt.pivot_label == base.pivot_label
    assert alt.terms == base.terms == {(0, -1): 1, (1, -1): 1}
    assert alt.formal_success is base.formal_success is True
