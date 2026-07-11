# Perf.12: "auto" rref_backend — heuristic selection tests (dict vs numba).
#
# select_rref_backend is a pure function, so the whole decision table is tested
# here without numba installed (availability is passed in explicitly). The one
# end-to-end numba test at the bottom skips cleanly when numba is missing.
# Selection never changes results — equivalence is enforced in
# test_rref_backend.py / test_rref_numba_backend.py; this file only checks the
# routing and its diagnostics.
from __future__ import annotations

import random

import pytest

from parametric_ibp_lf_reducer.sparse_rref import (
    AUTO_RREF_BACKEND,
    AUTO_RREF_THRESHOLDS,
    DEFAULT_RREF_BACKEND,
    NUMBA_RREF_BACKEND,
    RREF_BACKEND_CHOICES,
    RREF_BACKENDS,
    BackendUnavailable,
    rref_backend_available,
    rref_mod_p,
    select_rref_backend,
)

P = 2_147_483_629  # < 2**31: passes the int64 guard
BIG_PRIME = 2_305_843_009_213_693_951  # >= 2**31: fails the int64 guard

# Sizes that clear every auto gate (thresholds are intentionally conservative).
BIG = dict(
    n_rows=AUTO_RREF_THRESHOLDS["min_rows"],
    n_cols=AUTO_RREF_THRESHOLDS["min_cols"],
    initial_nnz=AUTO_RREF_THRESHOLDS["min_nnz"],
)


def _select(requested, *, numba=True, prime=P, **overrides):
    kw = {**BIG, **overrides}
    return select_rref_backend(requested, prime=prime, numba_available=numba, **kw)


# ---------------------------------------------------------------- choice list


def test_auto_is_a_choice_but_not_a_concrete_backend():
    assert AUTO_RREF_BACKEND == "auto"
    assert AUTO_RREF_BACKEND in RREF_BACKEND_CHOICES
    assert AUTO_RREF_BACKEND not in RREF_BACKENDS
    assert set(RREF_BACKEND_CHOICES) == {*RREF_BACKENDS, AUTO_RREF_BACKEND}
    # "auto" is always available: it degrades to dict, never errors.
    assert rref_backend_available(AUTO_RREF_BACKEND)


# ---------------------------------------------------------------- explicit requests


@pytest.mark.parametrize("backend", ["dict", "int_sparse_experimental"])
@pytest.mark.parametrize("numba", [True, False])
def test_explicit_backends_pass_through(backend, numba):
    assert _select(backend, numba=numba) == (backend, "explicit request")


def test_none_means_default_backend():
    assert _select(None, numba=True) == (DEFAULT_RREF_BACKEND, "explicit request")


def test_explicit_numba_with_numba_passes_through():
    assert _select(NUMBA_RREF_BACKEND, numba=True) == (NUMBA_RREF_BACKEND, "explicit request")


def test_explicit_numba_without_numba_is_loud():
    # Explicit requests are never silently substituted (Perf.12 contract).
    with pytest.raises(BackendUnavailable, match="speed"):
        _select(NUMBA_RREF_BACKEND, numba=False)


def test_unknown_backend_rejected():
    with pytest.raises(ValueError, match="unknown rref_backend"):
        _select("numpy")


# ---------------------------------------------------------------- auto decision table


def test_auto_picks_numba_when_everything_clears():
    chosen, reason = _select("auto", numba=True)
    assert chosen == NUMBA_RREF_BACKEND
    assert reason.startswith("auto:") and "clears thresholds" in reason


def test_auto_without_numba_falls_back_to_dict():
    chosen, reason = _select("auto", numba=False)
    assert chosen == DEFAULT_RREF_BACKEND
    assert "numba unavailable" in reason


def test_auto_big_prime_falls_back_to_dict():
    # int64 guard: products must fit; auto must not route into the downstream error.
    chosen, reason = _select("auto", numba=True, prime=BIG_PRIME)
    assert chosen == DEFAULT_RREF_BACKEND
    assert "2**31" in reason


@pytest.mark.parametrize("gate", ["n_rows", "n_cols", "initial_nnz"])
def test_auto_any_single_gate_below_threshold_means_dict(gate):
    chosen, reason = _select("auto", numba=True, **{gate: BIG[gate] - 1})
    assert chosen == DEFAULT_RREF_BACKEND
    assert "below thresholds" in reason


def test_auto_exactly_at_thresholds_means_numba():
    # Gates are inclusive (>=): the documented boundary belongs to numba.
    chosen, _ = _select("auto", numba=True)
    assert chosen == NUMBA_RREF_BACKEND


# ---------------------------------------------------------------- rref_mod_p integration


def _random_rows(rng, n_rows, cols, lo=2, hi=8):
    return [
        {c: rng.randrange(1, P) for c in rng.sample(cols, rng.randint(lo, min(hi, len(cols))))}
        for _ in range(n_rows)
    ]


def _label_cols(n):
    return [("T1", i, (i * 7 + 3) % 5) for i in range(n)]


def test_rref_mod_p_auto_small_matrix_matches_dict_and_reports_selection():
    rng = random.Random(20260711)
    cols = _label_cols(12)
    rows = _random_rows(rng, 30, cols)

    ref = rref_mod_p([dict(r) for r in rows], P, backend="dict")
    res = rref_mod_p([dict(r) for r in rows], P, backend="auto", collect_stats=True)

    # Tiny matrix: auto resolves to dict regardless of numba availability.
    assert res == ref
    stats = res.stats
    assert stats["requested_rref_backend"] == "auto"
    assert stats["selected_rref_backend"] == DEFAULT_RREF_BACKEND
    assert stats["backend"] == DEFAULT_RREF_BACKEND
    assert stats["backend_selection_reason"].startswith("auto:")
    assert stats["auto_thresholds_used"] == dict(AUTO_RREF_THRESHOLDS)
    assert isinstance(stats["numba_available"], bool)


def test_rref_mod_p_explicit_backend_reports_no_auto_thresholds():
    res = rref_mod_p([{0: 1, 1: 2}, {1: 3}], P, backend="dict", collect_stats=True)
    stats = res.stats
    assert stats["requested_rref_backend"] == "dict"
    assert stats["selected_rref_backend"] == "dict"
    assert stats["backend_selection_reason"] == "explicit request"
    assert stats["auto_thresholds_used"] is None


def test_rref_mod_p_rejects_auto_typo():
    with pytest.raises(ValueError, match="unknown rref_backend"):
        rref_mod_p([{0: 1}], P, backend="Auto")


# ---------------------------------------------------------------- end-to-end with numba


def test_rref_mod_p_auto_large_matrix_uses_numba_and_matches_dict():
    pytest.importorskip("numba")
    if not rref_backend_available(NUMBA_RREF_BACKEND):
        pytest.skip("numba import blocked")

    rng = random.Random(93)
    n_rows = AUTO_RREF_THRESHOLDS["min_rows"]
    cols = _label_cols(AUTO_RREF_THRESHOLDS["min_cols"])
    # ~8 nnz/row * 500 rows comfortably clears min_nnz=3000.
    rows = _random_rows(rng, n_rows, cols, lo=6, hi=10)

    ref = rref_mod_p([dict(r) for r in rows], P, backend="dict")
    res = rref_mod_p([dict(r) for r in rows], P, backend="auto", collect_stats=True)

    assert res == ref
    assert res.stats["selected_rref_backend"] == NUMBA_RREF_BACKEND
    assert "clears thresholds" in res.stats["backend_selection_reason"]
