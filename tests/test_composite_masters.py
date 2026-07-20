"""Tests for Method.3 composite locally-finite master feasibility.

Synthetic checks are hand-verified: for the one-variable family with
``G0 = 1 + x`` the labels A = (0, 1) -> (1+x)^-1 and B = (1, 0) -> x*(1+x)^-2
are both log-divergent as x -> oo with identical leading Laurent layers
(1/x - ... ), so the composite A - B = (1+x)^-2 is locally finite everywhere.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest
import sympy as sp

from parametric_ibp_lf_reducer import base_score, parse_family_text
from parametric_ibp_lf_reducer.composite_masters import (
    STATUS_FEASIBLE,
    STATUS_NO_COMPOSITE,
    CompositeCandidate,
    build_candidate_pool,
    composite_master_feasibility,
    feasibility_to_payload,
    leading_asymptotic_signature,
)

ONE_VAR = """
IBPInput = <|
  "Variables" -> {x}, "Parameters" -> {ep}, "Regulators" -> {ep},
  "Polynomials" -> <| "G0" -> 1 + x |>,
  "MonomialExponents" -> <| x -> 0 |>,
  "PolynomialExponents" -> <| "G0" -> -2 |>
|>
"""

LABEL_A = (0, 1)  # (1+x)^-1
LABEL_B = (1, 0)  # x*(1+x)^-2
RAY_INF = (-1,)  # x -> oo


@pytest.fixture(scope="module")
def one_var_family():
    return parse_family_text(ONE_VAR)


@pytest.fixture(scope="module")
def int2_family():
    path = (
        Path(__file__).resolve().parents[1]
        / "examples"
        / "external_int2_dimensionless_input.wl.txt"
    )
    return parse_family_text(path.read_text(encoding="utf-8"))


# Non-LF terms of the certified Method.2 corrected normal form (see
# notes/EXTERNAL_INT2_AUDIT.md and validation/external_int2_corrected_reduction.json).
INT2_NF_LABELS = [
    (-1, 0, 0, -1, -1, 0, 0),
    (-1, 0, 0, 0, -1, 0, -1),
    (0, 0, 0, -1, 0, 0, -1),
    (0, 0, 0, 0, -1, 0, 0),
    (0, 0, 0, 0, 0, -1, 0),
    (0, 0, 1, -1, 0, 0, -1),
]


# --- Phase A: leading asymptotic signatures ---------------------------------------------


def test_signature_matches_base_score_and_hand_series(one_var_family):
    y = sp.Symbol("y_x")
    s_a = leading_asymptotic_signature(one_var_family, LABEL_A, RAY_INF, order=3)
    s_b = leading_asymptotic_signature(one_var_family, LABEL_B, RAY_INF, order=3)
    for sig, label in ((s_a, LABEL_A), (s_b, LABEL_B)):
        assert sig.score == base_score(one_var_family, label, RAY_INF) == 0
    # (1+x)^-1 = 1/x - 1/x^2 + 1/x^3 - ...
    assert [sp.cancel(c) for c in s_a.coefficients] == [1 / y, -1 / y**2, y ** (-3)]
    # x*(1+x)^-2 = 1/x - 2/x^2 + 3/x^3 - ...
    assert [sp.cancel(c) for c in s_b.coefficients] == [1 / y, -2 / y**2, 3 / y**3]
    # identical leading layer -> cancellation is possible at depth 1
    assert sp.cancel(s_a.leading - s_b.leading) == 0


# --- candidate pool ---------------------------------------------------------------------


def test_candidate_pool_origins_dedup_and_determinism(one_var_family):
    pool = build_candidate_pool(
        one_var_family,
        [LABEL_A, LABEL_B],
        var_shift_axes=("x",),
        poly_shift_axes=("G0",),
        shift_depths=(-1,),
        numerator_vars=(),
        numerator_degree=0,
    )
    got = [(c.label, c.origin) for c in pool]
    assert got == [
        ((0, 1), "nf[0]"),
        ((-1, 1), "nf[0]+n_x-1"),
        ((0, 0), "nf[0]+m_G0-1"),
        ((1, 0), "nf[1]"),
        ((1, -1), "nf[1]+m_G0-1"),
    ]
    # dedup: (0, 0) reachable from both bases appears exactly once
    labels = [c.label for c in pool]
    assert len(labels) == len(set(labels))


# --- Phase B: feasibility ---------------------------------------------------------------


def test_two_candidate_cancellation_found(one_var_family):
    pool = [CompositeCandidate(LABEL_A, "manual[A]"), CompositeCandidate(LABEL_B, "manual[B]")]
    res = composite_master_feasibility(one_var_family, pool, RAY_INF)
    assert res.status == STATUS_FEASIBLE
    assert res.kernel_dimension == 1
    assert res.full_dimension == 1
    (vec,) = res.full_basis
    assert vec.fully_lf
    assert vec.failing_rays == () and vec.unknown_rays == ()
    nz = vec.nonzero(res.participants)
    assert {(lbl, sp.Integer(c)) for lbl, c in nz} == {
        (LABEL_A, sp.Integer(1)),
        (LABEL_B, sp.Integer(-1)),
    }


def test_single_divergent_candidate_has_trivial_kernel(one_var_family):
    res = composite_master_feasibility(
        one_var_family, [CompositeCandidate(LABEL_A, "manual[A]")], RAY_INF
    )
    assert res.status == STATUS_NO_COMPOSITE
    assert "primary-ray cancellation kernel is trivial" in res.notes


def test_all_lf_pool_reports_no_participants(one_var_family):
    res = composite_master_feasibility(
        one_var_family, [CompositeCandidate((0, 0), "manual[base]")], RAY_INF
    )
    assert res.status == STATUS_NO_COMPOSITE
    assert "no pool member is non-locally-finite on the primary ray" in res.notes


def test_feasibility_deterministic(one_var_family):
    pool = [CompositeCandidate(LABEL_A, "manual[A]"), CompositeCandidate(LABEL_B, "manual[B]")]
    first = feasibility_to_payload(composite_master_feasibility(one_var_family, pool, RAY_INF))
    second = feasibility_to_payload(composite_master_feasibility(one_var_family, pool, RAY_INF))
    assert first == second


# --- real External Int2 family ----------------------------------------------------------


def test_int2_leading_signatures_fast(int2_family):
    x5, x7, r = sp.symbols("x5 x7 r")
    s1 = leading_asymptotic_signature(int2_family, (0, 0, 0, 0, -1, 0, 0), (-1, 0, 0), order=2)
    s2 = leading_asymptotic_signature(int2_family, (0, 0, 0, 0, 0, -1, 0), (-1, 0, 0), order=2)
    assert s1.score == -1 and s2.score == -1
    assert sp.cancel(s1.leading - 1 / ((1 + x5) * (1 + x7) * (x7 + r * x5))) == 0
    assert sp.cancel(s2.leading - 1 / ((1 + x7) ** 2 * (x7 + r * x5))) == 0


@pytest.mark.skipif(
    os.environ.get("RUN_EXTERNAL_INT2") != "1",
    reason="heavy External Int2 composite feasibility; set RUN_EXTERNAL_INT2=1",
)
def test_int2_composite_feasibility_integration(int2_family):
    pool = build_candidate_pool(
        int2_family,
        INT2_NF_LABELS,
        var_shift_axes=("x2",),
        poly_shift_axes=("G0", "G3"),
        shift_depths=(-1, -2),
        numerator_vars=("x5", "x7"),
        numerator_degree=2,
    )
    assert len(pool) == 225
    res = composite_master_feasibility(int2_family, pool, (-1, 0, 0))
    assert res.status == STATUS_FEASIBLE
    assert len(res.participants) == 48
    assert res.kernel_dimension == 21
    assert res.full_dimension == 13
    assert all(vec.fully_lf for vec in res.full_basis)
    # the clean two-term difference composite must be in the refined basis
    expected = {
        ((-1, 0, 0, 0, -1, 0, 0), sp.Integer(1)),
        ((0, 0, 0, -1, -1, 0, 0), sp.Integer(-1)),
    }
    seen = [
        {(lbl, sp.cancel(c)) for lbl, c in vec.nonzero(res.participants)} for vec in res.full_basis
    ]
    assert expected in seen
