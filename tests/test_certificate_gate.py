"""Tests for the reducer's row-span certificate gate (Pass D4.5).

The D4.4 investigation showed that an on-grid holdout is only as independent as the grid: on a
degenerate (product-lattice) sample grid the dense interpolation can return a wrong function
that still "validates". The reducer therefore certifies the reconstructed relation at
independent off-sample points and, with ``require_certificate_for_success=True`` (the default),
``Success`` additionally requires ``certificate_status == "Passed"``. These tests pin:

* the fast synthetic regression of the product-grid false success — degenerate samples make
  interpolation return a wrong-but-holdout-passing function, and the certificate gate turns the
  old false ``Success`` into ``VerificationFailed``;
* a correct reduction still passes the gate (``Passed`` + ``Success``);
* uninformative certificates (all points rank-filtered or bad) are ``Insufficient`` -> failure;
* disabling the gate restores the old behaviour, explicitly and only on request.
"""

from __future__ import annotations

from fractions import Fraction

import pytest
import sympy as sp

from parametric_ibp_lf_reducer import (
    CERTIFICATE_FAILED,
    CERTIFICATE_INSUFFICIENT,
    CERTIFICATE_NOT_RUN,
    CERTIFICATE_PASSED,
    FAILURE_VERIFICATION_FAILED,
    STATUS_SUCCESS,
    ParamExpr,
    parse_family_text,
    reduce_rows_once,
)
from parametric_ibp_lf_reducer.row_generation import Row

PRIMES = [2_147_483_647, 2_147_483_629, 2_147_483_587]

ONE_VAR = """
IBPInput = <|
  "Variables" -> {x}, "Parameters" -> {ep}, "Regulators" -> {ep},
  "Polynomials" -> <| "G0" -> 1 + x |>,
  "MonomialExponents" -> <| x -> 0 |>,
  "PolynomialExponents" -> <| "G0" -> -2 |>
|>
"""

T, A, B = (0, 0), (1, 0), (2, 0)
LF = {T: False, A: True, B: True}


def _fam():
    return parse_family_text(ONE_VAR)


def _coeff_row(sympy_coeff: str) -> Row:
    """Row  J(0,0) - C(ep) J(1,0) = 0  with C(ep) given by a SymPy string."""
    ep = ("ep",)
    return Row("test", {}, {
        T: ParamExpr.from_int(1, ep),
        A: ParamExpr.from_sympy(sp.sympify(f"-({sympy_coeff})"), ep),
    })


def _samples(values):
    return [{"ep": Fraction(v)} for v in values]


# --- the product-grid false-success regression (fast synthetic analogue of D4.4) --------------
# True coefficient: C(ep) = 2 + (ep-3)(ep-4)(ep-5)(ep-6)(ep-7)(ep-8). On the degenerate sample
# grid ep=3..8 every value is exactly 2, so interpolation returns the WRONG constant 2 and the
# on-grid holdout (also on the vanishing locus) cannot notice — exactly the 6x6-lattice failure.
DEGENERATE_COEFF = "2 + (ep-3)*(ep-4)*(ep-5)*(ep-6)*(ep-7)*(ep-8)"
DEGENERATE_SAMPLES = _samples([3, 4, 5, 6, 7, 8])


def test_certificate_gate_catches_product_grid_false_interpolation():
    fam = _fam()
    res = reduce_rows_once(
        fam, T, [T, A], [_coeff_row(DEGENERATE_COEFF)], PRIMES, DEGENERATE_SAMPLES, lf_flags=LF
    )
    # the reconstruction is formally fine and every term is LF — but the relation is WRONG,
    # and the off-sample certificate is what catches it:
    assert res.status == FAILURE_VERIFICATION_FAILED
    assert res.success is False
    assert res.formal_success is True  # FormalSuccess stays honest in diagnostics
    assert res.all_locally_finite is not True  # never advertised as True on a failure
    assert all(t.locally_finite is True for t in res.terms)  # per-term flags stay truthful
    cert = res.diagnostics.extra["certificate"]
    assert cert["certificate_status"] == CERTIFICATE_FAILED
    assert cert["n_certificate_points_failed"] >= 1
    assert cert["first_nonzero_residual"]  # the surviving columns are reported
    assert cert["selected_rank"] == 1
    text = res.wolfram_style_text
    assert '"Status" -> "Failure"' in text
    assert '"Error" -> "VerificationFailed"' in text


def test_disabling_the_gate_restores_old_false_success():
    """Documents the pre-D4.5 behaviour: without the certificate the wrong constant is a formal
    Success. Opt-out must be explicit; nothing runs a certificate silently."""
    fam = _fam()
    res = reduce_rows_once(
        fam, T, [T, A], [_coeff_row(DEGENERATE_COEFF)], PRIMES, DEGENERATE_SAMPLES,
        lf_flags=LF, require_certificate_for_success=False,
    )
    assert res.status == STATUS_SUCCESS  # the (wrong) old behaviour, only on explicit request
    assert res.diagnostics.extra["certificate"]["certificate_status"] == CERTIFICATE_NOT_RUN


def test_explicit_points_reported_even_when_gate_disabled():
    fam = _fam()
    res = reduce_rows_once(
        fam, T, [T, A], [_coeff_row(DEGENERATE_COEFF)], PRIMES, DEGENERATE_SAMPLES,
        lf_flags=LF, require_certificate_for_success=False,
        certificate_points=[{"ep": Fraction(10)}],
    )
    # informative diagnostics: the certificate ran and failed, but the disabled gate lets the
    # (formal) Success through — the caller asked for exactly that.
    assert res.status == STATUS_SUCCESS
    assert res.diagnostics.extra["certificate"]["certificate_status"] == CERTIFICATE_FAILED


# --- correct reductions pass the gate ---------------------------------------------------------
def test_certificate_gate_passes_for_correct_reduction():
    fam = _fam()
    ep = sp.Symbol("ep")
    f = (2 * ep - 1) / (ep + 3)
    res = reduce_rows_once(
        fam, T, [T, A], [_coeff_row("(2*ep - 1)/(ep + 3)")], PRIMES, _samples(range(2, 12)),
        lf_flags=LF,
    )
    assert res.status == STATUS_SUCCESS
    assert sp.simplify(res.terms[0].coefficient - f) == 0
    cert = res.diagnostics.extra["certificate"]
    assert cert["certificate_status"] == CERTIFICATE_PASSED
    assert cert["n_certificate_points_passed"] >= 1
    assert cert["n_certificate_points_failed"] == 0
    assert any("certificate" in m for m in res.diagnostics.messages)


# --- uninformative certificates are Insufficient, never a silent pass -------------------------
def test_all_points_rank_filtered_is_insufficient():
    """A certificate point where the matrix drops rank is uninformative: it neither passes nor
    fails, and an all-filtered certificate cannot support Success."""
    fam = _fam()
    ep = ("ep",)
    row1 = Row("test", {}, {
        T: ParamExpr.from_int(1, ep),
        A: ParamExpr.from_int(-1, ep),
        B: ParamExpr.from_int(-1, ep),
    })
    row2 = Row("test", {}, {
        B: ParamExpr.from_sympy(sp.sympify("ep - 2"), ep),
        A: ParamExpr.from_sympy(sp.sympify("-(ep - 2)"), ep),
    })  # vanishes at ep=2 -> rank 1 there vs generic 2
    res = reduce_rows_once(
        fam, T, [T, A, B], [row1, row2], PRIMES[:2], _samples([3, 4, 5, 6, 7]),
        lf_flags=LF, certificate_points=[{"ep": Fraction(2)}],
    )
    assert res.status == FAILURE_VERIFICATION_FAILED
    cert = res.diagnostics.extra["certificate"]
    assert cert["certificate_status"] == CERTIFICATE_INSUFFICIENT
    assert cert["n_certificate_rank_filtered"] == 1
    assert cert["selected_rank"] == 2
    assert cert["certificate_rank_histogram"] == {1: 1}


def test_all_points_bad_specialization_is_insufficient():
    fam = _fam()
    res = reduce_rows_once(
        fam, T, [T, A], [_coeff_row("(2*ep - 1)/(ep - 20)")], PRIMES, _samples(range(2, 12)),
        lf_flags=LF, certificate_points=[{"ep": Fraction(20)}],  # row denominator pole
    )
    assert res.status == FAILURE_VERIFICATION_FAILED
    cert = res.diagnostics.extra["certificate"]
    assert cert["certificate_status"] == CERTIFICATE_INSUFFICIENT
    assert cert["n_certificate_bad_points"] == 1
    assert cert["n_certificate_points_passed"] == 0


# --- coefficients corrupted AFTER a validating reconstruction are still caught ----------------
def test_corrupted_coefficients_after_reconstruction_fail_certificate(monkeypatch):
    """Even when reconstruction itself validates, the gate re-checks the *final* coefficients:
    a post-reconstruction corruption (off-by-one on every C) must yield VerificationFailed."""
    import parametric_ibp_lf_reducer.reducer as reducer_mod

    real = reducer_mod.reconstruct_coefficients

    def corrupt(records, params, **kw):
        return {lab: expr + 1 for lab, expr in real(records, params, **kw).items()}

    monkeypatch.setattr(reducer_mod, "reconstruct_coefficients", corrupt)
    fam = _fam()
    res = reduce_rows_once(
        fam, T, [T, A], [_coeff_row("(2*ep - 1)/(ep + 3)")], PRIMES, _samples(range(2, 12)),
        lf_flags=LF,
    )
    assert res.status == FAILURE_VERIFICATION_FAILED
    cert = res.diagnostics.extra["certificate"]
    assert cert["certificate_status"] == CERTIFICATE_FAILED
    assert cert["first_nonzero_residual"]


# --- mixed points: rank-deficient certificate samples are skipped, not fatal ------------------
def _rank_drop_system():
    """row1: T = A + B always; row2: (ep-2)(B - A) = 0 — vanishes at ep=2 (rank 2 -> 1)."""
    ep = ("ep",)
    row1 = Row("test", {}, {
        T: ParamExpr.from_int(1, ep),
        A: ParamExpr.from_int(-1, ep),
        B: ParamExpr.from_int(-1, ep),
    })
    row2 = Row("test", {}, {
        B: ParamExpr.from_sympy(sp.sympify("ep - 2"), ep),
        A: ParamExpr.from_sympy(sp.sympify("-(ep - 2)"), ep),
    })
    return [row1, row2]


def test_rank_deficient_certificate_points_skipped_but_generic_ones_decide():
    """A deficient point among the certificate points is counted as filtered and does NOT block
    Success as long as >= min_certificate_points generic points pass."""
    fam = _fam()
    res = reduce_rows_once(
        fam, T, [T, A, B], _rank_drop_system(), PRIMES[:2], _samples([3, 4, 5, 6, 7]),
        lf_flags=LF,
        certificate_points=[{"ep": Fraction(2)}, {"ep": Fraction(10)}, {"ep": Fraction(11)}],
    )
    assert res.status == STATUS_SUCCESS  # T = 2*A, certified at the two generic points
    cert = res.diagnostics.extra["certificate"]
    assert cert["certificate_status"] == CERTIFICATE_PASSED
    assert cert["n_certificate_points_passed"] == 2
    assert cert["n_certificate_rank_filtered"] == 1
    assert cert["n_certificate_rank_exceeded"] == 0
    assert cert["certificate_rank_histogram"] == {1: 1, 2: 2}


# --- a certificate point with rank ABOVE selected_rank proves selected_rank was not generic ----
def test_certificate_point_with_higher_rank_fails_honestly():
    """If every reduction sample was rank-deficient, the reconstruction's selected_rank is not
    the generic rank — a certificate point revealing a HIGHER rank must fail the gate (a
    specialization's rank can never exceed the generic one)."""
    fam = _fam()
    ep = ("ep",)
    vanish = "(ep-3)*(ep-4)*(ep-5)*(ep-6)*(ep-7)*(ep-8)"  # zero at EVERY reduction sample below
    row1 = Row("test", {}, {
        T: ParamExpr.from_int(1, ep),
        A: ParamExpr.from_int(-1, ep),
        B: ParamExpr.from_int(-1, ep),
    })
    row2 = Row("test", {}, {
        B: ParamExpr.from_sympy(sp.sympify(vanish), ep),
        A: ParamExpr.from_sympy(sp.sympify(f"-({vanish})"), ep),
    })
    res = reduce_rows_once(
        fam, T, [T, A, B], [row1, row2], PRIMES[:2], _samples([3, 4, 5, 6, 7, 8]),
        lf_flags=LF, certificate_points=[{"ep": Fraction(10)}],
    )
    assert res.status == FAILURE_VERIFICATION_FAILED
    cert = res.diagnostics.extra["certificate"]
    assert cert["certificate_status"] == CERTIFICATE_FAILED
    assert cert["selected_rank"] == 1  # built exclusively on deficient samples
    assert cert["n_certificate_rank_exceeded"] == 1
    assert cert["n_certificate_points_failed"] == 0  # failed by rank, not by residual


# --- config contract ---------------------------------------------------------------------------
def test_unsupported_certificate_rank_policy_rejected():
    fam = _fam()
    with pytest.raises(ValueError):
        reduce_rows_once(
            fam, T, [T, A], [_coeff_row("2")], PRIMES, _samples([2, 3, 4, 5]),
            lf_flags=LF, certificate_rank_policy="max_rank",
        )


def test_certificate_primes_override_is_used():
    fam = _fam()
    res = reduce_rows_once(
        fam, T, [T, A], [_coeff_row("(2*ep - 1)/(ep + 3)")], PRIMES[:2], _samples(range(2, 12)),
        lf_flags=LF, certificate_points=[{"ep": Fraction(15)}],
        certificate_primes=[PRIMES[2]],  # a prime not used for the reduction records
    )
    assert res.status == STATUS_SUCCESS
    assert res.diagnostics.extra["certificate"]["certificate_status"] == CERTIFICATE_PASSED
