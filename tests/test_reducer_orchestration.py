"""Tests for the reducer orchestration MVP (Pass 2I.2).

The orchestration is exercised mostly through ``reduce_rows_once`` on tiny *synthetic* rows over a
small *generic* family (arbitrary names) so we test wiring + failure mapping without the heavy
row-generation layer. One ``reduce_family_once`` smoke test drives the full pipeline.
"""

from __future__ import annotations

from fractions import Fraction

import sympy as sp

from parametric_ibp_lf_reducer import (
    ALL_FAILURE_REASONS,
    FAILURE_INTERPOLATION_FAILED,
    FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE,
    FAILURE_TARGET_NOT_REDUCIBLE,
    STATUS_SUCCESS,
    ParamExpr,
    ReducerConfig,
    ReductionResult,
    Row,
    parse_family_text,
    reduce_family_once,
    reduce_rows_once,
)

PRIMES = [2_147_483_647, 2_147_483_629, 2_147_483_587]

GENERIC_FAMILY_TEXT = """
IBPInput = <|
  "Variables" -> {u, v},
  "Parameters" -> {ep},
  "Regulators" -> {ep},
  "Domain" -> "PositiveOrthant",
  "Polynomials" -> <| "P0" -> 1 + u, "P1" -> 1 + v |>,
  "MonomialExponents" -> <| u -> -1 - ep, v -> ep |>,
  "PolynomialExponents" -> <| "P0" -> -1 + ep, "P1" -> -2 - ep |>,
  "TargetMultiplier" -> 1
|>;
"""

# A differently-named family to prove nothing is hardcoded.
ALT_FAMILY_TEXT = GENERIC_FAMILY_TEXT.replace("u", "a").replace("v", "b").replace("P0", "Q0").replace("P1", "Q1")

T = (0, 0, 0, 0)  # target = base integrand
M = (0, 0, -1, 0)  # P0^-1  -> "1/P0"
M2 = (0, 0, 0, -1)  # P1^-1


def _family():
    return parse_family_text(GENERIC_FAMILY_TEXT)


def _samples(vals):
    return [{"ep": Fraction(v)} for v in vals]


def _row(params, terms):
    row = Row(kind="synthetic", provenance={})
    for label, expr in terms.items():
        row.add_term(label, ParamExpr.from_sympy(sp.sympify(expr), params))
    return row


def _coeff_of(res, label):
    return next(t for t in res.terms if t.label == label).coefficient


# --- success ---------------------------------------------------------------------------------
def test_tiny_synthetic_success_through_reduce_rows_once():
    fam = _family()
    ep = sp.Symbol("ep")
    # 1*J[T] + (ep+3)*J[M] = 0  =>  J[T] = -(ep+3) J[M]
    row = _row(fam.parameters, {T: 1, M: "ep + 3"})
    res = reduce_rows_once(
        fam, T, [T, M], [row], PRIMES, _samples([1, 2, 3, 4, 5]), lf_flags={T: True, M: True}
    )
    assert isinstance(res, ReductionResult)
    assert res.status == STATUS_SUCCESS
    assert res.success is True
    assert res.all_locally_finite is True
    assert sp.simplify(_coeff_of(res, M) - (-(ep + 3))) == 0
    text = res.wolfram_style_text
    assert '"Status" -> "Success"' in text
    assert "**" not in text


# --- failure mapping -------------------------------------------------------------------------
def test_target_not_pivot_is_failure():
    fam = _family()
    row = _row(fam.parameters, {M: 1, M2: "ep + 3"})  # target T appears in no row
    res = reduce_rows_once(
        fam, T, [T, M, M2], [row], PRIMES, _samples([1, 2, 3, 4, 5]),
        lf_flags={T: True, M: True, M2: True},
    )
    assert res.status == FAILURE_TARGET_NOT_REDUCIBLE
    assert res.success is False
    assert res.diagnostics.extra["n_target_not_pivot"] > 0
    assert res.diagnostics.extra["n_reduced_records"] == 0


def test_interpolation_failed_is_failure():
    fam = _family()
    row = _row(fam.parameters, {T: 1, M: "ep + 3"})
    # only two distinct sample points -> reconstruction cannot validate
    res = reduce_rows_once(
        fam, T, [T, M], [row], PRIMES, _samples([1, 2]), lf_flags={T: True, M: True}
    )
    assert res.status == FAILURE_INTERPOLATION_FAILED
    assert res.success is False


def test_non_lf_reconstructed_term_is_failure():
    fam = _family()
    row = _row(fam.parameters, {T: 1, M: "ep + 3"})
    res = reduce_rows_once(
        fam, T, [T, M], [row], PRIMES, _samples([1, 2, 3, 4, 5]), lf_flags={T: True, M: False}
    )
    assert res.status == FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE
    assert res.success is False
    assert res.all_locally_finite is False


def test_unknown_lf_reconstructed_term_is_failure():
    fam = _family()
    row = _row(fam.parameters, {T: 1, M: "ep + 3"})
    res = reduce_rows_once(
        fam, T, [T, M], [row], PRIMES, _samples([1, 2, 3, 4, 5]),
        lf_flags={T: True, M: "Unknown"},
    )
    assert res.status == FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE
    assert res.success is False
    assert res.all_locally_finite == "Unknown"


# --- bad specializations skipped, not patched ------------------------------------------------
def test_bad_samples_skipped_but_valid_samples_succeed():
    fam = _family()
    ep = sp.Symbol("ep")
    # coefficient (ep+3)/ep -> denominator vanishes at ep=0 (bad specialization at that sample)
    row = _row(fam.parameters, {T: 1, M: "(ep + 3)/ep"})
    res = reduce_rows_once(
        fam, T, [T, M], [row], PRIMES, _samples([0, 1, 2, 3, 4, 5, 6]),
        lf_flags={T: True, M: True},
    )
    assert res.success is True
    assert res.diagnostics.extra["n_bad_specializations"] == len(PRIMES)  # only ep=0, all primes
    assert sp.simplify(_coeff_of(res, M) - (-(ep + 3) / ep)) == 0


# --- diagnostics -----------------------------------------------------------------------------
def test_run_diagnostics_include_rows_records_and_skips():
    fam = _family()
    row = _row(fam.parameters, {T: 1, M: "ep + 3"})
    samples = _samples([1, 2, 3, 4, 5])
    res = reduce_rows_once(fam, T, [T, M], [row], PRIMES, samples, lf_flags={T: True, M: True})
    ex = res.diagnostics.extra
    assert ex["n_rows"] == 1
    assert ex["n_records"] == len(PRIMES) * len(samples)
    assert ex["n_reduced_records"] == len(PRIMES) * len(samples)
    assert ex["n_skipped_records"] == 0
    assert ex["n_bad_specializations"] == 0
    assert "row_diagnostics" in ex and "reconstruction_diagnostics" in ex


# --- genericity ------------------------------------------------------------------------------
def test_no_hardcoded_family_names():
    fam = parse_family_text(ALT_FAMILY_TEXT)  # variables a,b and polynomials Q0,Q1
    row = _row(fam.parameters, {T: 1, M: "ep + 3"})
    res = reduce_rows_once(
        fam, T, [T, M], [row], PRIMES, _samples([1, 2, 3, 4, 5]), lf_flags={T: True, M: True}
    )
    assert res.success is True
    assert {t.integrand_text for t in res.terms} == {"1/Q0"}  # reflects THIS family's names


# --- full pipeline smoke ---------------------------------------------------------------------
def test_reduce_family_once_runs_end_to_end_smoke():
    fam = _family()
    cfg = ReducerConfig(
        primes=PRIMES[:2],
        samples=_samples([1, 2, 3, 4]),
        label_box=((0, 0), (-1, 0)),  # n fixed at 0, m in {-1,0} per polynomial
        max_ibp_degree=1,
    )
    res = reduce_family_once(fam, T, cfg)
    assert isinstance(res, ReductionResult)
    assert res.status in ({STATUS_SUCCESS} | ALL_FAILURE_REASONS)  # honest typed outcome
    assert "n_rows" in res.diagnostics.extra
    assert "n_records" in res.diagnostics.extra
