"""Tests for the reduction result model + strict Success gate (Pass 2I.1).

Everything here uses a small *generic* family (arbitrary variable/polynomial names) so the gate
and the Wolfram-like export are exercised without hardcoding any validation case.
"""

from __future__ import annotations

import sympy as sp

from parametric_ibp_lf_reducer import (
    FAILURE_INTERPOLATION_FAILED,
    FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE,
    FAILURE_RESOURCE_LIMIT_REACHED,
    FAILURE_TARGET_NOT_REDUCIBLE,
    FAILURE_VERIFICATION_FAILED,
    STATUS_SUCCESS,
    ReductionResult,
    ReductionTerm,
    build_reduction_result_from_reconstruction,
    parse_family_text,
)

GENERIC_FAMILY_TEXT = """
IBPInput = <|
  "Variables" -> {w, z},
  "Parameters" -> {ep, r},
  "Regulators" -> {ep},
  "Domain" -> "PositiveOrthant",
  "Polynomials" -> <|
    "R0" -> 1 + w,
    "R1" -> 1 + z
  |>,
  "MonomialExponents" -> <|
    w -> -1 - ep,
    z -> ep
  |>,
  "PolynomialExponents" -> <|
    "R0" -> -1 + ep,
    "R1" -> -2 - ep
  |>,
  "TargetMultiplier" -> 1
|>;
"""

TARGET = (0, 0, 0, 0)
L1 = (0, 0, -1, 0)  # R0^-1  -> "1/R0"
L2 = (0, 0, 0, -1)  # R1^-1  -> "1/R1"


def _family():
    return parse_family_text(GENERIC_FAMILY_TEXT)


def _verified_kwargs(**over):
    """Reconstruction evidence for a genuinely-verified reduction (overridable per test)."""
    base = dict(
        reconstruction_verified=True,
        independent_validation_passed=True,
        formal_success=True,
    )
    base.update(over)
    return base


def _coeffs():
    ep, r = sp.symbols("ep r")
    return {L1: ep + r, L2: 2 * ep}


# --- Success path ----------------------------------------------------------------------------
def test_success_only_when_all_lf_and_reconstruction_verified():
    res = build_reduction_result_from_reconstruction(
        _family(), TARGET, _coeffs(), {L1: True, L2: True}, **_verified_kwargs()
    )
    assert isinstance(res, ReductionResult)
    assert res.status == STATUS_SUCCESS
    assert res.success is True
    assert res.all_locally_finite is True
    assert len(res.terms) == 2
    assert all(t.locally_finite is True for t in res.terms)
    # every RHS integrand renders with the family's own polynomial names, no ** anywhere
    text = res.wolfram_style_text
    assert '"Status" -> "Success"' in text
    assert '"Status" -> "Failure"' not in text
    assert '"Error"' not in text


# --- local-finiteness failures ---------------------------------------------------------------
def test_non_lf_term_is_failure_not_success():
    res = build_reduction_result_from_reconstruction(
        _family(), TARGET, _coeffs(), {L1: True, L2: False}, **_verified_kwargs()
    )
    assert res.status == FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE
    assert res.success is False
    assert res.all_locally_finite is False
    assert L2 in res.diagnostics.non_lf_terms


def test_unknown_lf_term_is_failure_not_success():
    res = build_reduction_result_from_reconstruction(
        _family(), TARGET, _coeffs(), {L1: True, L2: "Unknown"}, **_verified_kwargs()
    )
    assert res.status == FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE
    assert res.success is False
    assert res.all_locally_finite == "Unknown"
    assert L2 in res.diagnostics.unknown_lf_terms


def test_formal_success_with_non_lf_term_is_failure_but_formalsuccess_true():
    res = build_reduction_result_from_reconstruction(
        _family(), TARGET, _coeffs(), {L1: True, L2: False},
        **_verified_kwargs(formal_success=True),
    )
    assert res.success is False
    assert res.status == FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE  # precise internal reason
    assert res.formal_success is True  # a formal normal form was found...
    text = res.wolfram_style_text
    assert '"FormalSuccess" -> True' in text  # ...and it is reported, but status is a failure
    assert '"Status" -> "Failure"' in text  # coarse exported status
    assert '"Error" -> "NormalFormNotLocallyFinite"' in text  # concrete reason in Error


# --- reconstruction / verification / target failures -----------------------------------------
def test_interpolation_failed_is_failure():
    res = build_reduction_result_from_reconstruction(
        _family(), TARGET, None, {}, **_verified_kwargs(interpolation_failed=True)
    )
    assert res.status == FAILURE_INTERPOLATION_FAILED
    assert res.success is False


def test_missing_coefficients_maps_to_interpolation_failed():
    res = build_reduction_result_from_reconstruction(
        _family(), TARGET, None, {}, **_verified_kwargs()
    )
    assert res.status == FAILURE_INTERPOLATION_FAILED


def test_unverified_reconstruction_is_verification_failed():
    res = build_reduction_result_from_reconstruction(
        _family(), TARGET, _coeffs(), {L1: True, L2: True},
        reconstruction_verified=True, independent_validation_passed=False, formal_success=True,
    )
    assert res.status == FAILURE_VERIFICATION_FAILED
    assert res.success is False
    # even though the (unconfirmed) terms are all LF=True, a failure must not advertise True
    assert res.all_locally_finite == "Unknown"


def test_target_not_reducible_is_failure():
    res = build_reduction_result_from_reconstruction(
        _family(), TARGET, _coeffs(), {L1: True, L2: True},
        **_verified_kwargs(target_reducible=False),
    )
    assert res.status == FAILURE_TARGET_NOT_REDUCIBLE
    assert res.success is False


def test_resource_limit_reached_is_failure():
    res = build_reduction_result_from_reconstruction(
        _family(), TARGET, _coeffs(), {L1: True, L2: True},
        **_verified_kwargs(resource_limit_reached=True),
    )
    assert res.status == FAILURE_RESOURCE_LIMIT_REACHED


# --- zero-reduction handling -----------------------------------------------------------------
def test_explicit_zero_reduction_allowed_is_success():
    res = build_reduction_result_from_reconstruction(
        _family(), TARGET, {}, {}, allow_zero_reduction=True, **_verified_kwargs()
    )
    assert res.status == STATUS_SUCCESS
    assert res.terms == ()
    assert res.diagnostics.zero_reduction is True


def test_empty_reduction_without_zero_flag_is_failure():
    res = build_reduction_result_from_reconstruction(
        _family(), TARGET, {}, {}, **_verified_kwargs()
    )
    assert res.status == FAILURE_VERIFICATION_FAILED


# --- Wolfram-like output ---------------------------------------------------------------------
def test_wolfram_output_has_status_and_no_python_power():
    ep, r = sp.symbols("ep r")
    coeffs = {L1: ep**2 + r, L2: sp.Integer(1)}  # ep**2 must render as ep^2, never ep**2
    res = build_reduction_result_from_reconstruction(
        _family(), TARGET, coeffs, {L1: True, L2: True}, **_verified_kwargs()
    )
    text = res.wolfram_style_text
    assert '"Status"' in text
    assert "AllLocallyFinite" in text
    assert "Terms" in text
    assert "Diagnostics" in text
    assert "**" not in text  # no Python power operator leaks into the Wolfram text
    assert "ep^2" in text


def test_failure_output_includes_error_and_formalsuccess():
    res = build_reduction_result_from_reconstruction(
        _family(), TARGET, None, {}, **_verified_kwargs(interpolation_failed=True)
    )
    text = res.wolfram_style_text
    assert '"Status" -> "Failure"' in text
    assert '"Error" -> "InterpolationFailed"' in text  # concrete reason surfaced in Error
    assert '"FormalSuccess"' in text
    assert "**" not in text


def test_failure_export_contract_is_uniform():
    """Every failure exports Status->Failure + Error-><reason>, AllLocallyFinite never True."""
    res = build_reduction_result_from_reconstruction(
        _family(), TARGET, _coeffs(), {L1: True, L2: False}, **_verified_kwargs()
    )
    text = res.wolfram_style_text
    assert '"Status" -> "Failure"' in text
    assert '"Error" -> "NormalFormNotLocallyFinite"' in text
    assert res.exported_status == "Failure"
    assert res.failure_reason == FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE
    assert res.all_locally_finite in (False, "Unknown")
    assert "**" not in text


# --- genericity: no hardcoded variable/polynomial names --------------------------------------
def test_integrand_text_uses_family_specific_names():
    fam = _family()
    res = build_reduction_result_from_reconstruction(
        fam, TARGET, _coeffs(), {L1: True, L2: True}, **_verified_kwargs()
    )
    integrands = {t.integrand_text for t in res.terms}
    assert integrands == {"1/R0", "1/R1"}  # reflects THIS family's names, nothing hardcoded
    assert all(isinstance(t, ReductionTerm) for t in res.terms)
