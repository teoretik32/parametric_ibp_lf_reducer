"""External Int2 Method.13 gates: complete rays, exact surface criteria, contract provenance.

Phase A item 5 (ray regression), Phase B item 9 (coordinate facet criterion), Phase C item 13
(normal-flux synthetics), Phase 0 (revoked-artifact provenance keys). Fast exact checks only —
the heavy Level-0 re-audit lives in ``scripts/audit_surface_validation.py``.
"""
from __future__ import annotations

from pathlib import Path

import pytest
import sympy as sp

import parametric_ibp_lf_reducer.valuations as valuations
from parametric_ibp_lf_reducer import (
    ParamExpr,
    SparsePoly,
    SurfacePolicy,
    base_score,
    boundary_ray_data,
    complete_polyhedral_rays,
    coordinate_primitive_surface_free,
    heuristic_candidate_rays,
    is_locally_finite,
    parse_family_text,
    vector_field_surface_free,
    zero_label,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
INT2_INPUT = REPO_ROOT / "examples" / "external_int2_dimensionless_input.wl.txt"
RESULT = REPO_ROOT / "validation" / "external_int2_lf_result.m"

TARGET = (0, 0, 0, 0, 0, 0, 0)
L1 = (-1, 0, 0, -1, -1, 0, 0)
L2 = (-1, 0, 0, 0, -1, 0, -1)
L3 = (0, 0, 0, -1, 0, 0, -1)
L4 = (0, 0, 1, -1, 0, 0, -1)
MIXED_RAY = (0, -1, -1)

EXPECTED_INT2_RAYS = {
    (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1),
    (-1, 0, -1), (-1, 0, 1), (-1, 1, 0), (-1, 1, 1), (-1, 2, 1), (0, -1, -1),
    (0, 1, 1), (1, -2, -1), (1, -1, -1), (1, -1, 0), (1, 0, -1), (1, 0, 1),
}


@pytest.fixture(scope="module")
def int2():
    return parse_family_text(INT2_INPUT.read_text(encoding="utf-8"))


# --- Phase A: complete ray construction --------------------------------------------------------
def test_int2_complete_set_has_all_18_expected_rays(int2):
    rays, complete = complete_polyhedral_rays(int2)
    assert complete is True
    assert {r.direction for r in rays} == EXPECTED_INT2_RAYS
    assert MIXED_RAY in {r.direction for r in rays}


def test_old_heuristic_set_reproduces_the_false_positive(int2):
    """Score on every OLD candidate ray is positive for L4 — that is exactly how the false
    LF=True verdict happened; the joint-infinity ray was simply never tested."""
    old = [r.direction for r in heuristic_candidate_rays(int2)]
    assert MIXED_RAY not in old
    for d in old:
        assert sp.simplify(base_score(int2, L4, d)) > 0, d
    # the complete set catches it: strict score zero -> non-LF
    assert base_score(int2, L4, MIXED_RAY) == 0
    assert is_locally_finite(int2, L4) is False


def test_verified_verdicts_retained(int2):
    for lab in (L1, L2, L3):
        assert is_locally_finite(int2, lab) is True
    assert is_locally_finite(int2, TARGET) is False


def test_ray_set_is_deterministic(int2):
    a = [(r.direction, r.kind) for r in complete_polyhedral_rays(int2)[0]]
    b = [(r.direction, r.kind) for r in complete_polyhedral_rays(int2)[0]]
    assert a == b
    # primitive integer normalization: gcd of every direction is 1
    from math import gcd

    for d, _k in a:
        assert gcd(*[abs(x) for x in d[:2]] or [1]) >= 0
        g = 0
        for x in d:
            g = gcd(g, abs(x))
        assert g == 1


def test_budget_overflow_degrades_to_unknown_not_true(int2, monkeypatch):
    """Resource-budget overflow -> Unknown, never True; random rays cannot certify True."""
    monkeypatch.setattr(valuations, "RAY_ENUMERATION_BUDGET", 1)
    fresh = parse_family_text(INT2_INPUT.read_text(encoding="utf-8"))  # fresh family cache
    rays, complete = complete_polyhedral_rays(fresh)
    assert complete is False
    assert len(rays) == 2 * fresh.nvars  # coordinate rays only
    # L1 is genuinely LF, but with an incomplete enumeration the verdict must degrade
    assert is_locally_finite(fresh, L1) == "Unknown"
    # a witnessed divergence is still allowed to be False (coordinate rays catch the target)
    assert is_locally_finite(fresh, TARGET) is False


# --- Phase B: coordinate facet criterion -------------------------------------------------------
def test_method12r_regression_row_is_rejected(int2):
    """label [-1,0,0,-3,0,-3,0], variable x5, P = x5^2: mixed ray (1,-1,0) has facet score -1."""
    lab = (-1, 0, 0, -3, 0, -3, 0)
    for pol in (SurfacePolicy.limit("minus"), SurfacePolicy.chamber({"ep": "-3/5"})):
        assert coordinate_primitive_surface_free(int2, lab, 1, (0, 2, 0), policy=pol) is False


def test_component_local_pass_but_mixed_face_fail():
    """P = x with y^-3: both x-facet exponents pass componentwise, the (1,1) facet score is -1."""
    fam = parse_family_text("""
    IBPInput = <|
      "Variables" -> {x, y}, "Parameters" -> {ep}, "Regulators" -> {ep},
      "Polynomials" -> <| "G0" -> 1 + x + y |>,
      "MonomialExponents" -> <| x -> 0, y -> -3 |>,
      "PolynomialExponents" -> <| "G0" -> -2 |>
    |>
    """)
    lab = zero_label(2, 1)
    e, f = fam.exponent_at_label(lab)
    # component-local exponents of P=x at the x facets (the OLD criterion) both pass:
    exp_zero = 1 + e[0].to_sympy()  # 1 > 0
    exp_inf = 1 + e[0].to_sympy() + f[0].to_sympy() * 1  # 1 - 2 = -1 < 0
    assert exp_zero == 1 and exp_inf == -1
    # ... but the corrected transverse-face criterion rejects (facet integral diverges at y->0):
    assert coordinate_primitive_surface_free(fam, lab, 0, (1, 0)) is False


def test_ordinary_valid_coordinate_primitive_accepted():
    fam = parse_family_text("""
    IBPInput = <|
      "Variables" -> {x, y}, "Parameters" -> {ep}, "Regulators" -> {ep},
      "Polynomials" -> <| "G0" -> 1 + x + y |>,
      "MonomialExponents" -> <| x -> 1, y -> 0 |>,
      "PolynomialExponents" -> <| "G0" -> -3 |>
    |>
    """)
    assert coordinate_primitive_surface_free(fam, zero_label(2, 1), 0, (0, 0)) is True


def test_chamber_and_limit_use_the_same_complete_rays(int2):
    """Both policies must test the same complete ray set — only the sign reading differs."""
    dirs, complete = boundary_ray_data(int2)
    assert complete and set(dirs) == EXPECTED_INT2_RAYS
    # the 12R regression row is rejected under BOTH policies by the same mixed ray (see above);
    # and a chamber point cannot resurrect a row the limit policy rejects on an ep-free score.
    lab = (-1, 0, 0, -3, 0, -3, 0)
    for pol in (SurfacePolicy.limit("minus"), SurfacePolicy.limit("plus"),
                SurfacePolicy.chamber({"ep": "-3/5"}), SurfacePolicy.chamber({"ep": "-1/7"})):
        assert coordinate_primitive_surface_free(int2, lab, 1, (0, 2, 0), policy=pol) is False


def test_strict_zero_facet_score_fails():
    # 1-var: P=1 with e=0, G0^-2: facet score at x->0 is exactly 0 -> Failed (0*Inf boundary)
    fam = parse_family_text("""
    IBPInput = <|
      "Variables" -> {x}, "Parameters" -> {ep}, "Regulators" -> {ep},
      "Polynomials" -> <| "G0" -> 1 + x |>,
      "MonomialExponents" -> <| x -> 0 |>,
      "PolynomialExponents" -> <| "G0" -> -2 |>
    |>
    """)
    assert coordinate_primitive_surface_free(fam, zero_label(1, 1), 0, (0,)) is False


# --- Phase C: normal-flux criterion ------------------------------------------------------------
FLUX2 = """
IBPInput = <|
  "Variables" -> {x, y}, "Parameters" -> {ep}, "Regulators" -> {ep},
  "Polynomials" -> <| "G0" -> 1 + x + y |>,
  "MonomialExponents" -> <| x -> 0, y -> 0 |>,
  "PolynomialExponents" -> <| "G0" -> -2 |>
|>
"""


def _field(coeff_x: int, coeff_y: int):
    return [
        SparsePoly.monomial(2, (1, 0), ParamExpr.from_int(coeff_x, ())),
        SparsePoly.monomial(2, (0, 1), ParamExpr.from_int(coeff_y, ())),
    ]


def test_marginal_components_with_cancelling_normal_flux_are_valid():
    """Q = (x, -y): each component alone is marginal on (-1,-1) (score 0), but the normal flux
    N_d = d_x - d_y * 1 cancels EXACTLY on the diagonal rays -> the row is valid."""
    fam = parse_family_text(FLUX2)
    assert vector_field_surface_free(fam, zero_label(2, 1), _field(1, -1)) is True


def test_componentwise_valid_but_surviving_mixed_flux_is_invalid():
    """Q = (x, y): every component passes the coordinate rays, but N_(-1,-1) = -2 survives with
    facet score exactly 0 -> invalid."""
    fam = parse_family_text(FLUX2)
    field = _field(1, 1)
    coord_only = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    assert vector_field_surface_free(fam, zero_label(2, 1), field, rays=coord_only) is True
    assert vector_field_surface_free(fam, zero_label(2, 1), field) is False


def test_exact_vs_accidental_modular_cancellation():
    """Coefficients 1 and 2147483646 cancel modulo 2147483647 but NOT exactly — no cancellation
    may be granted, the flux survives and the row is rejected."""
    fam = parse_family_text(FLUX2)
    assert vector_field_surface_free(fam, zero_label(2, 1), _field(1, 2147483646)) is False
    # sanity: the exact cancellation of the same shape IS granted
    assert vector_field_surface_free(fam, zero_label(2, 1), _field(1, -1)) is True


def test_epsilon_dependent_exact_cancellation():
    """Cancellation must be exact in the COEFFICIENT RING (here: polynomials in ep)."""
    fam = parse_family_text(FLUX2)
    ep_pos = ParamExpr.from_sympy(sp.Symbol("ep"), ("ep",))
    ep_neg = ParamExpr.from_sympy(-sp.Symbol("ep"), ("ep",))
    field = [SparsePoly.monomial(2, (1, 0), ep_pos), SparsePoly.monomial(2, (0, 1), ep_neg)]
    assert vector_field_surface_free(fam, zero_label(2, 1), field) is True


# --- Phase 0: provenance of the revoked artifact ----------------------------------------------
def test_revoked_artifact_distinguishes_row_span_from_integral_identity():
    text = RESULT.read_text(encoding="utf-8")
    assert '"Status" -> "Revoked(Method.12R)"' in text
    assert '"AllLocallyFinite" -> False' in text
    assert '"FormalRowSpanCertificate" -> "Passed"' in text
    assert '"IntegralIdentityStatus" -> "Revoked"' in text
    assert '"SurfaceValidationStatus" -> "Failed"' in text
