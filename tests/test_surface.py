"""Tests for surface-free filters: coordinate = component-only; vector = toric flux."""

from __future__ import annotations

import sympy as sp

from parametric_ibp_lf_reducer import (
    ParamExpr,
    SparsePoly,
    coordinate_primitive_surface_free,
    parse_family_text,
    regulated_sign,
    vector_field_surface_free,
    zero_label,
)

# G0 = 1+x+y ; x boundaries are fine, but y (and the diagonal) diverge.
ASYMMETRIC = """
IBPInput = <|
  "Variables" -> {x, y}, "Parameters" -> {ep}, "Regulators" -> {ep},
  "Polynomials" -> <| "G0" -> 1 + x + y |>,
  "MonomialExponents" -> <| x -> 0, y -> -3 |>,
  "PolynomialExponents" -> <| "G0" -> -2 |>
|>
"""

FLUX_FAMILY = """
IBPInput = <|
  "Variables" -> {x, y}, "Parameters" -> {ep}, "Regulators" -> {ep},
  "Polynomials" -> <| "G0" -> 1 + x + y |>,
  "MonomialExponents" -> <| x -> 0, y -> 0 |>,
  "PolynomialExponents" -> <| "G0" -> -2 |>
|>
"""


def test_regulated_sign_numeric_and_epsilon_directions():
    assert regulated_sign(sp.Integer(3), []) == "pos"
    assert regulated_sign(sp.Integer(-2), []) == "neg"
    assert regulated_sign(sp.Integer(0), []) == "zero"
    # -ep: at eps->0^- it is positive; at eps->0^+ it is negative.
    assert regulated_sign(sp.sympify("-ep"), ["ep"], "minus") == "pos"
    assert regulated_sign(sp.sympify("-ep"), ["ep"], "plus") == "neg"
    assert regulated_sign(sp.sympify("ep"), ["ep"], "minus") == "neg"
    # a non-regulator symbol cannot be decided
    assert regulated_sign(sp.sympify("r"), ["ep"], "minus") == "unknown"


def test_coordinate_surface_checks_transverse_faces_not_component_only():
    """Method.13: the x-facet boundary term is the TRANSVERSE integral of P*F over y.

    For ASYMMETRIC (y^-3 G0^-2) the component-local x-exponents of P = x look fine, but the
    transverse y-integral of P*F diverges at y -> 0 on every x-slice: the mixed ray (1, 1) has
    facet score ``score(PF, d) - d_x = 0 - 1 = -1``. Pre-Method.13 this row was (wrongly)
    accepted; the corrected filter rejects it.
    """
    fam = parse_family_text(ASYMMETRIC)
    lab = zero_label(2, 1)
    xi = fam.variables.index("x")
    yi = fam.variables.index("y")
    assert coordinate_primitive_surface_free(fam, lab, xi, multiplier_exps=(1, 0)) is False
    # The y-check fails already component-locally (y^-3 at y -> 0).
    assert coordinate_primitive_surface_free(fam, lab, yi, multiplier_exps=(0, 0)) is False


def test_coordinate_surface_accepts_ordinary_valid_primitive():
    """An ordinary valid primitive stays accepted: enough decay on every facet ray."""
    fam = parse_family_text("""
    IBPInput = <|
      "Variables" -> {x, y}, "Parameters" -> {ep}, "Regulators" -> {ep},
      "Polynomials" -> <| "G0" -> 1 + x + y |>,
      "MonomialExponents" -> <| x -> 1, y -> 0 |>,
      "PolynomialExponents" -> <| "G0" -> -3 |>
    |>
    """)
    lab = zero_label(2, 1)
    xi = fam.variables.index("x")
    # x * G0^-3, P = 1: facet scores are 1 (+e_x), 2 (-e_x), 2 ((1,1)), 1 ((-1,-1)) — all > 0.
    assert coordinate_primitive_surface_free(fam, lab, xi, multiplier_exps=(0, 0)) is True


def test_vector_field_surface_uses_toric_flux():
    fam = parse_family_text(FLUX_FAMILY)
    lab = zero_label(2, 1)
    one = ParamExpr.one(())
    qx = SparsePoly.monomial(2, (1, 0), one)  # Q_x = x
    qy = SparsePoly.monomial(2, (0, 1), one)  # Q_y = y
    field = [qx, qy]

    # Restricted to the coordinate rays only, the flux is surface-free.
    coord_only = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    assert vector_field_surface_free(fam, lab, field, rays=coord_only) is True
    # Adding the diagonal toric ray (-1,-1) makes the flux fail -> the toric ray is used.
    assert vector_field_surface_free(fam, lab, field, rays=[(-1, -1)]) is False
    # The default (full toric candidate set) therefore is not surface-free.
    assert vector_field_surface_free(fam, lab, field) is not True


def test_tangent_field_tangency_sanity_case_11_1():
    # 11.1 sanity: G = 1+x+y, Q = (x*y, -x*y) is tangent to G (Q.grad G = 0), so a div(Q F) row
    # would not shift m. Here we only verify tangency via sparse-poly algebra (no row generation).
    g = SparsePoly.from_sympy(sp.sympify("1 + x + y"), ["x", "y"], [])
    dgx = g.derivative(0)  # d/dx (1+x+y) = 1
    dgy = g.derivative(1)  # d/dy (1+x+y) = 1
    one = ParamExpr.one(())
    qx = SparsePoly.monomial(2, (1, 1), one)  # x*y
    qy = SparsePoly.monomial(2, (1, 1), one.scale_int(-1))  # -x*y
    defect = qx.mul(dgx).add(qy.mul(dgy))  # Q . grad G
    assert defect.is_zero
