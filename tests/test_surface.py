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


def test_coordinate_surface_is_component_local_not_toric_overstrict():
    fam = parse_family_text(ASYMMETRIC)
    lab = zero_label(2, 1)
    xi = fam.variables.index("x")
    yi = fam.variables.index("y")
    # Coordinate IBP in x with multiplier P = x : surface-free at x=0 and x=inf.
    assert coordinate_primitive_surface_free(fam, lab, xi, multiplier_exps=(1, 0)) is True
    # The y-component genuinely diverges at y=0, so its own coordinate check fails -- but that
    # does NOT make the x-check fail (component-local, not toric-overstrict).
    assert coordinate_primitive_surface_free(fam, lab, yi, multiplier_exps=(0, 0)) is False


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
