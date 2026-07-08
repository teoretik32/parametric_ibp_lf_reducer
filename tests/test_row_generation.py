"""Tests for algebraic and coordinate-IBP row generation (Pass 2A)."""

from __future__ import annotations

import sympy as sp

from conftest import load_example
from parametric_ibp_lf_reducer import (
    ParamExpr,
    algebraic_row,
    coordinate_ibp_primitive_row,
    generate_algebraic_rows,
    generate_coordinate_ibp_rows,
    make_label,
    parse_family_text,
    render_row,
    zero_label,
)

PRIME = 2_147_483_647

ONE_VAR = """
IBPInput = <|
  "Variables" -> {x}, "Parameters" -> {ep}, "Regulators" -> {ep},
  "Polynomials" -> <| "G0" -> 1 + x |>,
  "MonomialExponents" -> <| x -> 0 |>,
  "PolynomialExponents" -> <| "G0" -> -2 |>
|>
"""


def _coeff_sympy(row, label):
    return sp.simplify(row.terms[label].to_sympy())


def test_algebraic_row_structure_one_var():
    fam = parse_family_text(ONE_VAR)
    row = algebraic_row(fam, zero_label(1, 1), 0)  # G0 = 1 + x
    # J(0,0) - 1*J(0,-1) - 1*J(1,-1) = 0
    assert set(row.terms.keys()) == {(0, 0), (0, -1), (1, -1)}
    assert _coeff_sympy(row, (0, 0)) == 1
    assert _coeff_sympy(row, (0, -1)) == -1
    assert _coeff_sympy(row, (1, -1)) == -1


def test_algebraic_rows_need_no_surface_check():
    fam = parse_family_text(ONE_VAR)
    result = generate_algebraic_rows(fam, [zero_label(1, 1)])
    assert len(result) == 1
    assert result.rejected == []  # exact identities are never surface-rejected


def test_coordinate_ibp_row_matches_hand_expansion():
    fam = parse_family_text(ONE_VAR)
    # P = x, i = x : d/dx(x * F) with e_x = 0, f_G0 = -2
    #   -> (p_i + e_i) J[0,0] + f_0 * b_i * c J[1,-1] = 1*J[0,0] - 2*J[1,-1]
    row = coordinate_ibp_primitive_row(fam, zero_label(1, 1), 0, (1,))
    assert set(row.terms.keys()) == {(0, 0), (1, -1)}
    assert _coeff_sympy(row, (0, 0)) == 1
    assert _coeff_sympy(row, (1, -1)) == -2


def test_coordinate_surface_filter_accepts_and_rejects():
    fam = parse_family_text(ONE_VAR)
    result = generate_coordinate_ibp_rows(fam, [zero_label(1, 1)], max_degree=2)
    # P=x is surface-free -> accepted; P=1 and P=x^2 are not -> rejected.
    assert len(result) >= 1
    accepted_ps = {row.provenance["P"] for row in result.rows}
    assert (1,) in accepted_ps
    reasons = {r.reason for r in result.rejected}
    assert "surface_not_free" in reasons
    rejected_ps = {r.provenance["P"] for r in result.rejected}
    assert (0,) in rejected_ps  # P=1: marginal at x=0 -> not surface-free
    assert (2,) in rejected_ps  # P=x^2: marginal at x=inf -> not surface-free


def test_rows_are_over_integer_labels_with_paramexpr_coeffs():
    fam = parse_family_text(ONE_VAR)
    result = generate_coordinate_ibp_rows(fam, [zero_label(1, 1)], max_degree=1)
    for row in result.rows:
        for label, coeff in row.terms.items():
            assert isinstance(label, tuple) and len(label) == fam.nvars + fam.npolys
            assert all(isinstance(x, int) for x in label)
            assert isinstance(coeff, ParamExpr)


def test_d4_algebraic_row_has_parametric_coefficients():
    fam = parse_family_text(load_example("d4_explicit_family.wl.txt"))
    g2 = fam.poly_names.index("G2")  # G2 = r*x1*x2 + r*x2*x3 + r*x1*x2*x4 + x3*x4
    row = algebraic_row(fam, zero_label(4, 3), g2)
    r = sp.Symbol("r")
    # anchor +1
    assert _coeff_sympy(row, zero_label(4, 3)) == 1
    # r*x1*x2 monomial -> label n=(1,1,0,0), m lowers G2 -> coeff -r
    lab_r = make_label([1, 1, 0, 0], [0, 0, -1])
    assert sp.simplify(row.terms[lab_r].to_sympy() - (-r)) == 0
    # x3*x4 monomial -> coeff -1
    lab_1 = make_label([0, 0, 1, 1], [0, 0, -1])
    assert _coeff_sympy(row, lab_1) == -1


def test_render_row_is_wolfram_like_text():
    fam = parse_family_text(ONE_VAR)
    row = coordinate_ibp_primitive_row(fam, zero_label(1, 1), 0, (1,))
    text = render_row(fam, row)
    assert "J[" in text and text.endswith("= 0")
    assert "**" not in text  # Wolfram-like, not Python


def test_coordinate_rows_respect_eps_direction_regulation():
    # G0 exponent -1-ep, x exponent 0. For P=x the infinity-boundary exponent is
    #   p + e_x + f*maxpow_x(G0) = 1 + 0 + (-1-ep) = -ep : marginal at eps=0.
    # Vanishing at x=inf needs it < 0, i.e. eps>0 -> surface-free for eps->0^+ only.
    fam = parse_family_text("""
    IBPInput = <|
      "Variables" -> {x}, "Parameters" -> {ep}, "Regulators" -> {ep},
      "Polynomials" -> <| "G0" -> 1 + x |>,
      "MonomialExponents" -> <| x -> 0 |>,
      "PolynomialExponents" -> <| "G0" -> -1 - ep |>
    |>
    """)
    minus = {row.provenance["P"] for row in generate_coordinate_ibp_rows(
        fam, [zero_label(1, 1)], 1, eps_direction="minus").rows}
    plus = {row.provenance["P"] for row in generate_coordinate_ibp_rows(
        fam, [zero_label(1, 1)], 1, eps_direction="plus").rows}
    assert (1,) in plus and (1,) not in minus  # direction flips this row's acceptance
    assert minus != plus
