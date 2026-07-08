"""Tests for tangent IBP row generation (Pass 2C): div(Q F), no m-shift, surface-gated."""

from __future__ import annotations

from collections import Counter

import sympy as sp

from parametric_ibp_lf_reducer import (
    ParamExpr,
    SparsePoly,
    TangentField,
    generate_tangent_fields,
    generate_tangent_ibp_rows,
    is_locally_finite,
    make_label,
    parse_family_text,
    tangent_ibp_primitive_row,
    verify_tangent,
    zero_label,
)


def _mono(nvars, exps, coeff_int, params=()):
    return SparsePoly.monomial(nvars, exps, ParamExpr.from_int(coeff_int, params))


def _xy_field(fam):
    """The 11.1 field Q = (x*y, -x*y), verified against the family's polynomials."""
    p = fam.parameters
    qx, qy = _mono(2, (1, 1), 1, p), _mono(2, (1, 1), -1, p)
    ok, hs = verify_tangent(fam, [qx, qy])
    assert ok
    return TangentField((qx, qy), hs, (2, 0))


def test_tangent_row_11_1_structural_coefficients():
    # G=1+x+y, x^ep y^(2ep) G^-1 : div(Q F) with Q=(xy,-xy) gives
    #   (e_x+1) J[n_x, n_y+1; m] - (e_y+1) J[n_x+1, n_y; m] = 0,  e_x=ep, e_y=2ep.
    fam = parse_family_text("""
    IBPInput = <|
      "Variables" -> {x, y}, "Parameters" -> {ep}, "Regulators" -> {ep},
      "Polynomials" -> <| "G" -> 1 + x + y |>,
      "MonomialExponents" -> <| x -> ep, y -> 2*ep |>,
      "PolynomialExponents" -> <| "G" -> -1 |>
    |>
    """)
    row = tangent_ibp_primitive_row(fam, zero_label(2, 1), _xy_field(fam))
    ep = sp.Symbol("ep")
    assert set(row.terms.keys()) == {(0, 1, 0), (1, 0, 0)}
    assert sp.simplify(row.terms[(0, 1, 0)].to_sympy() - (ep + 1)) == 0
    assert sp.simplify(row.terms[(1, 0, 0)].to_sympy() - (-(2 * ep + 1))) == 0


def test_tangent_row_has_no_m_shift():
    fam = parse_family_text("""
    IBPInput = <|
      "Variables" -> {x, y}, "Parameters" -> {ep}, "Regulators" -> {ep},
      "Polynomials" -> <| "G" -> 1 + x + y |>,
      "MonomialExponents" -> <| x -> ep, y -> 2*ep |>,
      "PolynomialExponents" -> <| "G" -> -1 |>
    |>
    """)
    field = _xy_field(fam)
    for source_m in ((0,), (-2,), (3,)):
        source = make_label([0, 0], list(source_m))
        row = tangent_ibp_primitive_row(fam, source, field)
        assert not row.is_trivial()
        assert all(label[fam.nvars :] == source_m for label in row.terms)


def test_tangent_row_h_nonzero_adds_f0_times_h():
    # Q = (x+y+1, 0) is tangent to G=1+x+y with H=(1). The H-term contributes f_0 * H_0 to the
    # source-label coefficient with no m-shift. Comparing two labels whose f_0 differ by 1
    # isolates that contribution.
    fam = parse_family_text("""
    IBPInput = <|
      "Variables" -> {x, y}, "Parameters" -> {}, "Regulators" -> {},
      "Polynomials" -> <| "G" -> 1 + x + y |>,
      "MonomialExponents" -> <| x -> 0, y -> 0 |>,
      "PolynomialExponents" -> <| "G" -> -2 |>
    |>
    """)
    qx = SparsePoly.from_sympy(sp.sympify("x + y + 1"), fam.variables, fam.parameters)
    qy = SparsePoly.zero(2, ())
    ok, hs = verify_tangent(fam, [qx, qy])
    assert ok and not hs[0].is_zero  # H_0 = 1, genuinely nonzero
    field = TangentField((qx, qy), hs, (1, 0))

    row_a = tangent_ibp_primitive_row(fam, make_label([0, 0], [0]), field)  # f_0 = -2
    row_b = tangent_ibp_primitive_row(fam, make_label([0, 0], [-1]), field)  # f_0 = -3
    ca = row_a.terms[make_label([0, 0], [0])].to_sympy()
    cb = row_b.terms[make_label([0, 0], [-1])].to_sympy()
    # difference is exactly f_0(a) - f_0(b) = (-2) - (-3) = 1 = H_0 contribution, no m-shift.
    assert sp.simplify(ca - cb - 1) == 0
    assert all(label[fam.nvars :] == (0,) for label in row_a.terms)


def test_unverified_field_is_rejected_not_turned_into_row():
    fam = parse_family_text("""
    IBPInput = <|
      "Variables" -> {x, y}, "Parameters" -> {}, "Regulators" -> {},
      "Polynomials" -> <| "G" -> 1 + x + y |>,
      "MonomialExponents" -> <| x -> 0, y -> 0 |>,
      "PolynomialExponents" -> <| "G" -> -2 |>
    |>
    """)
    # Q=(x,0) is NOT tangent (Q.grad G = x); attach a bogus H=0.
    bad = TangentField((_mono(2, (1, 0), 1), SparsePoly.zero(2, ())), (SparsePoly.zero(2, ()),), (1, 0))
    assert not bad.is_tangent(fam)
    result = generate_tangent_ibp_rows(fam, [zero_label(2, 1)], [bad])
    assert len(result.rows) == 0  # never turned into a row
    assert any(r.reason == "field_not_tangent" for r in result.rejected)


def test_surface_gate_accepts_and_rejects_tangent_rows():
    # x^0 y^0 G^-6 : of the 7 tangent fields, exactly one passes the toric-flux surface gate.
    fam = parse_family_text("""
    IBPInput = <|
      "Variables" -> {x, y}, "Parameters" -> {}, "Regulators" -> {},
      "Polynomials" -> <| "G" -> 1 + x + y |>,
      "MonomialExponents" -> <| x -> 0, y -> 0 |>,
      "PolynomialExponents" -> <| "G" -> -6 |>
    |>
    """)
    fields = generate_tangent_fields(fam, [(1, 0), (2, 0)])
    result = generate_tangent_ibp_rows(fam, [zero_label(2, 1)], fields)
    assert len(result.rows) >= 1  # gate lets valid rows through
    reasons = Counter(r.reason for r in result.rejected)
    assert reasons["surface_not_free"] >= 1  # gate rejects flux-divergent rows
    assert len(result.rows) + len(result.rejected) == len(fields)


def test_tangent_rows_do_not_drop_non_lf_source_labels():
    # The primitive row builder must not skip a label just because it is not locally finite.
    fam = parse_family_text("""
    IBPInput = <|
      "Variables" -> {x, y}, "Parameters" -> {}, "Regulators" -> {},
      "Polynomials" -> <| "G" -> 1 + x + y |>,
      "MonomialExponents" -> <| x -> 0, y -> 0 |>,
      "PolynomialExponents" -> <| "G" -> -1 |>
    |>
    """)
    source = zero_label(2, 1)
    assert is_locally_finite(fam, source) is False  # base target is not LF
    row = tangent_ibp_primitive_row(fam, source, _xy_field(fam))
    assert not row.is_trivial()  # row still generated for the non-LF source
