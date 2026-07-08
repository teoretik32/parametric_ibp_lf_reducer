"""Row-generation hardening: arbitrary N (!=4), M>=2, non-x variable and non-G polynomial names.

These guard against any hidden hardcoding of dimensions or names in the row-generation layer,
and assert that no Success path exists yet.
"""

from __future__ import annotations

import sympy as sp

import pytest

import parametric_ibp_lf_reducer as pkg
from parametric_ibp_lf_reducer import (
    ParamExpr,
    SparsePoly,
    algebraic_row,
    generate_algebraic_rows,
    generate_coordinate_ibp_rows,
    make_label,
    parse_family_text,
    render_row,
    zero_label,
)

# N=3 (not 4), M=3 (>=2), variables u,v,w (not x_i), polynomials A,B,C (not G0/G1/H),
# parametric coefficient s (not r).
UVW_ABC = """
IBPInput = <|
  "Variables" -> {u, v, w},
  "Parameters" -> {ep, s},
  "Regulators" -> {ep},
  "Polynomials" -> <|
    "A" -> 1 + u + v,
    "B" -> 1 + w,
    "C" -> s*u*v + w
  |>,
  "MonomialExponents" -> <| u -> -ep, v -> 0, w -> 0 |>,
  "PolynomialExponents" -> <| "A" -> 2*ep, "B" -> 1 + 3*ep, "C" -> -2 - ep |>,
  "Assumptions" -> {s > 0}
|>
"""


def test_arbitrary_dimensions_and_names_parse():
    fam = parse_family_text(UVW_ABC)
    assert fam.variables == ("u", "v", "w")  # not x1..x4
    assert fam.poly_names == ("A", "B", "C")  # not G0/G1/H
    assert fam.nvars == 3 and fam.npolys == 3  # N != 4, M >= 2


def test_algebraic_row_custom_names_and_label_dimension():
    fam = parse_family_text(UVW_ABC)
    c_idx = fam.poly_names.index("C")  # C = s*u*v + w
    row = algebraic_row(fam, zero_label(3, 3), c_idx)
    # Every label has length N+M = 6, integer entries.
    for label in row.terms:
        assert len(label) == 6 and all(isinstance(x, int) for x in label)
    # s*u*v -> n=(1,1,0), m lowers C -> coeff -s
    lab_s = make_label([1, 1, 0], [0, 0, -1])
    s = sp.Symbol("s")
    assert sp.simplify(row.terms[lab_s].to_sympy() - (-s)) == 0
    # w -> n=(0,0,1), coeff -1
    lab_w = make_label([0, 0, 1], [0, 0, -1])
    assert sp.simplify(row.terms[lab_w].to_sympy() - (-1)) == 0
    assert sp.simplify(row.terms[zero_label(3, 3)].to_sympy() - 1) == 0


def test_algebraic_rows_one_per_polynomial_each_lowers_distinct_m():
    fam = parse_family_text(UVW_ABC)
    result = generate_algebraic_rows(fam, [zero_label(3, 3)])
    assert len(result) == fam.npolys == 3  # M >= 2, one row per polynomial
    assert result.rejected == []
    lowered = set()
    for row in result.rows:
        poly_idx = row.provenance["poly"]
        # some term has m lowered by exactly 1 in this polynomial's slot
        assert any(label[fam.nvars + poly_idx] == -1 for label in row.terms)
        lowered.add(poly_idx)
    assert lowered == {0, 1, 2}


def test_coordinate_rows_custom_dimension_and_m_shift():
    fam = parse_family_text(UVW_ABC)
    result = generate_coordinate_ibp_rows(fam, [zero_label(3, 3)], max_degree=1)
    assert len(result) >= 1  # surface filter admits some rows
    for row in result.rows:
        for label, coeff in row.terms.items():
            assert len(label) == 6
            assert isinstance(coeff, ParamExpr)
    # coordinate rows carry m-lowering terms (division by some polynomial G_l)
    assert any(
        any(label[fam.nvars + j] == -1 for j in range(fam.npolys))
        for row in result.rows
        for label in row.terms
    )


def test_render_uses_custom_names_not_default_ones():
    fam = parse_family_text(UVW_ABC)
    c_idx = fam.poly_names.index("C")
    text = render_row(fam, algebraic_row(fam, zero_label(3, 3), c_idx))
    assert "u" in text and "C" in text
    assert "x1" not in text and "G0" not in text and "H" not in text


def test_two_variable_two_polynomial_family_ab_pq():
    # A different arbitrary shape: N=2, M=2, variables a,b, polynomials P,Q.
    fam = parse_family_text("""
    IBPInput = <|
      "Variables" -> {a, b},
      "Parameters" -> {ep},
      "Regulators" -> {ep},
      "Polynomials" -> <| "P" -> 1 + a, "Q" -> 1 + b |>,
      "MonomialExponents" -> <| a -> 0, b -> 0 |>,
      "PolynomialExponents" -> <| "P" -> -2, "Q" -> -2 |>
    |>
    """)
    assert fam.variables == ("a", "b") and fam.poly_names == ("P", "Q")
    result = generate_algebraic_rows(fam, [zero_label(2, 2)])
    assert len(result) == 2
    # algebraic row for Q = 1 + b : J(0,0;0,0) - J(0,0;0,-1) - J(0,1;0,-1) = 0
    q_row = next(r for r in result.rows if r.provenance["poly"] == fam.poly_names.index("Q"))
    assert make_label([0, 0], [0, -1]) in q_row.terms
    assert make_label([0, 1], [0, -1]) in q_row.terms


def test_rows_do_not_depend_on_labels_being_locally_finite():
    # Rows are formal algebraic/IBP identities: they must be generated regardless of whether the
    # anchor label is locally finite. Use the (non-LF) base target label.
    fam = parse_family_text(UVW_ABC)
    base = zero_label(3, 3)
    alg = generate_algebraic_rows(fam, [base])
    coord = generate_coordinate_ibp_rows(fam, [base], max_degree=1)
    assert len(alg) == 3
    assert len(coord) >= 1  # generated even though the base integrand is not locally finite


def test_unimplemented_layers_stay_forbidden():
    # The text-in/text-out API is now available (Pass 2I.3) ...
    assert callable(pkg.reduce_wolfram_style_input)
    assert callable(pkg.reduce_wolfram_style_input_to_text)
    # ... but full family specialization and the CLI remain unimplemented/forbidden.
    fam = parse_family_text(UVW_ABC)
    with pytest.raises(NotImplementedError):
        fam.specialize({"ep": 1, "s": 1}, 2_147_483_647)
    from parametric_ibp_lf_reducer.__main__ import main

    with pytest.raises(SystemExit) as exc:
        main()
    assert exc.value.code != 0  # never a Success (exit 0)


def test_sparse_poly_layer_arbitrary_large_n():
    # Sanity that the underlying layer is not pinned to N=4 either.
    variables = [f"z{i}" for i in range(1, 10)]  # N=9
    poly = SparsePoly.from_sympy(sp.sympify("z1*z9 + z5**2 + 1"), variables, [])
    assert poly.nvars == 9
    support = set(poly.support())
    assert (1, 0, 0, 0, 0, 0, 0, 0, 1) in support  # z1*z9
    assert (0, 0, 0, 0, 2, 0, 0, 0, 0) in support  # z5^2
    assert poly.degree_in(4) == 2  # z5 is index 4
