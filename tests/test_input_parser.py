"""Parser tests for the explicit Wolfram-like association format (spec §3.1).

Uses the canonical D=4 example and the notebook-addendum fixtures. No reduction here — only
that variables/parameters/regulators and polynomial supports are preserved, and that ambiguous
/ integrand-only inputs are refused conservatively.
"""

from __future__ import annotations

import sympy as sp

import pytest

from conftest import load_example
from parametric_ibp_lf_reducer import (
    ParserNeedsExplicitFamily,
    parse_family_text,
    parse_mathematica_association,
    try_factor_integrand,
)
from parametric_ibp_lf_reducer.input_parser import MAssoc, MList


def _support_monomials(fam, poly_name):
    return set(fam.polynomial(poly_name).support())


def test_parse_canonical_d4_explicit_family():
    fam = parse_family_text(load_example("d4_explicit_family.wl.txt"))
    assert fam.variables == ("x1", "x2", "x3", "x4")
    assert fam.parameters == ("ep", "r")
    assert fam.regulators == ("ep",)
    assert fam.domain == "PositiveOrthant"
    assert fam.poly_names == ("G0", "G1", "G2")
    assert fam.nvars == 4 and fam.npolys == 3

    # G0 = 1 + x1 + x2 + x3 (x4 absent)
    assert _support_monomials(fam, "G0") == {
        (0, 0, 0, 0),
        (1, 0, 0, 0),
        (0, 1, 0, 0),
        (0, 0, 1, 0),
    }
    # G2 = r*x1*x2 + r*x2*x3 + r*x1*x2*x4 + x3*x4
    assert (1, 1, 0, 0) in _support_monomials(fam, "G2")
    assert (1, 1, 0, 1) in _support_monomials(fam, "G2")

    # Exponents preserved as ParamExpr (checked semantically via sympy).
    ep, r = sp.symbols("ep r")
    x4_index = fam.variables.index("x4")
    assert sp.simplify(fam.monomial_exponents[x4_index].to_sympy() - (-ep)) == 0
    g2_index = fam.poly_names.index("G2")
    assert sp.simplify(fam.polynomial_exponents[g2_index].to_sympy() - (-2 - ep)) == 0


def test_parse_association_ast_shapes():
    raw = parse_mathematica_association(load_example("d4_explicit_family.wl.txt"))
    assert isinstance(raw, MAssoc)
    assert "Variables" in raw and isinstance(raw["Variables"], MList)
    assert "Polynomials" in raw and isinstance(raw["Polynomials"], MAssoc)


@pytest.mark.parametrize(
    "fname,exp_vars,exp_params",
    [
        ("notebook_example1_d4_alt_explicit_family.wl.txt", ("x1", "x2", "x3", "x4"), ("ep", "r")),
        ("notebook_example2_n3_five_term_explicit_family.wl.txt", ("x", "y", "z"), ("ep", "r")),
        ("id4example2_candidate_explicit_family.wl.txt", ("x2", "x4", "x5", "x7"), ("ep",)),
        ("id3example3_x4_squared_candidate_family.wl.txt", ("x4", "x7", "x8"), ("ep",)),
    ],
)
def test_notebook_examples_parse_and_preserve_names(fname, exp_vars, exp_params):
    fam = parse_family_text(load_example(fname))
    assert fam.variables == exp_vars
    assert fam.parameters == exp_params
    assert fam.regulators == ("ep",)
    # Every declared polynomial parsed into a non-empty SparsePoly over the same var count.
    for name in fam.poly_names:
        poly = fam.polynomial(name)
        assert poly.nvars == len(exp_vars)
        assert not poly.is_zero


def test_notebook_example2_g3_has_parametric_coefficient():
    # G3 = r*x + y + z + y*z ; the r*x monomial must carry an r-dependent coefficient.
    fam = parse_family_text(load_example("notebook_example2_n3_five_term_explicit_family.wl.txt"))
    g3 = fam.polynomial("G3")
    x_index = fam.variables.index("x")
    rx_monomial = tuple(1 if i == x_index else 0 for i in range(fam.nvars))
    assert rx_monomial in g3.terms
    r = sp.Symbol("r")
    assert sp.simplify(g3.terms[rx_monomial].to_sympy() - r) == 0


def test_starred_example3_polynomial_support_parser_fixture():
    # Example 3* is known-value-only (no LF decomposition); used here ONLY as a parser fixture,
    # not as a family for main-spec 11.4.
    fam = parse_family_text(load_example("id4example2_candidate_explicit_family.wl.txt"))
    # G1 = 1 + x7 + x2*x7 + x4*x7 + x2*x4*x7 ; check the degree-3 mixed monomial x2*x4*x7.
    idx = {v: fam.variables.index(v) for v in ("x2", "x4", "x7")}
    mono = [0] * fam.nvars
    for v in ("x2", "x4", "x7"):
        mono[idx[v]] = 1
    assert tuple(mono) in fam.polynomial("G1").terms


def test_integrand_only_input_is_refused_conservatively():
    text = """
    IBPInput = <|
      "Variables" -> {x1, x2},
      "Parameters" -> {ep},
      "Regulators" -> {ep},
      "Integrand" -> x1^(-ep)*(1 + x1 + x2)^(2*ep)
    |>
    """
    raw = parse_mathematica_association(text)
    with pytest.raises(ParserNeedsExplicitFamily):
        # parse_explicit_family requires explicit polynomial family keys
        from parametric_ibp_lf_reducer import parse_explicit_family

        parse_explicit_family(raw)
    with pytest.raises(ParserNeedsExplicitFamily):
        try_factor_integrand(raw)


def test_undeclared_symbol_in_polynomial_is_rejected():
    text = """
    IBPInput = <|
      "Variables" -> {x1},
      "Parameters" -> {ep},
      "Regulators" -> {ep},
      "Polynomials" -> <| "G0" -> 1 + x1 + q |>,
      "MonomialExponents" -> <| x1 -> 0 |>,
      "PolynomialExponents" -> <| "G0" -> ep |>
    |>
    """
    raw = parse_mathematica_association(text)
    with pytest.raises(ValueError):
        from parametric_ibp_lf_reducer import parse_explicit_family

        parse_explicit_family(raw)
