"""Tests for SparsePoly: arbitrary N, degree>1 monomials, derivative, valuation, eval_mod_p."""

from __future__ import annotations

import sympy as sp

import pytest

from conftest import load_example
from parametric_ibp_lf_reducer import ParamExpr, SparsePoly, parse_family_text

PRIME = 2_147_483_647


def _poly(expr, variables, parameters=()):
    return SparsePoly.from_sympy(sp.sympify(expr.replace("^", "**")), variables, parameters)


def test_arbitrary_nvars_construction_and_support():
    variables = [f"x{i}" for i in range(1, 8)]  # N=7, no hardcoded count in the module
    p = _poly("x1*x3 + x7^3 + 5", variables)
    assert p.nvars == 7
    supp = set(p.support())
    assert (1, 0, 1, 0, 0, 0, 0) in supp
    assert (0, 0, 0, 0, 0, 0, 3) in supp  # x7^3
    assert (0, 0, 0, 0, 0, 0, 0) in supp  # constant term
    assert p.total_degree() == 3
    assert p.degree_in(6) == 3  # x7


def test_degree_two_monomial_x4_squared_starred_example4():
    # Example 4* is known-value-only (no LF decomposition); used here ONLY as a sparse-poly
    # fixture to exercise the x4^2 monomial, not as a family for main-spec 11.5.
    fam = parse_family_text(load_example("id3example3_x4_squared_candidate_family.wl.txt"))
    g3 = fam.polynomial("G3")  # x4 + x4^2 + x8 + x4*x8 + x4*x7*x8
    x4 = fam.variables.index("x4")
    x4_squared = tuple(2 if i == x4 else 0 for i in range(fam.nvars))
    assert x4_squared in g3.terms, "x4^2 monomial must be present in the sparse support"
    assert g3.degree_in(x4) == 2
    # mixed degree-3 monomial x4*x7*x8 exists
    idx = {v: fam.variables.index(v) for v in ("x4", "x7", "x8")}
    mixed = tuple(1 if i in idx.values() else 0 for i in range(fam.nvars))
    assert mixed in g3.terms
    assert g3.total_degree() == 3


def test_derivative_power_rule():
    p = _poly("x^3 + 2*x", ["x"])
    dp = p.derivative(0)  # 3*x^2 + 2
    assert dp.terms.keys() == {(2,), (0,)}
    # coefficient of x^2 is 3, constant is 2
    assert dp.terms[(2,)].eval_mod_p({}, PRIME) == 3
    assert dp.terms[(0,)].eval_mod_p({}, PRIME) == 2


def test_derivative_of_missing_variable_is_zero():
    p = _poly("1 + y", ["x", "y"])
    assert p.derivative(0).is_zero  # d/dx (1 + y) = 0


def test_valuation_min_convention():
    p = _poly("1 + x1 + x2", ["x1", "x2"])
    assert p.valuation((1, 1)) == 0  # min(0, 1, 1)
    assert p.valuation((-1, -1)) == -1  # min(0, -1, -1)
    assert p.valuation((2, 3)) == 0
    p2 = _poly("x1*x2 + x1^2", ["x1", "x2"])
    assert p2.valuation((1, 1)) == 2  # min(1+1, 2)


def test_add_mul_pow_binomial():
    x = _poly("x", ["x"])
    one = SparsePoly.one(1, ())
    onep_x = one.add(x)  # 1 + x
    sq = onep_x.pow_small(2)  # 1 + 2x + x^2
    assert sq.terms.keys() == {(0,), (1,), (2,)}
    assert sq.terms[(1,)].eval_mod_p({}, PRIME) == 2
    assert onep_x.pow_small(0) == SparsePoly.one(1, ())
    # mul consistency: (1+x)*(1+x) == (1+x)^2
    assert onep_x.mul(onep_x) == sq


def test_eval_mod_p_specializes_parametric_coefficients():
    # G2 = r*x1*x2 + x3*x4 with r a parameter -> at r=3 the x1*x2 coefficient is 3 mod p.
    variables = ["x1", "x2", "x3", "x4"]
    g = _poly("r*x1*x2 + x3*x4", variables, ("r",))
    spec = g.eval_mod_p({"r": 3}, PRIME)
    assert spec[(1, 1, 0, 0)] == 3
    assert spec[(0, 0, 1, 1)] == 1


def test_eval_mod_p_bad_point_returns_none():
    g = SparsePoly.monomial(1, (1,), ParamExpr.from_sympy(sp.sympify("1/ep"), ("ep",)))
    assert g.eval_mod_p({"ep": 0}, PRIME) is None


def test_scalar_mul_and_monomial_mul():
    variables = ["x", "y"]
    p = _poly("1 + x", variables)
    two = ParamExpr.from_int(2, ())
    assert p.scalar_mul(two).terms[(0, 0)].eval_mod_p({}, PRIME) == 2
    shifted = p.monomial_mul((0, 1))  # multiply by y -> y + x*y
    assert set(shifted.support()) == {(0, 1), (1, 1)}


def test_valuation_zero_polynomial_raises():
    with pytest.raises(ValueError):
        SparsePoly.zero(2, ()).valuation((1, 1))
