"""Tests for ParamExpr: exact parsing + finite-field evaluation cross-checked against SymPy."""

from __future__ import annotations

from fractions import Fraction

import sympy as sp

import pytest

from conftest import load_validation
from parametric_ibp_lf_reducer import ParamExpr

PRIME = 2_147_483_647  # 2**31 - 1, a Mersenne prime


def _sympy_value_mod_p(expr_str, subs, prime):
    """Reference: exact SymPy evaluation reduced modulo ``prime`` (or None if denom vanishes)."""
    expr = sp.sympify(expr_str.replace("^", "**"))
    val = sp.Rational(expr.subs(subs))
    num, den = int(val.p), int(val.q)
    if den % prime == 0:
        return None
    return (num % prime) * pow(den % prime, prime - 2, prime) % prime


def test_eval_matches_sympy_simple():
    params = ("ep", "r")
    pe = ParamExpr.from_sympy(sp.sympify("(2*ep - 1)/ep**2"), params)
    sample = {"ep": Fraction(3, 5), "r": Fraction(7, 4)}
    got = pe.eval_mod_p(sample, PRIME)
    expected = _sympy_value_mod_p("(2*ep - 1)/ep^2", {sp.Symbol("ep"): Fraction(3, 5)}, PRIME)
    assert got == expected


def test_denominator_zero_returns_none():
    pe = ParamExpr.from_sympy(sp.sympify("1/ep"), ("ep",))
    assert pe.eval_mod_p({"ep": 0}, PRIME) is None


def test_zero_and_constant():
    params = ("ep", "r")
    assert ParamExpr.zero(params).is_zero
    assert ParamExpr.from_sympy(sp.Integer(0), params).is_zero
    one = ParamExpr.one(params)
    assert one.eval_mod_p({"ep": 4, "r": 9}, PRIME) == 1


def test_arithmetic_consistent_with_eval():
    params = ("ep", "r")
    a = ParamExpr.from_sympy(sp.sympify("(ep + 2*r)/(ep - 1)"), params)
    b = ParamExpr.from_sympy(sp.sympify("r/(ep**2)"), params)
    sample = {"ep": Fraction(5, 3), "r": Fraction(2, 7)}
    av = a.eval_mod_p(sample, PRIME)
    bv = b.eval_mod_p(sample, PRIME)
    assert (a + b).eval_mod_p(sample, PRIME) == (av + bv) % PRIME
    assert (a - b).eval_mod_p(sample, PRIME) == (av - bv) % PRIME
    assert (a * b).eval_mod_p(sample, PRIME) == (av * bv) % PRIME
    assert a.scale_int(-3).eval_mod_p(sample, PRIME) == (-3 * av) % PRIME


def test_rejects_non_parameter_symbols():
    with pytest.raises(ValueError):
        ParamExpr.from_sympy(sp.sympify("ep + x1"), ("ep",))


@pytest.mark.parametrize(
    "vfile",
    [
        "expected_d4_coefficients.json",
        "notebook_example1_d4_alt_expected.json",
        "notebook_example2_n3_five_term_expected.json",
    ],
)
def test_validation_coefficients_parse_and_eval(vfile):
    """Every validation coefficient parses to a ParamExpr and evaluates consistently mod p."""
    data = load_validation(vfile)
    # Coefficients live under "coefficients" as a dict or a list depending on the file.
    coeffs = data["coefficients"]
    coeff_strings = list(coeffs.values()) if isinstance(coeffs, dict) else list(coeffs)
    params = ("ep", "r")
    sample = {"ep": Fraction(4, 9), "r": Fraction(11, 7)}
    subs = {sp.Symbol("ep"): Fraction(4, 9), sp.Symbol("r"): Fraction(11, 7)}
    for cs in coeff_strings:
        pe = ParamExpr.from_sympy(sp.sympify(cs.replace("^", "**")), params)
        got = pe.eval_mod_p(sample, PRIME)
        expected = _sympy_value_mod_p(cs, subs, PRIME)
        assert got == expected, f"mismatch for coefficient {cs!r}"


def test_d4_rational_check_point_evaluates():
    """The documented rational check point ep=-3/4, r=7/5 must be evaluable (nonzero denom)."""
    data = load_validation("expected_d4_coefficients.json")
    pt = data["rational_check_point"]
    ep = Fraction(pt["ep"])
    r = Fraction(pt["r"])
    sample = {"ep": ep, "r": r}
    subs = {sp.Symbol("ep"): ep, sp.Symbol("r"): r}
    for cs in data["coefficients"].values():
        pe = ParamExpr.from_sympy(sp.sympify(cs.replace("^", "**")), ("ep", "r"))
        got = pe.eval_mod_p(sample, PRIME)
        assert got is not None
        assert got == _sympy_value_mod_p(cs, subs, PRIME)
