"""Tests for Wolfram-like text export: '^' not '**', exact rationals, faithful round-trip."""

from __future__ import annotations

import sympy as sp

from conftest import load_validation
from parametric_ibp_lf_reducer import (
    IntegrandFactor,
    ParamExpr,
    coeff_to_wolfram_text,
    integrand_to_wolfram_text,
    label_to_wolfram_text,
)


def _reparse_wolfram(text):
    return sp.sympify(text.replace("^", "**"))


def test_rational_prints_as_fraction_not_decimal():
    assert coeff_to_wolfram_text(sp.Rational(3, 7)) == "3/7"
    text = coeff_to_wolfram_text(sp.Rational(3, 7))
    assert "0.4" not in text and "**" not in text


def test_coefficient_export_uses_caret_and_roundtrips():
    data = load_validation("expected_d4_coefficients.json")
    for cs in data["coefficients"].values():
        original = sp.sympify(cs.replace("^", "**"))
        pe = ParamExpr.from_sympy(original, ("ep", "r"))
        text = coeff_to_wolfram_text(pe)
        assert "**" not in text, f"Python power operator leaked into {text!r}"
        # Faithful representation: reparsing the Wolfram-like text recovers the coefficient.
        assert sp.simplify(_reparse_wolfram(text) - original) == 0


def test_coefficient_export_is_factorized():
    # A coefficient with an obvious factorization should print in product form.
    pe = ParamExpr.from_sympy(sp.sympify("(2*ep - 1)/ep**2"), ("ep", "r"))
    text = coeff_to_wolfram_text(pe)
    assert text.count("/") == 1
    assert "ep^2" in text
    assert "**" not in text


def test_label_to_wolfram_text_format():
    assert label_to_wolfram_text((0, 0, 1, -2, -1)) == "{0,0,1,-2,-1}"


def test_integrand_to_wolfram_text_matches_factor():
    factor = IntegrandFactor(
        variables=("x1", "x2", "x3", "x4"),
        poly_names=("G0", "G1", "G2"),
        monomial_powers=(0, 1, 1, 0),
        poly_powers=(-2, -1, 0),
    )
    assert integrand_to_wolfram_text(factor) == "x2*x3/(G0^2*G1)"
    assert "**" not in integrand_to_wolfram_text(factor)
