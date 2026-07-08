"""Tests for label-projection methods on ParametricFamily (Pass 1B)."""

from __future__ import annotations

import sympy as sp

import pytest

from conftest import load_example
from parametric_ibp_lf_reducer import IntegrandFactor, make_label, parse_family_text

PRIME = 2_147_483_647


def test_label_to_factor_and_wolfram_text_d4_master():
    fam = parse_family_text(load_example("d4_explicit_family.wl.txt"))
    # Master M1 relative factor: x2*x3/(G0^2*G1)  ->  n=(0,1,1,0), m=(-2,-1,0).
    label = make_label([0, 1, 1, 0], [-2, -1, 0])
    factor = fam.label_to_factor(label)
    assert isinstance(factor, IntegrandFactor)
    assert factor.monomial_powers == (0, 1, 1, 0)
    assert factor.poly_powers == (-2, -1, 0)
    assert fam.label_to_wolfram_text(label) == "x2*x3/(G0^2*G1)"


def test_integrand_factor_single_denominator_no_parens():
    # x4*x7*x8/G3 : single denominator factor must not be parenthesised.
    factor = IntegrandFactor(
        variables=("x4", "x7", "x8"),
        poly_names=("G3",),
        monomial_powers=(1, 1, 1),
        poly_powers=(-1,),
    )
    assert factor.to_wolfram_text() == "x4*x7*x8/G3"


def test_integrand_factor_pure_numerator_and_pure_denominator():
    variables = ("x1", "x2")
    polys = ("G0",)
    assert IntegrandFactor(variables, polys, (2, 0), (0,)).to_wolfram_text() == "x1^2"
    assert IntegrandFactor(variables, polys, (0, 0), (-1,)).to_wolfram_text() == "1/G0"


def test_exponent_at_label_shifts_base_exponents():
    fam = parse_family_text(load_example("d4_explicit_family.wl.txt"))
    # Base: x4 exponent -ep, G2 exponent -2-ep. Shift n_x4=+1, m_G2=+1.
    label = make_label([0, 0, 0, 1], [0, 0, 1])
    e, f = fam.exponent_at_label(label)
    ep = sp.Symbol("ep")
    x4 = fam.variables.index("x4")
    g2 = fam.poly_names.index("G2")
    assert sp.simplify(e[x4].to_sympy() - (1 - ep)) == 0  # (-ep) + 1
    assert sp.simplify(f[g2].to_sympy() - (-1 - ep)) == 0  # (-2-ep) + 1


def test_specialize_polynomials_is_real_partial():
    fam = parse_family_text(load_example("d4_explicit_family.wl.txt"))
    spec = fam.specialize_polynomials({"ep": 3, "r": 5}, PRIME)
    assert spec is not None
    # G2 = r*x1*x2 + ... ; at r=5 the x1*x2 coefficient is 5 mod p.
    x = {v: fam.variables.index(v) for v in fam.variables}
    x1x2 = tuple(1 if i in (x["x1"], x["x2"]) else 0 for i in range(fam.nvars))
    assert spec["G2"][x1x2] == 5


def test_specialize_full_is_honest_stub():
    fam = parse_family_text(load_example("d4_explicit_family.wl.txt"))
    with pytest.raises(NotImplementedError):
        fam.specialize({"ep": 3, "r": 5}, PRIME)
