"""Wolfram-like textual export helpers.

This is a *text* format only (spec §5.11): powers use ``^`` (never Python ``**``), rational
numbers print as ``p/q`` (never decimals), coefficients are factorized, integrands and labels
render in Wolfram-like list/fraction syntax. SymPy is used only for the final
factorization/pretty-printing here, never in a hot loop.
"""

from __future__ import annotations

import sympy as sp

from .coefficients import ParamExpr
from .family import IntegrandFactor
from .labels import Label


def sympy_to_wolfram_text(expr) -> str:
    """Render a SymPy expression as Wolfram-like text (``^`` for powers, exact rationals)."""
    text = sp.sstr(sp.sympify(expr))
    return text.replace("**", "^")


def coeff_to_wolfram_text(coeff, factor: bool = True) -> str:
    """Render a coefficient (``ParamExpr`` or SymPy expr) as factorized Wolfram-like text."""
    expr = coeff.to_sympy() if isinstance(coeff, ParamExpr) else sp.sympify(coeff)
    if factor:
        expr = sp.factor(expr)
    return sympy_to_wolfram_text(expr)


def integrand_to_wolfram_text(factor: IntegrandFactor) -> str:
    """Render a relative integrand factor, e.g. ``x2*x3/(G0^2*G1)``."""
    return factor.to_wolfram_text()


def label_to_wolfram_text(label: Label) -> str:
    """Render a label as a Wolfram-like integer list, e.g. ``{0,0,1,-2,-1}``."""
    return "{" + ",".join(str(int(x)) for x in label) + "}"
