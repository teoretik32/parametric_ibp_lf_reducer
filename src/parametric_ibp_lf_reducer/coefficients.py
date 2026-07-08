"""Parametric rational coefficient expressions with fast finite-field evaluation.

A :class:`ParamExpr` is an exact rational function of a fixed, ordered tuple of external
parameters (e.g. ``("ep", "r")``). It is built once from a SymPy expression during the
setup/parse phase; after that, all *hot-loop* work goes through :meth:`ParamExpr.eval_mod_p`,
which performs only integer arithmetic modulo a prime. No SymPy is used at evaluation time.

Internal representation: numerator and denominator are each stored as a sorted tuple of
``(monomial_exponents, integer_coefficient)`` pairs, where ``monomial_exponents`` is a tuple
of non-negative integer exponents over ``params``. Coefficients are normalized to integers
(common denominator cleared, integer content removed, denominator sign made positive), so the
representation is canonical up to a multivariate polynomial gcd (which we deliberately do not
compute here — semantic equality is available via :meth:`ParamExpr.equals`).
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd
from typing import Union

import sympy as sp

Number = Union[int, Fraction]
_Term = tuple[tuple[int, ...], Fraction]
_Terms = tuple[_Term, ...]


class _ModularZeroDivision(Exception):
    """Raised internally when a denominator vanishes modulo the working prime."""


def _lcm(a: int, b: int) -> int:
    return a // gcd(a, b) * b if a and b else (a or b)


def _clean(terms) -> list[_Term]:
    """Combine like monomials, drop zeros, return list (unsorted)."""
    acc: dict[tuple[int, ...], Fraction] = {}
    for monom, coeff in terms:
        c = Fraction(coeff)
        if c == 0:
            continue
        acc[monom] = acc.get(monom, Fraction(0)) + c
    return [(m, c) for m, c in acc.items() if c != 0]


def _canonicalize(nparams: int, num, den) -> tuple[_Terms, _Terms]:
    zmon = (0,) * nparams
    num = _clean(num)
    den = _clean(den)
    if not den:
        raise ZeroDivisionError("ParamExpr denominator is identically zero")
    if not num:
        return ((), ((zmon, Fraction(1)),))
    # Clear all fraction denominators so every coefficient becomes integer-valued.
    scale = 1
    for _, c in num + den:
        scale = _lcm(scale, c.denominator)
    numi = [(m, int(c * scale)) for m, c in num]
    deni = [(m, int(c * scale)) for m, c in den]
    # Remove common integer content across numerator and denominator.
    content = 0
    for _, c in numi + deni:
        content = gcd(content, c)
    if content == 0:
        content = 1
    numi = [(m, c // content) for m, c in numi]
    deni = [(m, c // content) for m, c in deni]
    # Make the leading denominator coefficient positive for a canonical sign.
    lead = sorted(deni)[-1][1]
    if lead < 0:
        numi = [(m, -c) for m, c in numi]
        deni = [(m, -c) for m, c in deni]
    num_t = tuple(sorted((m, Fraction(c)) for m, c in numi))
    den_t = tuple(sorted((m, Fraction(c)) for m, c in deni))
    return num_t, den_t


@dataclass(frozen=True)
class ParamExpr:
    """Exact rational function of the ordered parameter tuple ``params``."""

    params: tuple[str, ...]
    num: _Terms
    den: _Terms

    # ---- construction -------------------------------------------------------
    @classmethod
    def _make(cls, params: tuple[str, ...], num, den) -> "ParamExpr":
        n, d = _canonicalize(len(params), num, den)
        return cls(tuple(params), n, d)

    @classmethod
    def zero(cls, params) -> "ParamExpr":
        return cls._make(tuple(params), (), (((0,) * len(tuple(params)), Fraction(1)),))

    @classmethod
    def from_int(cls, k: int, params) -> "ParamExpr":
        params = tuple(params)
        zmon = (0,) * len(params)
        return cls._make(params, ((zmon, Fraction(k)),), ((zmon, Fraction(1)),))

    @classmethod
    def one(cls, params) -> "ParamExpr":
        return cls.from_int(1, params)

    @classmethod
    def from_sympy(cls, expr, params) -> "ParamExpr":
        """Build from a SymPy expression whose free symbols are a subset of ``params``."""
        params = tuple(str(p) for p in params)
        psyms = [sp.Symbol(p) for p in params]
        expr = sp.sympify(expr)
        extra = expr.free_symbols - set(psyms)
        if extra:
            raise ValueError(
                f"coefficient/exponent expression has non-parameter symbols {sorted(map(str, extra))}; "
                f"declared parameters are {list(params)}"
            )
        num_expr, den_expr = sp.fraction(sp.together(expr))
        num = _sympy_poly_terms(num_expr, psyms)
        den = _sympy_poly_terms(den_expr, psyms)
        if not den:
            raise ZeroDivisionError("coefficient expression has zero denominator")
        return cls._make(params, num, den)

    # ---- arithmetic (setup phase only; not used in hot loop) ----------------
    def _check(self, other: "ParamExpr") -> None:
        if self.params != other.params:
            raise ValueError(f"parameter mismatch: {self.params} vs {other.params}")

    def __add__(self, other: "ParamExpr") -> "ParamExpr":
        self._check(other)
        num = _padd(_pmul(self.num, other.den), _pmul(other.num, self.den))
        den = _pmul(self.den, other.den)
        return ParamExpr._make(self.params, num, den)

    def __sub__(self, other: "ParamExpr") -> "ParamExpr":
        return self + other.scale_int(-1)

    def __mul__(self, other: "ParamExpr") -> "ParamExpr":
        self._check(other)
        num = _pmul(self.num, other.num)
        den = _pmul(self.den, other.den)
        return ParamExpr._make(self.params, num, den)

    def scale_int(self, k: int) -> "ParamExpr":
        if k == 0:
            return ParamExpr.zero(self.params)
        num = tuple((m, c * k) for m, c in self.num)
        return ParamExpr._make(self.params, num, self.den)

    # ---- queries ------------------------------------------------------------
    @property
    def is_zero(self) -> bool:
        return len(self.num) == 0

    def equals(self, other: "ParamExpr") -> bool:
        """Semantic equality via SymPy (setup phase / tests only)."""
        return sp.simplify(self.to_sympy() - other.to_sympy()) == 0

    def to_sympy(self):
        psyms = [sp.Symbol(p) for p in self.params]

        def build(terms):
            e = sp.Integer(0)
            for monom, coeff in terms:
                term = sp.Rational(coeff.numerator, coeff.denominator)
                for s, ex in zip(psyms, monom):
                    if ex:
                        term *= s**ex
                e += term
            return e

        return build(self.num) / build(self.den)

    # ---- hot-loop evaluation ------------------------------------------------
    def eval_mod_p(self, sample: dict[str, Number], prime: int) -> int | None:
        """Evaluate modulo ``prime`` at ``sample`` (name -> int/Fraction).

        Returns an integer in ``[0, prime)`` or ``None`` when a denominator vanishes modulo
        ``prime`` (i.e. this ``(sample, prime)`` is a bad point that callers must skip).
        """
        try:
            vals = [_residue(sample[p], prime) for p in self.params]
            nv = _eval_terms(self.num, vals, prime)
            dv = _eval_terms(self.den, vals, prime)
        except _ModularZeroDivision:
            return None
        if dv % prime == 0:
            return None
        return (nv * pow(dv, prime - 2, prime)) % prime


# ---- module-level helpers (integer / Fraction only in the hot path) ---------
def _residue(value: Number, prime: int) -> int:
    if isinstance(value, int):
        return value % prime
    fr = Fraction(value)
    d = fr.denominator % prime
    if d == 0:
        raise _ModularZeroDivision
    return (fr.numerator % prime) * pow(d, prime - 2, prime) % prime


def _eval_terms(terms: _Terms, vals: list[int], prime: int) -> int:
    total = 0
    for monom, coeff in terms:
        d = coeff.denominator % prime
        if d == 0:
            raise _ModularZeroDivision
        c = (coeff.numerator % prime) * pow(d, prime - 2, prime) % prime
        for ex, v in zip(monom, vals):
            if ex:
                c = c * pow(v, ex, prime) % prime
        total = (total + c) % prime
    return total


def _padd(a: _Terms, b: _Terms) -> list[_Term]:
    acc: dict[tuple[int, ...], Fraction] = {}
    for monom, coeff in list(a) + list(b):
        acc[monom] = acc.get(monom, Fraction(0)) + coeff
    return [(m, c) for m, c in acc.items() if c != 0]


def _pmul(a: _Terms, b: _Terms) -> list[_Term]:
    acc: dict[tuple[int, ...], Fraction] = {}
    for m1, c1 in a:
        for m2, c2 in b:
            m = tuple(x + y for x, y in zip(m1, m2))
            acc[m] = acc.get(m, Fraction(0)) + c1 * c2
    return [(m, c) for m, c in acc.items() if c != 0]


def _sympy_poly_terms(expr, psyms) -> list[_Term]:
    """Expand a SymPy polynomial in ``psyms`` into ``[(monom, Fraction), ...]``."""
    expr = sp.expand(expr)
    if expr == 0:
        return []
    if not psyms:
        return [((), Fraction(sp.Rational(expr).p, sp.Rational(expr).q))]
    poly = sp.Poly(expr, *psyms)
    out: list[_Term] = []
    for monom, coeff in poly.terms():
        r = sp.Rational(coeff)
        out.append((tuple(int(e) for e in monom), Fraction(int(r.p), int(r.q))))
    return out
