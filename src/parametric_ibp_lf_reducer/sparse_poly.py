"""Sparse multivariate polynomials in the integration variables ``x_1..x_N``.

A :class:`SparsePoly` maps an integer exponent tuple (length ``nvars``) to a
:class:`~parametric_ibp_lf_reducer.coefficients.ParamExpr` coefficient (a rational function of
the external parameters). This is the canonical internal format for the ``G_l`` polynomials.

Design constraints (from spec §5.2):
- no expanded SymPy in the hot loop — arithmetic here is pure Python over ``ParamExpr``;
- monomials stored canonically as dict keys (sorted on demand);
- arbitrary monomial degrees (``x_i^2``, ``x_i^k``) with no special-casing;
- parametric coefficients specialize quickly to ``int mod p`` via :meth:`eval_mod_p`.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import sympy as sp

from .coefficients import ParamExpr

ExponentTuple = tuple[int, ...]


@dataclass
class SparsePoly:
    """Polynomial in ``nvars`` integration variables with ``ParamExpr`` coefficients."""

    nvars: int
    params: tuple[str, ...]
    terms: dict[ExponentTuple, ParamExpr] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.params = tuple(self.params)
        # Drop identically-zero coefficients so the representation stays canonical.
        self.terms = {
            m: c for m, c in self.terms.items() if not c.is_zero and len(m) == self.nvars
        }
        for m in self.terms:
            if len(m) != self.nvars:  # pragma: no cover - guarded above
                raise ValueError(f"exponent tuple {m} length != nvars={self.nvars}")

    # ---- constructors -------------------------------------------------------
    @classmethod
    def zero(cls, nvars: int, params) -> "SparsePoly":
        return cls(nvars, tuple(params), {})

    @classmethod
    def constant(cls, nvars: int, coeff: ParamExpr) -> "SparsePoly":
        return cls(nvars, coeff.params, {(0,) * nvars: coeff})

    @classmethod
    def one(cls, nvars: int, params) -> "SparsePoly":
        return cls.constant(nvars, ParamExpr.one(params))

    @classmethod
    def monomial(cls, nvars: int, exps: ExponentTuple, coeff: ParamExpr) -> "SparsePoly":
        exps = tuple(int(e) for e in exps)
        if len(exps) != nvars:
            raise ValueError(f"exponent tuple {exps} length != nvars={nvars}")
        return cls(nvars, coeff.params, {exps: coeff})

    @classmethod
    def from_sympy(cls, expr, variables, parameters) -> "SparsePoly":
        """Build from a SymPy polynomial in ``variables`` with coefficients over ``parameters``."""
        variables = [str(v) for v in variables]
        parameters = tuple(str(p) for p in parameters)
        vsyms = [sp.Symbol(v) for v in variables]
        psyms = {sp.Symbol(p) for p in parameters}
        expr = sp.expand(sp.sympify(expr))
        extra = expr.free_symbols - set(vsyms) - psyms
        if extra:
            raise ValueError(
                f"polynomial has undeclared symbols {sorted(map(str, extra))}; "
                f"variables={variables}, parameters={list(parameters)}"
            )
        terms: dict[ExponentTuple, ParamExpr] = {}
        if expr == 0:
            return cls(len(variables), parameters, {})
        poly = sp.Poly(expr, *vsyms)
        for monom, coeff in poly.terms():
            terms[tuple(int(e) for e in monom)] = ParamExpr.from_sympy(coeff, parameters)
        return cls(len(variables), parameters, terms)

    # ---- arithmetic ---------------------------------------------------------
    def _check(self, other: "SparsePoly") -> None:
        if self.nvars != other.nvars or self.params != other.params:
            raise ValueError("SparsePoly nvars/params mismatch")

    def add(self, other: "SparsePoly") -> "SparsePoly":
        self._check(other)
        out = dict(self.terms)
        for m, c in other.terms.items():
            out[m] = out[m] + c if m in out else c
        return SparsePoly(self.nvars, self.params, out)

    def scalar_mul(self, coeff: ParamExpr) -> "SparsePoly":
        if coeff.is_zero:
            return SparsePoly.zero(self.nvars, self.params)
        return SparsePoly(self.nvars, self.params, {m: c * coeff for m, c in self.terms.items()})

    def mul(self, other: "SparsePoly") -> "SparsePoly":
        self._check(other)
        out: dict[ExponentTuple, ParamExpr] = {}
        for m1, c1 in self.terms.items():
            for m2, c2 in other.terms.items():
                m = tuple(a + b for a, b in zip(m1, m2))
                prod = c1 * c2
                out[m] = out[m] + prod if m in out else prod
        return SparsePoly(self.nvars, self.params, out)

    def pow_small(self, k: int) -> "SparsePoly":
        if k < 0:
            raise ValueError("pow_small requires a non-negative exponent")
        result = SparsePoly.one(self.nvars, self.params)
        base = self
        e = k
        while e > 0:  # exponentiation by squaring
            if e & 1:
                result = result.mul(base)
            e >>= 1
            if e:
                base = base.mul(base)
        return result

    def monomial_mul(self, exps: ExponentTuple, coeff: ParamExpr | None = None) -> "SparsePoly":
        exps = tuple(int(e) for e in exps)
        if len(exps) != self.nvars:
            raise ValueError(f"exponent tuple {exps} length != nvars={self.nvars}")
        c = coeff if coeff is not None else ParamExpr.one(self.params)
        out: dict[ExponentTuple, ParamExpr] = {}
        for m, cm in self.terms.items():
            out[tuple(a + b for a, b in zip(m, exps))] = cm * c
        return SparsePoly(self.nvars, self.params, out)

    def derivative(self, var_index: int) -> "SparsePoly":
        if not 0 <= var_index < self.nvars:
            raise IndexError(f"var_index {var_index} out of range for nvars={self.nvars}")
        out: dict[ExponentTuple, ParamExpr] = {}
        for m, c in self.terms.items():
            e = m[var_index]
            if e == 0:
                continue
            new_m = m[:var_index] + (e - 1,) + m[var_index + 1 :]
            nc = c.scale_int(e)
            out[new_m] = out[new_m] + nc if new_m in out else nc
        return SparsePoly(self.nvars, self.params, out)

    # ---- queries ------------------------------------------------------------
    @property
    def is_zero(self) -> bool:
        return len(self.terms) == 0

    def support(self) -> list[ExponentTuple]:
        return sorted(self.terms.keys())

    def total_degree(self) -> int:
        if not self.terms:
            return -1
        return max(sum(m) for m in self.terms)

    def degree_in(self, var_index: int) -> int:
        if not self.terms:
            return -1
        return max(m[var_index] for m in self.terms)

    def valuation(self, ray) -> int:
        """Tropical valuation ``min_{a in support} (a . ray)`` (min-convention)."""
        if not self.terms:
            raise ValueError("valuation of the zero polynomial is undefined")
        ray = tuple(int(r) for r in ray)
        if len(ray) != self.nvars:
            raise ValueError(f"ray length {len(ray)} != nvars={self.nvars}")
        return min(sum(a * r for a, r in zip(m, ray)) for m in self.terms)

    # ---- hot-loop specialization -------------------------------------------
    def eval_mod_p(self, sample: dict, prime: int) -> dict[ExponentTuple, int] | None:
        """Specialize parameters modulo ``prime`` -> ``{monomial: int}`` (dropping zeros).

        Returns ``None`` if any coefficient denominator vanishes modulo ``prime`` (bad point).
        """
        out: dict[ExponentTuple, int] = {}
        for m, c in self.terms.items():
            v = c.eval_mod_p(sample, prime)
            if v is None:
                return None
            if v % prime != 0:
                out[m] = v % prime
        return out

    def to_sympy(self, variables):
        vsyms = [sp.Symbol(str(v)) for v in variables]
        e = sp.Integer(0)
        for m, c in self.terms.items():
            term = c.to_sympy()
            for s, ex in zip(vsyms, m):
                if ex:
                    term *= s**ex
            e += term
        return e

    def __eq__(self, other) -> bool:
        if not isinstance(other, SparsePoly):
            return NotImplemented
        return (
            self.nvars == other.nvars
            and self.params == other.params
            and self.terms == other.terms
        )
