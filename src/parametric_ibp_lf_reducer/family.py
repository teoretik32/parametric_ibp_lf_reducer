"""Container for a parsed parametric integral family plus label-projection helpers.

Mathematical meaning of the stored data (spec §2):

    F_base = prod_i x_i^(monomial_exponents[i]) * prod_l G_l^(polynomial_exponents[l])
             * target_multiplier
    F_label = F_base * prod_i x_i^(n_i) * prod_l G_l^(m_l),   label = (n_1..n_N, m_1..m_M)

Pass 1B adds the "which family member does a label denote" projections (``label_to_factor``,
``exponent_at_label``) and a Wolfram-like relative-factor rendering. It does NOT add
valuations, surface tests, local-finiteness decisions, or the finite-field specialization used
by row generation — those are later passes. ``specialize`` is therefore an honest stub.
"""

from __future__ import annotations

from dataclasses import dataclass

from .coefficients import ParamExpr
from .labels import Label, split_label
from .sparse_poly import SparsePoly


@dataclass(frozen=True)
class IntegrandFactor:
    """The relative factor ``prod_i x_i^(n_i) * prod_l G_l^(m_l)`` denoted by a label.

    Stores integer powers only (the shift relative to the base integrand). Renders itself to
    Wolfram-like text (``^`` for powers, denominator grouped in parentheses when needed), e.g.
    ``x2*x3/(G0^2*G1)`` or ``x4*x7*x8/G3``.
    """

    variables: tuple[str, ...]
    poly_names: tuple[str, ...]
    monomial_powers: tuple[int, ...]  # n_i, aligned with variables
    poly_powers: tuple[int, ...]  # m_l, aligned with poly_names

    @staticmethod
    def _term(name: str, power: int) -> str:
        return name if power == 1 else f"{name}^{power}"

    def to_wolfram_text(self) -> str:
        num: list[str] = []
        den: list[str] = []
        for name, p in zip(self.variables, self.monomial_powers):
            if p > 0:
                num.append(self._term(name, p))
            elif p < 0:
                den.append(self._term(name, -p))
        for name, p in zip(self.poly_names, self.poly_powers):
            if p > 0:
                num.append(self._term(name, p))
            elif p < 0:
                den.append(self._term(name, -p))
        num_s = "*".join(num) if num else "1"
        if not den:
            return num_s
        den_s = "*".join(den)
        if len(den) > 1:
            den_s = f"({den_s})"
        return f"{num_s}/{den_s}"


@dataclass
class ParametricFamily:
    variables: tuple[str, ...]
    parameters: tuple[str, ...]
    regulators: tuple[str, ...]
    domain: str
    poly_names: tuple[str, ...]
    polynomials: dict[str, SparsePoly]
    monomial_exponents: tuple[ParamExpr, ...]  # aligned with ``variables``
    polynomial_exponents: tuple[ParamExpr, ...]  # aligned with ``poly_names``
    target_multiplier: ParamExpr
    assumptions: tuple[str, ...] = ()
    options: dict | None = None

    def __post_init__(self) -> None:
        if self.options is None:
            self.options = {}
        if len(self.monomial_exponents) != len(self.variables):
            raise ValueError("monomial_exponents must align with variables")
        if len(self.polynomial_exponents) != len(self.poly_names):
            raise ValueError("polynomial_exponents must align with poly_names")
        if set(self.poly_names) != set(self.polynomials):
            raise ValueError("poly_names must match polynomials keys")

    # ---- sizes / accessors --------------------------------------------------
    @property
    def nvars(self) -> int:
        return len(self.variables)

    @property
    def npolys(self) -> int:
        return len(self.poly_names)

    def polynomial(self, name: str) -> SparsePoly:
        return self.polynomials[name]

    def ordered_polynomials(self) -> list[SparsePoly]:
        return [self.polynomials[name] for name in self.poly_names]

    def _split(self, label: Label) -> tuple[Label, Label]:
        return split_label(label, self.nvars, self.npolys)

    # ---- label projections (Pass 1B) ---------------------------------------
    def label_to_factor(self, label: Label) -> IntegrandFactor:
        """Return the relative integrand factor (integer powers) denoted by ``label``."""
        n, m = self._split(label)
        return IntegrandFactor(self.variables, self.poly_names, n, m)

    def label_to_wolfram_text(self, label: Label) -> str:
        return self.label_to_factor(label).to_wolfram_text()

    def exponent_at_label(self, label: Label) -> tuple[tuple[ParamExpr, ...], tuple[ParamExpr, ...]]:
        """Return ``(e_i, f_l)`` = full parametric exponents at ``label``.

        ``e_i = monomial_exponents[i] + n_i`` and ``f_l = polynomial_exponents[l] + m_l``.
        """
        n, m = self._split(label)
        e = tuple(
            self.monomial_exponents[i] + ParamExpr.from_int(n[i], self.parameters)
            for i in range(self.nvars)
        )
        f = tuple(
            self.polynomial_exponents[j] + ParamExpr.from_int(m[j], self.parameters)
            for j in range(self.npolys)
        )
        return e, f

    # ---- finite-field specialization ---------------------------------------
    def specialize_polynomials(self, sample: dict, prime: int) -> dict[str, dict] | None:
        """Specialize only the ``G_l`` polynomials' parameters modulo ``prime``.

        Returns ``{name: {monomial: int}}`` or ``None`` if the point is bad for some polynomial
        (a coefficient denominator vanishes mod ``prime``). This is a genuine partial building
        block; it does NOT specialize exponents (those feed row generation, a later pass).
        """
        out: dict[str, dict] = {}
        for name in self.poly_names:
            spec = self.polynomials[name].eval_mod_p(sample, prime)
            if spec is None:
                return None
            out[name] = spec
        return out

    def specialize(self, sample: dict, prime: int):
        """Full finite-field specialization of the family (exponents + polynomials).

        Not implemented in Pass 1B: this feeds row generation / the modular linear-algebra
        layer, which arrive later. Raising is deliberate — no fake ``FamilyModP`` is returned.
        Use :meth:`specialize_polynomials` for the currently-available partial specialization.
        """
        raise NotImplementedError(
            "full family specialization (exponents + polynomials -> FamilyModP) arrives with "
            "the modular row-generation layer; use specialize_polynomials(...) for now"
        )
