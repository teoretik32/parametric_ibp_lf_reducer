"""Surface-free filters for IBP rows (spec §7, method review §4).

An IBP identity ``0 = integral of a total derivative`` is only valid *for the integral* when the
primitive/flux contributes nothing on the boundary of the domain. This module decides that,
conservatively, in the regulated region ``epsilon -> 0^-`` (or ``0^+``):

- ``coordinate_primitive_surface_free`` — for a coordinate primitive ``d/dx_i (P F)`` it checks
  ONLY the two boundaries of that component, ``x_i = 0`` and ``x_i = infinity``. It deliberately
  does NOT demand vanishing along every toric ray (that would be over-strict and drop valid
  rows — spec §7.1).
- ``vector_field_surface_free`` — for a vector/tangent primitive ``div(Q F)`` it checks the
  normal flux across toric boundary rays.

Row generation itself is NOT here (later pass). Whenever assumptions are insufficient to decide
a sign, these return ``"Unknown"`` rather than a possibly-wrong ``True``.

Both filters accept an explicit :class:`SurfacePolicy` (regulated-limit reading vs an exact
rational chamber point); the default is the historical ``epsilon -> 0^-`` limit reading.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction

import sympy as sp

from .family import ParametricFamily
from .labels import Label
from .sparse_poly import SparsePoly
from .valuations import compute_candidate_rays, score_from_exponents


def regulated_sign(expr, regulators, direction: str = "minus") -> str:
    """Sign of ``expr`` in the limit ``epsilon -> 0^-`` (``minus``) or ``0^+`` (``plus``).

    Returns ``"pos"``, ``"neg"``, ``"zero"`` (marginal to first order) or ``"unknown"`` (a
    non-regulator parameter prevents a decision). Only the leading behaviour at the regulator's
    limit is used; if the value at ``epsilon = 0`` is nonzero it decides directly, otherwise the
    first-order coefficient and the chosen direction break the tie.
    """
    expr = sp.expand(sp.sympify(expr))
    reg_syms = {sp.Symbol(r) for r in regulators}
    if expr.free_symbols - reg_syms:
        return "unknown"
    if not reg_syms:
        if expr.is_positive:
            return "pos"
        if expr.is_negative:
            return "neg"
        return "zero"
    # Single primary regulator (MVP): use its limit.
    eps = next(iter(reg_syms))
    val0 = expr.subs(eps, 0)
    if val0.is_positive:
        return "pos"
    if val0.is_negative:
        return "neg"
    slope = sp.diff(expr, eps).subs(eps, 0)
    if slope == 0:
        return "zero"
    if direction == "minus":  # epsilon < 0, so sign(value) = sign(-slope)
        return "pos" if slope.is_negative else "neg"
    return "pos" if slope.is_positive else "neg"


@dataclass(frozen=True)
class SurfacePolicy:
    """Explicit sign policy for the surface-free filters (External Int2 Method.11, Phase A).

    Modes:

    - ``limit`` — signs are read from the leading behaviour in the regulated limit
      ``epsilon -> 0^-`` (``direction="minus"``) or ``0^+`` (``"plus"``) via
      :func:`regulated_sign`.  This is the package default (``minus``) and encodes the
      analytic-continuation convention used throughout: every IBP identity is proved in the
      regulated region on the chosen side of ``epsilon = 0`` where boundary terms vanish, and
      the resulting relation is then continued analytically in the regulators.
    - ``chamber`` — signs are decided by *exact rational substitution* of the regulators at an
      explicit chamber point, e.g. ``SurfacePolicy.chamber({"ep": Fraction(-3, 5)})``.  The
      identity set is the one valid in the convergence chamber containing that point; values of
      the reduction at other regulator values again follow by analytic continuation from that
      chamber.  Strictness: only a strictly positive / strictly negative exact value decides a
      sign; an exact zero reports ``"zero"`` and any leftover free symbol (non-regulator
      parameter, or a regulator missing from the point) reports ``"unknown"`` — both cause the
      row to be rejected, never a possibly-wrong ``True``.

    Instances are immutable.  Chamber values are exact :class:`~fractions.Fraction`; floats are
    rejected at construction (no floating-point sign decisions anywhere).
    """

    mode: str
    direction: str = "minus"
    point_items: tuple[tuple[str, Fraction], ...] = ()

    def __post_init__(self):
        if self.mode not in ("limit", "chamber"):
            raise ValueError(f"unknown SurfacePolicy mode {self.mode!r}")
        if self.direction not in ("minus", "plus"):
            raise ValueError(f"unknown limit direction {self.direction!r}")
        if self.mode == "chamber" and not self.point_items:
            raise ValueError("chamber policy requires at least one regulator value")
        if self.mode == "limit" and self.point_items:
            raise ValueError("limit policy takes no chamber point")

    @classmethod
    def limit(cls, direction: str = "minus") -> SurfacePolicy:
        """Regulated-limit policy (package default is ``direction="minus"``)."""
        return cls(mode="limit", direction=direction)

    @classmethod
    def chamber(cls, point: Mapping) -> SurfacePolicy:
        """Exact-rational chamber-point policy; values must be Fraction/int/string (no floats)."""
        items = []
        for name in sorted(point):
            value = point[name]
            if isinstance(value, float):
                raise TypeError(
                    f"chamber value for {name!r} must be exact (Fraction/int/str), got float"
                )
            items.append((str(name), Fraction(value)))
        return cls(mode="chamber", point_items=tuple(items))

    @property
    def point(self) -> dict[str, Fraction]:
        """Chamber point as a fresh dict (empty for limit mode)."""
        return dict(self.point_items)

    def describe(self) -> dict:
        """JSON-safe diagnostics record: mode plus direction or exact point values."""
        if self.mode == "limit":
            return {"mode": "limit", "direction": self.direction}
        return {"mode": "chamber", "point": {k: str(v) for k, v in self.point_items}}

    def sign_at(self, expr, regulators) -> str:
        """Sign of ``expr`` under this policy: ``"pos" | "neg" | "zero" | "unknown"``."""
        if self.mode == "limit":
            return regulated_sign(expr, regulators, self.direction)
        expr = sp.expand(sp.sympify(expr))
        subs = {
            sp.Symbol(name): sp.Rational(val.numerator, val.denominator)
            for name, val in self.point_items
        }
        val = sp.expand(expr.subs(subs))
        if val.free_symbols:
            return "unknown"
        if val.is_positive:
            return "pos"
        if val.is_negative:
            return "neg"
        if val.is_zero:
            return "zero"
        return "unknown"


def _resolve_policy(policy: SurfacePolicy | None, eps_direction: str) -> SurfacePolicy:
    return policy if policy is not None else SurfacePolicy.limit(eps_direction)


def _label_exps_symbolic(family: ParametricFamily, label: Label):
    e, f = family.exponent_at_label(label)
    return [pe.to_sympy() for pe in e], [pe.to_sympy() for pe in f]


def coordinate_primitive_surface_free(
    family: ParametricFamily,
    label: Label,
    var_index: int,
    multiplier_exps=None,
    eps_direction: str = "minus",
    policy: SurfacePolicy | None = None,
):
    """Is the coordinate primitive ``P * F_label`` surface-free at ``x_i = 0`` and ``x_i = inf``?

    ``multiplier_exps`` is the monomial ``P = prod_k x_k^(p_k)`` (a tuple of length ``nvars``;
    ``None`` means ``P = 1``). Only the ``x_i`` component matters at this component's boundaries:

        exp at x_i -> 0   = p_i + e_i + sum_l f_l * min_power_i(G_l)   must be > 0
        exp at x_i -> inf = p_i + e_i + sum_l f_l * max_power_i(G_l)   must be < 0

    Returns ``True`` / ``False`` / ``"Unknown"``. This is intentionally component-local: it does
    NOT require vanishing along mixed toric rays.
    """
    if not 0 <= var_index < family.nvars:
        raise IndexError(f"var_index {var_index} out of range")
    e_syms, f_syms = _label_exps_symbolic(family, label)
    p_i = 0 if multiplier_exps is None else int(multiplier_exps[var_index])
    unit_i = tuple(1 if k == var_index else 0 for k in range(family.nvars))

    exp_zero = p_i + e_syms[var_index]
    exp_inf = p_i + e_syms[var_index]
    for j, name in enumerate(family.poly_names):
        poly = family.polynomials[name]
        exp_zero += f_syms[j] * poly.valuation(unit_i)  # min power of x_i in G_l
        exp_inf += f_syms[j] * poly.degree_in(var_index)  # max power of x_i in G_l

    pol = _resolve_policy(policy, eps_direction)
    s_zero = pol.sign_at(exp_zero, family.regulators)
    s_inf = pol.sign_at(exp_inf, family.regulators)
    if s_zero == "pos" and s_inf == "neg":
        return True
    if s_zero == "unknown" or s_inf == "unknown":
        return "Unknown"
    return False


def vector_field_surface_free(
    family: ParametricFamily,
    label: Label,
    vector_field: list[SparsePoly],
    eps_direction: str = "minus",
    rays=None,
    policy: SurfacePolicy | None = None,
):
    """Is ``div(Q F_label)`` surface-free, i.e. does its normal flux vanish on all toric rays?

    ``vector_field`` is ``Q = (Q_1, ..., Q_N)`` as one :class:`SparsePoly` per variable. For each
    toric ray ``rho`` and each monomial ``c`` of each ``Q_i`` the flux term ``x^c / x_i * F`` must
    be boundary-suppressed (positive scaling score in the regulated region). Any non-positive or
    marginal contribution fails the row; undecidable ones give ``"Unknown"``.

    Unlike the coordinate check, this uses the full set of toric candidate rays (spec §7.2).
    """
    if len(vector_field) != family.nvars:
        raise ValueError(f"vector_field must have {family.nvars} components")
    e_syms, f_syms = _label_exps_symbolic(family, label)
    directions = (
        list(rays)
        if rays is not None
        else [ray.direction for ray in compute_candidate_rays(family)]
    )
    pol = _resolve_policy(policy, eps_direction)
    saw_unknown = False
    for direction in directions:
        for i in range(family.nvars):
            qi = vector_field[i]
            if qi.is_zero:
                continue
            for c in qi.support():
                e_shift = list(e_syms)
                for k in range(family.nvars):
                    e_shift[k] = e_shift[k] + c[k]
                e_shift[i] = e_shift[i] - 1
                score = score_from_exponents(e_shift, f_syms, family, direction)
                sign = pol.sign_at(score, family.regulators)
                if sign in ("neg", "zero"):
                    return False
                if sign == "unknown":
                    saw_unknown = True
    return "Unknown" if saw_unknown else True
