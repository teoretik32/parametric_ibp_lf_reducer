"""Surface-free filters for IBP rows (spec §7; corrected by External Int2 Method.13).

An IBP identity ``0 = integral of a total derivative`` is only valid *for the integral* when the
primitive/flux contributes nothing on the boundary of the domain. This module decides that with
the exact toric criteria derived in ``docs/TORIC_SURFACE_VALIDATION.md``:

- ``coordinate_primitive_surface_free`` — for a coordinate primitive ``d/dx_i (P F)`` the
  boundary term at the ``x_i`` facets is an ``(N-1)``-dimensional *transverse integral*, not a
  pointwise limit. It vanishes iff for EVERY boundary ray ``d`` of the complete polyhedral set
  with ``d_i != 0`` the facet score ``score(P F, d) - d_i`` is strictly positive, where ``score``
  is the full-measure scaling score. On the pure coordinate rays ``+-e_i`` this reduces exactly
  to the historical component-local exponent test; the mixed rays are the new content. The
  pre-Method.13 component-local acceptance admitted invalid rows (the External Int2 Method.12R
  counterexample) and is gone.
- ``vector_field_surface_free`` — for a vector/tangent primitive ``div(Q F)`` the flux across
  the face of ray ``d`` is governed by the *normal component* ``N_d = sum_i d_i Q_i / x_i``.
  ``N_d`` is assembled with exact coefficient arithmetic (terms from different components may
  cancel EXACTLY); every surviving monomial ``x^m`` must satisfy ``score(x^m F, d) > 0``.
  Components with ``d_i = 0`` are tangential to the face and contribute nothing.

Strictness (Method.13 items 8): exact arithmetic only; an exactly-zero score fails the row; an
unresolved symbolic sign gives ``"Unknown"``; an incomplete ray enumeration (budget) can never
produce ``True`` — the verdict degrades to ``"Unknown"``.

Both filters accept an explicit :class:`SurfacePolicy` (regulated-limit reading vs an exact
rational chamber point); the default is the historical ``epsilon -> 0^-`` limit reading. Both
policies use the SAME complete ray set — only the sign reading differs.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction

import sympy as sp

from .family import ParametricFamily
from .labels import Label
from .sparse_poly import SparsePoly
from .valuations import boundary_ray_data, score_from_exponents


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


# --- per-family exact score machinery (Method.13) ----------------------------------------------
# The scaling score along ray ``d`` of ``x^p * F_label`` decomposes as
#
#     score = S0(d) + sum_k d_k*(n_k + p_k) + sum_l val_d(G_l)*m_l
#
# where ``S0(d)`` is the (possibly regulator-dependent) score of the ZERO label and everything
# else is exact integer arithmetic. ``S0(d)`` and the per-ray polynomial valuations are cached on
# the family; sign decisions are memoized per (policy, ray, integer shift), keyed additionally by
# the identity of the module-level ``regulated_sign`` so a diagnostic monkeypatch of that function
# (Method.9-style) is never served a stale verdict.


def _surface_cache(family: ParametricFamily) -> dict:
    cache = family.__dict__.get("_surface_filter_cache")
    if cache is None:
        dirs, complete = boundary_ray_data(family)
        zero = tuple([0] * (family.nvars + family.npolys))
        e0, f0 = _label_exps_symbolic(family, zero)
        vals = {
            d: tuple(family.polynomials[name].valuation(d) for name in family.poly_names)
            for d in dirs
        }
        base = {d: sp.expand(score_from_exponents(e0, f0, family, d)) for d in dirs}
        cache = {
            "dirs": dirs,
            "complete": complete,
            "vals": vals,
            "base": base,
            "sign_memo": {},
        }
        object.__setattr__(family, "_surface_filter_cache", cache)
    return cache


def _direction_vals(family: ParametricFamily, cache: dict, d) -> tuple:
    vals = cache["vals"].get(d)
    if vals is None:
        vals = tuple(family.polynomials[name].valuation(d) for name in family.poly_names)
        cache["vals"][d] = vals
    return vals


def _base_expr(family: ParametricFamily, cache: dict, d):
    expr = cache["base"].get(d)
    if expr is None:
        zero = tuple([0] * (family.nvars + family.npolys))
        e0, f0 = _label_exps_symbolic(family, zero)
        expr = sp.expand(score_from_exponents(e0, f0, family, d))
        cache["base"][d] = expr
    return expr


def _score_sign(family: ParametricFamily, cache: dict, pol: SurfacePolicy, d, shift: int) -> str:
    """Sign of ``S0(d) + shift`` under ``pol`` — memoized, monkeypatch-safe (see above)."""
    patch_id = id(regulated_sign) if pol.mode == "limit" else 0
    key = (pol.mode, pol.direction, pol.point_items, patch_id, d, shift)
    memo = cache["sign_memo"]
    sign = memo.get(key)
    if sign is None:
        sign = pol.sign_at(_base_expr(family, cache, d) + shift, family.regulators)
        memo[key] = sign
    return sign


def coordinate_primitive_surface_free(
    family: ParametricFamily,
    label: Label,
    var_index: int,
    multiplier_exps=None,
    eps_direction: str = "minus",
    policy: SurfacePolicy | None = None,
):
    """Is the coordinate primitive ``P * F_label`` surface-free at the ``x_i`` facets?

    ``multiplier_exps`` is the monomial ``P = prod_k x_k^(p_k)`` (a tuple of length ``nvars``;
    ``None`` means ``P = 1``). The boundary term of ``integral d/dx_i (P F)`` at ``x_i -> 0`` /
    ``x_i -> infinity`` is the ``(N-1)``-dimensional transverse integral of ``P F`` over the
    facet, so pointwise decay in ``x_i`` alone is NOT sufficient (Method.12R counterexample).
    Exact criterion (derived in ``docs/TORIC_SURFACE_VALIDATION.md``): for every boundary ray
    ``d`` of the complete polyhedral set with ``d_i != 0``,

        surface_score(d, i) = score(P F, d) - d_i  must be strictly > 0,

    where ``score`` is the full-measure scaling score. On ``+-e_i`` this is exactly the
    historical ``exp_zero > 0`` / ``exp_inf < 0`` component test; the mixed rays additionally
    demand convergence-with-decay of the transverse facet integral.

    Returns ``True`` / ``False`` / ``"Unknown"``. Exactly-zero scores fail (marginal boundary
    term); undecidable symbolic signs and an over-budget (incomplete) ray enumeration give
    ``"Unknown"`` — never a possibly-wrong ``True``.
    """
    if not 0 <= var_index < family.nvars:
        raise IndexError(f"var_index {var_index} out of range")
    pol = _resolve_policy(policy, eps_direction)
    cache = _surface_cache(family)
    nvars = family.nvars
    n = label[:nvars]
    m = label[nvars:]
    p = (0,) * nvars if multiplier_exps is None else tuple(int(x) for x in multiplier_exps)

    saw_unknown = False
    for d in cache["dirs"]:
        di = d[var_index]
        if di == 0:
            continue  # the ray stays inside the x_i facet family's transverse directions
        vals = cache["vals"][d]
        shift = (
            sum(dk * (nk + pk) for dk, nk, pk in zip(d, n, p))
            + sum(v * ml for v, ml in zip(vals, m))
            - di
        )
        sign = _score_sign(family, cache, pol, d, shift)
        if sign in ("neg", "zero"):
            return False
        if sign == "unknown":
            saw_unknown = True
    if saw_unknown or not cache["complete"]:
        return "Unknown"
    return True


def vector_field_surface_free(
    family: ParametricFamily,
    label: Label,
    vector_field: list[SparsePoly],
    eps_direction: str = "minus",
    rays=None,
    policy: SurfacePolicy | None = None,
):
    """Is ``div(Q F_label)`` surface-free, i.e. does its normal flux vanish on all toric faces?

    ``vector_field`` is ``Q = (Q_1, ..., Q_N)`` as one :class:`SparsePoly` per variable. For each
    boundary ray ``d`` the flux across the corresponding face is governed by the exact *normal
    component* (``docs/TORIC_SURFACE_VALIDATION.md``)

        N_d(x) = sum_{i : d_i != 0} d_i * Q_i(x) / x_i,

    assembled with exact :class:`ParamExpr` coefficient arithmetic so that terms contributed by
    different components may cancel EXACTLY (and only exactly — no modular or floating
    cancellation). Every surviving monomial ``x^m`` of ``N_d`` must be boundary-suppressed:
    ``score(x^m F, d) > 0`` with the full-measure scaling score. A component with ``d_i = 0`` is
    tangential to the face and imposes no condition on that ray. A row is therefore NOT rejected
    merely because one component looks marginal in isolation, and NOT accepted merely because
    each component separately passes a per-term heuristic.

    ``rays=None`` uses the family's complete polyhedral boundary set (with the completeness flag
    honored: over budget -> ``"Unknown"``); an explicit ``rays`` list restricts the test to those
    directions (caller-owned scope, e.g. diagnostics).

    Returns ``True`` / ``False`` / ``"Unknown"`` with the same strictness rules as the
    coordinate filter: zero scores fail, undecidable signs give ``"Unknown"``.
    """
    if len(vector_field) != family.nvars:
        raise ValueError(f"vector_field must have {family.nvars} components")
    pol = _resolve_policy(policy, eps_direction)
    cache = _surface_cache(family)
    nvars = family.nvars
    n = label[:nvars]
    m = label[nvars:]
    if rays is None:
        directions = cache["dirs"]
        complete = cache["complete"]
    else:
        directions = [tuple(int(x) for x in d) for d in rays]
        complete = True  # caller-supplied scope is taken as authoritative for that scope

    saw_unknown = False
    for d in directions:
        vals = _direction_vals(family, cache, d)
        label_shift = sum(dk * nk for dk, nk in zip(d, n)) + sum(
            v * ml for v, ml in zip(vals, m)
        )
        # exact normal component N_d = sum_i d_i Q_i / x_i
        terms: dict[tuple, object] = {}
        for i in range(nvars):
            if d[i] == 0:
                continue
            qi = vector_field[i]
            if qi.is_zero:
                continue
            for c, coeff in qi.terms.items():
                mono = tuple(c[k] - (1 if k == i else 0) for k in range(nvars))
                contrib = coeff.scale_int(d[i])
                terms[mono] = terms[mono] + contrib if mono in terms else contrib
        for mono, coeff in terms.items():
            if coeff.is_zero:
                continue  # exact cancellation between components: no flux from this monomial
            shift = label_shift + sum(dk * mk for dk, mk in zip(d, mono))
            sign = _score_sign(family, cache, pol, d, shift)
            if sign in ("neg", "zero"):
                return False
            if sign == "unknown":
                saw_unknown = True
    if saw_unknown or not complete:
        return "Unknown"
    return True
