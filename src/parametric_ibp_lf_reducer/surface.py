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
"""

from __future__ import annotations

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


def _label_exps_symbolic(family: ParametricFamily, label: Label):
    e, f = family.exponent_at_label(label)
    return [pe.to_sympy() for pe in e], [pe.to_sympy() for pe in f]


def coordinate_primitive_surface_free(
    family: ParametricFamily,
    label: Label,
    var_index: int,
    multiplier_exps=None,
    eps_direction: str = "minus",
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

    s_zero = regulated_sign(exp_zero, family.regulators, eps_direction)
    s_inf = regulated_sign(exp_inf, family.regulators, eps_direction)
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
                sign = regulated_sign(score, family.regulators, eps_direction)
                if sign in ("neg", "zero"):
                    return False
                if sign == "unknown":
                    saw_unknown = True
    return "Unknown" if saw_unknown else True
