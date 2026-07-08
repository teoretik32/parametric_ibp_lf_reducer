"""Tangent (logarithmic) vector fields via a finite-degree syzygy ansatz (spec §5.7, §3).

A polynomial vector field ``Q = (Q_1, ..., Q_N)`` is *tangent* to the divisor ``G_1...G_M = 0``
when, for every polynomial ``G_l``, there is a polynomial ``H_l`` with

    Q . grad G_l = sum_i Q_i d_i G_l = H_l G_l.

Such fields generate ``div(Q F)`` IBP relations that do NOT shift the ``m`` indices, because the
``1/G_l`` from the logarithmic derivative cancels. This module only *finds and verifies* the
fields; turning them into rows is a later pass.

Solver (MVP, SymPy backend, setup phase only — never a hot loop): for a degree block
``(d_Q, d_H)`` write each ``Q_i`` (degree <= d_Q) and ``H_l`` (degree <= d_H) with unknown
coefficients, impose ``sum_i Q_i d_iG_l - H_l G_l == 0`` as a polynomial identity in ``x``,
and take the nullspace of the resulting homogeneous linear system over the field ``Q(params)``.
A field is accepted only if it is tangent to *all* ``G_l``; the zero field is dropped and
fields proportional (by a parameter-only scalar) are deduplicated.

SymPy MVP limitations (documented, see notes/assumptions.md A21): the nullspace is computed for
*generic* parameter values; specializations where a leading coefficient vanishes are not treated
separately. Only single-block dense ansatz is used (no Singular/Sage syzygy backend yet).
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from dataclasses import dataclass

import sympy as sp

from .coefficients import ParamExpr
from .family import ParametricFamily
from .sparse_poly import SparsePoly


@dataclass(frozen=True)
class TangentField:
    """A vector field ``Q`` with its multipliers ``H_l`` (``Q.grad G_l = H_l G_l``)."""

    components: tuple[SparsePoly, ...]  # Q_i, aligned with variables
    multipliers: tuple[SparsePoly, ...]  # H_l, aligned with poly_names
    degree_block: tuple[int, int]

    @property
    def is_zero(self) -> bool:
        return all(q.is_zero for q in self.components)

    def tangency_defect(self, family: ParametricFamily) -> tuple[SparsePoly, ...]:
        """Return ``(sum_i Q_i d_iG_l - H_l G_l)`` per polynomial; all zero iff tangent."""
        minus_one = ParamExpr.from_int(-1, family.parameters)
        defects = []
        for li, name in enumerate(family.poly_names):
            gl = family.polynomials[name]
            acc = SparsePoly.zero(family.nvars, family.parameters)
            for i in range(family.nvars):
                acc = acc.add(self.components[i].mul(gl.derivative(i)))
            acc = acc.add(self.multipliers[li].mul(gl).scalar_mul(minus_one))
            defects.append(acc)
        return tuple(defects)

    def is_tangent(self, family: ParametricFamily) -> bool:
        return all(d.is_zero for d in self.tangency_defect(family))


def _monomials(nvars: int, max_degree: int) -> list[tuple[int, ...]]:
    def rec(pos: int, remaining: int) -> Iterator[tuple[int, ...]]:
        if pos == nvars - 1:
            for e in range(remaining + 1):
                yield (e,)
            return
        for e in range(remaining + 1):
            for rest in rec(pos + 1, remaining - e):
                yield (e, *rest)

    return list(rec(0, max_degree))


def _mono_expr(xsyms, monom):
    expr = sp.Integer(1)
    for x, e in zip(xsyms, monom):
        if e:
            expr *= x**e
    return expr


def _coeff_dict(q_sympy, xsyms):
    """Map ``(component_index, x-monomial) -> coefficient`` for a field's SymPy components."""
    out = {}
    for i, expr in enumerate(q_sympy):
        e = sp.expand(expr)
        if e == 0:
            continue
        for monom, coeff in sp.Poly(e, *xsyms).terms():
            out[(i, tuple(monom))] = coeff
    return out


def _proportional(qa, qb, xsyms) -> bool:
    """True if the two SymPy fields are proportional by a parameter-only scalar."""
    da, db = _coeff_dict(qa, xsyms), _coeff_dict(qb, xsyms)
    keys = list(set(da) | set(db))
    va = [da.get(k, sp.Integer(0)) for k in keys]
    vb = [db.get(k, sp.Integer(0)) for k in keys]
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            if sp.simplify(va[i] * vb[j] - va[j] * vb[i]) != 0:
                return False
    return True


def verify_tangent(
    family: ParametricFamily, q_polys: Sequence[SparsePoly]
) -> tuple[bool, tuple[SparsePoly, ...] | None]:
    """Check whether ``q_polys`` is tangent to every ``G_l``; return ``(ok, H)`` if so.

    ``ok`` is ``True`` iff ``Q.grad G_l`` is divisible by ``G_l`` (as polynomials in ``x``, with
    coefficients rational in the parameters) for all ``l``; then ``H_l = Q.grad G_l / G_l``.
    """
    if len(q_polys) != family.nvars:
        raise ValueError(f"expected {family.nvars} field components, got {len(q_polys)}")
    variables = family.variables
    xsyms = [sp.Symbol(v) for v in variables]
    hs: list[SparsePoly] = []
    for name in family.poly_names:
        gl = family.polynomials[name]
        acc = SparsePoly.zero(family.nvars, family.parameters)
        for i in range(family.nvars):
            acc = acc.add(q_polys[i].mul(gl.derivative(i)))
        ratio = sp.cancel(acc.to_sympy(variables) / gl.to_sympy(variables))
        if sp.denom(ratio).free_symbols & set(xsyms):
            return False, None  # x remains in the denominator -> not divisible by G_l
        hs.append(SparsePoly.from_sympy(ratio, variables, family.parameters))
    return True, tuple(hs)


def generate_tangent_fields(
    family: ParametricFamily,
    degree_blocks: Sequence[tuple[int, int]],
    dedup: bool = True,
) -> list[TangentField]:
    """Find tangent vector fields for the given ``(d_Q, d_H)`` degree blocks.

    Returns fields that are tangent to ALL polynomials, with the zero field dropped and
    parameter-scalar-proportional duplicates removed.
    """
    variables = family.variables
    params = family.parameters
    nvars = family.nvars
    xsyms = [sp.Symbol(v) for v in variables]
    g_sym = [family.polynomials[name].to_sympy(variables) for name in family.poly_names]
    dg_sym = [[sp.diff(g, x) for x in xsyms] for g in g_sym]

    raw_fields = []  # (Q_sympy, H_sympy, block)
    for d_q, d_h in degree_blocks:
        q_mons = _monomials(nvars, d_q)
        h_mons = _monomials(nvars, d_h)
        q_expr, q_syms = [], []
        for i in range(nvars):
            expr = sp.Integer(0)
            for m in q_mons:
                s = sp.Symbol(f"__q_{i}_{'_'.join(map(str, m))}")
                q_syms.append(s)
                expr += s * _mono_expr(xsyms, m)
            q_expr.append(expr)
        h_expr, h_syms = [], []
        for li in range(family.npolys):
            expr = sp.Integer(0)
            for m in h_mons:
                s = sp.Symbol(f"__h_{li}_{'_'.join(map(str, m))}")
                h_syms.append(s)
                expr += s * _mono_expr(xsyms, m)
            h_expr.append(expr)

        unknowns = q_syms + h_syms
        eqs = []
        for li in range(family.npolys):
            residual = sum(q_expr[i] * dg_sym[li][i] for i in range(nvars)) - h_expr[li] * g_sym[li]
            for _, coeff in sp.Poly(sp.expand(residual), *xsyms).terms():
                eqs.append(coeff)
        if not eqs:
            continue
        a_mat, _ = sp.linear_eq_to_matrix(eqs, unknowns)
        for vec in a_mat.nullspace():
            sub = {u: vec[k] for k, u in enumerate(unknowns)}
            q_sol = tuple(sp.expand(qi.subs(sub)) for qi in q_expr)
            h_sol = tuple(sp.expand(hi.subs(sub)) for hi in h_expr)
            raw_fields.append((q_sol, h_sol, (d_q, d_h)))

    # Drop the zero field.
    raw_fields = [f for f in raw_fields if any(qi != 0 for qi in f[0])]
    # Deduplicate parameter-scalar-proportional fields.
    if dedup:
        kept = []
        for f in raw_fields:
            if not any(_proportional(f[0], g[0], xsyms) for g in kept):
                kept.append(f)
        raw_fields = kept

    out: list[TangentField] = []
    for q_sol, h_sol, block in raw_fields:
        field = TangentField(
            components=tuple(SparsePoly.from_sympy(qi, variables, params) for qi in q_sol),
            multipliers=tuple(SparsePoly.from_sympy(hi, variables, params) for hi in h_sol),
            degree_block=block,
        )
        if field.is_tangent(family):  # guard: keep only verified fields
            out.append(field)
    return out
