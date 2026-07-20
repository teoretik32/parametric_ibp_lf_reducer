"""Method.3 (diagnostic only): composite locally-finite master feasibility.

Given a family whose certified normal form contains non-locally-finite terms, test
whether *linear combinations* of labels ("composite masters")

    M = sum_i c_i(params) * J(label_i)

can be locally finite even though individual members are not: along each boundary
ray the leading asymptotic layers of the members may cancel identically.

Mathematical setup (same conventions as :mod:`valuations`, spec 5.4):

* A ray is a primitive integer direction ``d``; the boundary is approached as
  ``x_i = y_i * t^(d_i)`` with ``t -> 0+`` (so ``d_i = -1`` means ``x_i -> oo``).
* Each polynomial splits as ``G_l = t^(v_l) * (A_l0 + A_l1*t + ...)`` with
  ``v_l = min_{a in supp} a . d`` (tropical valuation) and ``A_l0 != 0``.
* At ``epsilon = 0`` the label integrand is the rational function
  ``R = prod x_i^(e0_i) * prod G_l^(f0_l)``; substituting the ray chart gives the
  exact Laurent expansion ``R = t^P * (c_0 + c_1*t + ...)`` computed here by
  truncated binomial series (no floating point, no heuristics).
* The strict scaling score of :func:`valuations.base_score` is ``P + sum_i d_i``;
  layer ``j`` of a member with score ``s`` sits at absolute level ``s + j``.
  STRICT RULE (same as the certified LF gate): every level ``<= 0`` must vanish.

Generic-epsilon soundness: exponents are affine in the regulators with
label-independent slopes, so ``J(label) = Phi * R`` where ``Phi`` (the full
epsilon-slope part) is one and the same invertible ``t``-power times unit series
for every label.  Cancelling the Laurent layers of ``sum c_i R_i`` with
epsilon-free ``c_i`` therefore cancels the corresponding layers of
``sum c_i J(label_i)`` at generic epsilon as well.

Scope guard: this module is purely diagnostic.  It never touches the reducer
core, the certificate, or the LF gates; it only *reads* the family via the same
helpers the LF gate uses (``base_score``/``compute_candidate_rays``) and reports.
"""

from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from .family import ParametricFamily
from .labels import Label
from .valuations import (
    Ray,
    _random_directions,
    base_score,
    compute_candidate_rays,
    exponents_at_eps0,
)

Direction = tuple[int, ...]

STATUS_FEASIBLE = "FeasibleCompositeBasis"
STATUS_NO_COMPOSITE = "NoCompositeFoundWithinAnsatz"
STATUS_BAD_SPECIALIZATION = "BadSpecialization"


# --------------------------------------------------------------------------- #
# Phase A: exact asymptotic signatures along a ray                            #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class AsymptoticSignature:
    """Exact leading boundary behaviour of one label along one ray.

    ``t_power`` is the Laurent offset ``P`` of the epsilon=0 integrand in the ray
    chart, ``score = P + sum_i d_i`` (identical to :func:`valuations.base_score`),
    and ``coefficients[j]`` is the exact layer-``j`` coefficient — a rational
    function of the face variables (``d_i = 0``), the ray coordinates
    ``y_<var>`` (``d_i != 0``) and the kinematic parameters.
    """

    label: Label
    direction: Direction
    t_power: int
    score: int
    coefficients: tuple[sp.Expr, ...]

    @property
    def leading(self) -> sp.Expr:
        return self.coefficients[0]


def _direction_of(ray) -> Direction:
    return ray.direction if isinstance(ray, Ray) else tuple(int(x) for x in ray)


def _var_symbols(family: ParametricFamily, direction: Direction):
    """Per-variable chart symbols: face variables keep their name, scaled ones
    become ``y_<name>`` (the coordinate along the ray)."""
    syms = []
    for i, name in enumerate(family.variables):
        syms.append(sp.Symbol(name if direction[i] == 0 else f"y_{name}"))
    return tuple(syms)


def _poly_layers(family, name, direction: Direction, chart_syms, order: int):
    """``G = t^v * (B_0 + B_1 t + ... )`` in the ray chart; returns ``(v, [B_k])``."""
    poly = family.polynomials[name]
    v = poly.valuation(direction)
    layers = [sp.Integer(0) for _ in range(order)]
    for monom, coeff in poly.terms.items():
        k = sum(int(a) * int(d) for a, d in zip(monom, direction)) - v
        if k >= order:
            continue
        term = coeff.to_sympy()
        for i, a in enumerate(monom):
            if a:
                term *= chart_syms[i] ** int(a)
        layers[k] += term
    return v, layers


def _series_mul(a, b, order: int):
    out = []
    for k in range(order):
        acc = sp.Integer(0)
        for i in range(k + 1):
            if a[i] == 0 or b[k - i] == 0:
                continue
            acc += a[i] * b[k - i]
        out.append(acc)
    return out


def _binomial_series(z, f: int, order: int):
    """Truncated ``(1 + z1*t + z2*t^2 + ...)^f`` for integer ``f`` (any sign)."""
    out = [sp.Integer(1)] + [sp.Integer(0)] * (order - 1)
    zpow = [sp.Integer(1)] + [sp.Integer(0)] * (order - 1)
    zser = [sp.Integer(0), *z][:order]
    while len(zser) < order:
        zser.append(sp.Integer(0))
    for j in range(1, order):
        zpow = _series_mul(zpow, zser, order)
        if all(c == 0 for c in zpow):
            break
        b = sp.binomial(f, j)
        for k in range(order):
            if zpow[k] != 0:
                out[k] += b * zpow[k]
    return out


def _int_exponents(family: ParametricFamily, label: Label):
    e0, f0 = exponents_at_eps0(family, label)
    exps = []
    for v in (*e0, *f0):
        if not v.is_Integer:
            raise ValueError(
                "composite-master analysis needs integer exponents at epsilon=0; "
                f"got {v} for label {label!r}"
            )
        exps.append(int(v))
    n = family.nvars
    return exps[:n], exps[n:]


def _label_layer_cache(cache: dict, family, label: Label, direction: Direction, order: int):
    key = (label, direction)
    have_order, data = cache.get(key, (0, None))
    if have_order >= order:
        return data
    e0, f0 = _int_exponents(family, label)
    chart = _var_symbols(family, direction)
    power = 0
    prefactor = sp.Integer(1)
    series = [sp.Integer(1)] + [sp.Integer(0)] * (order - 1)
    for i, d in enumerate(direction):
        if d == 0:
            prefactor *= chart[i] ** e0[i]
        else:
            power += d * e0[i]
            prefactor *= chart[i] ** e0[i]
    for j, name in enumerate(family.poly_names):
        if f0[j] == 0:
            continue
        v, layers = _poly_layers(family, name, direction, chart, order)
        power += f0[j] * v
        a0 = layers[0]
        if a0 == 0:  # pragma: no cover - valuation guarantees a0 != 0
            raise RuntimeError(f"empty leading layer for {name} along {direction}")
        prefactor *= a0 ** f0[j]
        z = [sp.cancel(c / a0) for c in layers[1:]]
        series = _series_mul(series, _binomial_series(z, f0[j], order), order)
    coeffs = tuple(sp.cancel(prefactor * c) for c in series)
    data = (power, coeffs)
    cache[key] = (order, data)
    return data


def leading_asymptotic_signature(
    family: ParametricFamily,
    label: Label,
    ray,
    order: int = 1,
    _cache: dict | None = None,
) -> AsymptoticSignature:
    """Exact ``t -> 0`` Laurent data of the epsilon=0 integrand along ``ray``.

    Returns the leading power, the leading coefficient as a rational function of
    the remaining (face) variables / ray coordinates / parameters and, for
    ``order > 1``, the next coefficients.  The reported ``score`` is asserted to
    agree with :func:`valuations.base_score` (cross-check against the LF gate).
    """
    direction = _direction_of(ray)
    cache = {} if _cache is None else _cache
    power, coeffs = _label_layer_cache(cache, family, label, direction, order)
    score = power + sum(direction)
    ref = base_score(family, label, direction)
    if sp.Integer(score) != sp.nsimplify(ref):
        raise RuntimeError(
            f"asymptotic score {score} disagrees with base_score {ref} "
            f"for {label!r} along {direction}"
        )
    return AsymptoticSignature(label, direction, power, score, coeffs[:order])


# --------------------------------------------------------------------------- #
# Phase B: candidate pool and cancellation kernel                             #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class CompositeCandidate:
    """One pool member: a label plus a human-readable provenance tag."""

    label: Label
    origin: str


def build_candidate_pool(
    family: ParametricFamily,
    base_labels,
    *,
    var_shift_axes=(),
    poly_shift_axes=(),
    shift_depths=(-1, -2),
    numerator_vars=(),
    numerator_degree: int = 2,
):
    """Deterministic candidate pool around ``base_labels``.

    For every base label: the label itself, single-axis shifts of the listed
    variable/polynomial exponents by ``shift_depths``, each optionally multiplied
    by a numerator monomial in ``numerator_vars`` up to total degree
    ``numerator_degree``.  Duplicates keep the first (most primitive) origin.
    """
    nvars = family.nvars
    var_idx = {name: i for i, name in enumerate(family.variables)}
    poly_idx = {name: j for j, name in enumerate(family.poly_names)}
    shifts: list[tuple[tuple[int, ...], str]] = [(tuple([0] * len(base_labels[0])), "")]
    for name in var_shift_axes:
        i = var_idx[name]
        for delta in shift_depths:
            vec = [0] * len(base_labels[0])
            vec[i] = delta
            shifts.append((tuple(vec), f"+n_{name}{delta:+d}"))
    for name in poly_shift_axes:
        j = poly_idx[name]
        for delta in shift_depths:
            vec = [0] * len(base_labels[0])
            vec[nvars + j] = delta
            shifts.append((tuple(vec), f"+m_{name}{delta:+d}"))
    monomials: list[tuple[tuple[int, ...], str]] = []
    if numerator_vars:
        idxs = [var_idx[name] for name in numerator_vars]

        def emit(pos, left, acc):
            if pos == len(idxs):
                vec = [0] * len(base_labels[0])
                tag = []
                for i, a in zip(idxs, acc):
                    vec[i] = a
                    if a:
                        tag.append(f"*{family.variables[i]}" + (f"^{a}" if a > 1 else ""))
                monomials.append((tuple(vec), "".join(tag)))
                return
            for a in range(left + 1):
                emit(pos + 1, left - a, [*acc, a])

        emit(0, numerator_degree, [])
        monomials.sort(key=lambda mv: (sum(mv[0]), mv[0]))
    else:
        monomials = [(tuple([0] * len(base_labels[0])), "")]
    pool: list[CompositeCandidate] = []
    seen: set[Label] = set()
    for k, base in enumerate(base_labels):
        for svec, stag in shifts:
            for mvec, mtag in monomials:
                label = tuple(b + s + m for b, s, m in zip(base, svec, mvec))
                if label in seen:
                    continue
                seen.add(label)
                pool.append(CompositeCandidate(label, f"nf[{k}]{stag}{mtag}"))
    return tuple(pool)


def _level_rows(family, participants, scores, direction: Direction, cache: dict):
    """Linear equations (rows over the participants) forcing every Laurent level
    ``<= 0`` of ``sum c_i R_i`` along ``direction`` to vanish identically."""
    levels: dict[int, list[tuple[int, sp.Expr]]] = {}
    for i, (cand, s) in enumerate(zip(participants, scores)):
        if s > 0:
            continue
        order = -s + 1
        _, coeffs = _label_layer_cache(cache, family, cand.label, direction, order)
        for j in range(order):
            sigma = s + j
            if sigma > 0:
                continue
            if coeffs[j] != 0:
                levels.setdefault(sigma, []).append((i, coeffs[j]))
    rows: list[list[sp.Expr]] = []
    for sigma in sorted(levels):
        entries = levels[sigma]
        dens = [sp.denom(f) for _, f in entries]
        big = dens[0]
        for d in dens[1:]:
            big = sp.lcm(big, d)
        numers = [(i, sp.expand(sp.cancel(f * big))) for i, f in entries]
        gens = sorted(
            {s for _, n in numers for s in n.free_symbols},
            key=lambda s: s.name,
        )
        if not gens:
            row = [sp.Integer(0)] * len(participants)
            for i, n in numers:
                row[i] = n
            rows.append(row)
            continue
        table: dict[tuple, list[sp.Expr]] = {}
        for i, n in numers:
            poly = sp.Poly(n, *gens, domain="EX")
            for monom, coeff in poly.terms():
                row = table.setdefault(monom, [sp.Integer(0)] * len(participants))
                row[i] += coeff
        for monom in sorted(table):
            rows.append([sp.cancel(c) for c in table[monom]])
    return rows


def _rref(rows, ncols: int):
    mat = [list(r) for r in rows]
    pivots: list[int] = []
    prow = 0
    for col in range(ncols):
        piv = None
        for i in range(prow, len(mat)):
            if mat[i][col] != 0:
                piv = i
                break
        if piv is None:
            continue
        mat[prow], mat[piv] = mat[piv], mat[prow]
        pv = mat[prow][col]
        mat[prow] = [sp.cancel(e / pv) for e in mat[prow]]
        for i in range(len(mat)):
            if i != prow and mat[i][col] != 0:
                f = mat[i][col]
                mat[i] = [sp.cancel(a - f * b) for a, b in zip(mat[i], mat[prow])]
        pivots.append(col)
        prow += 1
        if prow == len(mat):
            break
    return mat[:prow], pivots


def _kernel_basis(rows, ncols: int):
    rref_rows, pivots = _rref(rows, ncols)
    free = [c for c in range(ncols) if c not in pivots]
    basis = []
    for fc in free:
        vec = [sp.Integer(0)] * ncols
        vec[fc] = sp.Integer(1)
        for i, pc in enumerate(pivots):
            vec[pc] = sp.cancel(-rref_rows[i][fc])
        basis.append(_normalize_vector(vec))
    return basis, len(pivots)


def _normalize_vector(vec):
    dens = [sp.denom(v) for v in vec if v != 0]
    if not dens:
        return tuple(vec)
    big = dens[0]
    for d in dens[1:]:
        big = sp.lcm(big, d)
    scaled = [sp.expand(sp.cancel(v * big)) for v in vec]
    content = None
    for v in scaled:
        if v == 0:
            continue
        c = sp.factor_list(v)[0] if v.free_symbols else v
        content = c if content is None else sp.gcd(content, c)
    if content not in (None, 0):
        scaled = [sp.cancel(v / content) for v in scaled]
    for v in scaled:
        if v != 0:
            if (
                sp.LC(sp.Poly(v, *sorted(v.free_symbols, key=lambda s: s.name))) < 0
                if v.free_symbols
                else v < 0
            ):
                scaled = [sp.expand(-u) for u in scaled]
            break
    return tuple(scaled)


# --------------------------------------------------------------------------- #
# Phase B/C: feasibility over all rays                                        #
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class CompositeVector:
    """One composite candidate ``M = sum_i c_i * J(participants[i].label)``."""

    coefficients: tuple[sp.Expr, ...]
    fully_lf: bool
    failing_rays: tuple[Direction, ...]
    unknown_rays: tuple[Direction, ...]

    def nonzero(self, participants):
        return tuple((participants[i].label, c) for i, c in enumerate(self.coefficients) if c != 0)


@dataclass(frozen=True)
class CompositeFeasibilityResult:
    status: str
    primary_direction: Direction
    pool_size: int
    participants: tuple[CompositeCandidate, ...]
    participant_scores: tuple[int, ...]
    kernel_dimension: int
    primary_basis: tuple[CompositeVector, ...]
    full_dimension: int
    full_basis: tuple[CompositeVector, ...]
    witness_rays: tuple[Direction, ...]
    checked_rays: int
    notes: tuple[str, ...]


def _check_vector(family, participants, vec, directions, cache):
    failing: list[Direction] = []
    unknown: list[Direction] = []
    for direction in directions:
        try:
            scores = [
                int(sp.nsimplify(base_score(family, p.label, direction)))
                for p, c in zip(participants, vec)
                if c != 0
            ]
        except (TypeError, ValueError):
            unknown.append(direction)
            continue
        active = [(p, c) for p, c in zip(participants, vec) if c != 0]
        if not active or min(scores) > 0:
            continue
        levels: dict[int, sp.Expr] = {}
        for (p, c), s in zip(active, scores):
            if s > 0:
                continue
            order = -s + 1
            _, coeffs = _label_layer_cache(cache, family, p.label, direction, order)
            for j in range(order):
                sigma = s + j
                if sigma > 0:
                    continue
                levels[sigma] = levels.get(sigma, sp.Integer(0)) + c * coeffs[j]
        if any(sp.cancel(v) != 0 for v in levels.values()):
            failing.append(direction)
    return tuple(failing), tuple(unknown)


def composite_master_feasibility(
    family: ParametricFamily,
    pool,
    primary_ray,
    *,
    random_trials: int = 64,
    seed: int = 20260706,
) -> CompositeFeasibilityResult:
    """Solve for composite masters whose bad Laurent layers cancel.

    Solves the identical-vanishing conditions on ``primary_ray`` for the pool
    members that are non-locally-finite there ("participants"), then verifies
    every kernel vector on all candidate rays plus the deterministic random
    safety net, and finally refines the kernel to the subspace that is locally
    finite on *every* checked ray.  Coefficients live in the field of rational
    functions of the non-regulator parameters; a fixed-sample rank cross-check
    guards against special-locus artifacts (``BadSpecialization``).
    """
    primary = _direction_of(primary_ray)
    cache: dict = {}
    notes: list[str] = []

    scored = []
    for cand in pool:
        s = base_score(family, cand.label, primary)
        try:
            s_int = int(sp.nsimplify(s))
        except (TypeError, ValueError):
            return CompositeFeasibilityResult(
                STATUS_BAD_SPECIALIZATION,
                primary,
                len(pool),
                (),
                (),
                0,
                (),
                0,
                (),
                (),
                0,
                (f"symbolic score {s} on primary ray for {cand.label!r}",),
            )
        scored.append((cand, s_int))
    participants = tuple(c for c, s in scored if s <= 0)
    scores = tuple(s for _, s in scored if s <= 0)
    if not participants:
        return CompositeFeasibilityResult(
            STATUS_NO_COMPOSITE,
            primary,
            len(pool),
            (),
            (),
            0,
            (),
            0,
            (),
            (),
            0,
            ("no pool member is non-locally-finite on the primary ray",),
        )

    rows = _level_rows(family, participants, scores, primary, cache)
    basis, rank = _kernel_basis(rows, len(participants))

    # -- BadSpecialization guard: fixed generic samples must reproduce the rank.
    params = [sp.Symbol(p) for p in family.parameters if p not in family.regulators]
    samples = [
        {p: sp.Rational(3 + 4 * k, 7 + 2 * k) for k, p in enumerate(params)},
        {p: sp.Rational(11 + 2 * k, 5 + 4 * k) for k, p in enumerate(params)},
    ]
    if params and rows:
        sample_ranks = []
        for sub in samples:
            m = sp.Matrix([[sp.cancel(e.subs(sub)) for e in row] for row in rows])
            sample_ranks.append(m.rank())
        if any(rk != rank for rk in sample_ranks):
            notes.append(f"rank mismatch: symbolic {rank}, samples {sample_ranks}")
            return CompositeFeasibilityResult(
                STATUS_BAD_SPECIALIZATION,
                primary,
                len(pool),
                participants,
                scores,
                len(basis),
                (),
                0,
                (),
                (),
                0,
                tuple(notes),
            )

    directions: list[Direction] = []
    seen = set()
    for ray in compute_candidate_rays(family):
        if ray.direction not in seen:
            seen.add(ray.direction)
            directions.append(ray.direction)
    for d in _random_directions(family.nvars, random_trials, seed):
        if d not in seen:
            seen.add(d)
            directions.append(d)
    other = [d for d in directions if d != primary]

    primary_vectors = []
    witness: list[Direction] = []
    for vec in basis:
        failing, unk = _check_vector(family, participants, vec, other, cache)
        for d in failing:
            if d not in witness:
                witness.append(d)
        primary_vectors.append(CompositeVector(tuple(vec), not failing and not unk, failing, unk))

    # -- refine the kernel to the subspace locally finite on every checked ray.
    full_vectors: list[CompositeVector] = []
    if basis:
        cons: list[list[sp.Expr]] = []
        for d in other:
            sub_scores = []
            ok = True
            for p in participants:
                try:
                    sub_scores.append(int(sp.nsimplify(base_score(family, p.label, d))))
                except (TypeError, ValueError):
                    ok = False
                    break
            if not ok:
                notes.append(f"symbolic score on ray {d}; excluded from refinement")
                continue
            if min(sub_scores) > 0:
                continue
            for row in _level_rows(family, participants, tuple(sub_scores), d, cache):
                cons.append(
                    [sp.cancel(sum(v[i] * row[i] for i in range(len(participants)))) for v in basis]
                )
        kernel2, _ = (
            _kernel_basis(cons, len(basis))
            if cons
            else (
                [
                    tuple(sp.Integer(1) if i == k else sp.Integer(0) for i in range(len(basis)))
                    for k in range(len(basis))
                ],
                0,
            )
        )
        for avec in kernel2:
            vec = _normalize_vector(
                [
                    sp.cancel(sum(avec[k] * basis[k][i] for k in range(len(basis))))
                    for i in range(len(participants))
                ]
            )
            failing, unk = _check_vector(family, participants, vec, other, cache)
            if failing:  # pragma: no cover - refinement should have removed these
                raise RuntimeError(f"refined vector still fails on {failing}")
            full_vectors.append(CompositeVector(vec, not unk, failing, unk))

    fully = [v for v in full_vectors if v.fully_lf]
    status = STATUS_FEASIBLE if fully else STATUS_NO_COMPOSITE
    if not basis:
        notes.append("primary-ray cancellation kernel is trivial")
    elif not fully:
        notes.append("primary-ray kernel exists but no member is LF on all rays")
    return CompositeFeasibilityResult(
        status,
        primary,
        len(pool),
        participants,
        scores,
        len(basis),
        tuple(primary_vectors),
        len(fully),
        tuple(fully),
        tuple(witness),
        len(directions),
        tuple(notes),
    )


def feasibility_to_payload(result: CompositeFeasibilityResult) -> dict:
    """JSON-ready summary (strings for all symbolic data; fully deterministic)."""

    def vec_payload(v: CompositeVector) -> dict:
        return {
            "coefficients": [str(c) for c in v.coefficients],
            "nonzero": {str(list(lab)): str(c) for lab, c in v.nonzero(result.participants)},
            "fully_lf": v.fully_lf,
            "failing_rays": [list(d) for d in v.failing_rays],
            "unknown_rays": [list(d) for d in v.unknown_rays],
        }

    return {
        "status": result.status,
        "primary_direction": list(result.primary_direction),
        "pool_size": result.pool_size,
        "participants": [
            {"label": list(c.label), "origin": c.origin, "score": s}
            for c, s in zip(result.participants, result.participant_scores)
        ],
        "kernel_dimension": result.kernel_dimension,
        "primary_basis": [vec_payload(v) for v in result.primary_basis],
        "full_dimension": result.full_dimension,
        "full_basis": [vec_payload(v) for v in result.full_basis],
        "witness_rays": [list(d) for d in result.witness_rays],
        "checked_rays": result.checked_rays,
        "notes": list(result.notes),
        "reducer_core": {"modified": False, "certificate_gates": "untouched"},
    }
