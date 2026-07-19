"""Boundary rays, scaling scores, and the strict local-finiteness decision (spec §4.4, §5.4, §8).

Local finiteness is a property of the integral **at ``epsilon = 0``**. Along a boundary ray
``rho`` (a one-parameter degeneration ``x_i ~ lambda^(rho_i)``, ``lambda -> 0``) the integrand
times the measure scales as ``lambda^(kappa_rho - 1)``, where

    base_score(label, rho) = sum_i rho_i * (e_i + 1) + sum_l f_l * val_rho(G_l)

with ``e_i = a_i + n_i``, ``f_l = b_l + m_l`` evaluated at ``epsilon = 0`` and
``val_rho(G_l) = min_{a in supp G_l} (a . rho)``. The integral converges at that boundary iff
``base_score > 0``. Local finiteness requires this for *every* relevant ray.

STRICT RULE: ``base_score == 0`` at ``epsilon = 0`` is NOT locally finite, even if the
``epsilon`` term would regulate the integral for one sign of ``epsilon``. We never look at the
sign of the regulator here.

Ray candidates (MVP): coordinate zero rays ``+e_i``, infinity rays ``-e_i``, and simple
Newton-support rays (support monomials of each ``G_l`` and per-polynomial variable diagonals).
An adaptive random-ray safety net catches missed divergences; anything the deterministic test
cannot settle (symbolic exponents, possible bulk/interior zeros without positivity assumptions)
yields ``"Unknown"`` — never ``True``.
"""

from __future__ import annotations

import random
import re
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from math import gcd

import sympy as sp

from .family import ParametricFamily
from .labels import Label

Direction = tuple[int, ...]


@dataclass(frozen=True)
class Ray:
    """A primitive integer boundary direction with a provenance tag."""

    direction: Direction
    kind: str  # "coord0" | "coordInf" | "newton" | "mixed" | "random"


def _primitive(vec) -> Direction | None:
    g = 0
    for x in vec:
        g = gcd(g, abs(int(x)))
    if g == 0:
        return None
    return tuple(int(x) // g for x in vec)


def compute_candidate_rays(family: ParametricFamily) -> list[Ray]:
    """Deterministic candidate boundary rays (coordinate + simple Newton-support)."""
    n = family.nvars
    rays: list[Ray] = []
    seen: set[Direction] = set()

    def add(direction, kind: str) -> None:
        p = _primitive(direction)
        if p is None or p in seen:
            return
        seen.add(p)
        rays.append(Ray(p, kind))

    for i in range(n):
        add(tuple(1 if k == i else 0 for k in range(n)), "coord0")
        add(tuple(-1 if k == i else 0 for k in range(n)), "coordInf")
    for name in family.poly_names:
        poly = family.polynomials[name]
        occurring = [0] * n
        for c in poly.support():
            add(c, "newton")
            add(tuple(-x for x in c), "newton")
            for k in range(n):
                if c[k] != 0:
                    occurring[k] = 1
        add(tuple(occurring), "mixed")
        add(tuple(-x for x in occurring), "mixed")
    return rays


def valuation_poly(poly, ray) -> int:
    """Tropical valuation ``min_{a in supp} (a . ray)`` (delegates to SparsePoly.valuation)."""
    direction = ray.direction if isinstance(ray, Ray) else tuple(ray)
    return poly.valuation(direction)


def exponents_at_eps0(family: ParametricFamily, label: Label):
    """Return ``(e_i, f_l)`` as SymPy expressions with all regulators set to zero."""
    e, f = family.exponent_at_label(label)
    subs = {sp.Symbol(r): 0 for r in family.regulators}
    e0 = [sp.simplify(pe.to_sympy().subs(subs)) for pe in e]
    f0 = [sp.simplify(pe.to_sympy().subs(subs)) for pe in f]
    for v in (*e0, *f0):
        if v.has(sp.zoo, sp.oo, -sp.oo, sp.nan):
            raise ValueError("exponent is singular at epsilon=0")
    return e0, f0


def score_from_exponents(e_syms, f_syms, family: ParametricFamily, direction) -> sp.Expr:
    """``sum_i d_i (e_i + 1) + sum_l f_l * val_direction(G_l)`` as a SymPy expression."""
    total = sp.Integer(0)
    for i, ev in enumerate(e_syms):
        total += int(direction[i]) * (ev + 1)
    for j, name in enumerate(family.poly_names):
        total += f_syms[j] * family.polynomials[name].valuation(direction)
    return sp.expand(total)


def base_score(family: ParametricFamily, label: Label, ray) -> sp.Expr:
    """Scaling score along ``ray`` at ``epsilon = 0`` (see module docstring)."""
    e0, f0 = exponents_at_eps0(family, label)
    direction = ray.direction if isinstance(ray, Ray) else tuple(ray)
    return score_from_exponents(e0, f0, family, direction)


def _positive_symbols(assumptions) -> set[str]:
    out: set[str] = set()
    for a in assumptions:
        m = re.match(r"^\s*([A-Za-z]\w*)\s*>\s*0\s*$", str(a))
        if m:
            out.add(m.group(1))
    return out


def _classify(score: sp.Expr, positive_symbols: set[str]) -> str:
    s = sp.simplify(score)
    if s.free_symbols:
        repl = {
            sym: sp.Symbol(str(sym), positive=True)
            for sym in s.free_symbols
            if str(sym) in positive_symbols
        }
        s = sp.simplify(s.subs(repl))
    if s.is_positive:
        return "pos"
    if s.is_nonpositive:  # includes exactly zero -> STRICT RULE: not locally finite
        return "nonpos"
    return "unknown"


@lru_cache(maxsize=None)  # Perf.2: deterministic in (nvars, trials, seed, bound) — build once
def _random_directions(nvars: int, trials: int, seed: int, bound: int = 3) -> tuple[Direction, ...]:
    rng = random.Random(seed)
    out: list[Direction] = []
    seen: set[Direction] = set()
    attempts = 0
    while len(out) < trials and attempts < trials * 20:
        attempts += 1
        p = _primitive(tuple(rng.randint(-bound, bound) for _ in range(nvars)))
        if p is None or p in seen:
            continue
        seen.add(p)
        out.append(p)
    return tuple(out)


# --- Perf.2: per-family caches -----------------------------------------------------------------
# Everything cached here is a deterministic function of the (effectively immutable) family, so a
# cache hit is bit-for-bit identical to recomputation. The cache is stashed on the family
# instance itself (families are plain attribute holders — never hashed, compared or serialized
# through ``__dict__`` by the reducer), so its lifetime matches the family's.


def _family_cache(family: ParametricFamily) -> dict:
    cache = family.__dict__.get("_valuations_cache")
    if cache is None:
        possyms = _positive_symbols(family.assumptions)
        cache = {
            "rays": tuple(compute_candidate_rays(family)),
            "positive_symbols": possyms,
            "pos_subs": {sp.Symbol(s): sp.Symbol(s, positive=True) for s in possyms},
            "dir_vals": {},  # Direction -> per-polynomial tropical valuations
            "poly_pos": None,  # poly name -> all coefficients provably positive (lazy)
            "lf_memo": {},  # (label, random_trials, seed) -> True | False | "Unknown"
        }
        object.__setattr__(family, "_valuations_cache", cache)
    return cache


def _poly_valuations(family: ParametricFamily, cache: dict, direction: Direction):
    """``tuple(val_direction(G_l))`` in ``poly_names`` order — label-independent, cached."""
    vals = cache["dir_vals"].get(direction)
    if vals is None:
        vals = tuple(
            family.polynomials[name].valuation(direction) for name in family.poly_names
        )
        cache["dir_vals"][direction] = vals
    return vals


def _as_fractions(values) -> list[Fraction] | None:
    """Exact ``Fraction`` copies of SymPy rationals; ``None`` if any value is not Rational."""
    out: list[Fraction] = []
    for v in values:
        if v.is_Rational:
            out.append(Fraction(v.p, v.q))
        else:
            return None
    return out


def _classify_at_direction(
    e0, f0, e_frac, f_frac, family: ParametricFamily, cache: dict,
    direction: Direction, positive_symbols: set[str],
) -> str:
    """Classify ``base_score`` along ``direction`` using cached polynomial valuations.

    When all exponents are rational (``e_frac``/``f_frac`` given), the score is an exact
    rational number and classification is its strict sign — ``> 0`` is ``"pos"``, anything
    else (including exactly 0, STRICT RULE) is ``"nonpos"``; identical to what
    ``_classify(score_from_exponents(...))`` returns on a Rational. Otherwise the original
    symbolic path is used, only with the per-direction valuations shared across labels.
    """
    vals = _poly_valuations(family, cache, direction)
    if e_frac is not None:
        total = Fraction(0)
        for d, ev in zip(direction, e_frac):
            if d:
                total += d * (ev + 1)
        for fv, val in zip(f_frac, vals):
            if fv and val:
                total += fv * val
        return "pos" if total > 0 else "nonpos"
    total = sp.Integer(0)
    for i, ev in enumerate(e0):
        total += int(direction[i]) * (ev + 1)
    for j, fv in enumerate(f0):
        total += fv * vals[j]
    return _classify(sp.expand(total), positive_symbols)


def _bulk_safe(family: ParametricFamily, f0, positive_symbols: set[str]) -> bool:
    """True if every polynomial appearing as a denominator has provably positive coefficients.

    A polynomial with all-positive coefficients is strictly positive on ``R_+^N`` (a sum of
    positive monomials), so it has no interior zero and the boundary-ray test is sufficient.
    Otherwise a bulk/interior singularity cannot be ruled out without stronger assumptions.
    """
    subsyms = {sp.Symbol(s): sp.Symbol(s, positive=True) for s in positive_symbols}
    for j, name in enumerate(family.poly_names):
        fj = f0[j]
        # Only denominators (possibly-negative exponent) can create a blow-up.
        is_denominator = bool(fj.is_negative) or (bool(fj.free_symbols) and not fj.is_nonnegative)
        if not is_denominator:
            continue
        for coeff in family.polynomials[name].terms.values():
            c = coeff.to_sympy().subs(subsyms)
            if not c.is_positive:
                return False
    return True


def _bulk_safe_cached(family: ParametricFamily, f0, cache: dict) -> bool:
    """Same decision as :func:`_bulk_safe`; per-polynomial positivity is label-independent."""
    poly_pos = cache["poly_pos"]
    if poly_pos is None:
        subsyms = cache["pos_subs"]
        poly_pos = {
            name: all(
                bool(coeff.to_sympy().subs(subsyms).is_positive)
                for coeff in family.polynomials[name].terms.values()
            )
            for name in family.poly_names
        }
        cache["poly_pos"] = poly_pos
    for j, name in enumerate(family.poly_names):
        fj = f0[j]
        # Only denominators (possibly-negative exponent) can create a blow-up.
        is_denominator = bool(fj.is_negative) or (bool(fj.free_symbols) and not fj.is_nonnegative)
        if is_denominator and not poly_pos[name]:
            return False
    return True


def is_locally_finite(
    family: ParametricFamily,
    label: Label,
    random_trials: int = 64,
    seed: int = 20260706,
):
    """Decide local finiteness at ``epsilon = 0``. Returns ``True``, ``False`` or ``"Unknown"``.

    ``True`` requires: every candidate ray (and random safety-net ray) has ``base_score > 0``,
    all exponents are numeric at ``epsilon = 0``, and every denominator polynomial is provably
    positive (no bulk singularity). Any confirmed ``base_score <= 0`` gives ``False``. Anything
    undecidable gives ``"Unknown"`` (never ``True``).

    Perf.2: verdicts are memoized per ``(label, random_trials, seed)`` on the family, and the
    label-independent pieces (candidate rays, random safety-net directions, per-direction
    polynomial valuations, denominator positivity) are computed once per family. The decision
    itself is unchanged and deterministic.
    """
    cache = _family_cache(family)
    key = (label, random_trials, seed)
    memo = cache["lf_memo"]
    if key not in memo:
        memo[key] = _is_locally_finite_impl(family, label, random_trials, seed, cache)
    return memo[key]


def _is_locally_finite_impl(
    family: ParametricFamily, label: Label, random_trials: int, seed: int, cache: dict
):
    try:
        e0, f0 = exponents_at_eps0(family, label)
    except ValueError:
        return "Unknown"
    positive_symbols = cache["positive_symbols"]
    numeric = all(not v.free_symbols for v in (*e0, *f0))
    # Rational fast path (exact): available whenever every exponent is a SymPy Rational.
    e_frac = _as_fractions(e0)
    f_frac = _as_fractions(f0) if e_frac is not None else None
    if f_frac is None:
        e_frac = None

    saw_unknown = False
    for ray in cache["rays"]:
        cls = _classify_at_direction(
            e0, f0, e_frac, f_frac, family, cache, ray.direction, positive_symbols
        )
        if cls == "nonpos":
            return False
        if cls == "unknown":
            saw_unknown = True

    if numeric and not saw_unknown:
        for direction in _random_directions(family.nvars, random_trials, seed):
            cls = _classify_at_direction(
                e0, f0, e_frac, f_frac, family, cache, direction, positive_symbols
            )
            if cls == "nonpos":
                return False
            if cls == "unknown":  # pragma: no cover - numeric scores are always decidable
                saw_unknown = True

    if saw_unknown:
        return "Unknown"
    if not _bulk_safe_cached(family, f0, cache):
        return "Unknown"
    return True


# --- Method.1 (External Int2): explainable local-finiteness audit ------------------------------
# Append-only diagnostics layer. Everything below REUSES the primitives above (candidate rays,
# ``_classify_at_direction``, ``score_from_exponents``, per-family caches) and never mutates
# reducer state: the per-family cache is only read/extended through the same memoized helpers.
# The verdict is DELEGATED to :func:`is_locally_finite`; the per-ray table is recomputed with
# the same primitives (no short-circuit) and cross-checked, so the report can never silently
# disagree with the decision procedure actually used by the reducer.


@dataclass(frozen=True)
class RayVerdict:
    """Sign of ``base_score`` along one boundary ray, with the score itself for the report."""

    ray: Ray
    score: object  # sp.Expr | None (None when exponents are singular at epsilon=0)
    classification: str  # "pos" | "nonpos" | "unknown"
    detail: str = ""


@dataclass(frozen=True)
class ShiftRecommendation:
    """Effect of one unit label shift on ``base_score`` along every failing ray."""

    shift: tuple[int, ...]  # length nvars+npolys, exactly one entry is +1 or -1
    deltas_on_failing: tuple[int, ...]  # base_score delta per failing ray (same order)
    improves_all: bool  # strictly positive delta on EVERY failing ray (and some ray fails)


@dataclass(frozen=True)
class LocalFinitenessReport:
    """Explainable local-finiteness verdict for one label (Method.1 directional audit)."""

    label: Label
    verdict: object  # True | False | "Unknown" — always equals is_locally_finite(...)
    rays: tuple[RayVerdict, ...]
    failing_rays: tuple[RayVerdict, ...]  # every "nonpos" ray (no short-circuit)
    unknown_rays: tuple[RayVerdict, ...]
    shift_deltas: tuple[ShiftRecommendation, ...]  # all 2*(nvars+npolys) unit shifts
    recommended_shifts: tuple[ShiftRecommendation, ...]  # improves_all, best total delta first
    bulk_safe: bool | None  # None when exponents are singular at epsilon=0
    notes: tuple[str, ...]


def _unit_shift_records(
    family: ParametricFamily, cache: dict, failing: tuple[RayVerdict, ...]
) -> tuple[tuple[ShiftRecommendation, ...], tuple[ShiftRecommendation, ...]]:
    """All unit-shift deltas over the failing rays + the strictly-improving subset.

    A shift of ``n_i`` by ``d`` changes ``base_score`` along ray ``rho`` by ``d * rho_i``; a
    shift of ``m_l`` by ``d`` changes it by ``d * val_rho(G_l)`` (both integers). Deterministic
    order: axes in label order with ``+1`` before ``-1``; recommendations sorted by total delta
    (descending) with the enumeration order as tie-break (stable sort).
    """
    nvars, npolys = family.nvars, len(family.poly_names)
    records: list[ShiftRecommendation] = []
    for k in range(nvars + npolys):
        for d in (1, -1):
            shift = tuple(d if i == k else 0 for i in range(nvars + npolys))
            deltas = []
            for rv in failing:
                direction = rv.ray.direction
                if k < nvars:
                    delta = d * direction[k]
                else:
                    delta = d * _poly_valuations(family, cache, direction)[k - nvars]
                deltas.append(int(delta))
            improves = bool(failing) and all(x > 0 for x in deltas)
            records.append(ShiftRecommendation(shift, tuple(deltas), improves))
    recommended = tuple(
        sorted(
            (s for s in records if s.improves_all),
            key=lambda s: -sum(s.deltas_on_failing),
        )
    )
    return tuple(records), recommended


def explain_local_finiteness(
    family: ParametricFamily,
    label: Label,
    random_trials: int = 64,
    seed: int = 20260706,
) -> LocalFinitenessReport:
    """Explain the strict local-finiteness verdict for ``label`` ray by ray.

    The verdict itself is delegated to :func:`is_locally_finite` (same memo, same decision);
    this function additionally reports the score and sign of every candidate ray, every failing
    random safety-net ray, and the unit label shifts that would strictly improve the score on
    all failing rays. A mismatch between the reconstructed verdict and the delegated one raises
    ``RuntimeError`` (it would indicate a bug, not a data condition).
    """
    verdict = is_locally_finite(family, label, random_trials, seed)
    cache = _family_cache(family)
    positive_symbols = cache["positive_symbols"]
    notes: list[str] = []

    try:
        e0, f0 = exponents_at_eps0(family, label)
    except ValueError:
        notes.append("exponents are singular at epsilon=0; every ray is undecidable")
        rays = tuple(
            RayVerdict(ray, None, "unknown", "exponent singular at epsilon=0")
            for ray in cache["rays"]
        )
        if verdict != "Unknown":  # pragma: no cover - impl returns "Unknown" on ValueError
            raise RuntimeError(
                f"explain_local_finiteness disagrees with is_locally_finite for {label!r}: "
                f"'Unknown' != {verdict!r}"
            )
        return LocalFinitenessReport(
            label, verdict, rays, (), rays, (), (), None, tuple(notes)
        )

    numeric = all(not v.free_symbols for v in (*e0, *f0))
    e_frac = _as_fractions(e0)
    f_frac = _as_fractions(f0) if e_frac is not None else None
    if f_frac is None:
        e_frac = None

    ray_verdicts: list[RayVerdict] = []
    saw_unknown = False
    for ray in cache["rays"]:
        score = sp.simplify(score_from_exponents(e0, f0, family, ray.direction))
        cls = _classify_at_direction(
            e0, f0, e_frac, f_frac, family, cache, ray.direction, positive_symbols
        )
        detail = ""
        if cls == "nonpos" and score == 0:
            detail = "score == 0 (STRICT RULE: not locally finite)"
        if cls == "unknown":
            saw_unknown = True
        ray_verdicts.append(RayVerdict(ray, score, cls, detail))

    # Random safety net under exactly the same gating as ``_is_locally_finite_impl``; only
    # failing random rays are added to the report (the net is a witness generator, not a table).
    if numeric and not saw_unknown:
        for direction in _random_directions(family.nvars, random_trials, seed):
            cls = _classify_at_direction(
                e0, f0, e_frac, f_frac, family, cache, direction, positive_symbols
            )
            if cls == "nonpos":
                score = sp.simplify(score_from_exponents(e0, f0, family, direction))
                ray_verdicts.append(
                    RayVerdict(Ray(direction, "random"), score, cls, "random safety-net ray")
                )

    failing = tuple(rv for rv in ray_verdicts if rv.classification == "nonpos")
    unknown = tuple(rv for rv in ray_verdicts if rv.classification == "unknown")
    bulk = _bulk_safe_cached(family, f0, cache)

    if failing:
        explain_verdict: object = False
    elif unknown:
        explain_verdict = "Unknown"
        notes.append("no failing ray, but at least one ray sign is undecidable")
    elif not bulk:
        explain_verdict = "Unknown"
        notes.append(
            "all rays positive, but a denominator polynomial is not provably positive "
            "(possible bulk singularity)"
        )
    else:
        explain_verdict = True
    if explain_verdict != verdict:
        raise RuntimeError(
            f"explain_local_finiteness disagrees with is_locally_finite for {label!r}: "
            f"{explain_verdict!r} != {verdict!r}"
        )
    if failing and all(rv.ray.kind == "random" for rv in failing):
        notes.append("failure witnessed only by the random safety net (no candidate ray fails)")

    shift_deltas, recommended = _unit_shift_records(family, cache, failing)
    return LocalFinitenessReport(
        label,
        verdict,
        tuple(ray_verdicts),
        failing,
        unknown,
        shift_deltas,
        recommended,
        bulk,
        tuple(notes),
    )


def report_to_payload(report: LocalFinitenessReport) -> dict:
    """JSON-safe dict for a report (deterministic; SymPy scores rendered as strings)."""

    def ray_payload(rv: RayVerdict) -> dict:
        return {
            "direction": list(rv.ray.direction),
            "kind": rv.ray.kind,
            "score": None if rv.score is None else str(rv.score),
            "classification": rv.classification,
            "detail": rv.detail,
        }

    def shift_payload(s: ShiftRecommendation) -> dict:
        return {
            "shift": list(s.shift),
            "deltas_on_failing": list(s.deltas_on_failing),
            "improves_all": s.improves_all,
        }

    return {
        "label": list(report.label),
        "verdict": report.verdict if isinstance(report.verdict, bool) else str(report.verdict),
        "rays": [ray_payload(rv) for rv in report.rays],
        "failing_rays": [ray_payload(rv) for rv in report.failing_rays],
        "unknown_rays": [ray_payload(rv) for rv in report.unknown_rays],
        "shift_deltas": [shift_payload(s) for s in report.shift_deltas],
        "recommended_shifts": [shift_payload(s) for s in report.recommended_shifts],
        "bulk_safe": report.bulk_safe,
        "notes": list(report.notes),
    }
