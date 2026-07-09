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
