"""Finite-numerator locally-finite basis search (single-integrand semantics).

Direction correction after Method.3 (see ``docs/FINITE_NUMERATOR_BASIS_DESIGN.md``):
composite masters mixing different denominator sectors are sums of separately
divergent integrals and are NOT locally finite basis integrals — they were never
usable as HyperInt input. This module searches instead for INDIVIDUALLY locally
finite masters

    M = F_sector(x) * N(x),      N(x) = sum_{|alpha|_1 <= d} c_alpha * x^alpha,

where the COMPLETE numerator-decorated integrand must pass the strict LF gate on
every relevant ray by itself. No cancellation after integration is ever used; a
candidate is accepted only on a full ``is_locally_finite = True`` verdict for the
whole integrand, and it is exported as ONE integrand plus its exact defining
expansion into monomial labels.

Mathematical basis (proved in the design doc, machine-checked here):

* Lemma 1 (graded lowest layer): along any ray the lowest populated Laurent level
  of ``F_S * N`` has coefficient ``c_0 * N_init``, a product of nonzero elements
  of an integral domain. Under the STRICT RULE (every level <= 0 must vanish
  identically) the complete integrand is therefore LF on a ray iff EVERY monomial
  piece ``J(S + (alpha|0))`` is individually LF on that ray. The
  leading-cancellation solver below consequently must return exactly the span of
  the "clearing" monomials; the search *asserts* this (``lemma_consistent``)
  instead of assuming it.
* Lemma 2 (ray-sign obstruction): if every failing ray direction of a sector is
  componentwise <= 0, then ``min_alpha alpha . d <= 0`` for every polynomial
  numerator at every degree, so no polynomial numerator can cure the sector
  (``impossible_any_degree``).

Scope guard: purely diagnostic, read-only reuse of the same primitives as the LF
gate (``base_score`` / ``compute_candidate_rays`` / ``is_locally_finite``) and of
the exact-layer machinery of :mod:`composite_masters`. The reducer core, the
certificate and the LF gates are untouched.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import sympy as sp

from .coefficients import ParamExpr
from .composite_masters import _kernel_basis, _level_rows
from .family import IntegrandFactor, ParametricFamily
from .labels import Label
from .row_generation import monomials_up_to
from .sparse_poly import SparsePoly
from .valuations import (
    _classify,
    _family_cache,
    base_score,
    compute_candidate_rays,
    is_locally_finite,
)

Direction = tuple[int, ...]

STATUS_ALREADY_LF = "SectorAlreadyLF"
STATUS_FOUND = "FiniteNumeratorCandidatesFound"
STATUS_IMPOSSIBLE = "NumeratorCureImpossibleAnyDegree"
STATUS_NONE = "NoFiniteNumeratorWithinDegree"


@dataclass(frozen=True)
class _Member:
    """Minimal participant shim for :func:`composite_masters._level_rows`."""

    label: Label


# --------------------------------------------------------------------------- #
# Ansatz and labels


def numerator_ansatz(nvars: int, degree: int) -> list[tuple[int, ...]]:
    """All monomial exponents ``alpha >= 0`` with ``|alpha|_1 <= degree`` (incl. 1)."""
    return sorted(monomials_up_to(nvars, degree))


def decorated_label(sector: Label, alpha: tuple[int, ...]) -> Label:
    """Label of the monomial piece ``x^alpha * F_sector`` (n-part shifted by alpha)."""
    k = len(alpha)
    return tuple(sector[i] + alpha[i] for i in range(k)) + tuple(sector[k:])


# --------------------------------------------------------------------------- #
# Ray diagnostics


def sector_ray_table(family: ParametricFamily, sector: Label):
    """Per candidate ray: ``(direction, score, classification)`` for the bare sector."""
    cache = _family_cache(family)
    pos = cache["positive_symbols"]
    out = []
    for ray in compute_candidate_rays(family):
        s = base_score(family, sector, ray)
        out.append((tuple(ray.direction), s, _classify(s, pos)))
    return out


def failing_directions(family: ParametricFamily, sector: Label) -> list[Direction]:
    """Candidate-ray directions where the bare sector is not strictly positive."""
    return [d for d, _s, cls in sector_ray_table(family, sector) if cls != "pos"]


def impossible_any_degree(family: ParametricFamily, sector: Label) -> bool:
    """Lemma 2 flag: every failing ray componentwise <= 0 => no numerator (any
    degree) can cure the sector."""
    fails = failing_directions(family, sector)
    return bool(fails) and all(all(c <= 0 for c in d) for d in fails)


# --------------------------------------------------------------------------- #
# Honest leading-cancellation solver (Phase 2 step 5)


def leading_cancellation_kernel(
    family: ParametricFamily,
    sector: Label,
    alphas: list[tuple[int, ...]],
):
    """Impose vanishing of every non-integrable Laurent level of the COMPLETE
    integrand ``F_sector * sum c_alpha x^alpha`` on every candidate ray; solve the
    linear conditions for ``c_alpha``.

    Returns ``(kernel_basis, n_conditions)`` where each kernel vector is a tuple
    of exact coefficients aligned with ``alphas``. By Lemma 1 the kernel equals
    the span of the clearing monomials; callers cross-check this.
    """
    members = [_Member(decorated_label(sector, a)) for a in alphas]
    layer_cache: dict = {}
    rows: list[list[sp.Expr]] = []
    for ray in compute_candidate_rays(family):
        direction = tuple(ray.direction)
        scores: list[int] = []
        for m in members:
            s = sp.nsimplify(base_score(family, m.label, direction))
            if not s.is_Integer:
                raise ValueError(
                    f"non-integer ep=0 score {s} at {m.label} on ray {direction}; "
                    "not supported by the finite-numerator MVP"
                )
            scores.append(int(s))
        rows.extend(_level_rows(family, members, scores, direction, layer_cache))
    kernel, _rank = _kernel_basis(rows, len(members))
    return kernel, len(rows)


def score_clearing_alphas(
    family: ParametricFamily, sector: Label, alphas: list[tuple[int, ...]]
) -> list[tuple[int, ...]]:
    """Alphas whose monomial piece has strictly positive score on every candidate ray."""
    cache = _family_cache(family)
    pos = cache["positive_symbols"]
    out = []
    for a in alphas:
        lab = decorated_label(sector, a)
        ok = True
        for ray in compute_candidate_rays(family):
            if _classify(base_score(family, lab, ray), pos) != "pos":
                ok = False
                break
        if ok:
            out.append(a)
    return out


def gate_clearing_alphas(
    family: ParametricFamily,
    sector: Label,
    alphas: list[tuple[int, ...]],
    random_trials: int = 64,
    seed: int = 20260706,
) -> list[tuple[int, ...]]:
    """Alphas whose monomial piece passes the FULL LF gate (``True`` only)."""
    return [
        a
        for a in alphas
        if is_locally_finite(family, decorated_label(sector, a), random_trials, seed) is True
    ]


# --------------------------------------------------------------------------- #
# Full-integrand LF verdict


def full_integrand_lf(
    family: ParametricFamily,
    sector: Label,
    numerator: SparsePoly,
    random_trials: int = 64,
    seed: int = 20260706,
):
    """Strict-gate verdict for the COMPLETE integrand ``F_sector * N``.

    By Lemma 1 this equals the conjunction of the full gate over the monomial
    pieces; an independent min-valuation cross-check on the candidate rays guards
    the implementation (must never disagree for a ``True`` verdict).
    """
    if numerator.is_zero:
        raise ValueError("zero numerator has no LF verdict")
    verdict: object = True
    for alpha in numerator.support():
        v = is_locally_finite(family, decorated_label(sector, alpha), random_trials, seed)
        if v is False:
            verdict = False
            break
        if v != True:  # noqa: E712  -- verdict may be the string "Unknown"
            verdict = "Unknown"
    if verdict is True:
        # Independent cross-check: minimal score of the decorated integrand.
        cache = _family_cache(family)
        pos = cache["positive_symbols"]
        for ray in compute_candidate_rays(family):
            mins = min(
                sp.nsimplify(base_score(family, decorated_label(sector, a), ray))
                for a in numerator.support()
            )
            if _classify(mins, pos) != "pos":  # pragma: no cover - Lemma 1 guard
                raise RuntimeError(
                    f"Lemma 1 violation: gate True but min score {mins} on ray "
                    f"{tuple(ray.direction)} for sector {sector}"
                )
    return verdict


# --------------------------------------------------------------------------- #
# Exported object


@dataclass(frozen=True)
class FiniteNumeratorIntegral:
    """ONE numerator-decorated integrand ``F_sector(x) * N(x)`` with its verdict.

    ``defining_expansion`` is the exact identity
    ``M - sum_alpha c_alpha * J(sector + (alpha|0)) = 0`` used to embed the
    master into the ordinary IBP row system.
    """

    variables: tuple[str, ...]
    poly_names: tuple[str, ...]
    sector: Label
    numerator: SparsePoly
    lf_verdict: object

    @property
    def degree(self) -> int:
        return max(sum(a) for a in self.numerator.support())

    def defining_expansion(self) -> list[tuple[Label, str]]:
        return [
            (decorated_label(self.sector, alpha), str(sp.nsimplify(coeff.to_sympy())))
            for alpha, coeff in sorted(self.numerator.terms.items())
        ]

    def numerator_text(self) -> str:
        parts = []
        for alpha, coeff in sorted(self.numerator.terms.items()):
            factor = IntegrandFactor(
                self.variables, self.poly_names, alpha, (0,) * len(self.poly_names)
            ).to_wolfram_text()
            c = sp.nsimplify(coeff.to_sympy())
            if c == 1:
                parts.append(factor)
            elif factor == "1":
                parts.append(f"({c})")
            else:
                parts.append(f"({c})*{factor}")
        return " + ".join(parts)

    def to_wolfram_text(self) -> str:
        """The single HyperInt-ready integrand ``(N) * F_sector``."""
        k = len(self.variables)
        factor = IntegrandFactor(
            self.variables, self.poly_names, self.sector[:k], self.sector[k:]
        ).to_wolfram_text()
        return f"({self.numerator_text()})*({factor})"

    def payload(self) -> dict:
        return {
            "sector": list(self.sector),
            "numerator": self.numerator_text(),
            "degree": self.degree,
            "lf_verdict": str(self.lf_verdict),
            "integrand": self.to_wolfram_text(),
            "defining_expansion": [
                {"label": list(lab), "coeff": c} for lab, c in self.defining_expansion()
            ],
        }


# --------------------------------------------------------------------------- #
# Per-sector search and overall report


@dataclass
class SectorSearchReport:
    sector: Label
    degree: int
    bare_lf: object
    failing_rays: list[tuple[Direction, str]]
    ansatz_size: int
    kernel_dim: int
    n_conditions: int
    clearing: list[tuple[int, ...]]
    lemma_consistent: bool
    impossible_any_degree: bool
    status: str
    candidates: list[FiniteNumeratorIntegral] = field(default_factory=list)

    def payload(self) -> dict:
        return {
            "sector": list(self.sector),
            "degree": self.degree,
            "bare_lf": str(self.bare_lf),
            "failing_rays": [
                {"direction": list(d), "score": s} for d, s in self.failing_rays
            ],
            "ansatz_size": self.ansatz_size,
            "kernel_dim": self.kernel_dim,
            "n_conditions": self.n_conditions,
            "clearing_monomials": [list(a) for a in self.clearing],
            "lemma_consistent": self.lemma_consistent,
            "numerator_cure_impossible_any_degree": self.impossible_any_degree,
            "status": self.status,
            "candidates": [c.payload() for c in self.candidates],
        }


def _one_coeff(params) -> ParamExpr:
    return ParamExpr.from_int(1, tuple(params))


def search_sector(
    family: ParametricFamily,
    sector: Label,
    degree: int,
    random_trials: int = 64,
    seed: int = 20260706,
) -> SectorSearchReport:
    """Honest finite-numerator search for one denominator sector.

    Builds the full ansatz, solves the leading-cancellation conditions of the
    complete integrand, cross-checks the kernel against Lemma 1, and accepts
    candidates only on a full-integrand ``is_locally_finite = True`` verdict.
    """
    bare = is_locally_finite(family, sector, random_trials, seed)
    table = sector_ray_table(family, sector)
    fails = [(d, str(s)) for d, s, cls in table if cls != "pos"]
    alphas = numerator_ansatz(family.nvars, degree)
    kernel, n_conditions = leading_cancellation_kernel(family, sector, alphas)
    score_clear = score_clearing_alphas(family, sector, alphas)
    clear_idx = {alphas.index(a) for a in score_clear}
    lemma_consistent = len(kernel) == len(score_clear) and all(
        all(v == 0 for j, v in enumerate(vec) if j not in clear_idx) for vec in kernel
    )
    gate_clear = gate_clearing_alphas(family, sector, alphas, random_trials, seed)

    candidates: list[FiniteNumeratorIntegral] = []
    params = family.parameters
    for a in gate_clear:
        num = SparsePoly.monomial(family.nvars, a, _one_coeff(params))
        v = full_integrand_lf(family, sector, num, random_trials, seed)
        assert v is True  # gate-clearing monomials are LF by construction
        candidates.append(
            FiniteNumeratorIntegral(
                family.variables, family.poly_names, sector, num, v
            )
        )
    nontrivial = [a for a in gate_clear if sum(a) > 0]
    if len(nontrivial) > 1:
        # Representative genuinely polynomial candidate (all clearing monomials).
        num = SparsePoly(
            family.nvars,
            tuple(params),
            {a: _one_coeff(params) for a in nontrivial},
        )
        v = full_integrand_lf(family, sector, num, random_trials, seed)
        if v is True:
            candidates.append(
                FiniteNumeratorIntegral(
                    family.variables, family.poly_names, sector, num, v
                )
            )

    if bare is True:
        status = STATUS_ALREADY_LF
    elif any(c.degree > 0 for c in candidates):
        status = STATUS_FOUND
    elif impossible_any_degree(family, sector):
        status = STATUS_IMPOSSIBLE
    else:
        status = STATUS_NONE

    return SectorSearchReport(
        sector=sector,
        degree=degree,
        bare_lf=bare,
        failing_rays=fails,
        ansatz_size=len(alphas),
        kernel_dim=len(kernel),
        n_conditions=n_conditions,
        clearing=score_clear,
        lemma_consistent=lemma_consistent,
        impossible_any_degree=impossible_any_degree(family, sector),
        status=status,
        candidates=candidates,
    )


def finite_numerator_search(
    family: ParametricFamily,
    sectors: dict[str, Label],
    degrees: tuple[int, ...] = (1, 2),
    random_trials: int = 64,
    seed: int = 20260706,
) -> dict:
    """Scan sectors x degrees; return a JSON-ready payload with honest statuses."""
    reports = []
    for name, sector in sectors.items():
        for d in degrees:
            rep = search_sector(family, sector, d, random_trials, seed)
            entry = rep.payload()
            entry["sector_name"] = name
            reports.append(entry)
    statuses = sorted({r["status"] for r in reports})
    found = any(
        r["status"] == STATUS_FOUND and not r["numerator_cure_impossible_any_degree"]
        for r in reports
    )
    return {
        "method": "finite-numerator LF basis search (single-integrand semantics)",
        "degrees": list(degrees),
        "n_sectors": len(sectors),
        "reports": reports,
        "statuses": statuses,
        "new_lf_masters_found": found,
        "lemma_consistent_everywhere": all(r["lemma_consistent"] for r in reports),
    }
