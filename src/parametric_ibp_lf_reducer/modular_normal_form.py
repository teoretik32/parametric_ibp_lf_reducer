"""Matrix assembly + single-sample modular normal-form extraction (spec §6, §7).

At one ``(prime, sample)`` point this specializes the parametric rows to integer rows modulo
``prime`` (via :meth:`ParamExpr.eval_mod_p`), runs the sparse RREF with the ranking's column
order, and reads off the target's normal form as a combination of the free (master) labels.

Strict rules honoured here:
- a bad specialization (a coefficient denominator vanishing modulo ``prime``) rejects the whole
  sample — it is never patched or skipped silently;
- coefficients are integers modulo ``prime`` only — no floating point;
- output is deterministic (ranking + RREF are deterministic, terms are label-sorted);
- ``formal_success`` (target reduced to free labels at this point) is reported separately from
  physical success; local finiteness of the resulting terms is only *diagnosed* here. No
  ``Success`` is produced.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field

from .family import ParametricFamily
from .labels import Label
from .ranking import RankedLabels, rank_labels
from .row_generation import Row
from .sparse_rref import rref_mod_p
from .timing import StageTimings
from .valuations import is_locally_finite

STATUS_REDUCED = "Reduced"
STATUS_TARGET_NOT_REDUCIBLE = "TargetNotReducible"
STATUS_BAD_SPECIALIZATION = "BadSpecialization"
STATUS_EMPTY_SYSTEM = "EmptySystem"


class BadSpecialization(Exception):
    """Raised when a row coefficient has a vanishing denominator at ``(sample, prime)``."""


@dataclass
class NormalFormResult:
    status: str
    target_label: Label
    prime: int
    sample: dict
    formal_success: bool
    terms: dict = field(default_factory=dict)  # label -> coeff mod prime (target = sum terms)
    pivot_label: Label | None = None
    all_terms_lf: object = None  # True | False | "Unknown" | None
    non_lf_terms: list = field(default_factory=list)
    unknown_lf_terms: list = field(default_factory=list)
    nrows: int = 0
    rank: int = 0


def assemble_rows_mod_p(
    family: ParametricFamily, rows: Iterable[Row], sample: dict, prime: int
) -> list[dict]:
    """Specialize parametric rows to integer sparse rows modulo ``prime``.

    Raises :class:`BadSpecialization` if any coefficient's denominator vanishes modulo ``prime``.
    Zero coefficients (vanishing numerators) are simply dropped.
    """
    matrix: list[dict] = []
    for row in rows:
        specialized: dict = {}
        for label, coeff in row.terms.items():
            v = coeff.eval_mod_p(sample, prime)
            if v is None:
                raise BadSpecialization(
                    f"coefficient denominator vanishes mod {prime} at {sample}"
                )
            v %= prime
            if v:
                specialized[label] = v
        if specialized:
            matrix.append(specialized)
    return matrix


def modular_normal_form(
    family: ParametricFamily,
    rows: Iterable[Row],
    target_label: Label,
    sample: dict,
    prime: int,
    preferred_masters: Iterable[Label] = (),
    lf_map: dict | None = None,
    timings: StageTimings | None = None,
    ranking: RankedLabels | None = None,
) -> NormalFormResult:
    """Extract the target's normal form at one ``(prime, sample)`` point.

    ``timings`` (Perf.0) optionally accumulates per-stage wall-clock seconds; it never
    affects the result.

    ``ranking`` (Perf.1) optionally supplies a precomputed :class:`RankedLabels` covering (a
    superset of) the labels appearing in ``rows``; when given, ``rank_labels`` is not called
    here. The elimination order is total per-label (``(tier, -complexity, label)``), so
    filtering a superset ranking to the labels present at this point (done by ``rref_mod_p``)
    yields exactly the same pivot order as ranking the subset directly — results are identical.
    """
    t = timings if timings is not None else StageTimings()
    rows = list(rows)
    base = dict(status=STATUS_EMPTY_SYSTEM, target_label=target_label, prime=prime,
                sample=dict(sample), formal_success=False)
    if not rows:
        return NormalFormResult(**base)

    try:
        with t.stage("assemble_rows_mod_p"):
            matrix = assemble_rows_mod_p(family, rows, sample, prime)
    except BadSpecialization:
        return NormalFormResult(**{**base, "status": STATUS_BAD_SPECIALIZATION})

    if not matrix:
        return NormalFormResult(**base)

    labels = sorted({c for r in matrix for c in r})
    nrows = len(matrix)
    if target_label not in set(labels):
        return NormalFormResult(
            **{**base, "status": STATUS_TARGET_NOT_REDUCIBLE, "nrows": nrows}
        )

    if ranking is not None:
        ranked = ranking  # Perf.1: hoisted — built once per run, reused at every point
    else:
        with t.stage("ranking"):
            ranked = rank_labels(
                family, labels, target=target_label, preferred_masters=preferred_masters,
                lf_map=lf_map,
            )
    with t.stage("rref_mod_p"):
        res = rref_mod_p(matrix, prime, column_order=ranked.ordered)

    if target_label not in res.pivots:
        return NormalFormResult(
            **{**base, "status": STATUS_TARGET_NOT_REDUCIBLE, "nrows": nrows, "rank": res.rank}
        )

    with t.stage("extract_normal_form"):
        pivot_row = res.pivots[target_label]  # {target: 1, free cols...}
        # target + sum v*col = 0  =>  target = sum (-v)*col
        terms = {
            c: (prime - v) % prime
            for c, v in sorted(pivot_row.items())
            if c != target_label
        }

        non_lf: list = []
        unknown: list = []
        for c in terms:
            verdict = ranked.lf.get(c, None)
            if verdict is None:
                verdict = is_locally_finite(family, c)
            if verdict is False:
                non_lf.append(c)
            elif verdict is not True:
                unknown.append(c)
        if non_lf:
            all_lf: object = False
        elif unknown:
            all_lf = "Unknown"
        else:
            all_lf = True

    return NormalFormResult(
        status=STATUS_REDUCED,
        target_label=target_label,
        prime=prime,
        sample=dict(sample),
        formal_success=True,
        terms=terms,
        pivot_label=target_label,
        all_terms_lf=all_lf,
        non_lf_terms=non_lf,
        unknown_lf_terms=unknown,
        nrows=nrows,
        rank=res.rank,
    )
