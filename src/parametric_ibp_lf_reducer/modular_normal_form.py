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
                raise BadSpecialization(f"coefficient denominator vanishes mod {prime} at {sample}")
            v %= prime
            if v:
                specialized[label] = v
        if specialized:
            matrix.append(specialized)
    return matrix


def _extract_target_terms(
    family: ParametricFamily, ranked, pivot_row: dict, target_label: Label, prime: int
) -> tuple[dict, object, list, list]:
    """Read a target's normal form off its (fully reduced) RREF pivot row.

    Returns ``(terms, all_terms_lf, non_lf_terms, unknown_lf_terms)`` exactly as the strict
    single-target path always has; shared by the single- and multi-target extractions.
    """
    # target + sum v*col = 0  =>  target = sum (-v)*col
    terms = {c: (prime - v) % prime for c, v in sorted(pivot_row.items()) if c != target_label}

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
    return terms, all_lf, non_lf, unknown


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
    rref_backend: str | None = None,
) -> NormalFormResult:
    """Extract the target's normal form at one ``(prime, sample)`` point.

    ``rref_backend`` (Perf.11) selects the :func:`sparse_rref.rref_mod_p` implementation
    (``None`` means the default ``"dict"`` backend); every backend returns identical
    RREF results by construction, so this never affects the mathematical outcome.

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
    base = dict(
        status=STATUS_EMPTY_SYSTEM,
        target_label=target_label,
        prime=prime,
        sample=dict(sample),
        formal_success=False,
    )
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
        return NormalFormResult(**{**base, "status": STATUS_TARGET_NOT_REDUCIBLE, "nrows": nrows})

    if ranking is not None:
        ranked = ranking  # Perf.1: hoisted — built once per run, reused at every point
    else:
        with t.stage("ranking"):
            ranked = rank_labels(
                family,
                labels,
                target=target_label,
                preferred_masters=preferred_masters,
                lf_map=lf_map,
            )
    with t.stage("rref_mod_p"):
        res = rref_mod_p(matrix, prime, column_order=ranked.ordered, backend=rref_backend)

    if target_label not in res.pivots:
        return NormalFormResult(
            **{**base, "status": STATUS_TARGET_NOT_REDUCIBLE, "nrows": nrows, "rank": res.rank}
        )

    with t.stage("extract_normal_form"):
        pivot_row = res.pivots[target_label]  # {target: 1, free cols...}
        terms, all_lf, non_lf, unknown = _extract_target_terms(
            family, ranked, pivot_row, target_label, prime
        )

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


def modular_normal_forms_multi(
    family: ParametricFamily,
    rows: Iterable[Row],
    target_labels: Iterable[Label],
    sample: dict,
    prime: int,
    preferred_masters: Iterable[Label] = (),
    lf_map: dict | None = None,
    timings: StageTimings | None = None,
    ranking: RankedLabels | None = None,
    rref_backend: str | None = None,
) -> dict[Label, NormalFormResult]:
    """Perf.5: extract *several* targets' normal forms from ONE shared RREF at one point.

    ``rref_backend`` (Perf.11): see :func:`modular_normal_form` — backend selection only,
    identical results by construction.

    The expensive per-point work (``assemble_rows_mod_p`` + ``rref_mod_p``) does not depend on
    the target, so it is done once; every requested target is then read off the same reduced
    system. All targets sit in ranking tier 0 (see :func:`ranking.rank_labels` ``targets``), so
    each reducible target's pivot row is fully reduced against the other targets and its normal
    form contains free (master) columns only.

    Returns ``{target_label: NormalFormResult}`` in input target order. Per-target semantics are
    the strict single-target ones (same statuses, same LF diagnosis via
    :func:`_extract_target_terms`); a bad specialization rejects the whole point for every target
    (assembly is target-independent), never patched. For a single target this is bit-identical
    to :func:`modular_normal_form` (the tier-0 set is the same singleton).

    NB (honesty): for **two or more** targets the elimination order differs from any
    single-target run (the other targets are promoted to tier 0), so masters/coefficients are
    not guaranteed identical to per-target runs. Each result is still a genuine row-span
    relation and passes through the unchanged strict gates downstream.
    """
    targets = list(dict.fromkeys(target_labels))  # dedup, preserve order
    if not targets:
        raise ValueError("modular_normal_forms_multi requires at least one target label")
    t = timings if timings is not None else StageTimings()
    rows = list(rows)

    def _base(target_label: Label) -> dict:
        return dict(
            status=STATUS_EMPTY_SYSTEM,
            target_label=target_label,
            prime=prime,
            sample=dict(sample),
            formal_success=False,
        )

    if not rows:
        return {tgt: NormalFormResult(**_base(tgt)) for tgt in targets}

    try:
        with t.stage("assemble_rows_mod_p"):
            matrix = assemble_rows_mod_p(family, rows, sample, prime)
    except BadSpecialization:
        return {
            tgt: NormalFormResult(**{**_base(tgt), "status": STATUS_BAD_SPECIALIZATION})
            for tgt in targets
        }

    if not matrix:
        return {tgt: NormalFormResult(**_base(tgt)) for tgt in targets}

    labels = sorted({c for r in matrix for c in r})
    label_set = set(labels)
    nrows = len(matrix)

    if ranking is not None:
        ranked = ranking  # hoisted — built once per run (Perf.1), reused at every point
    else:
        with t.stage("ranking"):
            ranked = rank_labels(
                family,
                labels,
                targets=targets,
                preferred_masters=preferred_masters,
                lf_map=lf_map,
            )
    with t.stage("rref_mod_p"):
        res = rref_mod_p(matrix, prime, column_order=ranked.ordered, backend=rref_backend)

    out: dict[Label, NormalFormResult] = {}
    for tgt in targets:
        if tgt not in label_set:
            # mirrors the single-target pre-RREF check: nrows only, rank left at 0
            out[tgt] = NormalFormResult(
                **{**_base(tgt), "status": STATUS_TARGET_NOT_REDUCIBLE, "nrows": nrows}
            )
            continue
        if tgt not in res.pivots:
            out[tgt] = NormalFormResult(
                **{
                    **_base(tgt),
                    "status": STATUS_TARGET_NOT_REDUCIBLE,
                    "nrows": nrows,
                    "rank": res.rank,
                }
            )
            continue
        with t.stage("extract_normal_form"):
            terms, all_lf, non_lf, unknown = _extract_target_terms(
                family, ranked, res.pivots[tgt], tgt, prime
            )
        out[tgt] = NormalFormResult(
            status=STATUS_REDUCED,
            target_label=tgt,
            prime=prime,
            sample=dict(sample),
            formal_success=True,
            terms=terms,
            pivot_label=tgt,
            all_terms_lf=all_lf,
            non_lf_terms=non_lf,
            unknown_lf_terms=unknown,
            nrows=nrows,
            rank=res.rank,
        )
    return out
