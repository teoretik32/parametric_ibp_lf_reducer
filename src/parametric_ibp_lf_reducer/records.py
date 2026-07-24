"""Multi-sample modular normal-form record collection (Pass 2G.1).

Orchestrates :func:`modular_normal_form` over a grid of ``(prime, param-sample)`` points for a
*fixed* family + row system + target label, producing a flat, deterministic list of
:class:`NormalFormRecord`. This is the honest, serialization-friendly unit that feeds coefficient
reconstruction (:mod:`reconstruction`).

Strict rules honoured here:
- every ``(prime, sample)`` point is recorded — bad specializations and non-reducible points keep
  their real status (``coeffs`` empty), never silently dropped;
- output order is deterministic (samples outer, primes inner, both taken in input order);
- coefficients are integers modulo ``prime`` only (no floats);
- **no reconstruction and no ``Success`` here** — this pass only *collects* records. Turning them
  into exact rational coefficients (and the multivariate case) is a later pass.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, MutableMapping, Sequence
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from fractions import Fraction

from .family import ParametricFamily
from .labels import Label
from .modular_normal_form import (
    STATUS_BAD_SPECIALIZATION,
    STATUS_EMPTY_SYSTEM,
    STATUS_REDUCED,
    STATUS_TARGET_NOT_REDUCIBLE,
    NormalFormResult,
    modular_normal_form,
    modular_normal_forms_multi,
)
from .ranking import RankedLabels, rank_labels
from .row_generation import Row
from .timing import StageTimings


@dataclass
class NormalFormRecord:
    """One collected ``(prime, sample)`` normal-form point.

    ``coeffs`` maps master label -> coefficient modulo ``prime`` (only populated for a
    ``Reduced`` record; the target equals ``sum coeffs[label] * J[label]`` at this point).
    A missing label means an *exact zero* at this point (see reconstruction's union support),
    not a dropped value.
    """

    prime: int
    sample: dict
    target_label: Label
    status: str
    formal_success: bool
    coeffs: dict = field(default_factory=dict)  # label -> coeff mod prime
    support: tuple = ()  # labels present in coeffs (sorted)
    all_terms_lf: object = None  # True | False | "Unknown" | None
    non_lf_terms: tuple = ()
    unknown_lf_terms: tuple = ()
    rank: int = 0
    diagnostics: dict = field(default_factory=dict)


def record_from_result(result: NormalFormResult) -> NormalFormRecord:
    """Adapt a single-sample :class:`NormalFormResult` into a :class:`NormalFormRecord`."""
    coeffs = dict(result.terms)
    return NormalFormRecord(
        prime=result.prime,
        sample=dict(result.sample),
        target_label=result.target_label,
        status=result.status,
        formal_success=result.formal_success,
        coeffs=coeffs,
        support=tuple(sorted(coeffs)),
        all_terms_lf=result.all_terms_lf,
        non_lf_terms=tuple(result.non_lf_terms),
        unknown_lf_terms=tuple(result.unknown_lf_terms),
        rank=result.rank,
        diagnostics={
            "nrows": result.nrows,
            "pivot_label": result.pivot_label,
        },
    )


def sample_prime_key(sample: Mapping, prime: int) -> str:
    """Canonical, JSON-friendly identity of one ``(sample, prime)`` point (Perf.12 record cache).

    Parameter values are canonicalized through :class:`~fractions.Fraction` and names are
    sorted, so logically equal samples map to the same key regardless of dict order or value
    representation (``2`` vs ``Fraction(2, 1)``). The key deliberately excludes the row system:
    a cache file is only valid for the exact system it was produced from, and the *caller*
    owning the cache is responsible for scoping it (one cache file per system + target).
    """
    parts = ",".join(
        f"{name}={Fraction(value)}"
        for name, value in sorted(sample.items(), key=lambda kv: str(kv[0]))
    )
    return f"p={prime};{parts}"


# --- Perf.3: process-parallel point collection ------------------------------------------------
# The per-run inputs (family/rows/ranking/...) are installed once per worker process by the
# executor initializer, so each ``(sample, prime)`` task ships only the point itself.
_POINT_CTX: dict = {}


def _init_point_worker(
    family, rows, target_label, preferred_masters, lf_map, ranking, rref_backend=None
) -> None:  # pragma: no cover - runs inside worker processes
    """Executor initializer: stash the shared per-run inputs in this worker process."""
    _POINT_CTX["ctx"] = (family, rows, target_label, preferred_masters, lf_map, ranking)
    _POINT_CTX["rref_backend"] = rref_backend


def _run_point(task: tuple) -> NormalFormRecord:  # pragma: no cover - runs inside workers
    """Compute one ``(sample, prime)`` record in a worker process (math identical to serial)."""
    sample, prime = task
    family, rows, target_label, preferred_masters, lf_map, ranking = _POINT_CTX["ctx"]
    result = modular_normal_form(
        family,
        rows,
        target_label,
        dict(sample),
        prime,
        preferred_masters=preferred_masters,
        lf_map=lf_map,
        ranking=ranking,
        rref_backend=_POINT_CTX.get("rref_backend"),
    )
    return record_from_result(result)


def collect_normal_form_records(
    family: ParametricFamily,
    rows: Iterable[Row],
    target_label: Label,
    primes: Sequence[int],
    samples: Sequence[Mapping],
    preferred_masters: Iterable[Label] = (),
    lf_map: dict | None = None,
    timings: StageTimings | None = None,
    ranking: RankedLabels | None = None,
    jobs: int = 1,
    rref_backend: str | None = None,
    record_cache: MutableMapping | None = None,
    on_record: Callable[[str, NormalFormRecord], None] | None = None,
) -> list[NormalFormRecord]:
    """Run :func:`modular_normal_form` over ``samples x primes`` and collect every point.

    Perf.12 (record cache): when ``record_cache`` is given (a mutable mapping keyed by
    :func:`sample_prime_key`), points already present are **not** recomputed and are returned
    from the cache; newly computed records are inserted into it. ``on_record(key, record)`` is
    invoked once per *newly computed* record, in deterministic task order, as results become
    available — a caller can persist incrementally and survive interruption. With
    ``record_cache=None`` and ``on_record=None`` the behavior is bit-identical to the
    historical path. The cache is trusted: scoping it to the exact system/target is the
    caller's responsibility.

    ``rref_backend`` (Perf.11) selects the RREF implementation for every point (serial and
    parallel paths alike); all backends return identical results by construction.

    Iterates samples in the outer loop and primes in the inner loop, both in the given order, so
    the returned list is deterministic. Bad/absent-target points are recorded honestly (their
    status is preserved), never skipped — reconstruction decides what to consume.

    Perf.1: the label ranking is built **once** here (timing key ``ranking_once``) from the
    union of all row labels and reused at every ``(prime, sample)`` point — ``rank_labels`` is
    no longer called per record. The elimination order is total per-label, so filtering the
    once-built superset order to the labels present at a point (done inside ``rref_mod_p``)
    is identical to ranking that point's labels directly; results are bit-for-bit unchanged.
    A caller may also pass a precomputed ``ranking`` covering all row labels.

    Perf.3: ``jobs`` (int >= 1) selects how many worker *processes* compute the independent
    ``(prime, sample)`` points. ``jobs=1`` (default) is the exact serial path. For ``jobs>1``
    the shared inputs are pickled once per worker (executor initializer) and results come back
    via ``ProcessPoolExecutor.map``, which preserves task order — the returned list is
    bit-identical to the serial one (same records, same order). Math is untouched. Caveat:
    per-point stage keys (``assemble_rows_mod_p``, ``rref_mod_p``, ``extract_normal_form``)
    accumulate inside the workers and are *not* merged back, so they read ``0.0`` in the
    caller's ``timings`` when ``jobs>1``; the caller-side ``records_total`` stage is unaffected.
    """
    if not isinstance(jobs, int) or isinstance(jobs, bool) or jobs < 1:
        raise ValueError(f"jobs must be an int >= 1, got {jobs!r}")
    rows = list(rows)
    primes = list(primes)
    samples = list(samples)
    preferred_masters = tuple(preferred_masters)
    t = timings if timings is not None else StageTimings()
    if ranking is None:
        with t.stage("ranking_once"):
            labels = sorted({c for row in rows for c in row.terms} | {target_label})
            ranking = rank_labels(
                family,
                labels,
                target=target_label,
                preferred_masters=preferred_masters,
                lf_map=lf_map,
            )
    tasks = [(sample, prime) for sample in samples for prime in primes]
    if record_cache is None:
        todo = list(tasks)
    else:
        todo = [
            (sample, prime)
            for sample, prime in tasks
            if sample_prime_key(sample, prime) not in record_cache
        ]

    def _note(sample, prime, record: NormalFormRecord) -> None:
        key = sample_prime_key(sample, prime)
        if record_cache is not None:
            record_cache[key] = record
        if on_record is not None:
            on_record(key, record)

    computed: list[NormalFormRecord] = []
    if jobs > 1 and len(todo) > 1:
        max_workers = min(jobs, len(todo))
        chunksize = max(1, len(todo) // (max_workers * 4))
        with ProcessPoolExecutor(
            max_workers=max_workers,
            initializer=_init_point_worker,
            initargs=(family, rows, target_label, preferred_masters, lf_map, ranking, rref_backend),
        ) as pool:
            for (sample, prime), record in zip(
                todo, pool.map(_run_point, todo, chunksize=chunksize)
            ):
                _note(sample, prime, record)
                computed.append(record)
    else:
        for sample, prime in todo:
            result = modular_normal_form(
                family,
                rows,
                target_label,
                dict(sample),
                prime,
                preferred_masters=preferred_masters,
                lf_map=lf_map,
                timings=timings,
                ranking=ranking,
                rref_backend=rref_backend,
            )
            record = record_from_result(result)
            _note(sample, prime, record)
            computed.append(record)
    if record_cache is None:
        return computed  # todo == tasks: identical to the historical path
    return [record_cache[sample_prime_key(sample, prime)] for sample, prime in tasks]


# --- Perf.5: multi-target collection over ONE shared RREF per point ---------------------------
def _run_point_multi(task: tuple) -> tuple:  # pragma: no cover - runs inside workers
    """Compute one ``(sample, prime)`` multi-target point in a worker (math identical to serial)."""
    sample, prime = task
    family, rows, target_labels, preferred_masters, lf_map, ranking = _POINT_CTX["ctx"]
    results = modular_normal_forms_multi(
        family,
        rows,
        target_labels,
        dict(sample),
        prime,
        preferred_masters=preferred_masters,
        lf_map=lf_map,
        ranking=ranking,
        rref_backend=_POINT_CTX.get("rref_backend"),
    )
    return tuple(record_from_result(results[tgt]) for tgt in target_labels)


def collect_normal_form_records_multi(
    family: ParametricFamily,
    rows: Iterable[Row],
    target_labels: Sequence[Label],
    primes: Sequence[int],
    samples: Sequence[Mapping],
    preferred_masters: Iterable[Label] = (),
    lf_map: dict | None = None,
    timings: StageTimings | None = None,
    ranking: RankedLabels | None = None,
    jobs: int = 1,
    rref_backend: str | None = None,
) -> dict[Label, list[NormalFormRecord]]:
    """Perf.5: collect records for *several* targets from ONE RREF per ``(prime, sample)`` point.

    ``rref_backend`` (Perf.11): see :func:`collect_normal_form_records` — backend selection
    only, identical results by construction.

    Same iteration order and per-target record semantics as :func:`collect_normal_form_records`
    (samples outer, primes inner; every point recorded honestly), but the per-point
    assemble + RREF work is shared across all targets via
    :func:`modular_normal_forms_multi`. Returns ``{target: [records...]}`` with every list in
    the same deterministic point order.

    The hoisted ranking (Perf.1) is built once with **all** targets in tier 0
    (``rank_labels(..., targets=target_labels)``); a caller may pass a precomputed ``ranking``
    built the same way. For a single target this is bit-identical to the single-target
    collector; for several targets the elimination order (hence masters) may differ from
    per-target runs — see the :func:`modular_normal_forms_multi` honesty note.

    ``jobs`` behaves exactly as in :func:`collect_normal_form_records` (Perf.3):
    ``jobs=1`` is the exact serial path; ``jobs>1`` computes points in worker processes with
    order-preserving ``ProcessPoolExecutor.map`` — bit-identical records, with the same caveat
    that per-point stage timings accumulate inside workers and read ``0.0`` in the caller's
    ``timings``.
    """
    if not isinstance(jobs, int) or isinstance(jobs, bool) or jobs < 1:
        raise ValueError(f"jobs must be an int >= 1, got {jobs!r}")
    targets = list(dict.fromkeys(target_labels))
    if not targets:
        raise ValueError("collect_normal_form_records_multi requires at least one target label")
    rows = list(rows)
    primes = list(primes)
    samples = list(samples)
    preferred_masters = tuple(preferred_masters)
    t = timings if timings is not None else StageTimings()
    if ranking is None:
        with t.stage("ranking_once"):
            labels = sorted({c for row in rows for c in row.terms} | set(targets))
            ranking = rank_labels(
                family,
                labels,
                targets=targets,
                preferred_masters=preferred_masters,
                lf_map=lf_map,
            )
    tasks = [(sample, prime) for sample in samples for prime in primes]
    out: dict[Label, list[NormalFormRecord]] = {tgt: [] for tgt in targets}
    if jobs > 1 and len(tasks) > 1:
        max_workers = min(jobs, len(tasks))
        chunksize = max(1, len(tasks) // (max_workers * 4))
        with ProcessPoolExecutor(
            max_workers=max_workers,
            initializer=_init_point_worker,
            initargs=(
                family,
                rows,
                tuple(targets),
                preferred_masters,
                lf_map,
                ranking,
                rref_backend,
            ),
        ) as pool:
            for point_records in pool.map(_run_point_multi, tasks, chunksize=chunksize):
                for tgt, rec in zip(targets, point_records):
                    out[tgt].append(rec)
        return out
    for sample, prime in tasks:
        results = modular_normal_forms_multi(
            family,
            rows,
            targets,
            dict(sample),
            prime,
            preferred_masters=preferred_masters,
            lf_map=lf_map,
            timings=timings,
            ranking=ranking,
            rref_backend=rref_backend,
        )
        for tgt in targets:
            out[tgt].append(record_from_result(results[tgt]))
    return out


def summarize_records(records: Iterable[NormalFormRecord]) -> dict:
    """Count records by status (diagnostics helper; does not decide success)."""
    counts = {
        STATUS_REDUCED: 0,
        STATUS_TARGET_NOT_REDUCIBLE: 0,
        STATUS_BAD_SPECIALIZATION: 0,
        STATUS_EMPTY_SYSTEM: 0,
    }
    total = 0
    for rec in records:
        total += 1
        counts[rec.status] = counts.get(rec.status, 0) + 1
    return {"total": total, "by_status": counts, "reduced": counts[STATUS_REDUCED]}
