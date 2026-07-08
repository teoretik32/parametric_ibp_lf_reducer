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

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field

from .family import ParametricFamily
from .labels import Label
from .modular_normal_form import (
    STATUS_BAD_SPECIALIZATION,
    STATUS_EMPTY_SYSTEM,
    STATUS_REDUCED,
    STATUS_TARGET_NOT_REDUCIBLE,
    NormalFormResult,
    modular_normal_form,
)
from .row_generation import Row


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


def collect_normal_form_records(
    family: ParametricFamily,
    rows: Iterable[Row],
    target_label: Label,
    primes: Sequence[int],
    samples: Sequence[Mapping],
    preferred_masters: Iterable[Label] = (),
    lf_map: dict | None = None,
) -> list[NormalFormRecord]:
    """Run :func:`modular_normal_form` over ``samples x primes`` and collect every point.

    Iterates samples in the outer loop and primes in the inner loop, both in the given order, so
    the returned list is deterministic. Bad/absent-target points are recorded honestly (their
    status is preserved), never skipped — reconstruction decides what to consume.
    """
    rows = list(rows)
    primes = list(primes)
    samples = list(samples)
    records: list[NormalFormRecord] = []
    for sample in samples:
        for prime in primes:
            result = modular_normal_form(
                family,
                rows,
                target_label,
                dict(sample),
                prime,
                preferred_masters=preferred_masters,
                lf_map=lf_map,
            )
            records.append(record_from_result(result))
    return records


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
