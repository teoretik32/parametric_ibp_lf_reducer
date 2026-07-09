"""Reducer orchestration (MVP, single fixed pass) — spec §6.

Wires the existing layers into one non-adaptive reduction attempt:

    enumerate labels -> local-finiteness flags -> generate rows (algebraic + coordinate-IBP,
    optionally tangent-IBP) -> collect modular normal-form records over ``primes x samples`` ->
    reconstruct coefficient functions -> the strict Success gate in :mod:`result`.

This module does NOT search or expand label boxes, does NOT tune anything, and does NOT stamp
``Success`` itself — the only ``Success`` comes from
:func:`result.build_reduction_result_from_reconstruction`. ``reduce_family_once`` runs the full
pipeline for a parsed family; ``reduce_rows_once`` takes ready-made rows so the orchestration can
be tested on tiny synthetic systems without the heavy row-generation layer.

Failure mapping (honest, conservative):
- no reduced records / target never a pivot -> ``TargetNotReducible``;
- reconstruction/interpolation did not validate -> ``InterpolationFailed``;
- a reconstructed term is not locally finite -> ``NormalFormNotLocallyFinite``;
- a reconstructed term has ``"Unknown"`` local finiteness -> failure (never ``Success``);
- bad specializations are skipped and counted, never patched;
- rank-deficient reduced records (Pass D4.3) are excluded from reconstruction by the max-rank
  selection in :func:`reconstruction.select_records_for_reconstruction` — skipped + counted in
  ``record_selection`` diagnostics, never zero-filled into the value table;
- reconstructed relations are row-span certified at independent off-sample points (Pass D4.5,
  :func:`certificate.verify_reduction_relation_mod_p`): with
  ``require_certificate_for_success=True`` (the default) ``Success`` additionally requires
  ``certificate_status == "Passed"`` — a failed/uninformative certificate maps to
  ``VerificationFailed`` (``FormalSuccess`` may still be ``True`` in diagnostics). This is the
  guard against sampling-degeneracy false positives (an on-lattice holdout is not independent;
  see ``notes/D4_STATUS.md`` §D4.4 and assumption A30).
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from fractions import Fraction

from .certificate import verify_reduction_relation_mod_p
from .family import ParametricFamily
from .labels import Label, enumerate_box
from .modular_normal_form import (
    STATUS_BAD_SPECIALIZATION,
    STATUS_EMPTY_SYSTEM,
    STATUS_REDUCED,
    STATUS_TARGET_NOT_REDUCIBLE,
)
from .reconstruction import (
    InterpolationFailed,
    reconstruct_coefficients,
    select_records_for_reconstruction,
)
from .records import collect_normal_form_records
from .result import ReductionResult, build_reduction_result_from_reconstruction
from .row_generation import (
    Row,
    generate_algebraic_rows,
    generate_coordinate_ibp_rows,
    generate_tangent_ibp_rows,
)
from .tangent_fields import generate_tangent_fields
from .valuations import is_locally_finite


@dataclass
class ReducerConfig:
    """Inputs for one fixed reduction pass (no adaptivity)."""

    primes: Sequence[int] = ()
    samples: Sequence[Mapping] = ()
    labels: Sequence[Label] | None = None  # explicit column set...
    label_box: tuple | None = None  # ...or an (n_range, m_range) box for ``enumerate_box``
    max_ibp_degree: int = 1
    tangent_degree_blocks: Sequence[tuple[int, int]] | None = None
    min_valid_records: int = 1
    preferred_masters: Sequence[Label] = ()
    eps_direction: str = "minus"
    # row-span certification of the reconstructed relation (Pass D4.5 / Verify.1):
    certificate_points: Sequence[Mapping] = ()  # explicit off-sample points (auto if empty)
    certificate_primes: Sequence[int] = ()  # primes for certification (main ``primes`` if empty)
    require_certificate_for_success: bool = True  # Success needs certificate_status == "Passed"
    min_certificate_points: int = 1  # informative (rank-generic) passing points required
    certificate_rank_policy: str = "selected_rank"  # only supported policy (explicit contract)


@dataclass
class ReducerRunDiagnostics:
    """Per-run counters produced by the orchestration (never decides success)."""

    n_labels: int = 0
    n_rows: int = 0
    n_records: int = 0
    n_reduced_records: int = 0
    n_selected_records: int = 0  # reduced records surviving the rank-consistency filter
    n_bad_specializations: int = 0
    n_target_not_pivot: int = 0
    n_empty_system: int = 0
    n_skipped_records: int = 0
    row_diagnostics: dict = field(default_factory=dict)
    record_selection: dict = field(default_factory=dict)  # rank filter diagnostics (Pass D4.3)
    reconstruction_diagnostics: dict = field(default_factory=dict)
    certificate: dict = field(default_factory=dict)  # row-span certificate diagnostics (D4.5)


# --- label enumeration + row generation ------------------------------------------------------
def _enumerate_labels(family: ParametricFamily, config: ReducerConfig) -> list[Label]:
    if config.labels is not None:
        return list(config.labels)
    if config.label_box is not None:
        n_range, m_range = config.label_box
        return list(enumerate_box(family.nvars, family.npolys, n_range, m_range))
    raise ValueError("ReducerConfig requires either 'labels' or 'label_box'")


def _generate_rows(
    family: ParametricFamily, seed_labels: Sequence[Label], config: ReducerConfig
) -> tuple[list[Row], dict]:
    rows: list[Row] = []
    by_kind: dict[str, int] = {}
    rejected: dict[str, int] = {}

    def _tally(res, kind: str) -> None:
        by_kind[kind] = len(res.rows)
        rows.extend(res.rows)
        for rr in res.rejected:
            rejected[rr.reason] = rejected.get(rr.reason, 0) + 1

    _tally(generate_algebraic_rows(family, seed_labels), "algebraic")
    _tally(
        generate_coordinate_ibp_rows(
            family, seed_labels, config.max_ibp_degree, eps_direction=config.eps_direction
        ),
        "coordinate_ibp",
    )
    if config.tangent_degree_blocks:
        fields = generate_tangent_fields(family, list(config.tangent_degree_blocks))
        _tally(
            generate_tangent_ibp_rows(
                family, seed_labels, fields, eps_direction=config.eps_direction
            ),
            "tangent_ibp",
        )
    return rows, {"by_kind": by_kind, "rejected": rejected}


# --- row-span certificate step (Pass D4.5) ----------------------------------------------------
CERTIFICATE_PASSED = "Passed"
CERTIFICATE_FAILED = "Failed"
CERTIFICATE_INSUFFICIENT = "Insufficient"  # no informative rank-generic points at all
CERTIFICATE_NOT_RUN = "NotRun"

_AUTO_CERT_DENOMS = (13, 17, 19)  # odd denominators unlikely to collide with user grids


def _default_certificate_points(samples: Sequence[Mapping], n: int = 3) -> list[dict]:
    """Deterministic off-sample certificate points: strictly beyond the per-parameter maximum of
    the reduction samples (so no point coincides with a sample), with varied denominators so the
    points do not inherit the samples' lattice structure."""
    samples = [dict(s) for s in samples]
    if not samples:
        return []
    params = sorted({k for s in samples for k in s})
    maxima = {p: max(Fraction(s[p]) for s in samples if p in s) for p in params}
    points = []
    for k in range(n):
        points.append(
            {
                p: maxima[p] + Fraction(2 * k + 1, _AUTO_CERT_DENOMS[(k + j) % len(_AUTO_CERT_DENOMS)])
                for j, p in enumerate(params)
            }
        )
    return points


def _run_certificate_step(
    family: ParametricFamily,
    rows: Sequence[Row],
    target_label: Label,
    coeffs: dict,
    points: Sequence[Mapping],
    primes: Sequence[int],
    selected_rank: int | None,
    min_points: int,
    lhs_terms: Mapping | None = None,
) -> dict:
    """Certify the reconstructed relation at the given points; never stamps ``Success``.

    ``lhs_terms`` (label -> coefficient) optionally generalizes the certified left-hand side to
    a linear combination (see :func:`certificate.verify_reduction_relation_mod_p`); ``None``
    keeps the classic single-target LHS.

    A point is *informative* only when its matrix rank equals ``selected_rank`` (the rank the
    reconstruction was built on): rank-DEFICIENT points solve a smaller system and are counted
    as filtered, and points where a row/coefficient denominator vanishes are counted as bad —
    neither passes nor fails. A point whose rank EXCEEDS ``selected_rank`` is a hard failure:
    a specialization's rank can never exceed the generic rank, so the reconstruction's
    ``selected_rank`` was not generic and its coefficients cannot be trusted. Overall
    ``certificate_status``: ``"Failed"`` if any informative point has a nonzero residual or any
    point exceeds ``selected_rank``; ``"Passed"`` if at least ``min_points`` informative points
    pass and none fail; otherwise ``"Insufficient"``.
    """
    n_passed = n_failed = n_rank_filtered = n_rank_exceeded = n_bad = 0
    histogram: dict[int, int] = {}
    first_residual: tuple | None = None
    primes = list(primes)
    for k, point in enumerate(points):
        prime = primes[k % len(primes)]
        cert = verify_reduction_relation_mod_p(
            family, rows, target_label, coeffs, dict(point), prime, lhs_terms=lhs_terms
        )
        if cert.status not in ("InSpan", "NotInSpan"):
            n_bad += 1
            continue
        histogram[cert.rank] = histogram.get(cert.rank, 0) + 1
        if selected_rank is not None and cert.rank > selected_rank:
            n_rank_exceeded += 1
        elif selected_rank is not None and cert.rank < selected_rank:
            n_rank_filtered += 1
        elif cert.in_span:
            n_passed += 1
        else:
            n_failed += 1
            if first_residual is None:
                first_residual = tuple(sorted(cert.residual.items()))
    if n_failed or n_rank_exceeded:
        status = CERTIFICATE_FAILED
    elif n_passed >= min_points:
        status = CERTIFICATE_PASSED
    else:
        status = CERTIFICATE_INSUFFICIENT
    return {
        "certificate_status": status,
        "n_certificate_points": len(list(points)),
        "n_certificate_points_passed": n_passed,
        "n_certificate_points_failed": n_failed,
        "n_certificate_rank_filtered": n_rank_filtered,
        "n_certificate_rank_exceeded": n_rank_exceeded,
        "n_certificate_bad_points": n_bad,
        "selected_rank": selected_rank,
        "certificate_rank_histogram": dict(sorted(histogram.items())),
        "certificate_points_used": [dict(p) for p in points],
        "first_nonzero_residual": first_residual,
    }


# --- shared core -----------------------------------------------------------------------------
def _reduce_core(
    family: ParametricFamily,
    target_label: Label,
    labels: Sequence[Label],
    rows: Sequence[Row],
    primes: Sequence[int],
    samples: Sequence[Mapping],
    *,
    lf_map: dict,
    preferred_masters: Sequence[Label],
    min_valid_records: int,
    row_diagnostics: dict,
    certificate_points: Sequence[Mapping] = (),
    certificate_primes: Sequence[int] = (),
    require_certificate_for_success: bool = True,
    min_certificate_points: int = 1,
    certificate_rank_policy: str = "selected_rank",
) -> ReductionResult:
    if certificate_rank_policy != "selected_rank":
        raise ValueError(
            f"unsupported certificate_rank_policy {certificate_rank_policy!r} "
            "(only 'selected_rank' is implemented)"
        )
    run = ReducerRunDiagnostics(
        n_labels=len(labels), n_rows=len(rows), row_diagnostics=dict(row_diagnostics)
    )

    records = collect_normal_form_records(
        family, rows, target_label, list(primes), list(samples),
        preferred_masters=tuple(preferred_masters), lf_map=lf_map,
    )
    run.n_records = len(records)
    for rec in records:
        if rec.status == STATUS_REDUCED and rec.formal_success:
            run.n_reduced_records += 1
        elif rec.status == STATUS_BAD_SPECIALIZATION:
            run.n_bad_specializations += 1
        elif rec.status == STATUS_TARGET_NOT_REDUCIBLE:
            run.n_target_not_pivot += 1
        elif rec.status == STATUS_EMPTY_SYSTEM:
            run.n_empty_system += 1
    run.n_skipped_records = run.n_records - run.n_reduced_records

    formal_success = run.n_reduced_records > 0
    target_reducible = formal_success  # the target became a pivot in at least one point

    # Rank-consistency selection (Pass D4.3): only max-rank reduced records feed reconstruction;
    # rank-deficient specializations are skipped + counted, never union-0-filled into the table.
    selected, selection = select_records_for_reconstruction(records)
    run.n_selected_records = len(selected)
    run.record_selection = selection
    messages: list[str] = []
    if selection["n_rank_filtered_records"]:
        messages.append(
            "rank filter: kept {kept}/{valid} reduced records at rank {rank} "
            "(histogram {hist})".format(
                kept=len(selected),
                valid=selection["n_valid_records_before_rank_filter"],
                rank=selection["selected_rank"],
                hist=selection["rank_histogram"],
            )
        )

    coeffs: dict | None = None
    interpolation_failed = False
    verified = False
    if target_reducible:
        if run.n_selected_records < min_valid_records:
            interpolation_failed = True  # too few rank-consistent points to reconstruct
            messages.append(
                f"only {run.n_selected_records} rank-consistent records "
                f"(< min_valid_records={min_valid_records})"
            )
        else:
            try:
                coeffs = reconstruct_coefficients(selected, family.parameters)
                verified = True  # reconstruction validated on independent holdout points
            except InterpolationFailed:
                coeffs = None
                interpolation_failed = True

    # Row-span certification at independent off-sample points (Pass D4.5). An on-grid holdout is
    # only as independent as the grid itself (A30), so Success additionally requires the
    # reconstructed relation to be certified in the row span away from the sample grid.
    verification_failed = False
    if coeffs is not None:
        points = list(certificate_points)
        if not points and require_certificate_for_success:
            points = _default_certificate_points(samples, n=max(3, min_certificate_points))
        if points:
            run.certificate = _run_certificate_step(
                family, rows, target_label, coeffs, points,
                list(certificate_primes) or list(primes),
                selected_rank=selection["selected_rank"], min_points=min_certificate_points,
            )
        else:
            run.certificate = {"certificate_status": CERTIFICATE_NOT_RUN}
        cert_status = run.certificate["certificate_status"]
        if cert_status != CERTIFICATE_NOT_RUN:
            messages.append(
                "certificate: {status} (passed {p}/{n}, rank-filtered {rf}, "
                "rank-exceeded {re}, bad {b})".format(
                    status=cert_status,
                    p=run.certificate["n_certificate_points_passed"],
                    n=run.certificate["n_certificate_points"],
                    rf=run.certificate["n_certificate_rank_filtered"],
                    re=run.certificate["n_certificate_rank_exceeded"],
                    b=run.certificate["n_certificate_bad_points"],
                )
            )
            if run.certificate.get("first_nonzero_residual"):
                residual_labels = [lab for lab, _ in run.certificate["first_nonzero_residual"]]
                messages.append(f"certificate residual columns: {residual_labels[:8]}")
        if require_certificate_for_success and cert_status != CERTIFICATE_PASSED:
            verification_failed = True
    else:
        run.certificate = {"certificate_status": CERTIFICATE_NOT_RUN}

    lf_flags: dict = {}
    if coeffs is not None:
        for lab in coeffs:
            lf_flags[lab] = lf_map[lab] if lab in lf_map else is_locally_finite(family, lab)

    run.reconstruction_diagnostics = {
        "verified": verified,
        "n_coefficients": 0 if coeffs is None else len(coeffs),
    }

    result = build_reduction_result_from_reconstruction(
        family,
        target_label,
        coeffs,
        lf_flags,
        reconstruction_verified=verified,
        independent_validation_passed=verified,
        formal_success=formal_success,
        interpolation_failed=interpolation_failed,
        target_reducible=target_reducible,
        verification_failed=verification_failed,
        n_records=run.n_records,
        n_skipped_records=run.n_skipped_records,
        messages=messages,
    )
    _attach_run_diagnostics(result, run)
    return result


def _attach_run_diagnostics(result: ReductionResult, run: ReducerRunDiagnostics) -> None:
    result.diagnostics.extra.update(
        {
            "n_labels": run.n_labels,
            "n_rows": run.n_rows,
            "n_records": run.n_records,
            "n_reduced_records": run.n_reduced_records,
            "n_selected_records": run.n_selected_records,
            "n_bad_specializations": run.n_bad_specializations,
            "n_target_not_pivot": run.n_target_not_pivot,
            "n_empty_system": run.n_empty_system,
            "n_skipped_records": run.n_skipped_records,
            "row_diagnostics": run.row_diagnostics,
            "record_selection": run.record_selection,
            "reconstruction_diagnostics": run.reconstruction_diagnostics,
            "certificate": run.certificate,
        }
    )
    result.diagnostics.extra["run"] = run


# --- public MVP API --------------------------------------------------------------------------
def reduce_family_once(
    family: ParametricFamily, target_label: Label, config: ReducerConfig
) -> ReductionResult:
    """Run one fixed reduction pass for a parsed family and return a typed result."""
    labels = _enumerate_labels(family, config)
    lf_map = {lab: is_locally_finite(family, lab) for lab in labels}
    rows, row_diag = _generate_rows(family, labels, config)
    return _reduce_core(
        family,
        target_label,
        labels,
        rows,
        config.primes,
        config.samples,
        lf_map=lf_map,
        preferred_masters=config.preferred_masters,
        min_valid_records=config.min_valid_records,
        row_diagnostics=row_diag,
        certificate_points=config.certificate_points,
        certificate_primes=config.certificate_primes,
        require_certificate_for_success=config.require_certificate_for_success,
        min_certificate_points=config.min_certificate_points,
        certificate_rank_policy=config.certificate_rank_policy,
    )


def reduce_rows_once(
    family: ParametricFamily,
    target_label: Label,
    labels: Iterable[Label],
    rows: Iterable[Row],
    primes: Sequence[int],
    samples: Sequence[Mapping],
    lf_flags: dict | None = None,
    preferred_masters: Iterable[Label] | None = None,
    min_valid_records: int = 1,
    certificate_points: Sequence[Mapping] = (),
    certificate_primes: Sequence[int] = (),
    require_certificate_for_success: bool = True,
    min_certificate_points: int = 1,
    certificate_rank_policy: str = "selected_rank",
) -> ReductionResult:
    """Orchestrate a reduction over ready-made ``rows`` (for testing without row generation).

    ``lf_flags`` supplies local-finiteness verdicts for the labels; any label missing there
    (e.g. a master produced by the RREF) is evaluated with :func:`is_locally_finite`.
    """
    labels = list(labels)
    rows = list(rows)
    lf_map = dict(lf_flags) if lf_flags is not None else {
        lab: is_locally_finite(family, lab) for lab in labels
    }
    return _reduce_core(
        family,
        target_label,
        labels,
        rows,
        list(primes),
        list(samples),
        lf_map=lf_map,
        preferred_masters=tuple(preferred_masters or ()),
        min_valid_records=min_valid_records,
        row_diagnostics={"provided_rows": len(rows)},
        certificate_points=certificate_points,
        certificate_primes=certificate_primes,
        require_certificate_for_success=require_certificate_for_success,
        min_certificate_points=min_certificate_points,
        certificate_rank_policy=certificate_rank_policy,
    )
