"""Opt-in adaptive search over the fixed-pass reducer (Pass Adaptive.1).

``reduce_family_adaptive`` runs a **deterministic schedule** of :class:`SearchLevel`
configurations, calling the existing :func:`reduce_family_once` for each level and stopping at
the first *certified* ``Success``. It is controlled schedule expansion — **not** a universal
decision procedure for reducibility: exhausting the schedule proves nothing beyond "these
particular configurations did not certify a reduction".

Honesty contract:

* the fixed-pass reducer is called as-is — no math here, no new ``Success`` path; every level
  goes through the strict gate in :mod:`result` (reconstruction + row-span certificate +
  AllLocallyFinite), so adaptive mode can never return ``Success`` that a single fixed pass
  with the same config would not;
* when no level succeeds, the returned result is the **best partial failure** under a
  deterministic order (see ``_partial_key``), with the full per-level history attached at
  ``result.diagnostics.extra["adaptive"]``;
* resource limits (``max_levels`` / ``max_labels`` / ``max_rows`` / ``timeout_sec``) produce an
  honest typed failure or an honest early stop — never a silent downgrade. **No limit is
  hard-preemptive**: ``max_labels`` skips an oversized level *before* it starts, ``max_rows``
  is observed only *after* the offending level has already run to completion, and
  ``timeout_sec`` is checked **between** levels only (levels are atomic — a long level always
  runs to completion). ``timeout_sec`` is the only wall-clock knob, is disabled by default,
  and never changes the mathematical content of any completed level.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field, replace
from fractions import Fraction
from typing import Sequence

from .family import ParametricFamily
from .finite_field import generate_primes
from .labels import Label, _normalize_ranges
from .reducer import (
    CERTIFICATE_FAILED,
    CERTIFICATE_INSUFFICIENT,
    CERTIFICATE_NOT_RUN,
    CERTIFICATE_PASSED,
    ReducerConfig,
    reduce_family_once,
)
from .result import (
    FAILURE_INTERPOLATION_FAILED,
    FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE,
    FAILURE_PARSER_NEEDS_EXPLICIT_FAMILY,
    FAILURE_RESOURCE_LIMIT_REACHED,
    FAILURE_TARGET_NOT_REDUCIBLE,
    FAILURE_VERIFICATION_FAILED,
    ReductionResult,
    build_reduction_result_from_reconstruction,
)

# Keep in sync with ``api._DEFAULT_LABEL_BOX`` (n-shifts fixed at 0, each m-shift in {-1, 0}).
_DEFAULT_LABEL_BOX = ((0, 0), (-1, 0))

# Stop reasons reported in ``diagnostics.extra["adaptive"]["stop_reason"]``.
STOP_SUCCESS = "success"
STOP_LEVELS_EXHAUSTED = "levels_exhausted"
STOP_RESOURCE_LIMIT = "resource_limit"
STOP_NON_RETRYABLE = "non_retryable_failure"

# Escalating the configuration cannot help these failure kinds, so the loop stops immediately.
_NON_RETRYABLE = frozenset({FAILURE_RESOURCE_LIMIT_REACHED, FAILURE_PARSER_NEEDS_EXPLICIT_FAMILY})

# Certificate quality for the best-partial order: a *failed* certificate is worse than one that
# never ran (failed = the coefficients are demonstrably wrong in an independent point).
_CERT_QUALITY = {
    CERTIFICATE_PASSED: 3,
    CERTIFICATE_INSUFFICIENT: 2,
    CERTIFICATE_NOT_RUN: 1,
    CERTIFICATE_FAILED: 0,
}


# --- public dataclasses -------------------------------------------------------------------------
@dataclass(frozen=True)
class SearchLevel:
    """One deterministic configuration step. ``None`` fields inherit the base config.

    ``tangent_degree_blocks=()`` means "explicitly no tangent rows at this level" (``None``
    inherits the base config's blocks). ``extra_samples`` / ``extra_primes`` deterministically
    extend the base sample set / prime list (never replace user-provided values).
    """

    name: str = ""
    label_box: tuple | None = None  # (n_range, m_range) for enumerate_box
    labels: tuple[Label, ...] | None = None  # explicit column set (wins over label_box)
    max_ibp_degree: int | None = None
    tangent_degree_blocks: tuple[tuple[int, int], ...] | None = None
    min_valid_records: int | None = None
    extra_samples: int = 0  # extend via the deterministic scattered generator (A30-safe)
    extra_primes: int = 0  # extend with the next primes below 2^31 (descending)
    rref_backend: str | None = None  # None -> inherit (package default stays "dict")


@dataclass(frozen=True)
class AdaptiveSearchConfig:
    """Adaptive-loop policy. ``levels=None`` uses :func:`default_search_levels`.

    ``max_labels`` is a pre-flight gate (an oversized level is *not* run); ``max_rows`` is
    checked after a level completes and stops further escalation; ``timeout_sec`` is checked
    between levels only. All limits are honest: hitting one is reported exactly, and if nothing
    ran at all the result is a typed ``ResourceLimitReached`` failure.
    """

    levels: tuple[SearchLevel, ...] | None = None
    max_levels: int = 3
    max_labels: int | None = None
    max_rows: int | None = None
    timeout_sec: float | None = None


@dataclass(frozen=True)
class AdaptiveLevelReport:
    """JSON-safe summary of one scheduled level (attempted or pre-flight-skipped)."""

    level: int
    name: str
    label_box: tuple | None
    n_explicit_labels: int | None
    max_ibp_degree: int
    tangent_degree_blocks: tuple | None
    n_samples: int
    n_primes: int
    rref_backend: str | None
    ran: bool
    skipped_reason: str | None
    status: str | None
    failure_reason: str | None
    error: str | None  # short deterministic failure detail; None on success / not-ran
    certificate_status: str | None
    n_labels: int
    n_rows: int
    n_reduced_records: int
    selected_rank: int | None
    n_non_lf_terms: int
    n_unknown_lf_terms: int
    reconstruction_verified: bool
    recommendation: str | None
    wall_seconds: float  # observability only — never feeds control flow


@dataclass
class AdaptiveSearchDiagnostics:
    """Full adaptive-run history (attached as a plain dict to the returned result)."""

    n_levels_planned: int
    n_levels_run: int = 0
    stop_reason: str = STOP_LEVELS_EXHAUSTED
    resource_limit: dict | None = None
    best_level: int | None = None
    reports: tuple[AdaptiveLevelReport, ...] = field(default_factory=tuple)

    def to_json_dict(self) -> dict:
        from dataclasses import asdict

        return {
            "enabled": True,
            "n_levels_planned": self.n_levels_planned,
            "n_levels_run": self.n_levels_run,
            "stop_reason": self.stop_reason,
            "resource_limit": self.resource_limit,
            "best_level": self.best_level,
            "levels": [asdict(r) for r in self.reports],
        }


# --- deterministic default MVP schedule ----------------------------------------------------------
def _normalized_box(box, family: ParametricFamily):
    n_pairs = tuple(_normalize_ranges(box[0], family.nvars, "n"))
    m_pairs = tuple(_normalize_ranges(box[1], family.npolys, "m"))
    return n_pairs, m_pairs


def _expand_box(box, family: ParametricFamily, delta: int, n_mask: tuple[int, ...] | None = None):
    """Deepen every m-range by ``delta`` on the low side; masked n-axes widen symmetrically.

    ``n_mask`` (one 0/1 flag per n-axis) selects n-axes to widen to ``(lo-delta, hi+delta)``.
    Labels are unconstrained integer shift tuples (:func:`labels.enumerate_box` only requires
    ``lo <= hi``), so widening both directions is structurally valid; the caller bounds the
    multiplicative blow-up via the ``max_labels`` guard in :func:`default_search_levels`.
    ``n_mask=None`` keeps every n-range unchanged (the MVP default).
    """
    n_pairs, m_pairs = _normalized_box(box, family)
    if n_mask is not None and delta > 0:
        n_pairs = tuple(
            (lo - delta, hi + delta) if flag else (lo, hi)
            for (lo, hi), flag in zip(n_pairs, n_mask)
        )
    return (n_pairs, tuple((lo - delta, hi) for lo, hi in m_pairs))


def _count_box_labels(box, family: ParametricFamily) -> int:
    n_pairs, m_pairs = _normalized_box(box, family)
    total = 1
    for lo, hi in (*n_pairs, *m_pairs):
        total *= hi - lo + 1
    return total


def default_search_levels(
    family: ParametricFamily,
    config: ReducerConfig,
    *,
    rref_backend: str | None = None,
    expand_n: Sequence[int] | None = None,
    max_labels: int | None = None,
) -> tuple[SearchLevel, ...]:
    """The deterministic 3-level MVP schedule, derived only from the base config.

    * level 0 (``base``): the base label box, ``max_ibp_degree=1``, no tangent rows;
    * level 1 (``expand-1``): every m-range deepened by one, ``max_ibp_degree=2``,
      tangent blocks ``((1, 1),)``;
    * level 2 (``deep``): m-ranges deepened by two, tangent ``((1, 1), (2, 2))``, four extra
      scattered samples and two extra primes; ``rref_backend`` applies to this level only
      (recommended ``"auto"`` when the ``speed`` extra is installed; default inherits the
      base config, i.e. the package default ``dict`` backend).

    ``expand_n`` (opt-in) is a per-n-axis 0/1 mask: masked axes widen symmetrically by the
    level delta (``(lo-k, hi+k)`` at level *k*), unmasked axes stay frozen. Because n-expansion
    multiplies the box volume, ``expand_n`` **requires** ``max_labels`` — a build-time guard:
    every planned level's label count must stay within it (``ValueError`` otherwise). This
    guard is distinct from the runtime pre-flight gate ``AdaptiveSearchConfig.max_labels``
    (which *skips* oversized levels instead of refusing to build the schedule).

    An explicit ``config.labels`` list cannot be grown deterministically — provide explicit
    ``SearchLevel``s instead (``ValueError``).
    """
    if config.labels is not None:
        raise ValueError(
            "default_search_levels cannot grow an explicit `labels` list; "
            "pass AdaptiveSearchConfig(levels=...) with your own SearchLevels"
        )
    mask: tuple[int, ...] | None = None
    if expand_n is not None:
        mask = tuple(1 if bool(x) else 0 for x in expand_n)
        if len(mask) != family.nvars:
            raise ValueError(
                f"expand_n mask has {len(mask)} entries, family has {family.nvars} n-axes"
            )
        if not any(mask):
            raise ValueError(
                "expand_n selects no n-axes; omit expand_n instead of passing all zeros"
            )
        if max_labels is None:
            raise ValueError(
                "expand_n widens the label box multiplicatively; "
                "pass max_labels=... to bound the planned schedule"
            )
    base_box = config.label_box if config.label_box is not None else _DEFAULT_LABEL_BOX
    levels = (
        SearchLevel(
            name="base",
            label_box=_normalized_box(base_box, family),
            max_ibp_degree=1,
            tangent_degree_blocks=(),
        ),
        SearchLevel(
            name="expand-1",
            label_box=_expand_box(base_box, family, 1, mask),
            max_ibp_degree=2,
            tangent_degree_blocks=((1, 1),),
        ),
        SearchLevel(
            name="deep",
            label_box=_expand_box(base_box, family, 2, mask),
            max_ibp_degree=2,
            tangent_degree_blocks=((1, 1), (2, 2)),
            extra_samples=4,
            extra_primes=2,
            rref_backend=rref_backend,
        ),
    )
    if max_labels is not None:
        for lvl in levels:
            planned = _count_box_labels(lvl.label_box, family)
            if planned > max_labels:
                raise ValueError(
                    f"default_search_levels: level {lvl.name!r} plans {planned} labels "
                    f"> max_labels {max_labels}; shrink the base box or the expand_n mask, "
                    f"or raise max_labels"
                )
    return levels


# --- per-level config construction ---------------------------------------------------------------
def _point_key(pt) -> tuple:
    return tuple(sorted((str(k), Fraction(v)) for k, v in pt.items()))


def _extended_samples(family: ParametricFamily, base_samples: Sequence, extra: int) -> list[dict]:
    """Base samples plus ``extra`` fresh deterministic scattered points (no duplicates)."""
    from .api import default_scattered_samples  # local import: api imports this module

    out = [dict(pt) for pt in base_samples]
    seen = {_point_key(pt) for pt in out}
    pool = default_scattered_samples(family.parameters, len(out) + extra + 16)
    for pt in pool:
        key = _point_key(pt)
        if key in seen:
            continue
        out.append(pt)
        seen.add(key)
        if len(out) >= len(base_samples) + extra:
            break
    return out


def _extended_primes(base_primes: Sequence[int], extra: int) -> list[int]:
    """Base primes plus the next ``extra`` distinct primes below 2^31 (descending)."""
    out = [int(p) for p in base_primes]
    have = set(out)
    for p in generate_primes(len(have) + extra + 8):
        if p in have:
            continue
        out.append(p)
        have.add(p)
        if len(out) >= len(base_primes) + extra:
            break
    return out


def _level_config(family: ParametricFamily, base: ReducerConfig, lvl: SearchLevel) -> ReducerConfig:
    kwargs: dict = {}
    if lvl.labels is not None:
        kwargs["labels"] = tuple(lvl.labels)
        kwargs["label_box"] = None
    elif lvl.label_box is not None:
        kwargs["label_box"] = lvl.label_box
        kwargs["labels"] = None
    if lvl.max_ibp_degree is not None:
        kwargs["max_ibp_degree"] = lvl.max_ibp_degree
    if lvl.tangent_degree_blocks is not None:
        kwargs["tangent_degree_blocks"] = tuple(lvl.tangent_degree_blocks) or None
    if lvl.min_valid_records is not None:
        kwargs["min_valid_records"] = lvl.min_valid_records
    if lvl.rref_backend is not None:
        kwargs["rref_backend"] = lvl.rref_backend
    if lvl.extra_samples > 0:
        kwargs["samples"] = _extended_samples(family, base.samples, lvl.extra_samples)
    if lvl.extra_primes > 0:
        kwargs["primes"] = _extended_primes(base.primes, lvl.extra_primes)
    return replace(base, **kwargs)


def _planned_n_labels(family: ParametricFamily, cfg: ReducerConfig) -> int:
    if cfg.labels is not None:
        return len(tuple(cfg.labels))
    box = cfg.label_box if cfg.label_box is not None else _DEFAULT_LABEL_BOX
    return _count_box_labels(box, family)


# --- failure recommendations + best-partial order ------------------------------------------------
def _recommendation(result: ReductionResult) -> str | None:
    if result.success:
        return None
    reason = result.status
    if reason == FAILURE_TARGET_NOT_REDUCIBLE:
        return "expand the label box, raise max_ibp_degree, or add tangent degree blocks"
    if reason == FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE:
        bad = sorted(
            list(result.diagnostics.non_lf_terms) + list(result.diagnostics.unknown_lf_terms)
        )
        where = f" around the remaining non-LF/Unknown labels {bad}" if bad else ""
        return f"expand the label box{where} — a different reduction path may avoid them"
    if reason == FAILURE_INTERPOLATION_FAILED:
        return "increase scattered samples and/or primes"
    if reason == FAILURE_VERIFICATION_FAILED:
        return (
            "add independent certificate points and re-run; never accept the current coefficients"
        )
    if reason == FAILURE_RESOURCE_LIMIT_REACHED:
        return "raise the reported resource limit or provide a smaller explicit schedule"
    if reason == FAILURE_PARSER_NEEDS_EXPLICIT_FAMILY:
        return "provide explicit Polynomials/MonomialExponents/PolynomialExponents"
    return None


def _partial_key(result: ReductionResult, level_index: int) -> tuple:
    """Deterministic best-partial order (higher tuple = better).

    Success > target reduced > fewer non-LF terms > fewer Unknown terms > reconstruction
    verified > certificate quality > earlier level (deterministic tie-break).
    """
    d = result.diagnostics
    cert_status = (d.extra.get("certificate") or {}).get("certificate_status")
    return (
        1 if result.success else 0,
        0 if result.status == FAILURE_TARGET_NOT_REDUCIBLE else 1,
        -len(d.non_lf_terms),
        -len(d.unknown_lf_terms),
        1 if d.reconstruction_verified else 0,
        _CERT_QUALITY.get(cert_status, 0),
        -level_index,
    )


# --- per-level report -----------------------------------------------------------------------------
def _error_text(result: ReductionResult | None) -> str | None:
    """Short deterministic failure detail for the per-level report (``None`` on success).

    Full failed ``ReductionResult``s are deliberately not retained (the adaptive payload stays
    small and JSON-safe), so this preserves the human-readable detail — the attempt's
    diagnostic messages, falling back to the typed failure reason. Messages are deterministic
    text (never timings), truncated to keep the payload bounded.
    """
    if result is None or result.success:
        return None
    msgs = [m for m in result.diagnostics.messages if m]
    text = "; ".join(msgs) if msgs else (result.failure_reason or result.status)
    return text[:500] if text else None


def _report_for(
    level_index: int,
    lvl: SearchLevel,
    cfg: ReducerConfig,
    *,
    ran: bool,
    skipped_reason: str | None = None,
    result: ReductionResult | None = None,
    wall: float = 0.0,
) -> AdaptiveLevelReport:
    extra = result.diagnostics.extra if result is not None else {}
    cert = extra.get("certificate") or {}
    selection = extra.get("record_selection") or {}
    return AdaptiveLevelReport(
        level=level_index,
        name=lvl.name or f"level-{level_index}",
        label_box=cfg.label_box if cfg.labels is None else None,
        n_explicit_labels=len(tuple(cfg.labels)) if cfg.labels is not None else None,
        max_ibp_degree=cfg.max_ibp_degree,
        tangent_degree_blocks=(
            tuple(cfg.tangent_degree_blocks) if cfg.tangent_degree_blocks else None
        ),
        n_samples=len(list(cfg.samples)),
        n_primes=len(list(cfg.primes)),
        rref_backend=cfg.rref_backend,
        ran=ran,
        skipped_reason=skipped_reason,
        status=result.status if result is not None else None,
        failure_reason=(
            result.failure_reason if result is not None and not result.success else None
        ),
        error=_error_text(result),
        certificate_status=cert.get("certificate_status"),
        n_labels=int(extra.get("n_labels", 0)),
        n_rows=int(extra.get("n_rows", 0)),
        n_reduced_records=int(extra.get("n_reduced_records", 0)),
        selected_rank=selection.get("selected_rank"),
        n_non_lf_terms=len(result.diagnostics.non_lf_terms) if result is not None else 0,
        n_unknown_lf_terms=(len(result.diagnostics.unknown_lf_terms) if result is not None else 0),
        reconstruction_verified=(
            bool(result.diagnostics.reconstruction_verified) if result is not None else False
        ),
        recommendation=_recommendation(result) if result is not None else None,
        wall_seconds=round(wall, 6),
    )


# --- main loop -------------------------------------------------------------------------------------
def reduce_family_adaptive(
    family: ParametricFamily,
    target_label: Label,
    config: ReducerConfig,
    search: AdaptiveSearchConfig | None = None,
) -> ReductionResult:
    """Run the deterministic adaptive schedule; return Success or the best partial failure.

    Each level is one ordinary :func:`reduce_family_once` pass (strict gates intact). The full
    per-level history lands in ``result.diagnostics.extra["adaptive"]``. See the module
    docstring for the honesty contract and resource-limit semantics.
    """
    search = search if search is not None else AdaptiveSearchConfig()
    if search.max_levels < 1:
        raise ValueError("max_levels must be >= 1")
    schedule = (
        tuple(search.levels) if search.levels is not None else default_search_levels(family, config)
    )
    if not schedule:
        raise ValueError("adaptive schedule is empty")
    schedule = schedule[: search.max_levels]

    reports: list[AdaptiveLevelReport] = []
    attempts: list[tuple[int, ReductionResult]] = []
    stop_reason = STOP_LEVELS_EXHAUSTED
    resource_limit: dict | None = None
    t0 = time.perf_counter()

    for i, lvl in enumerate(schedule):
        if i > 0 and search.timeout_sec is not None:
            elapsed = time.perf_counter() - t0
            if elapsed >= search.timeout_sec:
                stop_reason = STOP_RESOURCE_LIMIT
                resource_limit = {
                    "kind": "timeout_sec",
                    "limit": search.timeout_sec,
                    "observed_sec": round(elapsed, 6),
                    "level": i,
                }
                break

        cfg = _level_config(family, config, lvl)

        planned = _planned_n_labels(family, cfg)
        if search.max_labels is not None and planned > search.max_labels:
            stop_reason = STOP_RESOURCE_LIMIT
            resource_limit = {
                "kind": "max_labels",
                "limit": search.max_labels,
                "observed": planned,
                "level": i,
            }
            reports.append(
                _report_for(
                    i,
                    lvl,
                    cfg,
                    ran=False,
                    skipped_reason=(
                        f"planned n_labels {planned} exceeds max_labels {search.max_labels}"
                    ),
                )
            )
            break

        t_level = time.perf_counter()
        res = reduce_family_once(family, target_label, cfg)
        wall = time.perf_counter() - t_level
        attempts.append((i, res))
        reports.append(_report_for(i, lvl, cfg, ran=True, result=res, wall=wall))

        if res.success:
            stop_reason = STOP_SUCCESS
            break
        if res.status in _NON_RETRYABLE:
            stop_reason = STOP_NON_RETRYABLE
            break
        n_rows = int(res.diagnostics.extra.get("n_rows", 0))
        if search.max_rows is not None and n_rows > search.max_rows:
            stop_reason = STOP_RESOURCE_LIMIT
            resource_limit = {
                "kind": "max_rows",
                "limit": search.max_rows,
                "observed": n_rows,
                "level": i,
            }
            break

    if attempts:
        best_level, chosen = max(attempts, key=lambda pair: _partial_key(pair[1], pair[0]))
    else:
        # Nothing ran at all (pre-flight resource stop): honest typed failure, no coefficients.
        best_level = None
        detail = resource_limit or {"kind": "empty_schedule"}
        chosen = build_reduction_result_from_reconstruction(
            family,
            target_label,
            {},
            {},
            resource_limit_reached=True,
            messages=(f"adaptive search stopped before any level ran: {detail}",),
        )

    diagnostics = AdaptiveSearchDiagnostics(
        n_levels_planned=len(schedule),
        n_levels_run=sum(1 for r in reports if r.ran),
        stop_reason=stop_reason,
        resource_limit=resource_limit,
        best_level=best_level,
        reports=tuple(reports),
    )
    chosen.diagnostics.extra["adaptive"] = diagnostics.to_json_dict()
    chosen.diagnostics.messages = (
        *chosen.diagnostics.messages,
        (
            f"adaptive search: {diagnostics.n_levels_run}/{diagnostics.n_levels_planned} "
            f"levels run, stop_reason={stop_reason}, best_level={best_level}"
        ),
    )
    return chosen
