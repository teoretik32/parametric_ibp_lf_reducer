"""LF-constrained reduction feasibility modulo a prime (Method.1, External Int2).

Question answered here: does the row system admit ANY reduction of the target whose support
lies entirely inside the *allowed* (locally finite) labels — independent of the normal-form
ranking used by the reducer?

Projection argument. The rows span relations ``sum_l r(l) J[l] = 0``. A reduction
``J[t] = sum_{l allowed} c_l J[l]`` exists iff ``e_t - sum c_l e_l`` lies in the row span,
i.e. iff, after DELETING every allowed column, the projected target unit vector lies in the
span of the projected rows. Deleting columns is a linear projection, so the test is exact
(mod ``prime`` at one ``sample``): ``Feasible`` here is a certificate that the ranking-driven
normal form *could* have produced an LF-only answer at this point; ``Obstructed`` means no
linear combination of these rows eliminates the target through allowed labels only — a
statement about this row system and this label box, never a global impossibility claim.

Everything is append-only: this module reuses ``assemble_rows_mod_p``/``rref_mod_p`` and the
certificate's pivot-reduction helper, and never touches reducer state.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from .certificate import _reduce_vector_by_pivots
from .labels import Label
from .modular_normal_form import BadSpecialization, assemble_rows_mod_p
from .row_generation import Row
from .sparse_rref import rref_mod_p

STATUS_FEASIBLE = "Feasible"
STATUS_OBSTRUCTED = "Obstructed"
STATUS_BAD_SPECIALIZATION = "BadSpecialization"


@dataclass(frozen=True)
class LFFeasibilityResult:
    """Outcome of the LF-constrained span test at one ``(sample, prime)`` point."""

    status: str  # Feasible | Obstructed | BadSpecialization
    prime: int
    sample: tuple  # sorted (name, value) pairs — deterministic, JSON-safe
    rank: int  # rank of the projected system (0 when BadSpecialization)
    nrows: int  # nonzero specialized rows
    n_projected_rows: int  # rows surviving the deletion of allowed columns
    n_allowed: int  # allowed (LF-True) labels in the box, target excluded
    n_forbidden: int  # forbidden columns actually present (incl. out-of-box support)
    residual_support: tuple  # canonical residual of e_target (Obstructed witness), sorted
    detail: str = ""


def _sample_key(sample: Mapping) -> tuple:
    return tuple(sorted((str(k), str(v)) for k, v in sample.items()))


def _partition_labels(
    labels: Sequence[Label], target_label: Label, lf_flags: Mapping[Label, object]
) -> set[Label]:
    """Allowed = LF-True box labels, target excluded; everything else is forbidden.

    ``"Unknown"`` verdicts, ``False`` verdicts, labels missing from ``lf_flags`` and any row
    support outside ``labels`` are all treated as forbidden (conservative direction: a missing
    verdict can only turn a true ``Feasible`` into ``Obstructed``, never the reverse).
    """
    return {
        lab
        for lab in labels
        if lab != target_label and lf_flags.get(lab) is True
    }


def feasibility_to_payload(result: LFFeasibilityResult) -> dict:
    """JSON-safe dict (deterministic key order; labels as lists)."""
    return {
        "status": result.status,
        "prime": result.prime,
        "sample": [list(pair) for pair in result.sample],
        "rank": result.rank,
        "nrows": result.nrows,
        "n_projected_rows": result.n_projected_rows,
        "n_allowed": result.n_allowed,
        "n_forbidden": result.n_forbidden,
        "residual_support": [list(lab) for lab in result.residual_support],
        "detail": result.detail,
    }


def lf_reduction_feasible_mod_p(
    rows: Sequence[Row],
    labels: Sequence[Label],
    target_label: Label,
    lf_flags: Mapping[Label, object],
    sample: Mapping,
    prime: int,
    column_order: Sequence[Label] | None = None,
) -> LFFeasibilityResult:
    """Span test: can ``target_label`` be reduced through LF-True labels only (mod ``prime``)?

    Independent of the normal-form ranking: the projected system is reduced with a plain RREF
    whose ``column_order`` (default: forbidden labels sorted, target last) only affects pivot
    choice, never the answer. ``BadSpecialization`` is reported as a status, not raised.
    """
    allowed = _partition_labels(labels, target_label, lf_flags)
    try:
        # ``assemble_rows_mod_p`` never touches its family argument (kept for signature
        # compatibility with the certificate call sites), so ``None`` is safe here.
        matrix = assemble_rows_mod_p(None, rows, sample, prime)
    except BadSpecialization as exc:
        return LFFeasibilityResult(
            status=STATUS_BAD_SPECIALIZATION,
            prime=prime,
            sample=_sample_key(sample),
            rank=0,
            nrows=0,
            n_projected_rows=0,
            n_allowed=len(allowed),
            n_forbidden=0,
            residual_support=(),
            detail=str(exc),
        )

    projected: list[dict] = []
    forbidden_present: set[Label] = set()
    for row in matrix:
        proj = {c: v for c, v in row.items() if c not in allowed}
        if proj:
            projected.append(proj)
            forbidden_present.update(proj)
    forbidden_present.discard(target_label)

    if column_order is None:
        column_order = [*sorted(forbidden_present), target_label]
    res = rref_mod_p(projected, prime, column_order=column_order)
    residual = _reduce_vector_by_pivots({target_label: 1}, res.pivots, prime)

    status = STATUS_FEASIBLE if not residual else STATUS_OBSTRUCTED
    return LFFeasibilityResult(
        status=status,
        prime=prime,
        sample=_sample_key(sample),
        rank=res.rank,
        nrows=len(matrix),
        n_projected_rows=len(projected),
        n_allowed=len(allowed),
        n_forbidden=len(forbidden_present),
        residual_support=tuple(sorted(residual)),
        detail="" if not residual else "target unit vector not in projected row span",
    )


# --- explicit coefficients (only meaningful after Feasible) ------------------------------------


def _tracked_elimination(
    pairs: list[tuple[dict, dict]], prime: int, column_order: Sequence[Label]
) -> dict:
    """Mutually reduced pivots over the projected part, tracking the full-row companion.

    Mirrors the historical dict elimination of ``sparse_rref`` (including full back-
    substitution) but carries, for every pivot row, the same linear combination applied to the
    *unprojected* rows. Deterministic in the input order and ``column_order``.
    """
    pos = {c: i for i, c in enumerate(column_order)}
    fallback = len(pos)

    def _axpy2(dst: tuple[dict, dict], src: tuple[dict, dict], factor: int) -> None:
        for d, s in zip(dst, src):
            for c, v in s.items():
                nv = (d.get(c, 0) + factor * v) % prime
                if nv:
                    d[c] = nv
                elif c in d:
                    del d[c]

    pivots: dict = {}  # col -> (proj_row, full_row), proj pivot coefficient == 1
    for proj0, full0 in pairs:
        proj = dict(proj0)
        full = dict(full0)
        for col, (pproj, pfull) in pivots.items():
            f = proj.get(col, 0)
            if f:
                _axpy2((proj, full), (pproj, pfull), -f)
        if not proj:
            continue  # relation among allowed labels only — carries no information here
        col = min(proj, key=lambda c: (pos.get(c, fallback), c))
        inv = pow(proj[col], prime - 2, prime)
        proj = {c: v * inv % prime for c, v in proj.items()}
        full = {c: v * inv % prime for c, v in full.items()}
        for pc, (pproj, pfull) in pivots.items():
            f = pproj.get(col, 0)
            if f:
                _axpy2((pproj, pfull), (proj, full), -f)
        pivots[col] = (proj, full)
    return pivots


def lf_reduction_coefficients_mod_p(
    rows: Sequence[Row],
    labels: Sequence[Label],
    target_label: Label,
    lf_flags: Mapping[Label, object],
    sample: Mapping,
    prime: int,
    column_order: Sequence[Label] | None = None,
) -> tuple[LFFeasibilityResult, dict]:
    """One explicit LF-supported reduction ``J[t] = sum coeffs[l] * J[l]`` (mod ``prime``).

    Returns ``(feasibility_result, coeffs)``; ``coeffs`` is empty unless ``Feasible``. The
    particular solution is the deterministic Gaussian one (same pivot priorities as the span
    test), so at generic points its support is a fixed set and its values are evaluations of a
    single rational function of the parameters — exactly what reconstruction needs.
    """
    allowed = _partition_labels(labels, target_label, lf_flags)
    try:
        matrix = assemble_rows_mod_p(None, rows, sample, prime)
    except BadSpecialization as exc:
        return (
            LFFeasibilityResult(
                status=STATUS_BAD_SPECIALIZATION,
                prime=prime,
                sample=_sample_key(sample),
                rank=0,
                nrows=0,
                n_projected_rows=0,
                n_allowed=len(allowed),
                n_forbidden=0,
                residual_support=(),
                detail=str(exc),
            ),
            {},
        )

    pairs: list[tuple[dict, dict]] = []
    forbidden_present: set[Label] = set()
    for row in matrix:
        proj = {c: v for c, v in row.items() if c not in allowed}
        if proj:
            pairs.append((proj, row))
            forbidden_present.update(proj)
    forbidden_present.discard(target_label)

    if column_order is None:
        column_order = [*sorted(forbidden_present), target_label]
    pivots = _tracked_elimination(pairs, prime, column_order)

    # Reduce (e_t, e_t) by the tracked pivots: the projected part becomes the canonical
    # residual; when it vanishes the full part equals ``e_t - y`` with ``y`` in the row span
    # and ``proj(y) = e_t`` — its allowed-column entries ARE the reduction coefficients.
    vproj = {target_label: 1}
    vfull = {target_label: 1}
    for col, (pproj, pfull) in pivots.items():
        f = vproj.get(col, 0)
        if f:
            for c, v in pproj.items():
                nv = (vproj.get(c, 0) - f * v) % prime
                if nv:
                    vproj[c] = nv
                elif c in vproj:
                    del vproj[c]
            for c, v in pfull.items():
                nv = (vfull.get(c, 0) - f * v) % prime
                if nv:
                    vfull[c] = nv
                elif c in vfull:
                    del vfull[c]

    feasible = not vproj
    coeffs = {}
    if feasible:
        coeffs = {c: v % prime for c, v in vfull.items() if c in allowed and v % prime}
    result = LFFeasibilityResult(
        status=STATUS_FEASIBLE if feasible else STATUS_OBSTRUCTED,
        prime=prime,
        sample=_sample_key(sample),
        rank=len(pivots),
        nrows=len(matrix),
        n_projected_rows=len(pairs),
        n_allowed=len(allowed),
        n_forbidden=len(forbidden_present),
        residual_support=tuple(sorted(vproj)),
        detail="" if feasible else "target unit vector not in projected row span",
    )
    return result, coeffs
