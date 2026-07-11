"""Minimal sparse modular RREF over GF(p) with a caller-supplied pivot-column order (spec §5.9).

Rows are sparse dicts ``{column: value}`` where a *column* is any hashable, orderable key (an
integer index, or a label tuple). The system is homogeneous (``A . J = 0``); this routine
reduces it to reduced row-echelon form, choosing pivots on the highest-priority columns first so
that low-priority columns (the intended masters) tend to remain free.

This is deliberately minimal: no matrix assembly from parametric rows and no target normal-form
extraction — those come in a later pass. Everything here is pure integer arithmetic modulo p.

Perf.7 additions (both optional, defaults preserve the historical behavior exactly):

- ``collect_stats=True`` attaches a JSON-safe counters dict to ``RREFResult.stats``
  (row/column/nnz/fill-in/inversion counters; no stdout output).
- ``backend="int_sparse_experimental"`` runs the *same* elimination on integer column ids
  (labels are mapped to ints before the pivot loop and mapped back afterwards), avoiding
  tuple-hashing overhead inside ``_axpy``/``get``. The default backend remains ``"dict"``;
  results are identical by construction (one bijective relabeling of columns).

Perf.10 addition (optional, opt-in only):

- ``backend="numba_int_array_experimental"`` runs the same pivot loop as one ``@njit``
  kernel over sorted parallel ``int64`` arrays (see :mod:`.sparse_rref_numba`). It reuses
  the exact int relabeling front-end above, requires ``prime < 2**31`` (int64 product
  safety) and requires numba (``pip install .[speed]``); importing this package never
  touches numba. If numba is missing, requesting the backend raises
  :class:`BackendUnavailable`. The default backend is still ``"dict"``.
"""

from __future__ import annotations

import importlib.util
import time
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from statistics import median

Column = object
Row = dict

DEFAULT_RREF_BACKEND = "dict"
NUMBA_RREF_BACKEND = "numba_int_array_experimental"
RREF_BACKENDS = ("dict", "int_sparse_experimental", NUMBA_RREF_BACKEND)

#: Primes accepted by the numba backend must satisfy ``prime < 2**31`` (int64 products).
_NUMBA_MAX_PRIME_EXCLUSIVE = 1 << 31


class BackendUnavailable(RuntimeError):
    """The requested rref backend cannot run in this environment (e.g. numba missing)."""


def rref_backend_available(backend: str) -> bool:
    """True if ``backend`` is registered *and* can actually run here (cheap, import-free)."""
    if backend not in RREF_BACKENDS:
        return False
    if backend == NUMBA_RREF_BACKEND:
        return importlib.util.find_spec("numba") is not None
    return True


def _axpy(target: dict, source: dict, factor: int, p: int) -> None:
    """target += factor * source  (mod p), pruning zero entries."""
    for c, v in source.items():
        nv = (target.get(c, 0) + factor * v) % p
        if nv:
            target[c] = nv
        elif c in target:
            del target[c]


def _normalize_row(row: dict, p: int) -> dict:
    return {c: v % p for c, v in row.items() if v % p != 0}


@dataclass
class RREFResult:
    prime: int
    pivots: dict  # pivot_column -> reduced row (pivot coefficient == 1)
    pivot_order: list  # pivot columns, in the order they were eliminated
    free_cols: list  # columns present but without a pivot
    all_cols: list
    stats: dict | None = field(default=None, compare=False)  # Perf.7: optional counters

    @property
    def rank(self) -> int:
        return len(self.pivots)

    def reduced_rows(self) -> list[dict]:
        return [self.pivots[c] for c in self.pivot_order]


def _eliminate(active: list[dict], order: Sequence, prime: int) -> tuple[dict, list, int]:
    """Core pivot loop, generic over the column key type (labels or int ids).

    Returns ``(pivots, pivot_order, n_inversions)``. This is byte-for-byte the historical
    algorithm; both backends funnel through it so verdicts cannot diverge.
    """
    pivots: dict = {}
    pivot_order: list = []
    inversions = 0
    for col in order:
        idx = next((i for i, r in enumerate(active) if r.get(col, 0) != 0), None)
        if idx is None:
            continue
        prow = active.pop(idx)
        inv = pow(prow[col], prime - 2, prime)
        inversions += 1
        prow = {c: (v * inv) % prime for c, v in prow.items()}
        # eliminate this column from every other active row...
        for r in active:
            f = r.get(col, 0)
            if f:
                _axpy(r, prow, -f, prime)
        # ...and from already-established pivot rows (keeps full RREF).
        for pc in pivots:
            f = pivots[pc].get(col, 0)
            if f:
                _axpy(pivots[pc], prow, -f, prime)
        pivots[col] = prow
        pivot_order.append(col)
    return pivots, pivot_order, inversions


def _row_nnz_stats(rows: Iterable[dict]) -> tuple[int, int, float]:
    """(total nnz, max row nnz, median row nnz) — JSON-safe scalars, cheap single pass."""
    sizes = [len(r) for r in rows]
    if not sizes:
        return 0, 0, 0.0
    return sum(sizes), max(sizes), float(median(sizes))


def rref_mod_p(
    rows: Iterable[dict],
    prime: int,
    column_order: Sequence | None = None,
    *,
    backend: str | None = None,
    collect_stats: bool = False,
) -> RREFResult:
    """Reduce ``rows`` to RREF modulo ``prime``, preferring pivots per ``column_order``.

    ``column_order`` lists columns from highest to lowest pivot priority (e.g. the ranking's
    elimination order). Columns not listed are pivoted last, in sorted order. Returns the pivot
    rows (fully reduced, pivot coefficient 1), the pivot/free column split, and the rank.

    ``backend`` (Perf.7) selects the elimination representation; ``None`` uses
    :data:`DEFAULT_RREF_BACKEND` (``"dict"``, the historical implementation).
    ``collect_stats`` (Perf.7) attaches a JSON-safe counters dict as ``result.stats``;
    it never affects the result and prints nothing.
    """
    chosen = DEFAULT_RREF_BACKEND if backend is None else backend
    if chosen not in RREF_BACKENDS:
        raise ValueError(f"unknown rref_backend {chosen!r}; expected one of {RREF_BACKENDS}")

    active = [r for r in (_normalize_row(row, prime) for row in rows) if r]
    all_cols = sorted({c for r in active for c in r})
    seen = set()
    order: list = []
    for c in column_order or ():
        if c in all_cols and c not in seen:
            order.append(c)
            seen.add(c)
    for c in all_cols:
        if c not in seen:
            order.append(c)

    stats: dict | None = None
    if collect_stats:
        nnz0, max0, med0 = _row_nnz_stats(active)
        stats = {
            "backend": chosen,
            "n_rows": len(active),
            "n_cols": len(all_cols),
            "nnz_initial": nnz0,
            "row_nnz_before_max": max0,
            "row_nnz_before_median": med0,
        }
        t0 = time.perf_counter()

    if chosen in ("int_sparse_experimental", NUMBA_RREF_BACKEND):
        # Bijective column relabeling: elimination order becomes 0..k-1, so ordering is
        # preserved exactly; all dict keys inside the pivot loop are small ints.
        col_to_id = {c: i for i, c in enumerate(order)}
        int_active = [{col_to_id[c]: v for c, v in r.items()} for r in active]
        if chosen == NUMBA_RREF_BACKEND:
            if prime >= _NUMBA_MAX_PRIME_EXCLUSIVE:
                raise ValueError(
                    f"backend {chosen!r} requires prime < 2**31 "
                    f"(int64 products must fit); got {prime}"
                )
            try:
                from . import sparse_rref_numba
            except ImportError as exc:
                raise BackendUnavailable(
                    f"backend {chosen!r} requires numba; install the 'speed' extra "
                    f"(pip install .[speed]) or pick another backend"
                ) from exc
            int_pivots, int_pivot_order, inversions = sparse_rref_numba.eliminate_int_rows(
                int_active, len(order), prime
            )
        else:
            int_pivots, int_pivot_order, inversions = _eliminate(
                int_active, range(len(order)), prime
            )
        id_to_col = order  # id i -> order[i]
        pivots = {
            id_to_col[pc]: {id_to_col[c]: v for c, v in prow.items()}
            for pc, prow in int_pivots.items()
        }
        pivot_order = [id_to_col[i] for i in int_pivot_order]
    else:
        pivots, pivot_order, inversions = _eliminate(active, order, prime)

    if stats is not None:
        # Forward elimination and back-substitution are interleaved in this implementation
        # (full RREF is maintained incrementally), so a single elimination time is reported.
        stats["elimination_time_s"] = time.perf_counter() - t0
        nnz1, max1, med1 = _row_nnz_stats(pivots.values())
        stats.update(
            rank=len(pivots),
            pivot_count=len(pivot_order),
            inversions=inversions,
            nnz_final=nnz1,
            row_nnz_after_max=max1,
            row_nnz_after_median=med1,
            fill_in_ratio=(nnz1 / nnz0) if nnz0 else None,
        )

    pivot_set = set(pivots)
    free_cols = [c for c in all_cols if c not in pivot_set]
    return RREFResult(
        prime=prime,
        pivots=pivots,
        pivot_order=pivot_order,
        free_cols=free_cols,
        all_cols=all_cols,
        stats=stats,
    )
