"""Minimal sparse modular RREF over GF(p) with a caller-supplied pivot-column order (spec §5.9).

Rows are sparse dicts ``{column: value}`` where a *column* is any hashable, orderable key (an
integer index, or a label tuple). The system is homogeneous (``A . J = 0``); this routine
reduces it to reduced row-echelon form, choosing pivots on the highest-priority columns first so
that low-priority columns (the intended masters) tend to remain free.

This is deliberately minimal: no matrix assembly from parametric rows and no target normal-form
extraction — those come in a later pass. Everything here is pure integer arithmetic modulo p.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass

Column = object
Row = dict


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

    @property
    def rank(self) -> int:
        return len(self.pivots)

    def reduced_rows(self) -> list[dict]:
        return [self.pivots[c] for c in self.pivot_order]


def rref_mod_p(
    rows: Iterable[dict],
    prime: int,
    column_order: Sequence | None = None,
) -> RREFResult:
    """Reduce ``rows`` to RREF modulo ``prime``, preferring pivots per ``column_order``.

    ``column_order`` lists columns from highest to lowest pivot priority (e.g. the ranking's
    elimination order). Columns not listed are pivoted last, in sorted order. Returns the pivot
    rows (fully reduced, pivot coefficient 1), the pivot/free column split, and the rank.
    """
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

    pivots: dict = {}
    pivot_order: list = []
    for col in order:
        idx = next((i for i, r in enumerate(active) if r.get(col, 0) != 0), None)
        if idx is None:
            continue
        prow = active.pop(idx)
        inv = pow(prow[col], prime - 2, prime)
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

    pivot_set = set(pivots)
    free_cols = [c for c in all_cols if c not in pivot_set]
    return RREFResult(
        prime=prime,
        pivots=pivots,
        pivot_order=pivot_order,
        free_cols=free_cols,
        all_cols=all_cols,
    )
