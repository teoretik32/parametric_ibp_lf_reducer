"""Perf.10: experimental Numba int-array backend for the sparse modular RREF kernel.

This module is imported *lazily* by :mod:`.sparse_rref` when the caller explicitly requests
``backend="numba_int_array_experimental"``; importing the package never touches numba. If numba
is not installed the import of this module raises ``ImportError`` and ``rref_mod_p`` converts
that into a clear :class:`~.sparse_rref.BackendUnavailable`.

Representation: each row is a pair of parallel C-contiguous ``int64`` arrays ``(cols, vals)``
sorted by column id, where column ids are the *elimination-order positions* ``0..k-1`` produced
by the same bijective relabeling the ``int_sparse_experimental`` backend uses. Values live in
``[1, p-1]``. The whole pivot loop runs inside one ``@njit`` kernel (no per-row Python calls);
pivot selection and the axpy/merge arithmetic are byte-for-byte the historical algorithm from
``sparse_rref._eliminate``:

- pivots are chosen on the lowest remaining column id (== highest ranking priority) from the
  *first* active row (in original row order) with a nonzero entry there;
- the pivot row is scaled to leading 1 (one modular inversion per pivot, counted identically);
- that column is then eliminated from every other row, active and established-pivot alike,
  so full RREF is maintained incrementally exactly as in the dict backend.

Overflow safety: all products are ``a * b`` with ``a, b in [0, p)``; the wrapper enforces
``p < 2**31`` so products fit comfortably in ``int64`` (< 2**62). No float arithmetic anywhere.

The first call in a fresh environment pays a one-time JIT compile (a few seconds);
``cache=True`` persists the compiled kernels to ``__pycache__`` for later processes.
"""

from __future__ import annotations

import numpy as np

try:
    from numba import njit
    from numba.typed import List as NumbaList
except ImportError as exc:  # pragma: no cover - exercised only when numba is absent
    raise ImportError(
        "numba is required for the 'numba_int_array_experimental' rref backend"
    ) from exc

#: Exclusive upper bound for primes accepted by this backend (int64 product safety).
MAX_PRIME_EXCLUSIVE = 1 << 31


@njit(cache=True)
def _bisect(cols, col):
    """Index of ``col`` in the sorted array ``cols``, or -1 if absent."""
    lo = 0
    hi = cols.shape[0]
    while lo < hi:
        mid = (lo + hi) // 2
        if cols[mid] < col:
            lo = mid + 1
        else:
            hi = mid
    if lo < cols.shape[0] and cols[lo] == col:
        return lo
    return -1


@njit(cache=True)
def _powmod(base, exp, p):
    """base**exp mod p by binary exponentiation (all operands < 2**31)."""
    result = 1
    b = base % p
    e = exp
    while e > 0:
        if e & 1:
            result = (result * b) % p
        b = (b * b) % p
        e >>= 1
    return result


@njit(cache=True)
def _scale(vals, factor, p):
    """New array ``vals * factor mod p`` (factor nonzero, so no pruning needed)."""
    out = np.empty_like(vals)
    for i in range(vals.shape[0]):
        out[i] = (vals[i] * factor) % p
    return out


@njit(cache=True)
def _axpy_merge(tc, tv, sc, sv, factor, p):
    """``target + factor * source`` (mod p) as a fresh sorted pair, zeros pruned.

    ``factor`` must already be reduced into ``[0, p)``.
    """
    nt = tc.shape[0]
    ns = sc.shape[0]
    oc = np.empty(nt + ns, np.int64)
    ov = np.empty(nt + ns, np.int64)
    i = 0
    j = 0
    k = 0
    while i < nt and j < ns:
        ci = tc[i]
        cj = sc[j]
        if ci < cj:
            oc[k] = ci
            ov[k] = tv[i]
            i += 1
            k += 1
        elif ci > cj:
            v = (factor * sv[j]) % p
            if v != 0:
                oc[k] = cj
                ov[k] = v
                k += 1
            j += 1
        else:
            v = (tv[i] + factor * sv[j]) % p
            if v != 0:
                oc[k] = ci
                ov[k] = v
                k += 1
            i += 1
            j += 1
    while i < nt:
        oc[k] = tc[i]
        ov[k] = tv[i]
        i += 1
        k += 1
    while j < ns:
        v = (factor * sv[j]) % p
        if v != 0:
            oc[k] = sc[j]
            ov[k] = v
            k += 1
        j += 1
    # .copy() releases the oversized backing buffer instead of keeping a view alive.
    return oc[:k].copy(), ov[:k].copy()


@njit(cache=True)
def _eliminate_kernel(rows_c, rows_v, n_rows, n_cols, p):
    """Full pivot loop over column ids ``0..n_cols-1``; mirrors ``sparse_rref._eliminate``.

    Rows are never popped; ``is_pivot`` marks promoted rows, which preserves the relative
    order of the remaining active rows exactly like ``list.pop`` does in the dict backend.
    Returns ``(pivot_cols_in_elimination_order, pivot_col_of_row, n_inversions)``.
    """
    is_pivot = np.zeros(n_rows, np.uint8)
    pivot_of_row = np.full(n_rows, -1, np.int64)
    pivot_seq = np.full(n_rows, -1, np.int64)
    npiv = 0
    inversions = 0
    for col in range(n_cols):
        pidx = -1
        pj = -1
        for i in range(n_rows):
            if is_pivot[i] == 0:
                j = _bisect(rows_c[i], col)
                if j >= 0:
                    pidx = i
                    pj = j
                    break
        if pidx < 0:
            continue
        inv = _powmod(rows_v[pidx][pj], p - 2, p)
        inversions += 1
        rows_v[pidx] = _scale(rows_v[pidx], inv, p)
        pc = rows_c[pidx]
        pv = rows_v[pidx]
        # Eliminate this column from every other row: still-active rows *and*
        # already-established pivot rows (keeps full RREF, as in the dict backend).
        for i in range(n_rows):
            if i == pidx:
                continue
            j = _bisect(rows_c[i], col)
            if j >= 0:
                factor = p - rows_v[i][j]  # == (-f) mod p, in [1, p-1]
                nc, nv = _axpy_merge(rows_c[i], rows_v[i], pc, pv, factor, p)
                rows_c[i] = nc
                rows_v[i] = nv
        is_pivot[pidx] = 1
        pivot_of_row[pidx] = col
        pivot_seq[npiv] = col
        npiv += 1
    return pivot_seq[:npiv].copy(), pivot_of_row, inversions


def eliminate_int_rows(int_rows: list[dict], n_cols: int, prime: int) -> tuple[dict, list, int]:
    """Drop-in replacement for ``_eliminate(int_rows, range(n_cols), prime)``.

    ``int_rows`` are sparse rows keyed by elimination-order column ids with values already
    normalized into ``[1, prime-1]`` (as produced by ``rref_mod_p``'s relabeling front-end).
    Returns ``(pivots, pivot_order, n_inversions)`` with plain Python ``int`` keys/values,
    identical to the dict backend's output by construction.
    """
    if prime >= MAX_PRIME_EXCLUSIVE:
        raise ValueError(
            f"numba rref backend requires prime < 2**31 (int64 product safety); got {prime}"
        )
    if not int_rows:
        return {}, [], 0

    rows_c = NumbaList()
    rows_v = NumbaList()
    for r in int_rows:
        cols = np.fromiter(sorted(r), np.int64, count=len(r))
        vals = np.empty(len(r), np.int64)
        for k in range(cols.shape[0]):
            vals[k] = r[int(cols[k])]
        rows_c.append(cols)
        rows_v.append(vals)

    pivot_seq, pivot_of_row, inversions = _eliminate_kernel(
        rows_c, rows_v, len(int_rows), n_cols, prime
    )

    by_col: dict = {}
    for i in range(len(int_rows)):
        col = int(pivot_of_row[i])
        if col >= 0:
            by_col[col] = {int(c): int(v) for c, v in zip(rows_c[i], rows_v[i])}
    pivot_order = [int(c) for c in pivot_seq]
    # Reproduce the dict backend's insertion order (pivots appear in elimination order).
    pivots = {col: by_col[col] for col in pivot_order}
    return pivots, pivot_order, int(inversions)
