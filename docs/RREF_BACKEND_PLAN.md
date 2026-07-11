# RREF backend plan (Perf.7, branch `perf/rref-backend-prototype`)

Pure symbolic mathematics / finite-field linear algebra. Scope: make the
modular RREF kernel (`src/parametric_ibp_lf_reducer/sparse_rref.py`) faster
without changing any mathematical result, gate, or default behavior.

## Current implementation (backend `"dict"`)

- Rows are Python dicts `{column: int}`, where a column is a **label tuple**
  (or int); arithmetic is plain Python ints mod p.
- Pivot loop (spec §5.9): walk the caller-supplied elimination order; for each
  pivot column, find the first active row containing it, invert the pivot
  (`pow(x, p-2, p)`), then `_axpy` it out of every other active row **and**
  every established pivot row (full RREF maintained incrementally — forward
  and backward passes are interleaved, not separable).
- Complexity is dominated by `_axpy` dict get/set/del with **tuple hashing**
  on every access, and per-entry Python-int `%` arithmetic.

## Why RREF still dominates corrected Example 4* after Perf.5/Perf.6

- Perf.5 (shared-RREF for linear-LHS multi-target) and Perf.6
  (certificate-point RREF reuse) eliminated all *redundant* RREFs: every
  remaining elimination is now computed exactly once.
- The measured Perf.6 run spends ~2900s in one `rref_mod_p` call
  (12360 rows, 972 labels, `selected_rank=9924`) — ~2/3 of the ~1h15m wall.
- Orchestration-level reuse is therefore exhausted; only a faster **kernel**
  (cheaper per-entry work, better data layout) can reduce this further.

## Candidate backends

| # | Candidate | Idea | Assessment |
|---|---|---|---|
| A | Optimized Python integer-column sparse dict (**prototyped here**) | Map label columns to int ids before the pivot loop; identical algorithm, int keys only; map back after | Lowest risk; wins only tuple-hash overhead; result identical by construction |
| B | Sorted-list / int-array sparse rows | Rows as parallel sorted arrays `(cols[], vals[])`; `_axpy` becomes a merge | Enables C-friendly layout; more invasive; merge cost vs dict must be measured |
| C | CSR-like representation | Whole-matrix CSR with row slicing | Poor fit for row-mutating elimination (constant re-packing); only useful as an input format to D/E/F/G |
| D | Dense / block fallback for high-density rows | Once a row's density passes a threshold, switch it to a dense array (numpy `int64` with explicit mod) | Real win if fill-in creates dense tails; needs fill-in counters (added in this pass) to justify |
| E | Numba over int-array rows | JIT the merge/axpy kernels after B | Explicitly out of scope this pass; requires B first |
| F | C++/Rust extension | Native kernel, same semantics | Later; only after B/D show the layout is right |
| G | CUDA | GPU elimination | Only after a compact (B/C-style) representation exists; not before |

## Risks

- Silent divergence from the certified path → mitigated: experimental backend
  funnels through the **same** `_eliminate` loop (column relabeling is
  bijective), default stays `"dict"`, equivalence tests compare full
  `RREFResult` structures.
- Pivot-order sensitivity: the elimination order defines which columns become
  masters; any backend must preserve it exactly (tested).
- Fill-in blowup with array representations (B/D) — counters
  (`fill_in_ratio`, row-nnz before/after) now measurable per call.
- Windows/spawn parallelism already shown negative (Perf.3); backends must not
  depend on multiprocessing.

## Proposed first prototype (this branch)

1. **Counters** in `rref_mod_p(collect_stats=True)`: n_rows, n_cols,
   nnz_initial/final, row-nnz max/median before/after, rank, pivot count,
   modular inversions, fill-in ratio, elimination time. JSON-safe, optional,
   no stdout.
2. **Backend A** as `rref_backend="int_sparse_experimental"` (default remains
   `"dict"`), plus equivalence tests and small benchmarks.
3. Decision gate: if A shows a real win on representative shapes, wire an
   opt-in flag upward; if not, record the negative result here and proceed to
   B (int-array rows) as the next candidate, since it is the prerequisite for
   D/E/F/G anyway.

## Measured results (this branch)

See `tests/test_rref_backend.py` (equivalence, 15 tests) and
`scripts/bench_rref_backends.py` (reproducible; synthetic label-tuple-keyed
matrices, shapes chosen to mimic the real runs; prime `2_147_483_629`; single
process, wall-clock seconds, this machine, 2026-07-11).

| Matrix (rows×cols, ~nnz/row) | dict | int_sparse_experimental | ratio |
|---|---|---|---|
| tiny 200×150, ~6 (rank 150) | 0.29s | 0.30s | 1.01x |
| medium 1000×800, ~8 (rank 799) | 51.56s | 46.55s | 0.90x |
| D4-like 3000×2400, ~8 (rank 2400) | 1495.48s | 1372.65s | 0.92x |

Results identical across backends on every shape (checked in-script and in the
test suite).

### Decision (gate from step 3)

Backend A is **equivalence-verified but only ~8–10% faster** at scale: tuple
hashing is a minor cost; the elimination is dominated by `_axpy` dict traffic
and Python-int modular arithmetic, which A does not change. Verdict:

- Keep `"dict"` as the default; keep `int_sparse_experimental` as an opt-in
  kernel-level flag (harmless, mildly positive, and its column-relabeling
  pass is the required front-end for every array-based candidate).
- **Do not** plumb a user-facing option upward for a <10% win.
- Next candidate: **B (sorted int-array rows)** — it is the prerequisite for
  the dense-fallback (D), Numba (E), and native (F/G) paths, which is where
  the order-of-magnitude headroom lives. The new counters (fill-in ratio,
  row-nnz profiles) should be captured from a real Example 4* run first to
  size D's density threshold.
