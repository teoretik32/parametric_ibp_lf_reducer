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

## Perf.8: real-matrix stats + candidate-B feasibility (this branch)

Real (non-synthetic) profile from the corrected Example 4* row system,
medium subset (`scripts/profile_rref_real_matrix.py`, output committed at
`validation/rref_real_matrix_profile.json`; prime `2_147_483_647`, sample
`ep=15/7`, this machine, 2026-07-11):

| Stat | value |
|---|---|
| shape / rank | 512×917, rank 512 (128 algebraic + 384 coordinate-IBP rows) |
| nnz initial | 2640 (row-nnz max 7, median 5.5) |
| nnz final | 27762 (row-nnz max 114, median 62) |
| fill-in ratio | **10.5x** |
| elimination | dict 0.81s, int_sparse_experimental 0.56s (**0.69x**) |

Notes: the int backend's relative win is larger on the real matrix than on
the synthetic shapes (0.69x vs 0.90–0.92x) — real label tuples are longer and
hash-heavier than the synthetic ones. Final pivot-row density is only
~7–12% of 917 columns: fill-in is real (10.5x) but the rows stay **sparse**,
so a dense fallback (D) has no obvious trigger at this size; re-measure on
the full 12360×972 system before revisiting D.

### Candidate B micro-experiment (merge-axpy vs dict-axpy)

Throwaway micro-benchmark (recorded here for the negative result): axpy of
one sparse row into another at the *measured real* nnz profile, prime
`2_147_483_647`, 2000 reps over 400 pre-built rows, pure Python. The
sorted-int-array merge variant allocates fresh `(cols[], vals[])` output
arrays per axpy (the natural functional form); the dict variant mutates.

| row nnz | dict-axpy | merge-axpy | ratio |
|---|---|---|---|
| 6 (initial median) | 0.006s | 0.009s | 1.51x |
| 62 (final median) | 0.043s | 0.083s | **1.91x** |
| 114 (final max) | 0.074s | 0.124s | 1.68x |

### Decision (Perf.8)

**Candidate B in pure Python is rejected.** The merge kernel is 1.5–1.9x
*slower* than dict axpy exactly in the nnz range the real elimination
inhabits — Python-level list traversal cannot beat the C-implemented dict
fast path. B's layout only pays off once the merge itself leaves Python
bytecode (Numba/E or native F/G), so:

- Default stays `"dict"`; `int_sparse_experimental` remains the opt-in flag
  and is confirmed useful on real matrices (0.69x here).
- B/C/D are **not** implemented in pure Python. The next credible step is
  E/F (JIT/native merge kernel over int-array rows), which is out of scope
  for this codebase's pure-Python constraint — parked unless that constraint
  is lifted.
- Kernel work on this branch is complete; remaining RREF headroom is
  documented, not actionable within current constraints.
