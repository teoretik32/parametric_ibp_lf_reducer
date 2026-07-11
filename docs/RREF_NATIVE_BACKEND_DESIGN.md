# RREF native-backend design (Perf.9 — design only, no implementation)

Pure symbolic mathematics / finite-field linear algebra. This document
designs the next serious RREF backend for
`src/parametric_ibp_lf_reducer/sparse_rref.py`. **Nothing is implemented in
this pass**: the default backend stays `"dict"`, no gate/certificate/LF logic
changes, no tag/release. Inputs: `docs/RREF_BACKEND_PLAN.md` (Perf.7/Perf.8
verdicts), `validation/rref_real_matrix_profile.json` (real medium profile),
`notes/PERF_STATUS.md`.

## 1. Current bottleneck summary

- **Where the time goes.** On the corrected Example 4* full system the single
  remaining `rref_mod_p` call (12360 rows × 972 labels,
  `selected_rank=9924`) costs ~2900 s — ~2/3 of the ~1h15m wall. Perf.5/6
  removed all *redundant* RREFs; Perf.7/8 showed orchestration reuse is
  exhausted. RREF is the dominant cost and only a faster kernel helps.
- **Why it is slow.** Rows are Python dicts `{column: int}`; `_axpy` does a
  Python-level loop of dict get/set/del plus Python-int `(x*y) % p` per
  entry. Two costs dominate: (a) interpreter overhead per nonzero, (b) boxed
  arbitrary-precision int arithmetic where 64-bit machine arithmetic would
  do (`p < 2^31` in all runs: `2_147_483_647`, `..._629`, `..._587`).
- **Measured real profile** (medium subset, 512×917, rank 512): nnz
  2640 → 27762, **fill-in 10.5x**, final row-nnz median 62 / max 114 — i.e.
  ~7–12% density: rows stay sparse even after fill-in, so a *sparse* kernel
  remains the right shape. `int_sparse_experimental` (int column ids, same
  dict kernel) gives only 0.69x; pure-Python sorted-array rows (Candidate B,
  Perf.8) were 1.5–1.9x **slower** — the headroom is not reachable from
  Python bytecode.

## 2–3. Candidate designs

Common front-end for all candidates: the existing
`int_sparse_experimental` relabeling (`col_to_id` over the elimination
`order`, elimination on `0..k-1`, bijective map back). It already preserves
pivot order exactly and is equivalence-tested — every candidate below sits
*behind* it, so column semantics never change.

Exact GF(p) arithmetic, shared by A–F: all primes used are < 2^31, so values
fit int64 and products of two reduced values fit in 63 bits
(`(p-1)^2 < 2^62`); `(a * b) % p` in int64 is exact with no overflow.
Inversion stays `pow(x, p-2, p)` in Python (512–9924 calls per run — never
hot). No floats anywhere; any backend that cannot guarantee 64-bit exactness
(notably CUDA fp paths) is excluded until it can.

### A. Numba int-array sparse rows

Rows as parallel `int64` numpy arrays `(cols[], vals[])`, sorted by col;
`@njit(cache=True)` kernels for scale-row and merge-axpy (the Perf.8 loser —
but compiled, the merge is branch-cheap C, not bytecode).

- Expected speedup: **10–50x on `_axpy`**; whole-call est. 5–20x (pivot
  search + Python loop overhead remains). Real-profile check is the point of
  the prototype.
- Effort: **small (1–2 days)** — two kernels + row converters + tests; the
  elimination *structure* (`_eliminate`) is unchanged.
- Packaging risk: **moderate-low** — optional dependency, binary wheels on
  PyPI, no compiler toolchain for users; risk = numba's Python-version lag
  and import cost (~1–3 s JIT warmup, amortized by `cache=True`).
- Windows risk: **low** — first-class numba wheels; no fork/spawn issues
  (single-process, unlike the Perf.3 negative).
- Correctness/testing: identical `RREFResult` vs `"dict"` on the full
  existing equivalence suite + random/rank-deficient/pivot-order matrices;
  kernels are pure functions on int64 arrays — property-testable in
  isolation against the dict `_axpy`.
- GF(p): int64 mul-mod as above; assert `p < 2**31` at backend selection.

### B. C++ extension (pybind11/nanobind)

Same array layout, kernel in C++.

- Speedup: similar to A (maybe 1.2–2x over A once the whole pivot loop moves
  native). Effort: **medium-large (1–2 weeks)** — build system, CI wheels.
- Packaging: **high risk** — users need wheels per platform/Python or MSVC.
- Windows: MSVC toolchain risk for source installs.
- Testing: same equivalence suite; plus ABI/overflow unit tests.
- GF(p): `uint64_t` mul-mod, exact; optionally Montgomery later.

### C. Rust extension (PyO3/maturin)

As B with a safer language and nicer packaging (`maturin` wheels).

- Speedup: = B. Effort: **medium-large**, plus a Rust toolchain in CI.
- Packaging: **medium-high**; Windows: good (MSVC target mature).
- Testing/GF(p): as B (`u64` mul-mod; `checked_*` in debug builds).

### D. Dense/block fallback for high-density sectors

Switch a row (or trailing block) to a dense `int64` numpy vector once its
density passes a threshold; axpy becomes vectorized `(t + f*s) % p`.

- Speedup: only where density is high. **Measured density is 7–12%** at
  final state on the real medium profile — *no trigger exists today*; must
  re-measure on the full 12360×972 system before this is credible.
- Effort: small-medium, but threshold tuning + two representations in one
  loop. Packaging/Windows: none (numpy only). GF(p): int64 vector mod —
  exact. Verdict: **deferred** pending full-system density data.

### E. Hybrid sparse/dense rows

D generalized: per-row representation chosen dynamically, sparse↔dense
promotion/demotion inside `_eliminate`.

- Strictly a superset of D's complexity with the same missing trigger;
  highest bug surface (representation transitions under mutation).
  **Rejected as a first step**; revisit only if D's trigger appears and D
  alone underperforms.

### F. CUDA (later only)

GPU elimination requires a compact array representation (A's layout) plus
batched pivoting; exact 64-bit modular arithmetic on GPU is doable
(integer ops), but transfers dominate at 12360×972 scale.

- **Explicitly out of scope** until a compact CPU representation (A) exists
  and is adopted; recorded here only as ordering (matches Perf.7 plan G).

## 4. Recommended first prototype: **A — Numba int-array sparse rows**

Chosen because it is the *smallest* prototype that can prove speed on the
real profile: no toolchain, no wheels to build, ~2 kernels, reuses the
existing relabeling front-end and equivalence suite verbatim, and directly
attacks both measured costs (interpreter overhead + boxed ints). B/C pay
their packaging cost only if A's *layout* wins but JIT overhead disappoints;
A produces exactly the evidence needed to justify them. D/E lack a measured
trigger; F is sequenced after A by construction.

Constraint note: this requires lifting the pure-Python constraint to
"pure Python + optional `numba` extra". The backend must import lazily and
degrade to `"dict"` with a clear error if numba is absent — the core package
stays pure-Python.

## 5. Prototype boundary (implementation pass, not this one)

- **New module** `src/parametric_ibp_lf_reducer/sparse_rref_numba.py`:
  `_rows_to_arrays(rows) -> list[(cols, vals)]`, `_arrays_to_rows(...)`,
  `_axpy_merge(tc, tv, sc, sv, factor, p) -> (cols, vals)` (njit),
  `_scale_row(cols, vals, inv, p)` (njit), and
  `_eliminate_arrays(active, n_cols, prime)` mirroring `_eliminate`
  line-for-line. No changes to `_eliminate` itself.
- **`sparse_rref.py`**: add `"numba_int_array_experimental"` to
  `RREF_BACKENDS`; selection raises a clean `ValueError` ("numba not
  installed") if the import fails. `DEFAULT_RREF_BACKEND = "dict"`
  **unchanged**; no caller anywhere passes the new backend by default.
- **Comparison**: `RREFResult` equality (pivots, pivot_order, free_cols,
  all_cols) against `"dict"` — same check the A-backend uses today; stats
  counters must also agree on nnz/rank/pivot_count.
- **Benchmark matrix sizes**: existing synthetic trio (200×150 ~6/row,
  1000×800 ~8/row, 3000×2400 ~8/row via `scripts/bench_rref_backends.py`)
  + the real medium profile 512×917 (`scripts/profile_rref_real_matrix.py`)
  + one stretch run on the full corrected Example 4* 12360×972 system.
- **Pass/fail threshold**: adopt as supported opt-in if **≥3x** vs `"dict"`
  on the real medium profile with byte-identical results; plumb a
  user-facing option upward only if the full-system projection is **≥5x**
  (~2900 s → ≤600 s). Below 3x: record negative result, fall back to B/C
  evaluation. JIT warmup excluded from the ratio but reported.

## 6. Test plan

1. **Equivalence with dict backend**: parametrize the existing
   `tests/test_rref_backend.py` suite over the new backend
   (`pytest.importorskip("numba")`), all 15 cases must be identical.
2. **Random matrices**: seeded random sparse systems across densities
   (2–20 nnz/row), several primes (`2_147_483_647/629/587`), shapes up to
   ~300×300 in CI — full `RREFResult` equality.
3. **Rank-deficient**: duplicated rows, scalar multiples, zero rows,
   `rank < n_rows` asserted identical, `free_cols` identical.
4. **pivot_order**: permuted `column_order` inputs; assert the chosen pivot
   sequence matches `"dict"` exactly (order-sensitivity is a known risk).
5. **modular_normal_form integration**: run the normal-form path with the
   backend forced via its kernel-level flag; assert identical
   `NormalFormRecord`s on the tiny + medium fixtures.
6. **Corrected Example 4\* medium profile**: extend
   `scripts/profile_rref_real_matrix.py` with the new backend column;
   committed JSON must show identical rank/nnz/fill-in and the measured
   ratio.
7. **No certificate/LF changes**: certificate, LF-gate, and reducer
   orchestration tests run untouched (they never select the backend);
   `git diff` of the implementation pass must not touch certificate/LF
   modules.

## 7. Non-goals of the prototype

No default change, no multi-threading, no Montgomery/Barrett tuning, no
dense fallback, no GPU, no removal of `int_sparse_experimental`, no new
required dependencies. No implementation code was written in this design
pass (and no measuring helper was needed — all stats above already exist in
`validation/rref_real_matrix_profile.json` and the Perf.7/8 plan).
