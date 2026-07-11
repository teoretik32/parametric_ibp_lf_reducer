# PERF status

## Environment footgun (FIXED 2026, pre-Perf.1)

- Symptom: `python -m parametric_ibp_lf_reducer` resolved to the old
  site-packages copy **v0.1.0** (pre-Perf.0: no `timings` in diagnostics,
  old CLI that accepted a bare input path without the `reduce` subcommand),
  silently shadowing the working tree.
- Fix: editable install from the repo root:
  `python -m pip install -e .`
  Now: Version **0.1.1**, Location `C:\Python313\Lib\site-packages`,
  Editable project location `B:\soft\math_scratch`; import resolves to
  `B:\Soft\math_scratch\src\parametric_ibp_lf_reducer\__init__.py`
  (verified from an unrelated cwd).

## Rule for ALL perf commands

Run perf/timing commands against the working tree only — either rely on the
**editable install** (default now) or force `PYTHONPATH=src` explicitly.
After any reinstall/venv change, re-verify with:

```
python -c "import parametric_ibp_lf_reducer; print(parametric_ibp_lf_reducer.__file__)"
```

Expected: path under `B:\soft\math_scratch\src\`.

## Baseline (pre-Perf.1)

- Stage timing snapshots (tiny + fast D4 deg1): see `notes/perf_timings.md`.
- Hotspot: `ranking` = 96.7–98.7% of `records_total` → Perf.1 target.
- Post-fix health: `python -m pytest` — 165 passed; `ruff check .` — clean.

## Perf.1 — DONE: ranking hoisted out of per-sample records

- Change: label ranking is computed **once per run** (`ranking_once` stage)
  and reused by every per-sample normal-form record, instead of being
  recomputed inside each record. Old per-record `ranking` key is kept in
  `STAGE_KEYS` for schema stability and now reads `0.0`.
- Measured (same configs as the baseline snapshots, editable install):
  - Tiny CLI (`examples/tiny_success_input.wl.txt`):
    `records_total` 0.064 s, `ranking_once` 0.030 s (was: ranking 1.28 s of
    1.32 s total, 96.7%). Status still `Success`.
  - D4 deg1 fast config (labels=T+M1..M5, 2 primes x 9 samples,
    `max_ibp_degree=1`): **1.13 s wall** (was 13.87 s, ~12x);
    `records_total` 0.892 s, `ranking_once` 0.719 s (80.6%),
    `row_generation_total` 0.151 s. Status unchanged:
    `NormalFormNotLocallyFinite`, same non-LF terms.
- Guarantees held: no math/result changes, LF/certificate gates untouched,
  monkeypatch tests (`tests/test_perf1_ranking_hoist.py`) pin the
  call-count contract (ranking computed exactly once).
- Health after Perf.1: `python -m pytest -q` — 168 passed (dot count,
  72+72+24, no failure markers); `ruff check .` — clean.
- Next candidates (from the post-hoist D4 profile): `ranking_once` itself
  (0.72 s, 80.6% of records pass) and `coordinate_rows` /
  `assemble_rows_mod_p` + `rref_mod_p` (~0.15 s + ~0.17 s combined).

## Perf.2 — DONE: LF valuation caching inside `ranking_once`

- Change: `valuations.py` now memoizes the per-family tropical
  precomputation (family-level cache keyed by identity) and per-label
  `is_locally_finite` verdicts, so the LF sweep in `ranking_once` no longer
  recomputes valuations for repeated labels. No changes to ranking order,
  LF semantics, or certificate gates; no parallelism.
- Verdict-equality harness (vs pre-Perf.2 baseline `lf_baseline.json`):
  D4 208-label row-union + tiny 81-label box — **0 mismatches** on both.
- Measured (exact Perf.1 D4 deg1 fast config: labels=T+M1..M5,
  2 primes x 9 samples, `max_ibp_degree=1`, editable install):
  **0.821 s wall** (was 1.13 s); `ranking_once` **0.417 s** (was 0.719 s);
  `records_total` 0.589 s (was 0.892 s); `row_generation_total` 0.159 s.
  Status unchanged: `NormalFormNotLocallyFinite`, certificate Passed 3/3.
  Standalone cold LF sweep of the 208-label union: 0.428 s (was ~0.72 s
  inside `ranking_once`).
- Health after Perf.2 (verified, not dot-counted): collect-only —
  **240 tests collected**; `python -m pytest` — **233 passed, 7 skipped**
  (skips = `RUN_D4_FULL`-gated heavy e2e + env-gated, unchanged);
  `ruff check .` — clean. NB: earlier "165/168 passed" figures in the
  Perf.0/Perf.1 sections were dot-count estimates from `tail`-truncated
  output; actual collection was already larger.
- Next candidates (post-Perf.2 D4 profile): `ranking_once` remainder
  (0.42 s, ~71% of records pass — the non-LF part of ranking), then
  `coordinate_rows` 0.156 s, `assemble_rows_mod_p` 0.090 s,
  `rref_mod_p` 0.080 s.

## Perf.3 — DONE (capability landed; parallel is a LOSS at current sizes)

- Change: `collect_normal_form_records(..., jobs=N)` /
  `ReducerConfig.jobs` / CLI `--jobs` — independent `(prime, sample)`
  record points computed in worker **processes** via
  `ProcessPoolExecutor`. `jobs=1` (default) is the exact serial path,
  bit-for-bit identical to pre-Perf.3. Math, ranking, LF/certificate
  gates untouched; ranking is hoisted (Perf.1) and computed in the
  parent, workers only do per-point modular normal-form work.
- Correctness: `tests/test_perf3_jobs_equality.py` — 12 passed
  (serial vs parallel result equality incl. statuses, records,
  diagnostics); ad-hoc heavy-config key comparison (status, labels,
  record counts, certificate) equal for jobs=2/4/8.
- Measured (editable install, Windows 10, spawn start method;
  D4 deg1, labels=T+M1..M5, `max_ibp_degree=1`, 2 primes):
  - fast 2x9 (18 points): jobs=1 **0.814 s** wall; jobs=4 **3.929 s**.
  - heavy 2x49 (98 points): jobs=1 **1.636 s**; jobs=2 **2.693 s**;
    jobs=4 **4.107 s**; jobs=8 **6.329 s** — monotonically WORSE.
  - Root cause (measured, not guessed): Windows `spawn` pool
    spinup + package/sympy import round-trip costs **0.68 s** (1
    worker) to **1.56 s** (8 workers), while the ENTIRE serial
    records pass on the heaviest representative config is only
    1.39 s (`assemble_rows_mod_p` 0.50 s + `rref_mod_p` 0.46 s of
    parallelizable work). Fixed startup can never amortize; the
    per-point task (~14 ms) is far too small.
- Verdict: keep `jobs=1` as the default (it already is). Do NOT
  recommend `--jobs` unless per-run record work reaches tens of
  seconds (e.g. many primes × dense sample grids on much larger
  families). The plumbing is correct and tested, so it is ready if
  such workloads appear; revisit with a persistent worker pool or
  `fork`-style start (non-Windows) before expecting wins.
- Health after Perf.3: `python -m pytest` — **245 passed, 7 skipped**
  (233 + 12 new equality tests; skips unchanged); `ruff check .` —
  clean; `ruff format --check` on the five Perf.3-touched files
  (`records.py`, `reducer.py`, `api.py`, `cli.py`,
  `tests/test_perf3_jobs_equality.py`) — already formatted.

## Perf.4 — DONE: heavy corrected Example 4* profile (script-level stages)

- Instrumented `scripts/run_example4_star_corrected.py` with script-level
  stage timings (`_timed` wrapper, stderr log, exported under
  `perf4_timings` in `validation/example4_star_corrected_diagnostics.json`).
  Non-invasive (no library changes), so the instrumentation stays in the
  script permanently.
- Heavy run (RUN_D4_FULL-class config, two LHS targets, rank 9924,
  12360 rows x 972 labels; editable install): **total runtime ~2h24m**,
  status Success, combined certificate Passed 5/5, results identical to
  the certified baseline.
- Stage profile (seconds):
  - `row_generation_shared` **56.7** (shared; `rows_generated_once=True`,
    `row_generation_total=0.0` inside both subruns — sharing works);
  - `target_zero_reduction` **3724.7**;
  - `target_x7_reduction` **3571.7**;
  - `combined_certificate` **1300.5**;
  - `rref_mod_p` total across targets **5631.8** — **~77%** of the two
    reductions; certificate work total **~2877** (1576.2 per-target
    + 1300.5 combined) — **~40%** of wall-clock;
  - everything else is noise: `ranking_once` 19.3 (cached on the second
    target: 0.1), `assemble_rows_mod_p` 65.4, `reconstruction` 0.4.
- Conclusion / next step: **Perf.5 should target multi-target /
  linear-LHS normal-form reuse** (the RREF/ranking work is recomputed
  per target and per certificate point on the same row system;
  `rank_labels(..., target=...)` is target-dependent, so reuse needs a
  real design, not a one-line change) — **not multiprocessing**
  (Perf.3 verdict stands; per-point tasks are too small on Windows
  spawn, and here the cost is a few huge RREFs, not many small ones).

## Perf.5 — DONE: multi-target / linear-LHS normal-form reuse

- Change: `collect_normal_form_records_multi` / `reduce_rows_multi` —
  ONE shared pipeline (ranking, `assemble_rows_mod_p`, `rref_mod_p`,
  record collection, certificate points) over the SAME row system for
  several targets at once; per target only: record selection,
  reconstruction, Success gate, certificate verdicts.
  `scripts/run_example4_star_corrected.py` now reduces both LHS targets
  via the shared multi pass (`multi_target_reduction` stage replaces the
  two per-target stages; `sharing` note exported in `perf4_timings`).
- Correctness: `tests/test_perf5_multi_target.py` — 15 passed
  (multi vs serial equality incl. statuses, records, coefficients,
  plural certificate relations; singular/plural parity). Full suite +
  ruff clean after the edits. Heavy-run results **identical** to the
  certified baseline: same 2 combined terms
  (`(47703*ep^3-521*ep^2-57*ep-1)/(3300*ep^2)` on `{1,1,0,-1,0,0,0}`,
  `(816*ep^3+881*ep^2+66*ep+1)/(3300*ep^2)` on `{1,1,0,0,0,-1,0}`),
  `all_locally_finite=True`, rank 9924 for both targets, combined
  certificate **Passed 5/5** (filtered 0, exceeded 0, bad 0).
- Measured (same heavy corrected Example 4* config as Perf.4,
  12360 rows x 972 labels, editable install): **wall ~1h22m**
  (07:10 first log → 08:32 done) vs Perf.4 **~2h24m** — **~1.75x**.
  Stage profile (s): `multi_target_reduction` **3545.6** (was
  3724.7 + 3571.7 = 7296.4 across two targets), of which
  `rref_mod_p` **2715.1** ONCE (was 5631.8 = ~77% of two reductions),
  `records_total` 2768.3, `ranking_once` 19.1,
  `assemble_rows_mod_p` 32.5, `certificate_total` 776.9;
  `combined_certificate` 1293.3 (unchanged, ~= Perf.4's 1300.5);
  `row_generation_shared` 56.4, `rows_generated_once=True`.
  Caveat: the stage-timings snapshot is shared across targets, so
  per-target timing attribution is intentionally NOT claimed.
- Remaining hotspots (post-Perf.5): the single big `rref_mod_p`
  (2715 s) and certificate work (776.9 shared-pass + 1293.3 combined
  ≈ 2070 s). Further wins would need a faster mod-p RREF kernel
  (bit-packing / numpy) or certificate-point reuse — new designs,
  not reshuffling.

## Perf.6 — DONE: certificate-point RREF reuse (combined certificate)

- Change: the combined certificate reuses RREFs already computed for
  the multi-pass certificate points instead of redoing them — an RREF
  cache keyed by certificate point is filled during the shared
  multi-pass certificate and handed to the combined-certificate step
  (`_run_certificate_step(..., rref_cache)` in `reducer.py`; cache-aware
  `verify_reduction_relation_mod_p` path in `certificate.py`). Points
  not present in the cache fall back to the exact previous behavior.
- Heavy-run log line confirms the reuse:
  `combined certificate: 5 points, reusing 3 RREF(s) from the
  multi-pass certificate (Perf.6)`.
- Correctness: certificate-gate + perf suites re-run
  (`tests/test_certificate_gate.py`, `tests/test_perf5_multi_target.py`,
  `tests/test_perf3_jobs_equality.py`) and the full suite + ruff —
  all clean. Heavy-run results **identical** to the certified
  baseline: same 2 combined terms
  (`(47703*ep^3-521*ep^2-57*ep-1)/(3300*ep^2)` on `{1,1,0,-1,0,0,0}`,
  `(816*ep^3+881*ep^2+66*ep+1)/(3300*ep^2)` on `{1,1,0,0,0,-1,0}`),
  `all_locally_finite=True`, rank 9924 for both targets, combined
  certificate **Passed 5/5** (filtered 0, exceeded 0, bad 0;
  rank histogram `{9924: 5}`, `first_nonzero_residual=null`).
- Measured (same heavy corrected Example 4* config, 12360 rows x
  972 labels, editable install): **wall ~1h15m** (~07:39 launch →
  08:54 done) vs Perf.5 ~1h22m. Stage profile (s):
  `combined_certificate` **518.7** (was 1293.3, **~2.5x** on the
  targeted stage); `multi_target_reduction` 3798.5 (was 3545.6 —
  run-to-run noise on the big RREF, of which `rref_mod_p` 2949.4,
  `records_total` 3004.6, `certificate_total` 793.6, `ranking_once`
  19.8, `assemble_rows_mod_p` 33.7); `row_generation_shared` 58.7,
  `rows_generated_once=True`. Net targeted-stage sum
  (`multi_target_reduction` + `combined_certificate`): 4317 s vs
  Perf.5's 4839 s.
- Remaining hotspots (post-Perf.6): the single shared `rref_mod_p`
  (~2900 s, now ~2/3 of wall) and the certificate points that still
  need fresh RREFs. Further wins need a faster mod-p RREF kernel
  (bit-packing / numpy) — a new design, not reshuffling.

## Perf.7 — DONE (branch `perf/rref-backend-prototype`): RREF backend A

- See `docs/RREF_BACKEND_PLAN.md` for the full plan, candidate table,
  and measured tables. Summary: `rref_mod_p` gained optional
  `collect_stats=True` counters (nnz before/after, row-nnz profiles,
  fill-in ratio, inversions, elimination time; JSON-safe) and an
  opt-in `rref_backend="int_sparse_experimental"` (label→int column
  relabeling, identical `_eliminate` loop; default stays `"dict"`).
- Equivalence: `tests/test_rref_backend.py` (15 tests) — results
  identical across backends on every shape.
- Synthetic bench (`scripts/bench_rref_backends.py`): 0.90–0.92x at
  scale — verified but only ~8–10%; NOT plumbed upward.

## Perf.8 — DONE (same branch): real-matrix stats; candidate B REJECTED

- Real profile (`scripts/profile_rref_real_matrix.py` →
  `validation/rref_real_matrix_profile.json`): corrected Example 4*
  medium subset, 512×917 rank 512, nnz 2640→27762 (**fill-in 10.5x**,
  final row-nnz median 62 / max 114 of 917 cols — still sparse);
  elimination dict 0.81 s vs int backend 0.56 s (**0.69x** — the int
  backend wins MORE on real label tuples than on synthetic).
- Candidate B (sorted int-array rows, pure Python) micro-experiment at
  the measured real nnz profile: merge-axpy is **1.5–1.9x SLOWER**
  than dict-axpy (worst at the final median nnz 62: 1.91x). Numbers
  and methodology recorded in `docs/RREF_BACKEND_PLAN.md`.
- Verdict: B/C/D rejected in pure Python; dict stays default;
  `int_sparse_experimental` stays opt-in (confirmed useful on real
  matrices). Remaining RREF headroom requires leaving Python bytecode
  (Numba/native) — parked under the project's pure-Python constraint.
- Health: full suite green, `ruff check .` clean, profile script
  `ruff format --check` clean.

## Perf.7 / Perf.8 — CLOSED (2026-07-11)

- Perf.7 (RREF profiling + collect_stats infra) and Perf.8 (backend
  prototype: int_sparse_experimental opt-in, Candidate B rejected) are
  **closed**. Merged to `main` via no-ff merge `3a70cef`; HANDOFF note
  `185f2ea`. No tag, no release; version unchanged.
- Branch `perf/rref-backend-prototype` **deleted** (local tip was
  `e93a40f`; no remote copy existed). Only `main` remains.
- Default RREF backend: `"dict"`. Further RREF work = design only
  (`docs/RREF_BACKEND_PLAN.md`).

## Perf.10 — Numba int-array RREF backend prototype (branch `perf/numba-rref-backend`)

- New module `src/parametric_ibp_lf_reducer/sparse_rref_numba.py`:
  opt-in backend `numba_int_array_experimental`, lazily imported by
  `sparse_rref.rref_mod_p` only when explicitly requested — importing
  the package never touches numba; if numba is absent the backend
  reports unavailable with a clear error. **Not** the default; LF and
  certificate gates untouched; mathematics byte-for-byte the dict
  pivot algorithm (same pivot choice, same elimination order, one
  inversion per pivot, `p < 2**31` guard for int64 product safety).
- Tests: `tests/test_rref_numba_backend.py` (skips cleanly without
  numba) — 7-seed random-label equivalence, dense high-fill-in,
  rank-deficient/dup/zero rows, partial column order, small primes +
  `2**31 - 1` boundary, `>= 2**31` rejection, stats parity + plain-int
  results, and end-to-end `modular_normal_form` parity via backend
  swap. Full suite green; `ruff check .` + `ruff format` clean.
- Measured (editable install, this machine):
  - Synthetic bench `scripts/bench_rref_backends.py --fast`:
    tiny 200x150 ~6/row: dict 0.28s → numba 0.03s (**10.07x**);
    medium 1000x800 ~8/row: dict 51.57s → numba 2.04s (**25.34x**)
    (int_sparse_experimental: 1.12–1.15x).
  - Real ranking matrix `scripts/profile_rref_real_matrix.py`
    (512x917, rank 512, fill-in 10.5x): dict 0.828s →
    int_sparse 0.503s → **numba 0.083s** (~10x vs dict, ~6x vs
    int_sparse). JSON refreshed in
    `validation/rref_real_matrix_profile.json`.
  - One-time JIT compile on first use in a fresh env (a few seconds);
    `cache=True` persists kernels to `__pycache__`.
- Caveat: numba becomes an *optional* extra only; identical-results
  equivalence is enforced by tests, not assumed.
- Merged `--no-ff` into `main` (`81ec173`) and pushed; branch deleted.

## Perf.11 — RREF backend selection surfaced (ReducerConfig/API/CLI, branch `perf/rref-backend-plumbing`)

- `ReducerConfig` gains `rref_backend` (default `"dict"`); the value is
  threaded through `records.py` (serial path *and* ProcessPool workers
  via the per-point context), `modular_normal_form.py`, `reducer.py`
  and `certificate.py`, so ranking, normal-form and certificate RREF
  calls all honour the selection. Selection only — pivot choice,
  elimination order and all LF/certificate gates are untouched;
  results are backend-identical by construction (and enforced by the
  Perf.10 parity suite).
- Validation up front: backend names are checked against
  `sparse_rref.RREF_BACKENDS` (`"dict"`, `"int_sparse_experimental"`,
  numba opt-in) before any work starts; unknown names fail fast with a
  clear error.
- API: the public entry points accept `rref_backend` and pass it into
  `ReducerConfig` unchanged.
- CLI: new `--rref-backend` flag (choices = `RREF_BACKENDS`, default
  `dict`); an unknown value is a usage error (`EXIT_USAGE`).
- Tests: `tests/test_cli.py` extended — accepted value runs end-to-end
  and exits 0, unknown value is rejected with a usage error;
  `tests/test_rref_numba_backend.py` unchanged and green (skips
  cleanly without numba). Full suite + `ruff check`/`ruff format`
  clean.
