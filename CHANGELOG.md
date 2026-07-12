# Changelog

## Unreleased

### Added
- **Opt-in adaptive search (Pass Adaptive.1)**: `reduce_family_adaptive` /
  `reduce_wolfram_style_input_adaptive` / `AdaptiveSearchConfig` / `SearchLevel`
  / `default_search_levels` plus CLI flags `--adaptive` and
  `--adaptive-max-levels`. Runs a deterministic escalation schedule of ordinary
  fixed passes (expand label box m-ranges / IBP degree / tangent blocks /
  samples / primes), stops at the first *certified* `Success`, otherwise
  returns the deterministically best partial failure with a full per-level
  history and failure-specific recommendations under
  `diagnostics.extra["adaptive"]`. Resource limits (`max_labels` pre-flight,
  `max_rows` post-level, `timeout_sec` between levels) surface as typed
  `ResourceLimitReached` data — never as fabricated success. Docs:
  `docs/ADAPTIVE_SEARCH.md` / `docs/ADAPTIVE_SEARCH.ru.md`.

### Unchanged
- Without `--adaptive` the CLI/API path is the previous single fixed pass;
  every adaptive level goes through the same certificate gate, reconstruction
  verification and `AllLocallyFinite` check; exhausting the schedule proves
  nothing about non-reducibility (bounded schedule, not a prover).

## v0.1.4 — 2026

### Added
- **Optional `numba_int_array_experimental` RREF backend** (Perf.7–Perf.10):
  int64-array mod-p elimination kernel behind the same pivot/verdict contract
  as the `dict` reference backend.
- **`rref_backend="auto"` heuristic selection** (Perf.12): per-matrix choice of
  dict vs Numba using conservative size/prime gates.
- **Backend selection via `ReducerConfig` / Python API / CLI** (Perf.11); new
  `--rref-backend` CLI flag (`dict` / `numba_int_array_experimental` / `auto`).
- **Backend-selection diagnostics**: `requested_rref_backend`,
  `selected_rref_backend`, `backend_selection_reason`, `numba_available`,
  `auto_thresholds_used`.
- **No-Numba-safe lazy import and fallback**: `auto` silently falls back to
  `dict` when Numba is missing; an *explicit* Numba request fails fast with a
  clear error and is never substituted.

### Performance
- Corrected Example 4\* full-pipeline wall time (full box: 972 labels,
  12360 rows, selected rank 9924): **3963.4s (`dict`) → 803.8s (explicit
  Numba) → 766.5s (`auto` → Numba), ~5.17×**.
- `rref_mod_p` hotspot: **3124.1s → 689.2s → 656.1s (~4.76×)**.
- Certified full-box validation (Perf.13): 36/36 records valid, combined
  result `Success`, `AllLocallyFinite=True`, certificate `Passed`, same two
  certified coefficients — identical across all three backends.

### Correctness / unchanged
- **Exact equality** of mathematical outputs across `dict`, explicit Numba,
  and `auto`; certificate and LF gates unchanged.
- **Default backend remains `dict`** — Numba/auto are strictly opt-in
  (`pip install -e ".[speed]"`).
- Numba backend requires `prime < 2^31`; auto thresholds (unchanged this
  release): `min_rows=500`, `min_cols=400`, `min_nnz=3000`.

## v0.1.3 — 2026

### Performance
- **Certificate-point RREF reuse** (Perf.6, `88016a7`): reuse
  already-computed RREFs for overlapping certificate points instead of
  recomputing them in the combined certificate stage.
- Corrected Example 4\* combined certificate stage improved
  ~1293.3s → ~518.7s.
- Wall time improved ~1h22m → ~1h15m in the measured run.

### Unchanged
- **No mathematical result change**: Status `Success`, `AllLocallyFinite`
  True, combined certificate **Passed 5/5**, same two coefficients.

### Known hotspots
- Remaining cost is dominated by the single large modular RREF kernel
  (`rref_mod_p`); further wins need a faster mod-p RREF kernel
  (bit-packing / numpy-based elimination), not orchestration changes.

## v0.1.2 — 2026

### Performance
- **Shared-RREF reuse for linear-LHS / multi-target normal forms** (Perf.5,
  `e60763b`): when the LHS is a linear combination of targets over one shared
  row system, per-target normal forms are computed from a single RREF instead
  of re-running the full pipeline per target.
- Corrected Example 4\* runtime improved from ~2h24m to ~1h22m.
- RREF work reduced from ~5631.8s to ~2715.1s total (`rref_mod_p`) in the
  corrected Example 4\* profile.

### Tests
- New/updated tests for multi-target LHS equality
  (`tests/test_perf5_multi_target.py`) and the corrected Example 4\* path;
  full suite green (260 passed, 7 skipped), ruff clean.

### Unchanged
- **No math-result change**: coefficients and certificate remain unchanged
  (Status `Success`, `AllLocallyFinite` True, combined certificate
  **Passed 5/5**, `selected_rank=9924`).

### Known hotspots
- Remaining cost is dominated by one single large modular RREF (~2715s) and
  the certificate RREFs; further wins need a faster RREF kernel or
  certificate-point reuse, not orchestration changes.

## v0.1.1 — 2026 (candidate, not pushed)

### Added
- **Corrected Example 4\*** (exploratory, known-value-only): fixed integrand
  multiplier `15*ep + 24*ep*x7`, handled by linearity as
  `15*ep*J[{0,0,0,0,0,0,0}] + 24*ep*J[{0,1,0,0,0,0,0}]` with orchestration in
  `scripts/run_example4_star_corrected.py` (core `src/` unchanged; `lhs_terms`
  handling is fully generic). Certified `Success` (certificate `Passed`,
  `selected_rank=9924`); artifacts:
  `validation/example4_star_corrected_result.m`,
  `validation/example4_star_corrected_diagnostics.json`; tests:
  `tests/test_example4_star_corrected.py`.
- `docs/USAGE.md` / `docs/USAGE.ru.md`: "Corrected Example 4*" subsection.

### Unchanged
- Certified baseline remains D4 only; Example 4* stays exploratory
  (no reference LF decomposition, no numeric cross-check without
  master-integral values).

## v0.1.0 — 2026 (initial release)

### Added
- **Certified D4 LF reduction**: full-config D4 run reduces the target to a
  locally finite combination with support
  `{(0,1,1,0,-3,-1,0), (0,1,1,0,-2,-1,0), (1,1,0,0,-2,-1,0)}` ({M1, M2, M3});
  deterministic fingerprint: `n_rows=2092`, `n_records=108`, `n_selected=102`,
  `rank_histogram={1995: 6, 2041: 102}`.
- **CLI**: `python -m parametric_ibp_lf_reducer reduce <input.wl.txt>` with
  `--out`, `--diagnostics-json`, `--max-ibp-degree`, `--min-valid-records`;
  stable exit codes (`EXIT_SUCCESS`/`EXIT_FAILURE`/`EXIT_USAGE`).
- **Python API**: `reduce_wolfram_style_input` / `reduce_wolfram_style_input_to_text`
  returning typed `ReducerRunResult` with structured `ReductionDiagnostics`.
- **Certificate gate**: independent exact-modular verification of the reduction
  (rank-filter accounting, per-record pass/fail, `CERTIFICATE_PASSED` /
  `CERTIFICATE_FAILED` / `CERTIFICATE_INSUFFICIENT` / `CERTIFICATE_NOT_RUN`);
  D4 release certificate: **Passed 3/3** (rank-filtered 0, rank-exceeded 0, bad 0).
- Wolfram-like text parser/renderer (explicit-family requirement), examples
  (tiny success + D4 heavy), `docs/USAGE.md`, `scripts/final_check.sh|.ps1`.

### Limitations
- Wolfram-like text I/O only; no Mathematica/Wolfram runtime dependency, but
  also no other CAS formats.
- Input documents must carry an explicit parametric family
  (`FAILURE_PARSER_NEEDS_EXPLICIT_FAMILY` otherwise); no family inference.
- Single certified end-to-end configuration (D4); other families/dimensions
  run through the same pipeline but without a curated acceptance baseline.
- No adaptive search, no protected masters, no forced 5-term basis; failures
  are reported honestly via `ALL_FAILURE_REASONS` codes rather than retried.
- Heavy D4 acceptance is opt-in (`RUN_D4_FULL=1`, ~25–30 min).
