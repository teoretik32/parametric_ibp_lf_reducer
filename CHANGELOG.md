# Changelog

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
