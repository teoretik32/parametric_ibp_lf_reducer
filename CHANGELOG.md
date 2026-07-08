# Changelog

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
