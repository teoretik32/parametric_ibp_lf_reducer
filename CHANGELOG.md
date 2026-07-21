# Changelog

## Unreleased

### Added
- **External Int1 (standalone example): certified LF reduction.** Input
  `examples/external_int1_corrected_input.wl.txt`, runner
  `scripts/run_external_int1_corrected.py`, artifacts
  `validation/external_int1_corrected_reduction.m` /
  `validation/external_int1_corrected_full_formula.m` /
  `validation/external_int1_corrected_diagnostics.json`, tests
  `tests/test_external_int1_corrected.py`. Numeric original-vs-RHS check:
  rel_diff ≈ 1.38e-35. Reducer core unchanged.
- **External Int1 Laurent-structure audit — PASSED through `ep^0`.** Standalone
  high-precision script `scripts/audit_external_int1_laurent.py` (`mp.dps = 45`,
  Cauchy-circle Taylor extraction + PSLQ identification in a weight-graded
  basis, max PSLQ residual ≈ 2.6e-41); per-order match with the target
  `1/ep^4 − (π²/12)/ep² − (43·ζ3/6)/ep − π⁴/180`. Report:
  `notes/EXTERNAL_INT1_LAURENT_AUDIT.md`; machine-readable:
  `validation/external_int1_laurent_audit.json`. The audit is a high-precision
  numeric validation, **not** a formal symbolic proof. Reducer core unchanged.
- **External Int2 (dimensionless): certified partial reduction — stable
  negative LF result.** New example
  `examples/external_int2_dimensionless_input.wl.txt`
  (`ExternalInt2Dimensionless`, vars `x2,x5,x7`, polys `G0..G3`), gated runner
  `scripts/run_external_int2.py` (`RUN_EXTERNAL_INT2=1`), fast tests
  `tests/test_external_int2.py`. Heavy run #4 (`base+boost-s48-p6-x1`,
  60030.6 s): certificate **Passed 3/3** (rank filter 531/540 at rank 22361,
  histogram `{19722: 9, 22361: 531}`), reconstruction verified, independent
  validation passed. Outcome: `NormalFormNotLocallyFinite` — 6-term
  decomposition with two genuinely non-locally-finite residual terms (`1/G1`,
  `-(ep+1)/ep * 1/G2`), reproduced identically across boosted configurations
  and a deepened label box. Audit: `notes/EXTERNAL_INT2_AUDIT.md`; artifacts:
  `validation/external_int2_result.m` /
  `validation/external_int2_full_formula.m` /
  `validation/external_int2_diagnostics.json`. Reducer core unchanged.
- **Method.1 for External Int2: directional LF audit + LF-constrained
  feasibility mod p (diagnostic-only).** New module
  `src/parametric_ibp_lf_reducer/lf_feasibility.py` plus
  `explain_local_finiteness` in `valuations.py`, gated runner
  `scripts/run_external_int2_method1.py` (`RUN_EXTERNAL_INT2_M1=1`), tests
  `tests/test_lf_feasibility.py`, `tests/test_valuations_explain.py`,
  `tests/test_external_int2_method1.py`. Runs: Level A (base box, 648 labels,
  ~117 s) and Level B (expand-1, 2048 labels, ~2197 s); target LF verdict
  **False** at both levels (23/30 failing rays, 0 unknown). Phase B:
  Obstructed 0/6 (A) vs Mixed 2/6 (B) — both feasible points sit at the
  non-generic sample `ep=3` with a rank drop (18422 vs 20963), treated as a
  special-locus artifact; generic samples stay Obstructed at both primes.
  Artifacts: `validation/external_int2_method1_levelA.json` / `_levelB.json`;
  audit section in `notes/EXTERNAL_INT2_AUDIT.md`. Certificate and LF gates
  untouched; reducer core unchanged.
- **Method.2 for External Int2: wrapper-level leading-pole audit + prefactor
  fix (all 8 checks passed at full precision).** Standalone script
  `scripts/audit_external_int2_leading_pole.py` (no reducer import; heavy rerun
  opt-in via `RUN_INT2_POLE_AUDIT=1`). Exact `x7` preintegration
  `(B^ep - A^ep)/(ep*(B - A))` reduces Int2 to a 1-D form; crossover boundary
  poles cancel exactly (`K1 == C_B`), giving the pure leading pole
  `J2(ep,r) = -2/(3*r*ep^2) + O(1/ep)` (a naive `-1/(2*r*ep^2)` is excluded by
  the numeric Laurent fit). Corrected external prefactor:
  `EXTERNAL_PREFACTOR_TEXT` now carries `Exp[2*ep*EulerGamma]` (mirrored in
  `validation/external_int2_full_formula.m` /
  `validation/external_int2_diagnostics.json`; pinned by
  `test_prefactor_text_matches_p2`). Full object: with `r = s/t`,
  `P2*J2 = -4/(s*t^2*ep^4) + O(1/ep^3)`, matching the source `AnsvInt2` leading
  pole; `AnsvInt2` stays metadata only
  (`examples/external_int2_source_reference.wl.txt`), never a reducer
  coefficient. Report: `notes/EXTERNAL_INT2_LEADING_POLE_AUDIT.md`; JSON:
  `validation/external_int2_leading_pole_audit.json`; tests:
  `tests/test_external_int2_leading_pole.py`. Reducer core unchanged.
- **Method.3 for External Int2: composite locally-finite master feasibility.**
  New module `src/parametric_ibp_lf_reducer/composite_masters.py`, runner
  `scripts/run_external_int2_method3.py`, tests
  `tests/test_composite_masters.py` (heavy integration gated by
  `RUN_EXTERNAL_INT2=1`), artifact
  `validation/external_int2_composite_feasibility.json`. Outcome:
  **`FeasibleCompositeBasis`** — from a deterministic 225-candidate pool, the
  48-participant primary-ray cancellation kernel (dim 21) refines on 69
  checked rays to a 13-dimensional fully-LF composite basis; interpretable
  examples `J(1/(x2*G1)) - J(1/(G0*G1))` and `J((1+x5)/G1) - J((1+x7)/G2)`.
  Statements are scoped to this pool and ray set; `BadSpecialization` rank
  guard; reducer core, certificates and LF gates untouched.
- **Finite-numerator LF basis search for External Int2 (honest negative).**
  New module `src/parametric_ibp_lf_reducer/finite_numerator.py`
  (single-integrand semantics: a candidate is ONE decorated integrand
  `N(x)*F_S`, accepted only on a full `is_locally_finite = True` verdict;
  Lemma 1 graded-lowest-layer kernel cross-check, Lemma 2
  `numerator_cure_impossible_any_degree`), design doc
  `docs/FINITE_NUMERATOR_BASIS_DESIGN.md`, runner
  `scripts/run_external_int2_finite_numerator.py`, artifact
  `validation/external_int2_finite_numerator.json`, tests
  `tests/test_finite_numerator.py` (incl. an offset-convention regression and
  the defining-rows → `lf_reduction_feasible_mod_p` bridge). Verdict over the
  six certified normal-form sectors plus the probe `1/(G1*G3)` at degrees
  0–2: `NoFiniteNumeratorBasisWithinAnsatz` — the remnants `1/G1`, `1/G2`,
  `1/(G1*G3)` fail only on componentwise `<= 0` rays (`x -> oo`), which
  polynomial numerators can only worsen, so the cure is impossible at ANY
  degree; the other four sectors are `SectorAlreadyLF`; feasibility stage
  honestly `SkippedNoCandidates`. Reducer core unchanged.
- **Method.4 for External Int2: same-dimension LF-basis completeness audit
  (obstruction confirmed stable).** Gated runner
  `scripts/run_external_int2_method4.py` (`RUN_EXTERNAL_INT2=1`), artifact
  `validation/external_int2_method4.json` (+ probe JSONs), tests
  `tests/test_external_int2_method4.py`. Enriching the Method.1-style row
  system (level `deep`: 5000 labels, 77379 baseline rows) with richer
  tangent-IBP blocks `(3,3)`/`(4,4)` — 46 new vector fields, 39715 genuinely
  new rows (117094 total, rank 49559 → 54990) — flips **no** verdict at
  3 samples × 2 primes: generic points stay `Obstructed` ("target unit vector
  not in projected row span"), the special point `ep=3, r=54/11` stays
  Feasible; `flipped=0`. The Method.1 obstruction is not a row-basis
  truncation artifact; the viable route remains the Method.3 composite basis
  change. Read-only diagnostics; reducer core, certificates and LF gates
  unchanged. Elapsed ~19413 s (background).

## v0.2.0 — 2026

Release theme: **controlled adaptive search over certified fixed-pass
reductions** (Adaptive.1 / Adaptive.1a / Adaptive.2). No new math or
performance features; no adaptive-policy changes; heavy certified baselines
(D4, corrected Example 4\*) were deliberately not rerun.

### Added
- **Opt-in adaptive search (Pass Adaptive.1)**: `reduce_family_adaptive` /
  `reduce_wolfram_style_input_adaptive` / `AdaptiveSearchConfig` / `SearchLevel`
  / `AdaptiveLevelReport` / `AdaptiveSearchDiagnostics` /
  `default_search_levels` plus CLI flags `--adaptive` and
  `--adaptive-max-levels`. Runs a deterministic escalation schedule of ordinary
  fixed passes (expand label box m-ranges / IBP degree / tangent blocks /
  samples / primes), stops at the first *certified* `Success`, otherwise
  returns the deterministically best partial failure with a full per-level
  history and failure-specific recommendations under
  `diagnostics.extra["adaptive"]`. Resource limits (`max_labels` pre-flight,
  `max_rows` post-level, `timeout_sec` between levels) surface as typed
  `ResourceLimitReached` data — never as fabricated success. Docs:
  `docs/ADAPTIVE_SEARCH.md` / `docs/ADAPTIVE_SEARCH.ru.md`.
- **Adaptive.1a hardening**: opt-in `expand_n` mask for `default_search_levels`
  (masked n-axes widen symmetrically per level; **requires** a build-time
  `max_labels` guard — every planned level must fit, `ValueError` otherwise,
  distinct from the runtime pre-flight skip); per-level reports gain a bounded
  deterministic `error` detail (attempt's diagnostic messages, ≤500 chars,
  `None` on success; full failed results are deliberately not retained); docs
  spell out that no resource limit is hard-preemptive (levels are atomic).
### Real-family validation (Adaptive.2)
- The default schedule, started
  from a deliberately shallow base box on the real Example 2 five-term explicit
  family, escalates once (level 0 honest `NormalFormNotLocallyFinite` with a
  passed certificate + "expand the label box" recommendation → level 1 certified
  `Success`) and reproduces exactly the notebook basis and coefficients. Tests:
  `tests/test_adaptive_real_family.py` (fast API case in the normal suite;
  CLI e2e medium case gated behind `RUN_ADAPTIVE_MEDIUM=1`, config carried via
  document `Options`). Docs transcripts in `docs/ADAPTIVE_SEARCH.md` / `.ru.md`.
  No adaptive policy changes were needed.

### Correctness / unchanged
- Without `--adaptive` the CLI/API path is the previous single fixed pass,
  byte-for-byte; every adaptive level calls the existing fixed certified
  reducer (no new reduction path) and goes through the same certificate gate,
  reconstruction verification and `AllLocallyFinite` check.
- A certificate `Passed` never overrides a failed LF gate (layered gates — see
  the level-0 transcript in `docs/ADAPTIVE_SEARCH.md`); exhausting the
  schedule proves nothing about non-reducibility (bounded schedule, not a
  prover).

### Limitations
- `timeout_sec` is checked **between** atomic levels; `max_labels` is a
  pre-flight skip and `max_rows` a post-level limit — no resource limit
  hard-preempts a running level.
- Fixed explicit configurations remain the recommendation for reproducible
  research runs; adaptive search is an exploration tool.

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
