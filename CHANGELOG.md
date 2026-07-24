# Changelog

## Unreleased

### Added
- **External Int2 Method.11 Phase B: certified-solve runner (scaffold).**
  New runner `scripts/run_external_int2_method11.py`: assembles the full
  Level-0 system (baseline blocks + the `(3, 3)` extra tangent block) under
  `SurfacePolicy.chamber({"ep": -3/5})` through the first-class core API (no
  monkey-patching; the extra block gets the policy forwarded explicitly), then
  runs one fixed certified pass via `reduce_rows_once` with
  `require_certificate_for_success=True` (no opt-out flag by design). Heavy
  work is gated on `--allow-heavy` (otherwise an assembly-only artifact with
  verdict `NotRun(assembly-only)`); a comparable Method.10 artifact with a
  different chamber row total is an integrity stop
  (`Aborted(rows-mismatch)`, exit 2, bypass only via
  `--ignore-rows-mismatch`). Runner defaults: 3 packaged 31-bit primes, 12
  scattered samples (the `ep = 3` special sample stays in — the
  rank-consistency record filter is the authoritative discard),
  `--min-valid-records 6`, `--min-certificate-points 2`, `--jobs`,
  `--rref-backend`. Artifact `validation/external_int2_method11.json` (scope
  note, policy diagnostics, per-kind row counts, crosscheck, full result
  payload with per-term LF flags); Wolfram-style reduction text exported only
  on certified Success. Tests `tests/test_external_int2_method11.py` (20).
  The heavy certified run has NOT been executed yet — no Method.11 claims.
- **External Int2 Method.11 Phase A: explicit `SurfacePolicy` API (core).**
  New immutable `SurfacePolicy` (`surface.py`, exported at package level) with
  constructors `SurfacePolicy.limit(direction)` (the historical
  regulated-limit reading; package default `minus`) and
  `SurfacePolicy.chamber({"ep": Fraction(-3, 5)})` (exact rational chamber
  point; floats rejected at construction). Chamber signs are decided by exact
  rational substitution only — strictly positive/negative decides, an exact
  zero stays `zero`, leftover symbols give `unknown` (both reject the row).
  Threaded as an optional `policy` argument through
  `coordinate_primitive_surface_free`, `vector_field_surface_free`,
  `generate_coordinate_ibp_rows`, `generate_tangent_ibp_rows`, and as
  `ReducerConfig.surface_policy` into `reducer._generate_rows`; row
  diagnostics now record the active policy (`surface_policy: {mode,
  direction|point}`). Default behaviour is unchanged (`None` ->
  `SurfacePolicy.limit(eps_direction)`); the analytic-continuation reading of
  both modes is documented in the class docstring. Tests
  `tests/test_surface_policy.py` (9).
- **External Int2 Method.10: chamber-policy controlled LF-feasibility solve
  (full Level-0 box).** New runner `scripts/run_external_int2_method10.py`
  (nothing heavy without `--allow-heavy`; `--max-solves`,
  `--primes`) — re-assembles the FULL Level-0 system (3072 labels,
  chamber point `ep = -3/5`, target `(0,0,0,0,0,0,0)`) under the exact
  convergence-chamber sign policy (same diagnostic `surface.regulated_sign`
  swap as Method.9, script-local only), screens the chamber-only rows
  against the recorded Method.7 witnesses and, gated on breaks, runs
  witness-mode RREF solves with the span test constrained to LF-True
  labels. Tests `tests/test_external_int2_method10.py` (15 tests, incl.
  recorded-artifact schema/gates and the literal "NOT an LF basis" scope
  disclaimer), artifact `validation/external_int2_method10.json`. Outcome
  (exit 0, 4388.1s): limit-policy rows 46737 (392.8s), chamber-policy rows
  49439 (1411.1s); `chamber_only=2702` (`coordinate_ibp` 1024 +
  `tangent_ibp` 1678), `limit_only=0` — the chamber policy strictly
  extends production acceptance on the full box. Screening: 1685 breaking
  chamber-only rows per witness (identical dedup sets; 10110 total over 6
  witnesses), `rerun_justified=True`. Solve: 6/6 `Feasible` at rank 25234
  (3 generic samples x 2 primes `2147483647`/`2147483629`; nrows=49439,
  n_projected_rows=45443, n_allowed=1754, n_forbidden=28224, empty
  residual support) — verdict **`Feasible(modular)`**, vs `Obstructed`
  under the production policy (Method.6/7). Scope: per-(sample, prime)
  span evidence only, NOT an LF basis — reconstruction, row-span
  certificate and per-term LF export checks (gated by the analytic Laurent
  oracle) are still required before any Success claim. Production sign
  policy unchanged; recorded Method.6–9 artifacts untouched.
- **External Int2: full analytic Laurent oracle through `ep^0`
  (independent of the reducer).** New standalone script
  `scripts/audit_external_int2_full_laurent.py` — implements the exact `x7`
  preintegration, the one-dimensional hypergeometric ODE, the
  epsilon-recursion for `S=rQ` and the compact Laurent coefficients derived
  in `notes/EXTERNAL_INT2_FULL_LAURENT_DERIVATION.md`. Three symbolic checks
  (`exact_x7_preintegration`, `ode_epsilon_recurrence_s0_s3`,
  `full_laurent_coefficients_through_ep0`) — all pass; optional `--numeric`
  HPL fingerprints (mpmath, `--dps >= 20`; recorded at dps=40 for
  `(s,t)=(0.75,1)` and `(1.7,0.8)`). Artifact
  `validation/external_int2_full_laurent_audit.json`, WL cross-check
  `validation/external_int2_full_laurent_result.m`, tests
  `tests/test_external_int2_full_laurent.py` (6 tests, incl. a guard that
  the script does not import reducer core). Source notebook preserved
  verbatim at `examples/source/ParametricInt_examples_4_ChatGPT_v2.nb`;
  `examples/external_int2_source_reference.wl.txt` rewritten as
  pure-integrand + prefactor reference (`ExternalInt2PureIntegrand`,
  `ExternalPrefactor2`, `r = s/t`). Transplant bundle
  `external_int2_full_analytic_patch/` (byte-identical copies of all new
  files) with RU prompts `PROMPT_1_INTEGRATE_ANALYTIC_ORACLE_RU.md`,
  `PROMPT_2_FIND_TRUE_LF_BASIS_RU.md` and `README_PATCH_RU.md`. Scope: the
  oracle does NOT prove or construct an LF basis — it is the acceptance
  gate for future LF-basis work (any certified LF decomposition of External
  Int2 must reproduce these coefficients). Reducer core, certificates and
  LF gates untouched.
- **External Int2 Method.9: surface-policy audit (limit vs chamber sign
  policy) + tangent-module Singular export (Part C).** New runner
  `scripts/run_external_int2_method9.py` — pushes TWO sign policies through
  the IDENTICAL library filter code (production `regulated_sign` in the
  limit `ep -> 0^-` vs exact rational signs at the convergence-chamber
  point `ep = -3/5`) and diffs accept/reject decisions; pairing-only
  against the recorded Method.7 witnesses, NO RREF, read-only. Tests
  `tests/test_external_int2_method9.py` (10 tests incl. limit-mode
  equivalence with the library generator accept-sets and independent SymPy
  tangency residuals), artifact `validation/external_int2_method9.json`.
  Outcome (192 seeds, 586 checks, 16.1s): `agree_keep=41`,
  `agree_reject=506`, `chamber_only=39`, `limit_only=0` — the chamber
  policy strictly extends acceptance; all 39 chamber-only rows are
  `tangent_ibp` and break the recorded witnesses (`n_breaks_total=180`
  across 6 witnesses). Confirmed minimal flip: seed `(-1,0,1,0,-1,-1,0)`,
  tangent field 2, score `-3*ep-1` (limit `neg` -> rejected; chamber `pos`
  -> accepted). Part C: `scripts/export_external_int2_tangent_module.py`
  renders `scripts/external_int2_tangent_module.sing` (full tangent module
  via `syz`; Singular not installed — prepared for offline execution,
  regeneration-consistency tested). Production sign policy unchanged;
  recorded Method.6/7/8 artifacts untouched; per-(seed, field), per-box
  statements only — no cure claim.
- **External Int2 Method.8: targeted Level-0 re-elimination with the
  `ibp_deg3` rows (the one Method.7-breaking family).** New runner
  `scripts/run_external_int2_method8.py` (nothing runs without `--phase`;
  `--phase rerun` is HEAVY — witness-mode RREFs on the enlarged
  merged+`ibp_deg3` Level-0 system at 3 generic samples x 2 primes — and is
  gated behind `--allow-heavy`), tests `tests/test_external_int2_method8.py`
  (tiny box; integration gated by `RUN_EXTERNAL_INT2_M8=1`), artifact
  `validation/external_int2_method8_reelim.json`. Per point the runner
  reports `Feasible` (obstruction CURED at that point) or `Witness` (still
  obstructed; the NEW dual witness must annihilate every included `ibp_deg3`
  row, plus rank/nullity deltas vs the recorded Method.7 baseline).
  Recorded Method.6/7 artifacts are never overwritten; reducer core,
  certificates and LF gates untouched; no Level 1/2 rerun; per-(sample,
  prime), per-box statements only. Outcome: still `Obstructed` at all 6
  generic points (3 samples x 2 primes, 62241 rows) —
  rank=37478 / nullity=6779 / support=3796 identically (Δ vs the Method.7
  baseline: +12861 / +3350 / +1339); the new dual witness annihilates every
  included `ibp_deg3` row (0 breaks). The one breaking family is confirmed
  necessary-but-not-sufficient at this box.
- **External Int2 Method Audit.1 (Phase A): all-row-support LF feasibility —
  seed-box defect candidate refuted.** Domain-feedback hypothesis: the Level 0
  `Obstructed` could be an artifact of restricting the allowed set to LF-True
  labels inside the seed box. New opt-in all-row-support mode in
  `src/parametric_ibp_lf_reducer/lf_feasibility.py` (production seed-box
  default unchanged), runner `scripts/run_external_int2_audit1.py` (nothing
  runs without `--phase`; HEAVY behind `--allow-heavy`), tests
  `tests/test_lf_feasibility_all_support.py`, artifact
  `validation/external_int2_audit1_allsupport.json`. At the generic point
  `ep=15/7, r=32/11`, p=2147483647 the allowed set widens 1754 -> 13531
  (11777 newly allowed out-of-box LF-True labels), yet both modes stay
  `Obstructed` with the same `residual_support` = the target unit label:
  the target unit vector is not in the projected row span either way.
  `defect_confirmed=false`; read-only diagnostics, no new rows, no
  re-elimination; per-(sample, prime), per-box statement only.
- **External Int2 Method.7: dual-witness stability + pairing-only candidate
  screening at the medium T2 box (Level 0).** New runner
  `scripts/run_external_int2_method7.py` (nothing runs without `--phase`;
  `--phase witness` is HEAVY and gated behind `--allow-heavy`; `--phase screen`
  is cheap, pairing-only, and never re-eliminates), tests
  `tests/test_external_int2_method7.py`, recorded artifacts
  `validation/external_int2_method7_witness.json` /
  `validation/external_int2_method7_screen.json`. Phase A: all 6 points
  (3 generic samples x primes 2147483647/2147483629) are `Witness` with
  identical rank 24617 / nullity 3429 / support 2457 and a single shared
  support pattern (pairwise Jaccard 1.0), consistent with the recorded
  `external_int2_t2_rankrepair_level0.json`. Phase B/C: every existing row
  family annihilates all 6 witnesses (0 breaks); candidate `ibp_deg3` rows
  break the witness (6972 of 15504 new rows) while `ray_multipliers` deg 4-6
  and `tangent_(5,5)` are pairing-inert — `rerun_justified=true` for
  `candidate_ibp_deg3` only (a breaking row is necessary, NOT sufficient, to
  cure the obstruction; no re-elimination is run). Read-only diagnostics:
  reducer core, certificates and LF gates untouched; per-(sample, prime),
  per-label-box statements only, no global impossibility claim.
- **External Int2 Method.6: reproducibility cleanup + dual LF-obstruction
  certificate.** New library module
  `src/parametric_ibp_lf_reducer/lf_obstruction_witness.py` (exported from
  `__init__`): for an `Obstructed` LF-span system it builds an explicit dual
  witness `w` in the RIGHT nullspace of the projected matrix
  (`<row, w> == 0` for every projected row, `w[target] == 1`), with exact
  per-point checks (`check_annihilation`, `check_target_unit`), deterministic
  construction, JSON payload serializers, and a row-pairing helper
  (`test_rows_against_obstruction_witness`) that flags candidate rows which
  break vs annihilate a stored witness. New runner
  `scripts/run_external_int2_t2_rankrepair.py` reproduces the T2 rank-repair
  Levels 0-2 (`--describe`; nothing runs without `--levels`; heavy Levels 1/2
  gated behind `--allow-heavy`; `--witness` / `--probe-rows` Phase C modes),
  writing `*_repro.json` / `*_witness_level{N}.json` /
  `*_rowprobe_level{N}.json` artifacts that never overwrite the recorded files.
  Tests `tests/test_lf_obstruction_witness.py` and
  `tests/test_external_int2_t2_rankrepair.py`. Retro entries: Method.5
  (`validation/external_int2_method5.json`) and the T2 rank-repair artifacts
  (`validation/external_int2_t2_rankrepair_level{0,1,2}.json`) — both generic
  `Obstructed`. **No change to LF/certificate semantics; reducer core, gates
  and existing behavior are untouched.** No global impossibility claim is made.
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

### Docs
- **External Int2 audit — Method.5 / T2 rank-repair / Method.6 sections.**
  `notes/EXTERNAL_INT2_AUDIT.md` gains retro-docs for the Method.5 label-box
  geometry audit and the T2 rank-repair Levels 0-2, plus a Method.6 section with
  the dual-witness math, determinism rules and the four State bullets.
  `notes/HANDOFF.md` pass #40. **Codimension-one phrasing corrected**:
  `residual_support == [target]` does not imply the quotient dimension is one
  (quotient dimension = nullity = projected cols − rank, may exceed 1); only the
  prose is corrected — the recorded validation JSON `purpose` strings are
  historical and intentionally left byte-identical. Explicit: **no change to
  LF/certificate semantics.**

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
