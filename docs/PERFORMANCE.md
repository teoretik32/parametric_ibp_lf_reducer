# Performance notes (Perf.0–Perf.6, v0.1.3)

Optimization work targeted the corrected Example 4\* heavy run
(`scripts/run_example4_star_corrected.py`) without changing any math results:
coefficients, certificate verdicts, and diagnostics JSON schema are identical
before/after every step. LF/certificate gates were never weakened.

## Perf.0–Perf.6 summary

| Step | What | Outcome |
|---|---|---|
| Perf.0 | Baseline profiling of the heavy corrected Example 4\* run; per-stage wall-clock breakdown | `rref_mod_p` dominates; ranking and record orchestration secondary |
| Perf.1 | Lightweight stage-timing instrumentation: monotonic per-stage timers exported via `diagnostics.extra["timings"]` (CLI JSON included); pure observability, no control-flow influence | Stable timing keys (`row_generation_total`, `assemble_rows_mod_p`, `rref_mod_p`, `extract_normal_form`, `reconstruction`, `certificate_total`, `ranking_once`, …); timings never affect Success/Failure |
| Perf.2 | LF-valuation caching inside `ranking_once` (per-family tropical precomputation memoized; verdicts identical) | Ranking stage ~26% faster; records path untouched |
| Perf.3 | Parallel per-`(prime, sample)` normal-form record collection (`jobs=` / `--jobs`, `ProcessPoolExecutor`) | **Honest negative result on Windows/spawn** — kept experimental, `jobs=1` remains the default (see below) |
| Perf.4 | Ranking hoisted out of the per-record loop: ranking is a pure function of the label set, so it is computed once per family instead of per `(prime, sample)` point | Removes redundant re-ranking; bit-identical records/ranking output |
| Perf.5 | **Shared-RREF multi-target / linear-LHS normal-form reuse** (commit `e60763b`): when the LHS is a linear combination of targets over one shared row system, all per-target normal forms are extracted from a single RREF instead of re-running the pipeline per target | Corrected Example 4\* wall time ~2h24m → ~1h22m; `rref_mod_p` total ~5631.8s → ~2715.1s; results bit-identical |
| Perf.6 | **Certificate-point RREF reuse** (commit `88016a7`): the combined certificate reuses RREFs already computed for overlapping certificate points (same `(prime, sample)` row systems) instead of recomputing them per point | `combined_certificate` stage ~1293.3s → ~518.7s (~2.5x); wall ~1h22m → ~1h15m; certificate verdicts identical (**Passed 5/5**) |

## Why certificate-point RREF reuse works (Perf.6)

- Combined-certificate points overlap with RREFs already computed earlier in
  the run: the reduction and per-target certificate steps eliminate the same
  `(prime, sample)` row systems, so the elimination result can be cached and
  reused instead of recomputed per certificate point.
- It is a pure caching/orchestration change: every certificate point is still
  verified against the same exact-modular gate with identical verdicts
  (combined certificate **Passed 5/5**, same rank histogram); no gate was
  weakened and no math changed.

## Remaining hotspot / why further gains need a faster kernel (post-Perf.6)

- After Perf.6 the single large modular RREF (`rref_mod_p`, ~2900s in the
  measured run) is ~2/3 of wall time; orchestration-level reuse is exhausted —
  every remaining RREF is computed exactly once.
- Further wins therefore require a faster mod-p RREF **kernel** itself
  (bit-packing / numpy-based elimination), i.e. a new implementation, not
  reshuffling of when RREFs are computed.

## Why multiprocessing records was negative (Perf.3, Windows/spawn)

- Windows has no `fork`; `ProcessPoolExecutor` uses **spawn**: every worker
  pays fresh interpreter start + package/SymPy import before the first task.
- Per-record tasks are small; pickling inputs/results across process
  boundaries added overhead comparable to the work itself.
- The measured heavy-run comparison showed no win (and regressions at higher
  worker counts), so `jobs=1` stayed the default and the flag is documented
  as experimental. Equality tests guarantee parallel output, when used,
  matches the serial path exactly.

## Why shared-RREF linear-LHS reuse is the current recommendation (Perf.5)

- The corrected Example 4\* LHS is linear in the targets
  (`15*ep*J[...] + 24*ep*J[...]`) over the **same** row system, so the
  expensive elimination work is target-independent: one RREF serves all
  targets, and per-target normal forms are read off the shared elimination.
- It is a pure orchestration change: no new math, no new approximations,
  no gate changes — certificates re-verify each target as before
  (combined certificate **Passed 5/5**).
- It attacks the dominant cost (`rref_mod_p`) directly, unlike
  parallelization, which only wrapped the same work in process overhead.

## Corrected Example 4\* fingerprint (must reproduce)

- 972 labels, 12360 rows
- `selected_rank = 9924`
- combined certificate: **Passed 5/5**
- Status `Success`, `AllLocallyFinite` True
- final two coefficients unchanged (bit-identical to the v0.1.1 certified
  baseline in `validation/example4_star_corrected_result.m`)

## Remaining hotspots (post-Perf.5)

- One single large modular RREF (~2715s) and the certificate RREFs
  (~2070s combined) dominate; further wins require a faster RREF kernel
  (bit-packing / numpy-based elimination) or certificate-point reuse,
  not further orchestration reshuffling.
