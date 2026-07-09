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
