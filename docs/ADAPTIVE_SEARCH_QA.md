# Adaptive search QA record (v0.2.0)

Status: **complete** — Adaptive.1 (feature), Adaptive.1a (hardening),
Adaptive.2 (real-family validation). No adaptive-policy changes were needed
after validation.

## Scope

Opt-in bounded adaptive search (`--adaptive` / `--adaptive-max-levels`,
`reduce_family_adaptive`, `api.reduce_wolfram_style_input_adaptive`) layered
over the existing certified fixed-pass reducer. Out of scope for v0.2.0: any
new math, any performance/kernel work, any change to the default
(non-adaptive) path.

## Design invariants (verified)

1. **Opt-in only.** Without `--adaptive` the CLI/API path is the previous
   single fixed pass, byte-for-byte (fixed-path regression in
   `tests/test_cli.py`; flag gating asserted in the adaptive suites —
   `--adaptive-max-levels` without `--adaptive` is a usage error).
2. **No new reduction path.** Every level calls the existing fixed certified
   reducer (`reduce_family_once`); `Success` still requires reconstruction
   verification, the independent row-span certificate **and** all RHS terms
   `LF=True`.
3. **Determinism.** `default_search_levels` builds a fixed 3-level schedule
   (base → `expand-1` → `expand-2`: label-box m-deepening, IBP degree, tangent
   blocks, extra samples/primes); level configs and the best-partial
   tie-break are deterministic; reports are reproducible.
4. **Honest resource semantics.** `max_labels` (pre-flight skip, plus the
   distinct build-time `ValueError` guard for `expand_n`), `max_rows`
   (post-level) and `timeout_sec` (between atomic levels) surface as typed
   `ResourceLimitReached` data — never as fabricated success; no limit
   hard-preempts a running level.
5. **Layered gates.** A certificate `Passed` does not override a failed LF
   gate: level 0 of the real-family run is an honest
   `NormalFormNotLocallyFinite` *with* a passed certificate.
6. **Bounded, not a prover.** Exhausting the schedule returns the
   deterministically best partial failure with per-level history and
   failure-specific recommendations under `diagnostics.extra["adaptive"]`;
   it proves nothing about non-reducibility.

## Real-family validation (Adaptive.2)

- Family: real Example 2 five-term explicit family (`I3exampl2`), deliberately
  shallow base label box.
- Level 0: honest `NormalFormNotLocallyFinite`, certificate `Passed`,
  recommendation "expand the label box".
- Level 1 (`expand-1`, 72 labels / 1116 rows): certified `Success`; the
  notebook basis and coefficients are reproduced exactly
  (e.g. `C[(0,1,0,0,-2,-2,0)] = -2 + 2/ep^2`).
- Full transcripts: `docs/ADAPTIVE_SEARCH.md` / `docs/ADAPTIVE_SEARCH.ru.md`.

## Test inventory

- `tests/test_adaptive_search.py` — unit + loop-policy suite: schedule
  construction, `expand_n` guard, resource limits, diagnostics shape,
  best-partial selection, CLI flag gating.
- `tests/test_adaptive_real_family.py` — fast real-family API case in the
  normal suite (~25 s) + CLI e2e medium case gated behind
  `RUN_ADAPTIVE_MEDIUM=1` (config carried entirely via document `Options`:
  `LabelBox`, `PreferredMasters`, `Samples`, `RREFBackend`).
- Release gate: full fast suite + `ruff check .` green (see `CHANGELOG.md`
  v0.2.0 and `docs/FINAL_QA.md` process).

## Deliberately NOT rerun for v0.2.0

- Heavy D4 acceptance (`RUN_D4_FULL=1`) and the corrected Example 4\* heavy
  benchmark: no math or kernel changes in this release, so the certified
  v0.1.4 baselines remain authoritative.

## Verdict

Adaptive search ships in v0.2.0 as a controlled, opt-in, fully diagnosed
escalation layer over certified fixed passes; fixed explicit configurations
remain the recommendation for reproducible research runs.
