# Adaptive search (opt-in, Pass Adaptive.1)

Adaptive search runs a **deterministic escalation schedule** of ordinary fixed reduction
passes and stops at the first *certified* `Success`. It is **controlled schedule expansion,
not a universal proof of reducibility**: exhausting the schedule only means "none of these
particular configurations certified a reduction" — it never proves that no reduction exists.

## Honesty contract

- Every level is one ordinary `reduce_family_once` pass. Adaptive mode adds **no new math and
  no new `Success` path**: reconstruction verification, the row-span certificate gate and the
  `AllLocallyFinite` check apply to every level exactly as in a single fixed pass.
- A passed certificate on a non-locally-finite or `Unknown` normal form is still a typed
  failure — adaptive mode never upgrades it.
- Adaptive mode is **opt-in**. Without `--adaptive` (CLI) or the adaptive API the behavior is
  the unchanged single fixed pass.
- The default RREF backend stays `dict`; a schedule level may explicitly select another
  backend (e.g. `auto`, recommended when the `speed` extra is installed).

## Python API

```python
from parametric_ibp_lf_reducer import (
    AdaptiveSearchConfig, SearchLevel,
    reduce_family_adaptive, reduce_wolfram_style_input_adaptive,
    default_search_levels,
)

result = reduce_wolfram_style_input_adaptive(input_text)          # default 3-level schedule
result = reduce_wolfram_style_input_adaptive(
    input_text,
    search=AdaptiveSearchConfig(
        levels=(
            SearchLevel(name="small", label_box=((0, 0), (-1, 0)), max_ibp_degree=1,
                        tangent_degree_blocks=()),
            SearchLevel(name="wide", label_box=((0, 0), (-2, 0)), max_ibp_degree=2,
                        tangent_degree_blocks=((1, 1),)),
        ),
        max_labels=20_000, max_rows=200_000, timeout_sec=600.0,
    ),
)
```

`reduce_family_adaptive(family, target_label, config, search=None)` is the underlying entry
point for already-parsed families. `SearchLevel` fields set to `None` inherit the base
config; `tangent_degree_blocks=()` means "explicitly no tangent rows at this level";
`extra_samples` / `extra_primes` deterministically extend (never replace) the base scattered
samples and prime list.

## Default schedule (`default_search_levels`)

Derived only from the base config (base box = the document/override `LabelBox`, or the
package default `((0, 0), (-1, 0))`):

| level | name | label box | `max_ibp_degree` | tangent blocks | extras |
|---|---|---|---|---|---|
| 0 | `base` | base box | 1 | none | — |
| 1 | `expand-1` | every m-range deepened by 1 | 2 | `((1, 1),)` | — |
| 2 | `deep` | every m-range deepened by 2 | 2 | `((1, 1), (2, 2))` | +4 samples, +2 primes; optional `rref_backend` (recommended `"auto"`) |

By default n-ranges never change. Opt-in `expand_n` (a per-n-axis 0/1 mask argument of
`default_search_levels`) widens masked n-axes symmetrically by the level delta (`(lo-k, hi+k)`
at level *k*). Because n-expansion multiplies the box volume, `expand_n` **requires**
`max_labels` — a build-time guard: every planned level's label count must stay within it
(`ValueError` otherwise). This guard is distinct from the runtime pre-flight
`AdaptiveSearchConfig.max_labels`, which *skips* oversized levels instead of refusing to build
the schedule.

An explicit `labels` list cannot be grown deterministically —
supply explicit `SearchLevel`s instead (`ValueError` otherwise).

## Result and history

The returned object is an ordinary `ReductionResult`:

- first certified `Success`, or
- the **best partial failure** when no level succeeds, ordered deterministically by:
  Success → target reduced → fewer non-LF terms → fewer `Unknown` terms → reconstruction
  verified → certificate quality (`Passed` > `Insufficient` > `NotRun` > `Failed`) → earlier
  level (tie-break).

The full history is attached at `result.diagnostics.extra["adaptive"]` (and in the CLI
`--diagnostics-json` payload under `"adaptive"`): per-level status, failure reason, a short
deterministic `error` detail (the attempt's diagnostic messages, truncated to 500 chars,
`None` on success — full failed results are deliberately not retained),
certificate status, row/label/record counters, a failure-specific **recommendation**
(e.g. `TargetNotReducible` → expand box/degree/tangent; `InterpolationFailed` → more
samples/primes; `VerificationFailed` → add independent points, **never** accept the current
coefficients), and observability-only wall-clock times.

## Resource limits

| limit | semantics |
|---|---|
| `max_levels` | truncates the schedule (default 3) |
| `max_labels` | **pre-flight**: an oversized level is *not* run; if nothing ran at all, the result is a typed `ResourceLimitReached` failure |
| `max_rows` | checked **after** a level completes; stops further escalation (the completed level still counts) |
| `timeout_sec` | checked **between** levels only; a running level is never aborted |

**No limit is hard-preemptive**: `max_labels` skips an oversized level *before* it starts,
`max_rows` is observed only *after* the offending level has already run to completion, and a
long level always runs to completion (levels are atomic).

Every limit hit is reported exactly (`kind`, `limit`, `observed`, `level`) in
`extra["adaptive"]["resource_limit"]`. `timeout_sec` is the only wall-clock knob, is disabled
by default, and never changes the mathematical content of any completed level.

## CLI

```
python -m parametric_ibp_lf_reducer reduce input.m --adaptive
python -m parametric_ibp_lf_reducer reduce input.wl.txt \
    --adaptive --adaptive-max-levels 3 --rref-backend auto \
    --out result.m --diagnostics-json diagnostics.json
```

`--adaptive-max-levels` requires `--adaptive`; without `--adaptive` the CLI path is byte-for-
byte the previous fixed single pass. Exit codes are unchanged: `0` only for a certified
`Success` (of some level), `1` for an honest failure, `2` for usage errors.

## Real-family validation (Adaptive.2)

The schedule is validated on a real explicit family — Example 2 (`I3exampl2`,
`Examples_for_IBP_parametric.nb`): `examples/notebook_example2_n3_five_term_explicit_family.wl.txt`,
whose known result is the 5-term basis in
`validation/notebook_example2_n3_five_term_expected.json`.

Starting from a deliberately shallow base box (`n = ((0,0),(0,1),(0,0))`,
`m = ((0,0),(-1,0),(-1,0),(0,0))`), 16 scattered samples, the default primes,
`PreferredMasters` set to the known basis and `rref_backend="auto"`, the **default schedule**
escalates and certifies with no hand-crafted levels (~25 s total, Windows, auto backend):

| level | config | size | outcome |
|---|---|---|---|
| 0 `base` | degree 1, no tangent rows | 8 labels / 64 rows | `NormalFormNotLocallyFinite`, certificate `Passed`, recommendation "expand the label box …" |
| 1 `expand-1` | m-ranges −1, degree 2, tangent `((1,1),)` | 72 labels / 1116 rows | **certified `Success`** |

Level 1 reproduces exactly the five expected masters with the notebook coefficients (e.g.
`C[(0,1,0,0,-2,-2,0)] = -2 + 2/ep^2`). Note the honest level-0 report: the certificate can
*pass* (the coefficients are consistent at an independent point) while the result still fails
the `AllLocallyFinite` gate — and the attached recommendation ("expand the label box") is
precisely what the next scheduled level does.

Tests (`tests/test_adaptive_real_family.py`):

* `test_default_schedule_certifies_real_five_term_family` — API, normal suite (~25 s);
* `test_cli_adaptive_certifies_real_family_medium` — the same run through the CLI, with the
  whole configuration carried in the document `Options` (`LabelBox`, `PreferredMasters`,
  `Samples`, `RREFBackend`); gated: set `RUN_ADAPTIVE_MEDIUM=1` to run.

```
python -m parametric_ibp_lf_reducer reduce ex2_adaptive.m --adaptive \
    --adaptive-max-levels 2 --diagnostics-json diag.json
```

Exhausting the schedule on a too-small box remains an honest failure: the full per-level
history and recommendations land in `diag.json["adaptive"]` (see the honesty contract above).
