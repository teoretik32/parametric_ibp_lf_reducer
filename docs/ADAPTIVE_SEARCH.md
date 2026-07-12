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

n-ranges never change. An explicit `labels` list cannot be grown deterministically —
supply explicit `SearchLevel`s instead (`ValueError` otherwise).

## Result and history

The returned object is an ordinary `ReductionResult`:

- first certified `Success`, or
- the **best partial failure** when no level succeeds, ordered deterministically by:
  Success → target reduced → fewer non-LF terms → fewer `Unknown` terms → reconstruction
  verified → certificate quality (`Passed` > `Insufficient` > `NotRun` > `Failed`) → earlier
  level (tie-break).

The full history is attached at `result.diagnostics.extra["adaptive"]` (and in the CLI
`--diagnostics-json` payload under `"adaptive"`): per-level status, failure reason,
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

Every limit hit is reported exactly (`kind`, `limit`, `observed`, `level`) in
`extra["adaptive"]["resource_limit"]`. `timeout_sec` is the only wall-clock knob, is disabled
by default, and never changes the mathematical content of any completed level.

## CLI

```
python -m parametric_ibp_lf_reducer reduce input.m --adaptive
python -m parametric_ibp_lf_reducer reduce input.m --adaptive --adaptive-max-levels 2 \
    --diagnostics-json diag.json
```

`--adaptive-max-levels` requires `--adaptive`; without `--adaptive` the CLI path is byte-for-
byte the previous fixed single pass. Exit codes are unchanged: `0` only for a certified
`Success` (of some level), `1` for an honest failure, `2` for usage errors.
