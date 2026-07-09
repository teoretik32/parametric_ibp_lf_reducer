# Perf timing snapshots (Perf.0.1 baseline)

Wall-clock stage timings from `diagnostics.extra["timings"]` (exported in CLI
`--diagnostics-json` under `diagnostics.timings`). Machine: Teoretik win32 box,
CPython 3.13, pure-Python path (no Wolfram runtime).

> Footgun: the venv has the pip-installed `parametric-ibp-lf-reducer v0.1.0`
> (pre-Perf.0). `python -m parametric_ibp_lf_reducer` resolves to site-packages,
> NOT `src/`. For perf runs always use `PYTHONPATH=src` (and note the src CLI
> requires the `reduce` subcommand; the installed v0.1.0 accepted a bare input path).

## Snapshot 1: tiny example (examples/tiny_success_input.wl.txt, CLI, Success)

```json
{
 "row_generation_total": 0.0250,
 "algebraic_rows": 0.0007,
 "coordinate_rows": 0.0243,
 "tangent_fields": 0.0,
 "tangent_rows": 0.0,
 "lf_flags": 0.0097,
 "records_total": 1.3223,
 "assemble_rows_mod_p": 0.0257,
 "ranking": 1.2783,
 "rref_mod_p": 0.0157,
 "extract_normal_form": 0.0002,
 "reconstruction": 0.1507,
 "certificate_total": 0.0046,
 "certificate_points_total": 0.0046
}
```

## Snapshot 2: fast D4 diagnostic (deg1, 9 samples x 2 primes, no RUN_D4_FULL)

Config mirrors `test_d4_reduce_family_once_current_config` (labels=T+M1..M5,
`max_ibp_degree=1`); honest outcome is `NormalFormNotLocallyFinite` at this
degree — timings are still representative of the hot path.

```json
{
 "row_generation_total": 0.1497,
 "algebraic_rows": 0.0025,
 "coordinate_rows": 0.1471,
 "tangent_fields": 0.0,
 "tangent_rows": 0.0,
 "lf_flags": 0.0394,
 "records_total": 13.8668,
 "assemble_rows_mod_p": 0.0888,
 "ranking": 13.6928,
 "rref_mod_p": 0.0818,
 "extract_normal_form": 0.0001,
 "reconstruction": 0.0039,
 "certificate_total": 0.0517,
 "certificate_points_total": 0.0517
}
```

## Hotspot analysis

- **`ranking` dominates everything**: 96.7% of `records_total` on tiny
  (1.28/1.32 s) and 98.7% on D4 deg1 (13.69/13.87 s). This is the Perf.1
  target; `rref_mod_p` / `assemble_rows_mod_p` are two orders of magnitude
  cheaper.
- Certificate step is negligible (< 0.06 s even on D4 deg1), and
  `certificate_points_total` accounts for essentially all of
  `certificate_total` — no hidden overhead in the wrapper.
- `reconstruction` is visible only on the Success path (tiny: 0.15 s); D4
  deg1 fails before it matters.
- Heavy D4 / Example4* runs deliberately NOT taken yet (per Perf.0.1 scope);
  extrapolation from the deg1 ranking share suggests they are ranking-bound
  as well.
