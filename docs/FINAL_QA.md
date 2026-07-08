# Final QA — v0.1.0 (Finalize.1)

Environment: Windows 10, Python ≥3.11, package `parametric-ibp-lf-reducer` 0.1.0.

## Gates

| Check | Command | Result |
|---|---|---|
| Fast test suite | `python -m pytest` | **216 passed, 6 skipped** (all 6 skips are `RUN_D4_FULL`-gated heavy D4 tests; they pass in the heavy run below → all 222 tests green) |
| Lint | `ruff check .` | **clean** (`All checks passed!`) |
| Heavy D4 acceptance | `RUN_D4_FULL=1 python -m pytest tests/test_d4_vertical.py tests/test_d4_cli_e2e.py -q` | **9 passed**, exit 0 (~25–30 min) |

## D4 release fingerprint (must reproduce exactly)

- `status=Success`, `all_lf=True`
- `n_rows=2092`
- `n_records=108` (`n_reduced=108`)
- `n_selected=102` (`n_bad_spec=0`, `n_target_not_pivot=0`)
- `rank_histogram={1995: 6, 2041: 102}`, selected rank `2041`, 6 records rank-filtered
- certificate: **Passed 3/3** (rank-filtered 0, rank-exceeded 0, bad 0)
- terms (target support): `{M1, M2, M3}` =
  `{(0,1,1,0,-3,-1,0), (0,1,1,0,-2,-1,0), (1,1,0,0,-2,-1,0)}`

## Known limitations

- Wolfram-like text exchange format only (no CAS runtime integration).
- Explicit parametric family required in input
  (`FAILURE_PARSER_NEEDS_EXPLICIT_FAMILY` otherwise).
- D4 is the only curated/certified acceptance configuration.
- No adaptive search / protected masters / forced 5-term basis by design;
  failures surface as `ALL_FAILURE_REASONS` diagnostic codes.
- Heavy acceptance is opt-in via `RUN_D4_FULL=1`.

## Deliverable packaging

Zip command (PowerShell, excludes caches/venv/build artifacts) — see
release notes; caches (`__pycache__`, `.pytest_cache`, `.ruff_cache`) are
not part of the deliverable.
