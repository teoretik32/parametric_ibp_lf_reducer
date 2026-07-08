# USAGE — parametric_ibp_lf_reducer

Pure-Python locally-finite parametric IBP reducer. Wolfram/Mathematica syntax appears
ONLY as a text exchange format (input/output documents); no Wolfram runtime is ever
invoked.

## Install & test

```bash
python -m venv .venv
source .venv/bin/activate          # PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -e '.[dev,speed]'
python -m pytest                   # full fast suite
ruff check .
```

## Tiny example (~1-2 s, full pipeline including certificate)

```bash
python -m parametric_ibp_lf_reducer reduce examples/tiny_success_input.wl.txt \
    --out result.m --diagnostics-json diagnostics.json
```

Expected: exit code 0, `"Status" -> "Success"`, one locally finite term.
Details: `examples/tiny_success_expected_notes.md`.

Equivalent API call:

```python
from parametric_ibp_lf_reducer import api
result = api.reduce_wolfram_style_input(open("examples/tiny_success_input.wl.txt").read())
print(result.status, result.wolfram_style_text)
```

## D4 heavy example (SLOW: ~10-15 minutes)

`examples/d4_cli_example_input.wl.txt` is the D4 validation family with the full
verified configuration embedded in its `"Options"` association (36 scattered
non-lattice sample points, label box, preferred masters M1..M5, off-sample
rank-generic certificate points). It is self-contained:

```bash
python -m parametric_ibp_lf_reducer reduce examples/d4_cli_example_input.wl.txt \
    --out d4_result.m --diagnostics-json d4_diagnostics.json
```

Expected (certificate-verified, recorded outcome): exit 0, `"Status" -> "Success"`,
`"AllLocallyFinite" -> True`, `certificate_status "Passed"`, and a **3-term** LF basis
with labels in {M1, M2, M3} = {(0,1,1,0,-2,-1,0), (1,1,0,0,-2,-1,0), (0,1,1,0,-3,-1,0)}.
See "Limitations" below on why this differs from the 5-term reference basis.

The matching opt-in test (same warning — ~10-15 min per module):

```bash
RUN_D4_FULL=1 python -m pytest tests/test_d4_vertical.py tests/test_d4_cli_e2e.py
```

(PowerShell: `$env:RUN_D4_FULL = "1"` first.)

## Exit codes

| code | meaning |
|------|---------|
| 0 | reduction reached `Status -> "Success"` (LF + certificate gates passed) |
| 1 | reduction ran but did not reach `Success` — honest typed failure; result text and JSON are still written; reason on stderr and in `"error"` |
| 2 | usage / I/O / malformed-document error — nothing was reduced |

## Statuses

| status | meaning |
|--------|---------|
| `Success` | all terms locally finite **and** the certificate gate passed |
| `Failure` (coarse exported `Status`; concrete reason in `Error` / JSON `status`) | honest typed failure: `TargetNotReducible`, `InterpolationFailed`, `NormalFormNotLocallyFinite`, `ResourceLimitReached` |
| `VerificationFailed` | a reduction was found formally, but independent verification (reconstruction check / row-span certificate) rejected it |
| `ParserNeedsExplicitFamily` | the input document lacks an explicit parametric family; nothing was reduced |

## diagnostics JSON fields

- `status` / `exported_status` / `success` / `error` — typed outcome; `status` is either
  `"Success"` or one of the `ALL_FAILURE_REASONS` constants.
- `certificate_status` — `"Passed"` / `"Failed"` / `"Insufficient"` / `"NotRun"`
  (row-span certificate gate is default-ON; `Success` implies `"Passed"`).
- `certificate` — scalar certificate counters (`n_certificate_points`,
  `n_certificate_points_passed/failed`, `selected_rank`, ...).
- `target_label`, `all_locally_finite`
- `terms[]` — `label`, `coefficient` (text, `^` powers), `integrand` (text),
  `locally_finite` (True/False/"Unknown").
- `diagnostics` — `formal_success`, `reconstruction_verified`,
  `independent_validation_passed`, `n_terms`, `non_lf_terms`, `unknown_lf_terms`,
  `n_records`, `n_skipped_records`, `zero_reduction`, `messages`.

## Limitations (Release.1)

- **Explicit-family input only.** A document carrying just a whole `Integrand` returns a
  typed `Failure/ParserNeedsExplicitFamily` — integrand auto-factorization is not guessed.
- **No adaptive search.** One fixed pass per invocation: the label box / degrees /
  samples come from the document `Options` (or defaults); nothing enlarges them on failure.
- **No Mathematica runtime dependency** — and therefore no symbolic cross-check against
  Wolfram; verification is modular row-span certification at rational points.
- **Dense multivariate interpolation limits.** Coefficient reconstruction uses dense
  degree search over scattered rational samples; high coefficient degree or many
  parameters needs more samples than the defaults, and product-lattice grids can
  validate wrong interpolants (use scattered points, as the examples do).
- **D4 reduces to a 3-term LF basis** {M1,M2,M3}, equivalent to (and certified against
  the row span containing) the 5-term reference basis M1..M5; the reference basis is
  deliberately NOT forced.

## Example 4* (exploratory, known-value-only)

`examples/example4_star_input.wl.txt` is a *known-value-only* example: the
ε-expansion of the integral itself is known
(`validation/notebook_star_example4_known_value_expansion.txt`), but **no**
reference LF decomposition exists for it. An exploratory reducer run returned a
certified `Success` (certificate `Passed` 3/3) with a 2-term all-locally-finite
reduction — artifacts: `validation/example4_star_result.m`,
`validation/example4_star_diagnostics.json`; details:
`notes/example4_star_exploratory.md`. This is **not** part of the certified
baseline (the single curated end-to-end configuration remains D4), and a numeric
cross-check against the known value is impossible without the master-integral
values.

## Release sanity check

```bash
python -m pytest
ruff check .
# optional heavy acceptance (~25-30 min total):
RUN_D4_FULL=1 python -m pytest tests/test_d4_vertical.py tests/test_d4_cli_e2e.py
```

or `scripts/final_check.ps1` / `scripts/final_check.sh` (pass `-Heavy` / `--heavy` to
include the D4 runs).

Note: `ruff format --check` still wants to reformat 22 files that predate the current
formatting config; this is pre-existing and intentionally not mixed with logic changes.
`ruff check .` is clean.
