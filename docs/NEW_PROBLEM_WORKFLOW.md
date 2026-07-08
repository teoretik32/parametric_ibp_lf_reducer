# NEW PROBLEM WORKFLOW — bringing a fresh integrand into the reducer

Practical recipe for setting up a new reduction problem with
`parametric_ibp_lf_reducer` (Release.1). Companion to `USAGE.md`
(install / CLI / exit codes) and `docs/01_clean_spec_ru.md` (theory).

Reminder (Release.1 limitation): **explicit-family input only** and **no adaptive
search** — one fixed pass per invocation. Everything below is about doing that
factorization and search enlargement *yourself*.

---

## 1. From integrand to explicit-family input

A document carrying just a whole `Integrand` returns the typed failure
`Failure/ParserNeedsExplicitFamily` — auto-factorization is deliberately not guessed.
You must split the integrand by hand into the family form

```
x^(nu + n) * Product_i  P_i(x; params)^(alpha_i + m_i)
```

Steps:

1. **Identify variables** — the integration variables (`"Variables"`); the domain is
   the positive orthant (`"Domain" -> "PositiveOrthant"`).
2. **Identify parameters and regulators** — every free symbol that is not an
   integration variable goes into `"Parameters"`; the subset acting as dimensional
   regulators (shifting exponents, e.g. `ep`) also goes into `"Regulators"`.
3. **Factor the integrand** into irreducible (or at least fixed) polynomial factors
   `P_i` with parametric exponents, times a pure monomial in the variables.
   - polynomial *bases* go into `"Polynomials"` (named, e.g. `"G0"`, `"G1"`, ...);
   - the per-variable monomial exponents (the `nu` part, may contain regulators)
     go into `"MonomialExponents"`;
   - the per-polynomial exponents (the `alpha` part) go into `"PolynomialExponents"`.
4. **`"TargetMultiplier"`** — extra rational multiplier of the target (usually `1`).
5. **`"Assumptions"`** — positivity/genericity assumptions on parameters, e.g. `{r > 0}`.
6. Wrap everything as `IBPInput = <| ... |>;` (see template in §5).

The *target* of the reduction is the label `(n, m) = (0, ..., 0)` relative to these
base exponents unless `"TargetLabel"` overrides it in `"Options"`.

Sanity check of your factorization: multiply everything back symbolically (by hand or
any CAS) and compare with the original integrand *before* running the reducer.

---

## 2. Choosing labels and the label box

A **label** is a pair of integer shift vectors `(n, m)`:
`n` shifts the variable monomial exponents (one entry per variable, in `"Variables"`
order), `m` shifts the polynomial exponents (one entry per polynomial, in
`"Polynomials"` order). Example for 4 variables + 3 polynomials:
`{0,1,1,0,-2,-1,0}` = `n = (0,1,1,0)`, `m = (-2,-1,0)`.

`"LabelBox"` bounds the search space: per-coordinate `{min, max}` ranges, first the
`n` block, then the `m` block:

```
"LabelBox" -> {{{0,1},{0,1},{0,1},{0,0}}, {{-4,0},{-1,0},{0,0}}}
```

Guidelines:

- `"LabelBox" -> Automatic` is fine for a first probe; switch to an explicit box once
  you know where the masters live.
- Keep the box **small first** — cost grows with box volume. Enlarge only along
  directions suggested by failures (§3).
- Negative `m` entries mean *raising* inverse powers of that polynomial (deeper
  denominators); this is where masters usually sit for IBP-style reductions.
- `"Labels"` (explicit list) pins the candidate set exactly; `"PreferredMasters"`
  biases the ranking toward a known basis without forcing it (the D4 example prefers
  M1..M5 yet certifies a 3-term basis — equivalent within the certified row span).
- Prefer labels whose shifted exponents stay **locally finite** on the orthant;
  non-LF candidates show up as `non_lf_terms` / `unknown_lf_terms` in diagnostics.

---

## 3. Reading failure statuses

`status` in the diagnostics JSON is `"Success"` or one of `ALL_FAILURE_REASONS`
(exit code 1 — honest typed failure; result text + JSON are still written):

| status | meaning | first fix to try |
|--------|---------|------------------|
| `ParserNeedsExplicitFamily` | document has only an `Integrand`, no explicit family | do §1 by hand |
| `TargetNotReducible` | target label not in the row span generated inside the current box/degrees | enlarge `LabelBox`, raise `MaxIBPDegree`, add `TangentDegrees` blocks |
| `InterpolationFailed` | dense coefficient reconstruction could not fit the scattered samples | add more (scattered!) `Samples`, check degrees aren't huge |
| `NormalFormNotLocallyFinite` | reduction exists but some term is not LF on the orthant | move/extend the box toward LF labels; adjust `PreferredMasters` |
| `VerificationFailed` | reconstructed coefficients fail re-checks at independent points | more `Samples`, more/different `Primes`, raise `MinValidRecords` |
| `ResourceLimitReached` | configured resource cap hit before completion | shrink box/degrees, or accept a longer run with raised limits |

Independently, `certificate_status` gates `Success` (row-span certificate is
default-ON): `"Passed"` / `"Failed"` / `"Insufficient"` / `"NotRun"`.
`"Insufficient"` means too few usable rank-generic points — add off-sample
`CertificatePoints`. `Success` implies `"Passed"`.

Exit codes: `0` Success, `1` typed failure (see table), `2` usage / malformed
document (nothing was reduced — fix the input file, not the search).

---

## 4. Expanding the search

There is **no adaptive search** — each knob is enlarged manually between runs, one
axis at a time so you can attribute the effect:

1. **`MaxIBPDegree`** (CLI: `--max-ibp-degree`) — row-generation degree cap; the
   cheapest lever against `TargetNotReducible`. Raise by 1 at a time.
2. **`TangentDegrees`** (alias of `TangentDegreeBlocks`) — extra tangent-field degree
   blocks, e.g. `{{1,1},{2,2}}`; adds qualitatively new rows, not just more of the same.
3. **`LabelBox`** — widen coordinate ranges where failures point (usually more
   negative `m`). Volume ↗ cost, so widen selectively.
4. **`Samples`** — scattered rational points for `Parameters`, e.g.
   `<| "ep" -> 15/7, "r" -> 14/3 |>`. Use **non-lattice, scattered** points —
   product-lattice grids can validate wrong interpolants. Add points when
   `InterpolationFailed` / `VerificationFailed`.
5. **`MinValidRecords`** (CLI: `--min-valid-records`) — reconstruction evidence
   floor; raise it together with `Samples` for confidence, lower only for probing.
6. **`Primes`** — large machine primes for modular passes (e.g. `2147483647`,
   `2147483629`, `2147483587`); add primes on suspected unlucky-prime behavior.
7. **`CertificatePoints`** — off-sample rank-generic rational points for the
   certificate gate; add when `certificate_status` is `"Insufficient"`.

Unknown/informational `Options` keys (e.g. `"SurfaceMode"`, `"OutputFormat"`) are
tolerated and passed through — misspelled *known* keys silently fall back to
defaults, so double-check spelling against the list above.

---

## 5. Input file template

```
IBPInput = <|
  "Variables"  -> {x1, x2},
  "Parameters" -> {ep, r},
  "Regulators" -> {ep},
  "Domain"     -> "PositiveOrthant",

  "Polynomials" -> <|
    "G0" -> 1 + x1 + x2,
    "G1" -> 1 + r*x1*x2
  |>,

  "MonomialExponents"   -> <| x1 -> -1 - ep, x2 -> ep |>,
  "PolynomialExponents" -> <| "G0" -> 2*ep, "G1" -> -2 - ep |>,

  "TargetMultiplier" -> 1,
  "Assumptions" -> {r > 0},

  "Options" -> <|
    "MaxIBPDegree"    -> 2,
    "TangentDegrees"  -> {{1,1}},
    "LabelBox"        -> Automatic,          (* or explicit {{n-ranges},{m-ranges}} *)
    (* "TargetLabel"      -> {0,0, 0,0}, *)
    (* "PreferredMasters" -> { {0,1,-2,0}, ... }, *)
    "Primes"          -> {2147483647, 2147483629, 2147483587},
    "MinValidRecords" -> 16,
    "Samples" -> {
      <| "ep" -> 2,    "r" -> 17/6 |>,
      <| "ep" -> 15/7, "r" -> 14/3 |>
      (* ... more scattered non-lattice points ... *)
    }
    (* "CertificatePoints" -> { <| "ep" -> ..., "r" -> ... |>, ... }, *)
  |>
|>;
```

Run:

```bash
python -m parametric_ibp_lf_reducer reduce my_input.wl.txt \
    --out result.m --diagnostics-json diag.json
```

Working references: `examples/tiny_success_input.wl.txt` (minimal, ~1-2 s) and
`examples/d4_cli_example_input.wl.txt` (full verified configuration, ~10-15 min).

---

## 6. Diagnostics checklist

After every run, in order:

- [ ] **Exit code**: `2` → fix the document (parse/usage), don't touch search knobs.
- [ ] **`status`**: `"Success"` or which of `ALL_FAILURE_REASONS` → §3 table.
- [ ] **`certificate_status`**: `"Passed"`? `"Insufficient"` → add `CertificatePoints`;
      `"NotRun"` → the run failed before the gate.
- [ ] **`certificate` counters**: `n_certificate_points_passed/failed`, `selected_rank`.
- [ ] **`all_locally_finite`** and per-term `locally_finite`; count
      `non_lf_terms` / `unknown_lf_terms` (should be 0 for a clean answer).
- [ ] **`diagnostics.n_records` vs `n_skipped_records`**: many skips → samples are
      degenerate/unlucky, replace or add scattered points.
- [ ] **`zero_reduction`**: target collapsed to 0 — suspicious unless expected;
      re-verify the §1 factorization.
- [ ] **`reconstruction_verified`** and **`independent_validation_passed`** are true.
- [ ] **`messages`**: read them — they name the failing stage.
- [ ] **Cross-check** (when a reference value exists): compare against known results
      as in `examples/notebook_star_example3_known_value_input.wl.txt`.
- [ ] Changed exactly **one** search knob since the last run? If not, rerun the
      bisection before trusting conclusions.
