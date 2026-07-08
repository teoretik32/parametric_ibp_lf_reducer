# D4 STATUS (canonical 11.3)

Tracking the D4 vertical validation. See `VERTICAL_AUDIT.md` for the full-pipeline audit and
`tests/test_d4_vertical.py` for the tests. Nothing here changes core behaviour or fabricates
`Success`.

## D4.1 — modular row-span certificate: **PASSED** ✅ (now over 3 points)
The expected relation `J[T] = sum_i C_i J[M_i]` (C_i = exact reference coefficients) is certified
**in the generated row span** at 3 safe `(ep, r, prime)` points — `(2/3,3)`, `(3/2,4)`, `(5/4,6)`.
Row system: box `n=[(0,1),(0,1),(0,1),(0,0)]`, `m=[(-4,0),(-1,0),(0,0)]`, coord IBP deg2, tangent
`[(1,1),(2,2)]` → **2092 rows** (algebraic 240, coordinate 1826, tangent 26; surface-rejected
1558). Reducing `J[T] − ΣCi·J[Mi]` by the RREF pivots gives residual `{}` at every point. So row
generation is sufficient; this is not a one-point accident.

## D4.2 — full-config `reduce_family_once` with preferred_masters: **InterpolationFailed**
Config: full box, deg2, tangent `[(1,1),(2,2)]`, 3 primes × 16 `(ep,r)` samples,
`min_valid_records=16`, `preferred_masters=[M1..M5]`. Observed (~342 s, opt-in
`RUN_D4_FULL=1`):

| field | value |
|---|---|
| status | `InterpolationFailed` (all_locally_finite `Unknown`) |
| terms | `[]` (reconstruction produced nothing) |
| n_rows | 2092 |
| n_records | 48 |
| n_reduced_records | 48 |
| n_bad_specializations | 0 |
| n_target_not_pivot | 0 |

### What actually happens (single-point probes, `modular_normal_form` + preferred_masters)
- At **generic** points the target reduces to a **stable, all-LF subset `{M1,M2,M3}`** (3 terms,
  every term `is_locally_finite=True`, zero non-master labels). Verified at `ep=3,r=2`, `ep=4,r=5`,
  `ep=7/2,r=9/2` (Test `test_d4_target_reduces_to_lf_subset_of_masters`).
- So **`preferred_masters` already works**: only locally-finite masters are left free; M4, M5 are
  reducible to M1,M2,M3 inside this row system, giving a *smaller* LF basis than the reference
  5-term one. That is a legitimate basis-cardinality difference, not a non-LF leak.
- But the **generic rank is 2041 while some sampled integer points are rank-deficient**
  (e.g. `ep=2,r=3` → rank 1995; `ep=3/2,r=4` → 2011). At a rank-deficient point the target's
  normal-form **support shrinks/shifts** (e.g. only `{M2}`), and `collect_value_table`'s union
  support 0-fills the missing masters. Those spurious zeros make each master's value sequence
  incoherent as a rational function → `interpolate_multivariate` fails its holdout → the whole
  reconstruction reports `InterpolationFailed`.

### First blocker to exact 11.3 acceptance
**Reconstruction / sample rank-consistency — NOT `preferred_masters`, ranking, or row generation.**
The reducer finds a valid LF reduction at every generic point, but does not (a) detect and drop
rank-deficient specializations before reconstruction, nor (b) fix a single free-master support
across samples. Concretely the next pass should filter records to the modal/maximal-rank support
(or majority free-set) so reconstruction sees a coherent value table.

### Is `protected_masters` needed?
**No — not for this symptom.** The masters are already kept free and are LF; the failure is a
reconstruction-input consistency problem. A protected/forced free-basis might later be *one* way to
pin the support to exactly `{M1..M5}` if we want to match the reference basis, but the immediate
blocker is rank-consistent sample selection. (Also open: confirm `{M1,M2,M3}` is a *correct*
reduction, i.e. the rows reducing M4,M5 are surface-valid, before treating the 3-term basis as the
answer. Out of D4.2 scope.)

## D4.3 — rank-consistency filtering before reconstruction: **Success (3-term LF basis)**
> **⚠️ D4.4 correction:** the 6×6-lattice run below passed the strict gate, but the D4.4
> row-span certificate later showed its *interpolated coefficients are wrong off-lattice*
> (product-grid degeneracy fooled the on-lattice holdout — see §D4.4). The rank-filtering fix
> itself stands (it is what made the run reach reconstruction at all); the sampling grid was
> the remaining flaw.
Generic fix (no D4 hardcode): `reconstruction.select_records_for_reconstruction(records,
rank_policy="max_rank")` — only `Reduced`+`formal_success` records at the **maximal observed RREF
rank** feed reconstruction (a specialization's rank can only drop below the generic rank, so
max-rank records = generic points; rank-deficient records solve a smaller system and their
shrunken support must not be union-0-filled). `collect_value_table`/`reconstruct_coefficients`
apply it by default (`rank_policy="all"` = old behaviour, tests/debug only); the reducer checks
`min_valid_records` against the **post-filter** count and exports selection diagnostics
(`record_selection`: rank_histogram, selected_rank, n_rank_filtered_records,
support_after_rank_filter) plus a message. See assumption A29.

Full-config rerun (36-sample 6×6 grid `ep,r ∈ {2..7}`, 3 primes, deg2, tangent
`[(1,1),(2,2)]`, `preferred_masters=[M1..M5]`, 795 s, opt-in `RUN_D4_FULL=1`):

| field | value |
|---|---|
| status | **`Success`** (via the strict gate; all_locally_finite `True`) |
| terms | `{M1, M2, M3}` — all `LocallyFinite=True` |
| n_records / n_reduced | 108 / 108 (0 bad-spec, 0 target-not-pivot) |
| rank_histogram | `{1995: 18, 2041: 90}` → 6 rank-deficient samples × 3 primes filtered |
| n_selected | 90 (min_valid_records=16 satisfied post-filter) |

So the D4.2 blocker (rank-poisoned value table) is fixed and the reducer produces its **first
honest D4 Success**: reconstruction validated on independent holdout points, every term LF.

### Honest caveats / what this Success is NOT yet
- The basis is the reducer's **3-term LF basis `{M1,M2,M3}`**, not the reference 5-term
  `M1..M5` — M4, M5 are reducible inside this row system, so the symbolic C1..C5 comparison
  branch does not apply. Exact 11.3-reference acceptance (equal basis + coefficients) would need
  the basis pinned to M1..M5 (e.g. a protected/forced-free policy — still NOT implemented).
- ~~Still open (inherited from D4.2): independently confirm the 3-term reduction against the
  reference decomposition~~ → **closed by D4.4** (row-span equivalence certificate below). What
  remains open is only the *analytic* validity of the rows themselves (surface checks are the
  guarantee; they were applied to every generated row).

## D4.4 — acceptance + equivalence certificate
Generic helper (no D4 hardcode, no ``Success`` stamp):
`certificate.verify_reduction_relation_mod_p(family, rows, target, terms, sample, prime)` →
`CertificateResult` — assembles the rows mod p, evaluates the claimed coefficients exactly
(SymPy / `ParamExpr` / `int` / `Fraction`), reduces the relation vector
`J[target] − Σ C_i·J[label_i]` by the RREF pivots and reports `InSpan`/`NotInSpan` with the
honest residual (`BadSpecialization`/`EmptySystem` reject a point, never patch it). Fast unit
tests: `tests/test_certificate.py`.

### The certificate immediately caught a real bug in the D4.3 "Success" ⚠️
Re-running the 6×6 **integer-lattice** full config and certifying the reducer's own
reconstructed relation at off-grid points: **NotInSpan** at `ep=2/3, r=3` (rank-generic 2041!),
while the reference relation IS in span there. Direct probe: the TRUE normal form at that point
is `{M1,M2,M3}` with values `(−5, 5, −2) mod p` — the interpolated coefficient functions return
different values. **Root cause: the 6×6 product lattice is degenerate for the dense degree
search** — `Π(ep−k), k=2..7` has total degree 6 = `max_deg` and vanishes on the whole lattice,
*including the holdout points* (they lie on the same lattice), so a wrong candidate passed
"independent" validation. The D4.3 result was gate-formally Success but its coefficients were
wrong off-lattice. Lesson (A30): holdout points on the sampling lattice are not independent;
sample grids for multivariate reconstruction must be scattered/non-degenerate, and row-span
certification at off-sample points is the reliable acceptance check.

### Fix (test-config only; no core change)
The heavy fixture now samples **35 scattered rational points** (`ep = 2 + k/7`,
`r = 2 + ((11k+5) mod 36)/6` — no product structure, no low-degree curve through the points)
plus the known rank-deficient `(2,3)` to keep the D4.3 filter exercised. Reducer-output
certificates use `REDUCER_CERT_POINTS` — probe-verified **rank-generic (2041)** off-sample
points (`(2/3,3)`, `(5/4,6)`, `(7/3,9/2)`); note `(3/2,4)` from the reference `CERT_POINTS` is
rank-deficient (2011) and is deliberately NOT used for reducer-output certification.

### Final rerun outcome (scattered grid, 2026-07-07, 933 s, opt-in): **ALL CERTIFIED** ✅

| field | value |
|---|---|
| status | **`Success`** (strict gate; all_locally_finite `True`) |
| terms | `{M1, M2, M3}` — all `LocallyFinite=True` |
| n_records / n_reduced / n_selected | 108 / 108 / 102 (0 bad-spec, 0 target-not-pivot) |
| rank_histogram | `{1995: 6, 2041: 102}` → deficient records filtered (incl. the planted `(2,3)`) |
| reducer-output certificate | **InSpan** at all 3 off-sample rank-generic points |
| reference certificate | **InSpan** at the same points ⟹ **equivalence certified** |
| heavy suite | 13 passed (`RUN_D4_FULL=1`, includes both certificates + M4/M5 diagnostic) |

### D4 acceptance statement
1. **Reference certificate** — `T = Σ Ci·Mi` (C1..C5) is in the generated row span
   (`test_d4_expected_relation_is_in_row_span_mod_p`, via the generic helper; runs by default,
   3 points) — **PASSED**.
2. **Reducer-output certificate** — the full-config result `T = Σ C_red_i·L_i` is in the SAME
   row span at the independent points (`test_d4_success_result_is_row_span_certified`, opt-in
   `RUN_D4_FULL=1`) — **PASSED** on the scattered grid.
3. **Equivalence** — both relations vanish modulo the same row span at the same points ⟹
   `Σ Ci·Mi − Σ C_red_i·L_i` is itself in the span: the decompositions are **equivalent modulo
   the generated IBP/algebraic rows** (`test_d4_reducer_relation_equivalent_to_reference`,
   opt-in) — **PASSED**. No coefficient-by-coefficient comparison — the bases differ.
4. **Why 3-term** — diagnostic `test_d4_m4_m5_reduce_to_smaller_basis` (runs by default,
   **PASSED**): with M1..M3 preferred free, **M4 and M5 each reduce to combinations of
   {M1,M2,M3}** in this row system. The smaller basis is a legitimate LF basis, not a defect.
   An exact 5-term reference basis would require an optional forced/protected-basis mode — not
   needed for mathematical Success.

## D4.5 / Pass Verify.1 — certificate gate inside the reducer
The A30 lesson is now part of the pipeline, not just the tests. `reduce_family_once` /
`reduce_rows_once` certify the reconstructed relation at independent off-sample points
(`certificate.verify_reduction_relation_mod_p`) **before** the Success gate:

- **Config:** `certificate_points` (explicit; auto-generated deterministic off-sample points
  beyond the per-parameter sample maximum when empty), `certificate_primes` (defaults to the
  reduction primes), `require_certificate_for_success=True` (default),
  `min_certificate_points=1`, `certificate_rank_policy="selected_rank"` (only supported policy;
  anything else → `ValueError`).
- **Point classification (Verify.1):** rank **<** `selected_rank` (the reconstruction's rank) →
  rank-filtered (uninformative); rank **>** `selected_rank` → **hard failure**
  (`n_certificate_rank_exceeded`: a specialization's rank can never exceed the generic rank, so
  the reconstruction's `selected_rank` was not generic); row/coefficient pole → bad; otherwise
  pass/fail by residual.
- **`certificate_status`:** `Failed` (any informative nonzero residual, or any rank-exceeded
  point) / `Passed` (≥ `min_certificate_points` informative passes, no failures) /
  `Insufficient` (no informative points) / `NotRun` (no reconstruction, or gate disabled with
  no points).
- **Gate:** with the default `require_certificate_for_success=True`, `Success` requires
  `certificate_status == "Passed"`; anything else → `Status="Failure"`,
  `Error="VerificationFailed"`. `FormalSuccess` stays honest in Diagnostics; per-term
  `LocallyFinite` flags stay truthful; exported `AllLocallyFinite` is never `True` on a failure
  (2I.1a contract). Diagnostics in `extra["certificate"]`: n points / passed / failed /
  rank-filtered / bad, `selected_rank`, `certificate_rank_histogram`,
  `first_nonzero_residual`, points used; plus messages.
- **Regression (fast, `tests/test_certificate_gate.py`):** the product-grid false success is
  reproduced synthetically — true coefficient `2 + Π(ep−k), k=3..8` on the degenerate grid
  `ep=3..8` interpolates to the wrong constant `2` and passes the on-grid holdout; the gate now
  returns `VerificationFailed` (old behaviour available only via explicit
  `require_certificate_for_success=False`). Plus (Verify.1): correct reduction →
  `Passed`+`Success`; post-reconstruction coefficient corruption → `Failed`; a deficient point
  among certificate points is skipped+counted while generic passes still yield `Success`;
  rank-exceeded point → `Failed` honestly; all-rank-filtered / all-bad points →
  `Insufficient` → failure; unsupported `certificate_rank_policy` → `ValueError`;
  `certificate_primes` override honoured; `VerificationFailed` exports
  `Status -> "Failure"`, `Error -> "VerificationFailed"`.
- **D4:** the heavy fixture passes the 3 probe-verified rank-generic `REDUCER_CERT_POINTS`
  explicitly; the recorded D4 `Success` is now a *certified* Success by construction.

## Tests (in `tests/test_d4_vertical.py`)
- `test_d4_expected_labels_are_lf_and_target_is_not_lf` (fast) — target non-LF, M1..M5 all LF.
- `test_d4_expected_relation_is_in_row_span_mod_p[ep,r,prime]` (integration ×3) — the reference
  certificate, via the generic D4.4 helper.
- `test_d4_target_reduces_to_lf_subset_of_masters[ep,r,prime]` (integration ×2) — generic normal
  form ⊆ {M1..M5}, all LF.
- `test_d4_row_counts_are_reported` (integration) — records the 2092-row breakdown.
- `test_d4_reduce_family_once_current_config` — small-config diagnostic (currently
  `NormalFormNotLocallyFinite`).
- `test_d4_reduce_family_once_full_config_preferred_masters` (integration, **opt-in
  `RUN_D4_FULL=1`**, ~13 min shared via the `d4_full_result` module fixture) — records the
  full-config outcome above; asserts the rank filter is wired + engaged; branches per status
  without *requiring* `Success` (observed: `Success`).
- `test_d4_success_result_is_row_span_certified` (integration, opt-in) — D4.4 acceptance:
  requires `Success`, all LF, terms ⊆ {M1..M5}, reducer relation certified in span ×3 points.
- `test_d4_reducer_relation_equivalent_to_reference` (integration, opt-in) — D4.4 equivalence:
  reference + reducer relations both vanish mod the same row span at the same points.
- `test_d4_m4_m5_reduce_to_smaller_basis[M4|M5]` (integration, runs by default) — M4, M5 reduce
  to {M1,M2,M3}: why the basis is 3-term.
- Rank-filtering unit/orchestration tests: `tests/test_rank_filtering.py` (Pass D4.3);
  certificate-helper unit tests: `tests/test_certificate.py` (Pass D4.4).
