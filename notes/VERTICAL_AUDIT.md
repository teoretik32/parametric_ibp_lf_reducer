# VERTICAL AUDIT — parser/family → … → result gate (after Pass 2I.2)

Date: 2026-07-06. Method: code read of the chain + a throwaway D4 dry-run
(`scratchpad/d4_dryrun.py`, not production). No production code changed. Key finding: the vertical
slice is **fully connected** and already reduces + reconstructs D4; the only blocker to a correct
D4 Success is that `reduce_family_once` leaves **non-LF masters free → LF gate fails**.

> **UPDATE (Pass D4.1, 2026-07-07):** the guess below that the *row system is incomplete* is
> **REFUTED** — see the D4.1 addendum at the bottom. With the full box + coord deg2 + tangent, the
> expected relation `J[T]=ΣCi·J[Mi]` is **certified in the generated row span** (modular certificate
> passes). The real blocker is **basis-selection** (ranking/label set), not row generation.

## 1. What is actually wired into `reduce_family_once`?
Chain inside `reduce_family_once(family, target, config)`:
`enumerate_box/labels` → `is_locally_finite` (LF flags) → `row_generation` (algebraic +
coordinate-IBP surface-filtered + optional tangent-IBP) → `collect_normal_form_records`
(**which internally runs `ranking` → `assemble_rows_mod_p` → `sparse_rref` → `modular_normal_form`
+ LF diagnosis**) → `reconstruct_coefficients` (uni/multivariate) → `result` strict gate.
Note: **parser/family is upstream, NOT inside** `reduce_family_once` (it takes an already-parsed
`ParametricFamily`; `parse_family_text` is called by the tests / future `reduce_wolfram_style_input`).
All other layers are genuinely connected — confirmed by the D4 dry-run running the whole path.

## 2. Which tests use real generated rows vs synthetic?
- **Real generated rows:** `test_row_generation*.py`, `test_tangent_rows.py`,
  `test_modular_normal_form.py`, `test_records.py` (row level), and the reducer-level
  `test_reduce_family_once_runs_end_to_end_smoke` (tiny generic family, asserts a valid typed
  result, not Success).
- **Synthetic rows (hand-built `Row`+`ParamExpr`, injected `lf_flags`):** all other
  `test_reducer_orchestration.py` cases (via `reduce_rows_once`).
So exactly **one reducer-level test drives real row generation**; the failure-mapping tests are synthetic.

## 3. Can canonical D4 (11.3) run as a dry-run now?
**Yes.** The dry-run parses the family, generates rows, reduces the target at 18 `(prime,sample)`
points, and **reconstructs multivariate `C(ep,r)` that passes holdout validation** — it reaches the
LF gate and returns a typed `NormalFormNotLocallyFinite` result. Nothing in the chain crashes.

## 4. First concrete cause D4 is not a Success
**The LF gate (`NormalFormNotLocallyFinite`).** With `labels=[T,M1..M5]`, `max_ibp_degree=1`, no
tangent, target reduces but the free masters are `{0,0,0,0,-1,0,0}` (=`1/G0`) and
`{0,1,0,0,-1,0,0}` (=`x2/G0`) — both **non-LF** — instead of M1..M5. Not parser / not label-box /
not target-not-pivot (it *is* a pivot) / not insufficient-records / not reconstruction. Root cause:
**row system too small to eliminate those non-LF intermediates**, so ranking (which correctly tiers
non-LF first) has no relation to pivot them out.

## 5. D4 labels for M1..M5 (N=4 vars, M=3 polys; label = (x1,x2,x3,x4,G0,G1,G2))
- target `T = (0,0,0,0, 0,0,0)`
- M1 `x2*x3/(G0^2*G1)` = `(0,1,1,0, -2,-1,0)`
- M2 `x1*x2/(G0^2*G1)` = `(1,1,0,0, -2,-1,0)`
- M3 `x2*x3/(G0^3*G1)` = `(0,1,1,0, -3,-1,0)`
- M4 `x1*x2/(G0^3*G1)` = `(1,1,0,0, -3,-1,0)`
- M5 `x2*x3/(G0^4*G1)` = `(0,1,1,0, -4,-1,0)`
LF check (dry-run): target → **False** (correct), M1..M5 → **all True**.

## 6. Row / surface counts on D4 at a small search config
- algebraic: **18** rows, 0 rejected.
- coordinate-IBP deg1: **80** rows, **32** `surface_not_free`. deg2: **239** rows, **81** `surface_not_free`.
- tangent fields `[(1,1)]`: **0** fields (G2 needs degree-2 syzygies); `[(1,1),(2,2)]`: **4** fields
  → **4** tangent rows, **20** `surface_not_free`.
- `reduce_family_once` (labels=6, deg1, no tangent): n_rows=**98**, n_records=**18** (2 primes × 9
  samples), n_reduced=**18**, target_not_pivot=**0**, bad_spec=**0** → LF gate fail (2 non-LF masters).

## 7. Hardcoded synthetic assumptions in reducer tests?
**Yes, by design (unit-level):** `reduce_rows_once` tests hand-build single-row systems with
explicit `ParamExpr` coefficients, **inject `lf_flags`** rather than computing them, and assume the
pivot outcome. They test *wiring + failure mapping*, not real LF/row-gen. This is fine but means
reducer coverage of the real row→rank→rref→LF path rests on the single smoke test (+ the layer tests).

## 8. Shortest path to the first real D4 attempt
The pipeline already attempts D4; the shortest path to a **correct** attempt (not new architecture):
1. `preferred_masters = [M1..M5]` so ranking keeps them free while eliminating non-LF intermediates;
2. `tangent_degree_blocks = [(1,1),(2,2)]` (the `(2,2)` block yields the 4 tangent/syzygy rows for G2);
3. `max_ibp_degree = 2`;
4. a `label_box` (or explicit label set) that also **includes the non-LF intermediates** (e.g. `1/G0`,
   `x2/G0`, and the G0-power / x1,x2,x3 shift chain down to M5) so there are pivots to eliminate them;
5. enough distinct `(ep,r)` samples for dense multivariate reconstruction (≥ ~16, ≥4 per axis).
This is a **config + modest search** change plus a Wolfram-in wrapper (Pass 2I.3), not a rewrite.

---

# ADDENDUM — Pass D4.1 (2026-07-07): targeted D4 vertical validation

Test-only pass (`tests/test_d4_vertical.py`, 4 tests). No production code changed. Purpose: split
"row-generation completeness" from "basis-selection".

## Modular row-span certificate — RESULT: **IN SPAN** ✅
`test_d4_expected_relation_is_in_row_span_mod_p` (`@pytest.mark.integration`): with box
`n=[(0,1),(0,1),(0,1),(0,0)]`, `m=[(-4,0),(-1,0),(0,0)]`, coord deg2, tangent `[(1,1),(2,2)]`, at
`ep=2/3, r=3`, `prime=2147483647` — assemble → RREF → reduce the vector
`J[T] − ΣCi·J[Mi]` (Ci from `expected_d4_coefficients.json`, taken mod p) by the RREF pivot rows →
**residual == {} (zero)**. So the expected 11.3 relation with the *exact expected* C1..C5 is a
genuine linear combination of the generated rows. **Row generation is sufficient at this config.**

## D4 row counts (box above, deg2, tangent[(1,1),(2,2)])
seeds=80; algebraic=240; coordinate_ibp=1826 (rejected 1264 surface); tangent_fields=4,
tangent_ibp=26 (rejected 294 surface); **total_rows=2092**.

## Current `reduce_family_once` outcome
Config `labels=[T,M1..M5]`, deg1, no tangent (fast diagnostic) → **`NormalFormNotLocallyFinite`**:
target reduces but frees the **non-LF** masters `1/G0` = `{0,0,0,0,-1,0,0}` and
`x2/G0` = `{0,1,0,0,-1,0,0}` instead of M1..M5. This is a **different LF basis / basis-selection**
outcome, NOT a row-span failure (the certificate above proves the M1..M5 relation exists).

## First blocker to exact 11.3 acceptance
**Basis-selection, not row generation.** The reducer's current label set + ranking do not force
M1..M5 as the free masters. Shortest fix (later pass, per §8 above): run the *full* box (so the
non-LF intermediates `1/G0`,`x2/G0` are present as pivot columns to eliminate), `max_ibp_degree=2`,
tangent `[(1,1),(2,2)]`, `preferred_masters=[M1..M5]`, and ≥16 `(ep,r)` samples — then the LF-tiered
ranking eliminates the non-LF intermediates and leaves M1..M5 free. Reconstruction of C1..C5 is
already proven feasible (multivariate path validated in Pass 2H/2H.1). Not attempted here.

Perf note: the certificate/counts tests build 2092 rows (~38 s together via a shared module fixture)
→ marked `@pytest.mark.integration` (still run by default; deselect with `-m "not integration"`).
