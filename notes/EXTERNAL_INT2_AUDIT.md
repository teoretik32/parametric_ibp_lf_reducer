# External Int2 (dimensionless) — reduction audit

Status: **stable negative result** — the target is *not* reducible to a locally
finite normal form within the explored adaptive label boxes. The partial
reduction itself is certified and independently validated; the obstruction is
mathematical, not numerical.

## Input

- Document: `examples/external_int2_dimensionless_input.wl.txt`
  (`ExternalInt2Dimensionless`).
- Family: vars `('x2', 'x5', 'x7')`, polys `('G0', 'G1', 'G2', 'G3')`,
  target label `(0,0,0,0,0,0,0)`.
- Runner: `scripts/run_external_int2.py` (gated by `RUN_EXTERNAL_INT2=1`).

## Heavy runs

Four heavy runs were performed. Early configurations failed with
`InterpolationFailed` at deeper adaptive levels; boosting scattered
samples/primes (numeric hardening) removed the interpolation failures and
exposed the underlying outcome, which is stable across independent
configurations:

- Run #3 (boosted samples/primes): `NormalFormNotLocallyFinite`, non-LF labels
  `{0,0,0,0,-1,0,0}` and `{0,0,0,0,0,-1,0}`.
- Run #4 (`base+boost-s48-p6-x1`, deepened label box + boosted
  samples/primes, wall time 60030.6 s): identical status and identical
  non-LF labels.

Run #4 quality gates (from `validation/external_int2_diagnostics.json`):

- Certificate: **Passed 3/3** (rank-filtered 0, rank-exceeded 0, bad 0);
  rank filter kept 531/540 reduced records at rank 22361
  (histogram `{19722: 9, 22361: 531}`); 2048 labels.
- `ReconstructionVerified -> True`, `IndependentValidationPassed -> True`,
  `FormalSuccess -> True`.
- Dominant cost: record generation (~59430 s of the 60030.6 s total);
  certificate ~466 s, reconstruction ~15 s.

The built-in recommendation ("expand the label box around the remaining
non-LF/Unknown labels — a different reduction path may avoid them") was
already implemented in run #4 and did not change the outcome.

## Certified result (`validation/external_int2_result.m`)

Six-term decomposition, `AllLocallyFinite -> False`:

| Integrand        | Coefficient                                   | Locally finite |
|------------------|-----------------------------------------------|----------------|
| `1/(x2*G0*G1)`   | `(ep - 1)/(6*ep*r)`                           | yes |
| `1/(x2*G1*G3)`   | `-(ep - 1)^2*(r + 1)/(6*ep^2*r)`              | yes |
| `1/(G0*G3)`      | `(ep - 1)*(2*ep^2*r + ep^2 - 1)/(6*ep^3*r)`   | yes |
| `1/G1`           | `1`                                           | **no** |
| `1/G2`           | `-(ep + 1)/ep`                                | **no** |
| `x7/(G0*G3)`     | `(ep - 1)*(ep + 1)*(2*ep - 1)/(6*ep^3*r)`     | yes |

Outputs: `validation/external_int2_result.m`,
`validation/external_int2_full_formula.m`,
`validation/external_int2_diagnostics.json`.

## Independent numeric check

Not applicable in the Int1 sense: with `Status -> Failure` there is no
certified ep-expansion to compare against quadrature. The reduction identity
itself is covered by the internal 3-prime certificate, the multi-prime
reconstruction check and the independent validation pass (all green).

## Conclusion / follow-ups

- The two residual terms `1/G1` and `-(ep+1)/ep * 1/G2` carry labels with a
  `-1` entry (positions 5 and 6) and are genuinely non-locally-finite in this
  basis; deeper boxes and harder numerics reproduce them exactly.
- Any future attempt should change the *method*, not the knobs: e.g. a
  different basis/ordering, sector decomposition of the residual terms, or an
  analytic treatment of `1/G1`, `1/G2` outside the IBP reduction.

## Method.1: directional LF audit + LF-constrained feasibility (mod p)

Script: `scripts/run_external_int2_method1.py`; outputs
`validation/external_int2_method1_levelA.json`,
`validation/external_int2_method1_levelB.json`.

| Level | box | labels | rows | LF True/False/Unknown | Phase B verdict | Phase C |
|---|---|---|---|---|---|---|
| A (base) | `[0,1]^3 x [-2,0]^4` | 648 | 5904 | 249/399/0 | Obstructed (0/6) | 0/6, support n/a |
| B (expand-1) | `[0,1]^3 x [-3,0]^4` | 2048 | 29104 | 1169/879/0 | Mixed (2/6) | 2/6, support_stable=True, common_support_size=67 |

- Target label `[0,0,0,0,0,0,0]`: LF verdict **False** at both levels
  (23/30 failing rays, 0 unknown) — consistent with the certified
  non-locally-finite finding for the residual terms.
- The 2 feasible points at Level B are both primes at the sample
  `ep=3, r=54/11` (rank 18422 vs 20963 at the generic samples). Since `ep=3`
  is a non-generic integer point with a visible rank drop, this feasibility
  is treated as a special-locus artifact, not evidence of a generic
  reduction path; the two generic samples stay Obstructed at both primes.
- Per the scope note in the JSON, all statements are per-(sample, prime)
  about this row system and label box only; certificate and LF gates are
  untouched.
- Final gate after the runs: full pytest suite green, `ruff check .` clean.
- Level A elapsed ~117 s; Level B elapsed ~2197 s (fork run in background).

## Method.3: composite locally-finite master feasibility

Module: `src/parametric_ibp_lf_reducer/composite_masters.py`; runner:
`scripts/run_external_int2_method3.py`; output:
`validation/external_int2_composite_feasibility.json`; tests:
`tests/test_composite_masters.py` (heavy integration gated by
`RUN_EXTERNAL_INT2=1`).

Question: instead of demanding that every *individual* master be locally
finite (Method.1 says that fails), do integer linear combinations
`M = sum_i c_i * J(label_i)` exist whose bad Laurent layers cancel
identically? Answer: **yes — `FeasibleCompositeBasis`**.

- Pool: 225 candidates from the 6 non-LF terms of the certified normal form
  (shifts along `n_x2`, `m_G0`, `m_G3` at depths −1/−2 plus `x5`/`x7`
  numerator insertions up to degree 2).
- On the primary ray `(-1,0,0)` (`x2 -> oo`): 48 participants,
  primary cancellation kernel dimension 21.
- After verification on all 69 checked rays (candidate rays + deterministic
  random safety net, 27 witness rays) the kernel refines to a
  **13-dimensional fully-LF composite basis** — every basis vector is locally
  finite on every checked ray, with exact symbolic layer cancellation.
- Interpretable examples: `J(1/(x2*G1)) - J(1/(G0*G1))` (2 terms) and
  `J((1+x5)/G1) - J((1+x7)/G2)` (4 terms) — the certified non-LF residuals
  `1/G1`, `1/G2` *do* combine into fully LF composites once numerator
  insertions are allowed.
- Scope guard: statements are about this pool and checked-ray set only;
  coefficients live in rational functions of `r` with a fixed-sample rank
  cross-check (`BadSpecialization` guard); reducer core, certificates and LF
  gates are untouched.
- Follow-up: rerun the reduction in a basis where the non-LF residuals are
  replaced by these composites (basis change at the normal-form level), aiming
  at `AllLocallyFinite -> True` for Int2.
- Elapsed ~31 s (pool build + feasibility, foreground).

## Finite-numerator LF basis search (task #37)

Module: `src/parametric_ibp_lf_reducer/finite_numerator.py`; design:
`docs/FINITE_NUMERATOR_BASIS_DESIGN.md`; runner:
`scripts/run_external_int2_finite_numerator.py`; output:
`validation/external_int2_finite_numerator.json`; tests:
`tests/test_finite_numerator.py`.

Question: does ONE numerator-decorated integrand `N(x) * F_S` (genuine
polynomial `N`, complete-integrand semantics — no cancellation after
integration is ever assumed) exist per remnant sector whose full integrand is
locally finite? Answer: **no — `NoFiniteNumeratorBasisWithinAnsatz`**
(degrees 0–2, seven sectors: the six certified normal-form sectors plus the
probe `1/(G1*G3)`).

- Labels are SHIFTS (offset convention: total exponent = base exponent +
  label shift against the Int2 base
  `x2^(1+ep) * G0^ep * G1^ep * G2^(-1-ep) * G3^(-1+ep)`).
- `1/(x2*G0*G1)`, `1/(x2*G1*G3)`, `1/(G0*G3)`, `x7/(G0*G3)`:
  `SectorAlreadyLF` (bare integrand already LF; a numerator cure is moot).
- `1/G1`, `1/G2`, `1/(G1*G3)`: `NumeratorCureImpossibleAnyDegree` — every
  failing ray is componentwise `<= 0` (`x -> oo` type; the probe's only bad
  ray is `(-1,0,0)`, `x2 -> oo`), and polynomial numerators only increase
  those layer scores (Lemma 2), so the impossibility holds for EVERY degree,
  not just 0–2. The leading-cancellation kernel is empty at all searched
  degrees.
- Lemma 1 (graded lowest layer) is consistent everywhere
  (`lemma_consistent_everywhere: true`): separately divergent pieces never
  combine into an accepted candidate.
- Feasibility stage: `SkippedNoCandidates` — nothing to bridge into
  `lf_reduction_feasible_mod_p`; the bridge itself (defining monomial
  expansions as rows, expansion labels marked allowed per Lemma 1) is
  unit-tested on the synthetic curable family.
- Scope: read-only diagnostics; reducer state, certificates and LF gates
  untouched. Elapsed ~10 s.
- Conclusion: numerator decoration alone cannot replace the remnants; the
  viable routes remain the Method.3 composite basis change (label shifts with
  positive polynomial components) or analytic treatment of `1/G1`, `1/G2`
  outside the reducer.
