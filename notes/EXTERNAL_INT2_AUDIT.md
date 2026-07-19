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
