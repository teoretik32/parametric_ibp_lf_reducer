# External Int2 Method.11c: offline reconstruction audit

> **Superseded by the Method.11c resolution.** The zero-fill interpretation of
> `ep=6,r=57/11` below was wrong: the point is a basis/pivot specialization, not a
> special zero, and the four generic LF coefficients *are* reconstructible from the
> same 38x4 cache (see `notes/EXTERNAL_INT2_METHOD11_RECONSTRUCTION_AUDIT.md` and
> `validation/external_int2_method11_reconstruction_audit.json`). The Phase B-D
> conclusions below ("not reconstructed", the Stage-5 schedule) are historical.
> The on-disk value table stays at schema `m11c-value-table/v1` as recorded history;
> regenerating with `build_value_table` now emits `m11c-value-table/v2`, where the
> deviating sample is never zero-filled and appears under
> `verification["excluded_from_table"]`.

Offline audit of the Method.11 record cache (38 samples x 4 primes); no new modular
records, no RREF/reducer solve. Artifacts: `external_int2_method11_value_table.json`,
`external_int2_method11_reconstruction_audit.json`.

## Phase A - value table
- samples: 38, records: 152, primes: [2147483579, 2147483587, 2147483629, 2147483647]
- selected rank: 26984; support size: 4
- special-zero: ['ep=6,r=57/11'] missing [['(0, 0, 1, -1, 0, 0, -1)']]
- verification: OK

## Phase B - dense degree diagnostics
- `(-1, 0, 0, -1, -1, 0, 0)`: first admissible dense cell {'num_deg': 3, 'den_deg': 6, 'unknowns': 38, 'status': 'underdetermined (needs 37 fit points, have 36)'}

  Through (6,6): 33 determined cells are inconsistent (nullspace dim 0 mod all 4 primes, which proves dim 0 over QQ -- the coefficient is not a rational function of those degrees); 8 cells ([(3, 3), (3, 4), (4, 3), (3, 5), (4, 4), (5, 3), (4, 5), (5, 4)]) have nullspace dim >= 2 because 36 scattered points with all-distinct ep and r are degenerate for these supports (no unique candidate, holdout cannot arbitrate); the remaining 8 cells need 37..55 fit points but only 36 exist. Hence no dense cell <= (6,6) can validate: the true rational form lies outside the dense (6,6) ladder reachable with 38 samples.

- `(-1, 0, 0, 0, -1, 0, -1)`: first admissible dense cell {'num_deg': 3, 'den_deg': 6, 'unknowns': 38, 'status': 'underdetermined (needs 37 fit points, have 36)'}
- `(0, 0, 0, -1, 0, 0, -1)`: first admissible dense cell {'num_deg': 3, 'den_deg': 6, 'unknowns': 38, 'status': 'underdetermined (needs 37 fit points, have 36)'}
- `(0, 0, 1, -1, 0, 0, -1)`: first admissible dense cell {'num_deg': 3, 'den_deg': 6, 'unknowns': 38, 'status': 'underdetermined (needs 37 fit points, have 36)'}
- Stage-4 artifact crosscheck: {'available': True, 'n_attempts': 49, 'mismatches': [], 'matches': True}

## Phase C - structured experiments
- A: no repeated-r fibers exist in the cached samples (all 38 r values are distinct); univariate interpolation in ep is impossible without new data
- B: no repeated-ep fibers exist in the cached samples (all 38 ep values are distinct); univariate interpolation in r is impossible without new data
- C `(-1, 0, 0, -1, -1, 0, 0)`: 610 tensor cells, 418 inconsistent, 122 ambiguous, validated: 0
- C `(-1, 0, 0, 0, -1, 0, -1)`: 610 tensor cells, 516 inconsistent, 48 ambiguous, validated: 0
- C `(0, 0, 0, -1, 0, 0, -1)`: 610 tensor cells, 578 inconsistent, 8 ambiguous, validated: 0
- C `(0, 0, 1, -1, 0, 0, -1)`: 610 tensor cells, 578 inconsistent, 8 ambiguous, validated: 0
- D `(-1, 0, 0, -1, -1, 0, 0)`: forms ['(3*(7*ep) + 0*(11*r) + 7)', '(0*(7*ep) + 1*(11*r) + 0)', '(0*(7*ep) + 2*(11*r) + 1)', '(3*(7*ep) + 3*(11*r) + -1)'] cover=True -> denominator forms found, but value*denominator is not a polynomial of total degree <= 7 (the maximum determinable from 36 fit points)
- D `(-1, 0, 0, 0, -1, 0, -1)`: forms ['(3*(7*ep) + 0*(11*r) + 7)', '(0*(7*ep) + 1*(11*r) + 0)', '(1*(7*ep) + 0*(11*r) + 0)', '(1*(7*ep) + 1*(11*r) + 0)'] cover=True -> denominator forms found, but value*denominator is not a polynomial of total degree <= 7 (the maximum determinable from 36 fit points)
- D `(0, 0, 0, -1, 0, 0, -1)`: forms ['(3*(7*ep) + 0*(11*r) + 7)', '(1*(7*ep) + 0*(11*r) + 0)', '(0*(7*ep) + 1*(11*r) + 0)', '(2*(7*ep) + 3*(11*r) + -2)'] cover=False -> no product of small linear forms in (7*ep, 11*r) covers all denominator primes; factor-aware ansatz refuted within the searched form family
- D `(0, 0, 1, -1, 0, 0, -1)`: forms ['(3*(7*ep) + 0*(11*r) + 7)', '(1*(7*ep) + 0*(11*r) + 0)', '(0*(7*ep) + 1*(11*r) + 0)'] cover=True -> denominator forms found, but value*denominator is not a polynomial of total degree <= 7 (the maximum determinable from 36 fit points)
- E: no validated shared-denominator cell within scanned degrees
- G: no tested linear chamber boundary yields exact low-degree rational fits on both sides for all four coefficients (32 boundaries, 17 usable splits, near-misses [])

## Reconstructed coefficients
- `(-1, 0, 0, -1, -1, 0, 0)`: not reconstructed
- `(-1, 0, 0, 0, -1, 0, -1)`: not reconstructed
- `(0, 0, 0, -1, 0, 0, -1)`: not reconstructed
- `(0, 0, 1, -1, 0, 0, -1)`: not reconstructed

## Phase D - minimal new-sample schedule
- viable ansatz costs: [{'ansatz': 'dense (3,6) for (-1, 0, 0, -1, -1, 0, 0)', 'equations_per_sample': 1, 'additional_samples': 1, 'total_samples': 39}, {'ansatz': 'dense (3,6) for (-1, 0, 0, 0, -1, 0, -1)', 'equations_per_sample': 1, 'additional_samples': 1, 'total_samples': 39}, {'ansatz': 'dense (3,6) for (0, 0, 0, -1, 0, 0, -1)', 'equations_per_sample': 1, 'additional_samples': 1, 'total_samples': 39}, {'ansatz': 'dense (3,6) for (0, 0, 1, -1, 0, 0, -1)', 'equations_per_sample': 1, 'additional_samples': 1, 'total_samples': 39}, {'ansatz': 'shared-denominator frontier cell (7,1) beyond current data', 'equations_per_sample': 4, 'additional_samples': 1, 'total_samples': 39}, {'ansatz': 'factor-aware poly degree 8 for (-1, 0, 0, -1, -1, 0, 0)', 'equations_per_sample': 1, 'additional_samples': 8, 'total_samples': 46}, {'ansatz': 'factor-aware poly degree 8 for (-1, 0, 0, 0, -1, 0, -1)', 'equations_per_sample': 1, 'additional_samples': 8, 'total_samples': 46}, {'ansatz': 'factor-aware poly degree 8 for (0, 0, 1, -1, 0, 0, -1)', 'equations_per_sample': 1, 'additional_samples': 8, 'total_samples': 46}]
- minimum: {'ansatz': 'dense (3,6) for (-1, 0, 0, -1, -1, 0, 0)', 'equations_per_sample': 1, 'additional_samples': 1, 'total_samples': 39}
- recommended points: [{'ep': '15/7', 'r': '61/11', 'reason': 'fixed-ep fiber ep=15/7', 'design_rank_after': 37}, {'ep': '15/7', 'r': '62/11', 'reason': 'fixed-ep fiber ep=15/7', 'design_rank_after': 38}, {'ep': '15/7', 'r': '63/11', 'reason': 'fixed-ep fiber ep=15/7', 'design_rank_after': 39}, {'ep': '15/7', 'r': '64/11', 'reason': 'fixed-ep fiber ep=15/7', 'design_rank_after': 40}, {'ep': '15/7', 'r': '65/11', 'reason': 'fixed-ep fiber ep=15/7', 'design_rank_after': 41}, {'ep': '18/7', 'r': '6', 'reason': 'fixed-ep fiber ep=18/7', 'design_rank_after': 42}, {'ep': '18/7', 'r': '67/11', 'reason': 'fixed-ep fiber ep=18/7', 'design_rank_after': 43}, {'ep': '18/7', 'r': '68/11', 'reason': 'fixed-ep fiber ep=18/7', 'design_rank_after': 44}]

Notes: candidate points avoid ep=3, ep=6, r=57/11 and all existing samples; fiber points additionally enable univariate cross-checks (phases C.A/C.B); greedy rank growth uses seeded random residues as generic-position proxies (seed 20260726); ordering is fully deterministic
