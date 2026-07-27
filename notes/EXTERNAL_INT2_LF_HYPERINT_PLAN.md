# External Int2 — HyperInt plan for the certified LF masters (Method.12A)

Scope: linear-reducibility audit + exact input preparation ONLY.  The certified
LF basis (4 labels, 4 coefficients), surface/LF/certificate semantics and all
Method.11 artifacts are unchanged.  No RREF, no modular records, no integration
was started.

## Inputs

- Family: `examples/external_int2_dimensionless_input.wl.txt`
  (variables x2, x5, x7; parameter r; domain = positive orthant).
- Certified reduction: `validation/external_int2_lf_full_formula.m`
  (`ReductionIdentity` — an identity AFTER integration, via
  `J[f_] := Inactive[Integrate][f, {x2,0,Infinity}, {x5,0,Infinity}, {x7,0,Infinity}]`;
  no pointwise integrand identity is implied).
- Certificate provenance: `CertificateProvenance` in
  `validation/external_int2_lf_result.m` (Passed, 4 points, 8 checks,
  generic rank 26984, convergence_chamber, ep=-3/5).
- Physical wrapper: `ExternalPrefactor2 = Exp[2*ep*EulerGamma]*t^(-3-ep)
  *Gamma[1-ep]*Gamma[-ep]^3*Gamma[ep]/(Gamma[-1-3*ep]*Gamma[-2*ep])`, r = s/t.

## Linear-reducibility audit (Fubini polynomial reduction)

Audit script: `scripts/audit_external_int2_lf_linear_reducibility.py`
Machine-readable report: `validation/external_int2_lf_linear_reducibility.json`

Initial polynomial set for every master: {x2, x5, x7, 1+x2, 1+x5, 1+x7,
1+x7+x2*x7+r*x2*x5} (all four G-polynomials carry nonzero generic exponents in
each master).  Per master, all 3! = 6 integration orders over (x2, x5, x7) were
enumerated independently; each master is audited on its own.

| Master | Label                     | Reducible | Valid orders | Chosen order | ep order |
|--------|---------------------------|-----------|--------------|--------------|----------|
| L1     | [-1, 0, 0, -1, -1, 0, 0]  | yes       | 6/6          | x2, x5, x7   | ep^2     |
| L2     | [-1, 0, 0, 0, -1, 0, -1]  | yes       | 6/6          | x2, x5, x7   | ep^3     |
| L3     | [0, 0, 0, -1, 0, 0, -1]   | yes       | 6/6          | x2, x5, x7   | ep^4     |
| L4     | [0, 0, 1, -1, 0, 0, -1]   | yes       | 6/6          | x2, x5, x7   | ep^4     |

Alphabet (identical upper bound for all four masters, letters in r):

```
x2, x5, x7, 1+x2, 1+x5, 1+x7, r*x5 - 1, r*x5 + x7, x7 - r,
1 + x7 + x2*x7 + r*x2*x5
```

Final-variable letters are linear (roots -1 and r), so the results are
r-parametric multiple polylogarithms — no algebraic extension needed.
No obstruction polynomial exists for any master or any order.

## Prepared HyperInt inputs (review gate — do NOT run unattended)

- `validation/hyperint/external_int2_lf_L1.mpl` — through ep^2
- `validation/hyperint/external_int2_lf_L2.mpl` — through ep^3
- `validation/hyperint/external_int2_lf_L3.mpl` — through ep^4
- `validation/hyperint/external_int2_lf_L4.mpl` — through ep^4

Each file defines the exact integrand, the epsilon truncation and the audited
integration order; the `hyperInt`/`fibrationBasis` calls are present but
commented out.  Long integration starts only after human review.

## Gate status

All four masters linearly reducible => inputs prepared; integration deferred.
Next phase (separate approval): run HyperInt per master, assemble
`ReductionIdentity` with `ExternalPrefactor2`, and cross-check the analytic
Laurent expansion against the certified rational coefficients C1..C4.
