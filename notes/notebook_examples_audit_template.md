# Notebook examples audit

For each notebook example, classify it before writing tests:

| Example | Has initial integrand | Has JIntegrandList | Has JCoefficientListt | Has known integral value/epsilon expansion | Allowed test role now |
|---|---:|---:|---:|---:|---|
| Example 1 / I4exampl1 | yes | yes | yes | not the key role | parser/coefficient now; optional known-decomposition e2e later |
| Example 2 / I3exampl2 | yes | yes | yes | not the key role | parser/coefficient now; optional known-decomposition e2e later |
| Example 3* / I4exampl3 | yes | no | no | yes | parser/sparse-poly only; future value-level check |
| Example 4* / tm1Int | yes | no | no | yes | parser/sparse-poly only; x4^2 fixture; future value-level check |

Rules:

- Do not map starred examples to main-spec 11.4/11.5 without separate evidence.
- Do not use main-spec 11.4/11.5 expected outputs as expected outputs for starred examples.
- Do not use Gamma/prefactor objects in starred examples as reducer coefficients.
- Any e2e reducer test must have an explicit LF basis and rational IBP coefficients or be marked xfail/pending.
