# Correction: starred notebook examples are value-level examples, not LF decomposition examples

Earlier wording may have incorrectly described `Example 3*` and `Example 4*` from `Examples_for_IBP_parametric.nb` as candidate input families for main-spec regressions 11.4 and 11.5.

Correct policy:

1. Starred examples are a separate type of example.
2. For starred examples, the notebook gives an integral and an analytic answer/epsilon expansion for that integral.
3. For starred examples, the notebook does not give a decomposition into locally finite integrands with rational IBP coefficients.
4. Therefore, starred examples must not be used as e2e reducer expected outputs.
5. Starred examples may be used for parser/sparse-poly fixtures and future value-level checks only.
6. Main-spec regressions 11.4 and 11.5 remain pending unless their actual input families and LF decomposition expectations are provided separately.

Implementation implication:

- Do not create `test_reduce_star_example3_to_id4example2_expected_one_term`.
- Do not create `test_reduce_star_example4_to_id3example3_basis_inclusion`.
- Do create parser tests that verify the starred inputs can be represented without hardcoding names or monomial degrees.
- Do create a sparse-poly test that verifies `x4^2` support from starred Example 4 input.
