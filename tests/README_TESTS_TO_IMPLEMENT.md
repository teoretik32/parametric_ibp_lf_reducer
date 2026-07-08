# Tests to implement

Build real pytest tests from `docs/03_validation_cases_ru.md` and `validation/expected_d4_coefficients.json`.

Minimum order:

1. Parser tests for explicit Wolfram-like association syntax.
2. Sparse polynomial tests: arbitrary `N`, monomial degrees > 1, derivative/valuation.
3. N=2 tangent-field sanity test: `G = 1+x+y`, `Q=(xy,-xy)` gives row without `m` shift.
4. Surface tests: coordinate IBP checks only `x_i=0,∞`; vector/tangent uses toric flux.
5. Local-finiteness tests at `epsilon=0`; never count epsilon-regulated convergence as LF.
6. D=4 expected LF basis/coefficient checks using independent samples.
