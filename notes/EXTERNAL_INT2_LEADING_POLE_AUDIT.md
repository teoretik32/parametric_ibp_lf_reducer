# External Int2 — Method.2 leading-pole audit (wrapper level)

Status: all 8 checks passed at full precision; wrapper prefactor corrected.
Artifacts: `scripts/audit_external_int2_leading_pole.py`,
`validation/external_int2_leading_pole_audit.json`,
`tests/test_external_int2_leading_pole.py`.

## Scope guard (what this audit does NOT touch)

- No import of `parametric_ibp_lf_reducer`, no RREF, no reducer-core change
  (`reducer_core` flags in the JSON are all `false`; pinned by fast tests).
- The certified pure-family reduction of Int2 (its rows, labels, coefficients,
  certificate) is untouched. Everything below is wrapper/reference metadata.
- `AnsvInt2` is NOT invented and never enters the reducer; its source text is
  stored as metadata only in `examples/external_int2_source_reference.wl.txt`
  (copied mechanically from the notebook). Only the leading pole
  `-4/(s*t^2*ep^4)` is compared here.

## Method.2 in one page

1. **Exact x7 preintegration.** For `A = 1 + r*x2*x5 > 0`, `B = 1 + x2 > 0`:
   `Integral[(1+x7)^(-1-ep)*(A+B*x7)^(-1+ep), {x7,0,Inf}] = (B^ep - A^ep)/(ep*(B-A))`
   (checked symbolically and numerically), so
   `J2(ep,r) = (1/ep)*Int[x2^ep*(1+x2)^ep*(1+x5)^ep*((1+x2)^ep-(1+r*x2*x5)^ep)/(1-r*x5)]`
   over `{x2,x5} > 0`; the integrand is regular at `x5 = 1/r`.
2. **Reduced 1-D form.** The x2 integral is a Gauss 2F1 (G&R 3.197.1):
   `J2 = (G1(ep)/ep)*Q(ep,r)`, `G1 = Gamma[1+ep]*Gamma[-1-3*ep]/Gamma[-2*ep]`,
   `Q = Int[(1+x5)^ep*(1 - 2F1(-ep,1+ep;-2*ep;1-r*x5))/(1-r*x5), {x5,0,Inf}]`.
3. **Pole bookkeeping (the subtle point).** Two boundary layers degenerate to
   `x5^(-1)` as `ep -> 0`: the tail `K1`-layer (`x5 -> Inf`) and the head
   `C_B`-layer (`x5 -> 0`, hypergeometric crossover). `K1 == C_B` *identically*,
   the crossover poles cancel exactly, and only the tail "1"-layer survives:
   `Q = 1/(r*ep) + O(1)`, hence

   `J2(ep,r) = -2/(3*r*ep^2) + O(1/ep)`.

   A naive tail-only count would give `-1/(2*r*ep^2)`; the numeric Laurent fit
   (`leading_pole_numeric`) confirms `c2*r = -2/3` and excludes `-1/2` cleanly.
4. **Corrected external prefactor.** The net `gamma_E` count of the Gamma ratio
   is exactly `-2*ep`, so the wrapper previously missing normalisation
   `Exp[2*ep*EulerGamma]` is required for EulerGamma-free Laurent data:

   `P2 = Exp[2*ep*EulerGamma]*t^(-3-ep)*Gamma[1-ep]*Gamma[-ep]^3*Gamma[ep]`
   `     /(Gamma[-1-3*ep]*Gamma[-2*ep]) = 6/(t^3*ep^2) + O(1/ep)`.
5. **Leading pole of the full object.** With `r = s/t`:
   `P2*J2 = (6/(t^3*ep^2))*(-2/(3*(s/t)*ep^2)) = -4/(s*t^2*ep^4) + O(1/ep^3)`,
   matching the supplied `AnsvInt2` leading pole exactly.

## What changed (wrapper/metadata only)

- `scripts/run_external_int2.py`: `EXTERNAL_PREFACTOR_TEXT` now carries
  `Exp[2*ep*EulerGamma]*...`; the mpmath spot check multiplies both sides by the
  same factor (comparison unaffected); docstring/artifact header updated.
- `validation/external_int2_full_formula.m` and
  `validation/external_int2_diagnostics.json`: mirrored prefactor text.
- `tests/test_external_int2.py::test_prefactor_text_matches_p2` pins the
  corrected `P2` (sympy, exact).

## Checks (full precision, `--fast` skips only the 2-D cross-check)

| check | result |
| --- | --- |
| `x7_identity_symbolic` / `x7_identity_numeric` | exact / passes |
| `hyp2f1_connection_formulas` | passes |
| `decomposition_consistency` (h, X, N, delta invariance) | passes |
| `prefactor_series` (lead 6, corrected EulerGamma-free, old not) | passes |
| `leading_pole_exact` (`K1 == C_B`, `Q`-pole, product) | exact |
| `leading_pole_numeric` (Laurent fit, `c2*r = -2/3`) | passes |
| `j2_2d_cross_check` (2-D quadrature vs 1-D decomposition) | passes |

Quadrature note: the 2-D cross-check needs `dps = 30`; at `dps = 20` nested
tanh-sinh accuracy caps near `5e-7` relative (first full run flagged exactly
this, probe-verified), while `dps = 30` agrees to `~5e-9`. Tolerance is `1e-7`.

## Reproduce

```
python scripts/audit_external_int2_leading_pole.py           # full, writes JSON
python scripts/audit_external_int2_leading_pole.py --fast    # bounded rerun
python -m pytest tests/test_external_int2_leading_pole.py -q # fast guards
```
