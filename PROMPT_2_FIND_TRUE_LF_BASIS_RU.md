# Prompt 2 для Claude Code / Fable: доведение софта до настоящего LF-базиса External Int2

```text
Pure symbolic mathematics / parametric IBP implementation audit.
No security/network code.

External Int2 full Laurent oracle is now integrated and green.
The exact value through ep^0 is an acceptance oracle, not a substitute for an LF basis.

Domain requirement from Leonid:
A genuine solution must express the target through integrals that are EACH individually locally finite at ep=0 and usable separately in HyperInt.
Do not accept cancellation between separately divergent integrals.
Do not change the definition of success.
Do not pursue quasi-finite dimension shifts or a different problem definition at this stage.

Current evidence:
- certified six-term relation exists but contains 1/G1 and 1/G2, LF=False;
- large box/tangent/coordinate expansions remained Obstructed in the current implementation;
- finite-polynomial-numerator search in tested problematic sectors was negative;
- dual witnesses certify obstruction only for the generated/projected row spaces;
- full analytic Laurent answer through ep^0 is known and must be reproduced by any final LF decomposition.

Goal:
find the first concrete discrepancy between our implementation and the intended parametric-IBP method, repair it generically, and obtain an actual LF-basis reduction if possible.

Do NOT launch another blind huge RREF before a cheap diagnostic shows new valid rows that break the obstruction witness.

Phase A — exact family/convention audit

1. Reconstruct the source integrand symbolically from ParametricFamily target label zero.
Compare line by line with the notebook after r=s/t and external-prefactor extraction.
Check:
- variables and order;
- G0..G3;
- base monomial exponents;
- base polynomial exponents;
- target multiplier/label;
- external Exp[2*ep*EulerGamma] Gamma prefactor.
Add a symbolic equality test.

Phase B — surface-policy audit

2. Derive the exact boundary condition used by surface.py for coordinate and tangent rows.
Separate:
- LF-master classification at ep=0;
- validity of an IBP identity in an open convergence chamber.

3. Add diagnostic-only exact-rational surface evaluation at ep=-3/5, r=1.
Do NOT change production semantics yet.

4. On a SMALL Int2 box classify candidate rows:
- valid both near ep=0 and at ep=-3/5;
- invalid both;
- rejected near ep=0 but strictly surface-suppressed at ep=-3/5.

5. Pair newly chamber-valid rows with the recorded dual obstruction witnesses.
If none has nonzero pairing, stop this branch.
If valid rows break witnesses, run only a small LF-feasibility solve before full reconstruction.

Phase C — full tangent/syzygy module completeness

6. Current bounded SymPy degree ansatz is not a proof of tangent-module completeness.
Create a Singular script/backend that computes, for the Int2 polynomials:
- the syzygy module Q.grad(G_l)=H_l*G_l for each G_l;
- the required module intersection;
- exact module generators.

7. If Singular is unavailable locally:
- generate an executable .sing file;
- document install/run commands;
- do not fake results.

8. Verify every returned generator exactly in Python/SparsePoly.
Compare its span and degrees to current tangent fields.
Generate only surface-valid rows and pair them with obstruction witnesses.

Phase D — positive-control differential diagnosis

9. Create docs/PARAMETRIC_IBP_POSITIVE_CONTROL.md listing the exact data required from Leonid for one successful 2mh case:
- source integrand/family;
- target;
- known LF masters;
- coefficients or reduction output.

10. When provided, run four checks independently:
- LF checker recognizes known masters;
- row generation without surface filter reaches them;
- production surface filter reaches them;
- full Singular module vs bounded ansatz.
This must identify whether the discrepancy is family transcription, LF geometry, surface policy, or tangent completeness.

Phase E — controlled solve

11. Only after Phase B or C produces valid new witness-breaking rows:
- run LF-feasibility at >=3 generic samples and >=2 primes;
- require generic Feasible, not rank-deficient special points;
- reconstruct rational coefficients;
- require row-span certificate Passed;
- require every exported RHS integral LF=True;
- compare the resulting LF-master epsilon expansions with the integrated Laurent oracle through ep^0.

12. If no implementation discrepancy is found:
return a scoped audit result, not global impossibility, and list the exact missing positive-control data.

Tests:
- target-integrand reconstruction;
- exact surface sign at rational chamber point;
- diagnostic surface mode cannot alter production results;
- exact syzygy identities;
- witness-breaking screening;
- no False/Unknown LF term can enter Success;
- analytic Laurent oracle remains green.

Run:
python -m pytest
ruff check .

No release/tag in this pass.
No commit until review.

Report max 40 lines:
- family convention verdict;
- surface-mode row counts and witness-breaking count;
- Singular/module verdict;
- positive-control status;
- first confirmed implementation discrepancy;
- LF-feasibility result if justified;
- whether a genuine HyperInt-ready LF basis was obtained;
- exact next action.
```
