# Parametric-IBP positive control: exact data request (Phase D)

Purpose. External Int2 is currently `Obstructed` at every recorded generic (sample, prime) point
under the production surface policy, while the chamber-policy audit (Method.9/Method.10) shows the
first policy-dependent lever. To decide WHICH layer is wrong — family transcription, LF geometry,
surface policy, or tangent-module completeness — we need one KNOWN-GOOD case ("positive control")
from the same method family: a 2mh (or comparable) integral where the intended parametric-IBP
method demonstrably produced an LF-basis reduction.

Scope guard. The positive control never changes the definition of success: every exported RHS
integral must be individually locally finite at `ep = 0` and usable separately in HyperInt.
Cancellation between separately divergent integrals is not acceptable. No quasi-finite dimension
shifts. The External Int2 Laurent oracle (through `ep^0`) remains the acceptance value check.

## 1. Data required from Leonid (exact list)

For ONE successful case, all four items, machine-readable (Wolfram-style text is fine):

1. **Source integrand / family** — exactly as fed to the working implementation:
   - integration variables and their order;
   - all base polynomials `G_l` (fully expanded or factored, either way unambiguous);
   - base monomial exponents and base polynomial exponents (as functions of `ep` and any
     kinematic parameters), including the sign/offset conventions;
   - kinematic parameters with their positivity assumptions and any rescaling already applied
     (e.g. `r = s/t`, external prefactor split).
2. **Target** — the label/multiplier of the integral that was reduced (in the same offset
   convention), plus the external prefactor if it was factored out before reduction.
3. **Known LF masters** — the explicit list of final basis integrands (labels or integrands)
   that the working implementation produced, each individually locally finite at `ep = 0`;
   confirmation that each was integrated separately in HyperInt.
4. **Coefficients or reduction output** — the rational-function coefficients `c_i(ep, params)`
   of the target in terms of the masters (or the raw reduction output from which they can be
   read), so equality can be verified at exact rational sample points mod p.

Optional but very valuable:

- the label box / degree bounds the working run used (IBP multiplier degree, tangent/logarithmic
  vector-field degrees);
- which IBP identities were admitted at the boundary (surface-term policy: limit at `ep -> 0`,
  chamber point, or symbolic condition);
- any preprocessing (sector decomposition, variable ordering, partial fractioning) applied
  before row generation.

## 2. The four independent checks we will run on it

Recorded as `validation/positive_control_*.json`, one artifact per check:

| # | Check | Question answered |
|---|-------|-------------------|
| 1 | `is_locally_finite` on every known master | Does our LF checker recognize the known-good masters as LF-True? (LF geometry layer) |
| 2 | Row generation WITHOUT surface filter, span test to the masters | Do our generators even produce the identities that reach the masters? (row-completeness layer) |
| 3 | Same span test with the PRODUCTION surface filter | Does the surface policy reject identities the working method used? (surface-policy layer) |
| 4 | Full Singular syzygy module vs bounded SymPy ansatz on the control polynomials | Is the tangent module truncated? (tangent-completeness layer) |

Decision table: the first check that fails on the positive control identifies the defective
layer; if all four pass on the control but Int2 stays obstructed, the discrepancy is in the
family transcription of Int2 itself (re-audit Phase A against the sender's exact conventions).

## 3. Status

- Requested: not yet sent (this document is the request text).
- Received: —
- Checks run: —

Until the positive-control data arrives, Int2 verdicts remain scoped per-(sample, prime),
per-label-box statements; no global impossibility claim is made or implied.
