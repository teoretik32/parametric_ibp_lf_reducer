# Finite-numerator locally-finite basis — design

Status: Phase 1 design for the External Int2 direction correction.

## 0. Correction of direction (what this replaces)

Method.3 (`composite_masters.py`) searched for *linear combinations of
different denominator sectors* whose divergent Laurent layers cancel between
members. Such composites are sums of separately divergent integrals: they are
finite only as a combination, every member still needs its own regulated
integral, and no member can be handed to HyperInt on its own. **That is not an
LF basis and is not used going forward.** Method.3 stays in the tree as a
diagnostic only.

The corrected goal: **individually locally finite basis integrals with
polynomial numerators**. Every exported object is ONE integrand

    F_label(x) * N(x),

with `N` a polynomial numerator, and the complete integrand must pass the
LF gate (`is_locally_finite = True`) on every relevant ray by itself. No
cancellation after integration, no cancellation between separately divergent
members is ever used.

## 1. Literature-oriented model review

Two established models for trading non-LF masters for finite objects:

1. **Finite numerator basis** (local numerators): decorate a divergent sector
   with a numerator polynomial chosen so the *integrand itself* becomes
   integrable — the numerator vanishes fast enough on the divergent locus.
   In parametric (Feynman/Lee–Pomeransky style) representations, numerator
   monomials `x^alpha` shift the monomial exponents `e_i -> e_i + alpha_i`,
   i.e. they move the label inside the same family. A polynomial numerator is
   a fixed linear combination of such shifted labels *sharing one denominator
   sector*.
2. **Raised propagator / quasi-finite basis** (von Manteuffel, Panzer,
   Schabinger, arXiv:1411.7392): shift denominator exponents (dots,
   `m_l -> m_l - 1`) and/or the spacetime dimension so that each basis
   integral has at worst an overall `1/ep` pole from the Gamma prefactor while
   the parametric integral itself converges. In our label language these are
   *m-part shifts*, changing `f_l * val_d(G_l)` contributions to the score.

This pass implements model 1 (numerators). Model 2 (raised propagators /
dimension shift) is the designated next candidate if model 1 fails — see §6.

## 2. Definitions

Let the family have variables `x = (x_1..x_k)`, polynomials `G_1..G_p`, and
labels `(n | m) in Z^k x Z^p` denoting the relative factor
`prod x_i^(n_i) * prod G_l^(m_l)` (`n_i > 0` numerator monomial powers,
`m_l < 0` denominator powers).

- **Denominator sector (sector label).** A label `S = (n0 | m)` whose m-part
  fixes the denominator `prod G_l^(m_l)`; the n-part `n0` is the base
  numerator offset (usually `0`). All members of one finite-numerator ansatz
  share the sector.
- **SparsePoly numerator.** `N(x) = sum_{alpha in A} c_alpha * x^alpha`,
  stored as a `SparsePoly` (integer exponent tuples `alpha >= 0`, coefficients
  `c_alpha` free of the integration variables; in this pass rational
  constants). `A = supp(N)`.
- **Numerator degree.** `deg N = max_{alpha in A} |alpha|_1`. The MVP scans
  `d = 1, 2`.
- **FiniteNumeratorIntegral.** The single integrand `F_S(x) * N(x)` together
  with: the sector `S`, the numerator `N`, the full-integrand LF verdict, and
  the defining expansion (below). Exported as ONE HyperInt-ready integrand.
- **Defining expansion into monomial labels.** Exact identity, by linearity
  of the integral sign:

      M_a  -  sum_alpha c_{a,alpha} * J(S + (alpha | 0))  =  0,

  a linear "defining row" connecting the new master `M_a` to ordinary
  monomial labels. These rows embed the new masters into the existing IBP row
  system without touching the reducer core.
- **Full-integrand LF verdict.** The verdict of the strict gate applied to
  the complete integrand `F_S * N` (not to its pieces): along every candidate
  ray (and random safety-net directions) every Laurent level `<= 0` of the
  ep=0 integrand must vanish identically, exponents must be numeric at ep=0,
  and denominators must be provably positive in the bulk (same rules as
  `valuations.is_locally_finite`, spec 5.4).

## 3. What the strict toric gate can see (Lemma 1)

Along a ray `d` the sector factor expands as `F_S = t^P (c_0 + c_1 t + ...)`
with `c_0 != 0` a rational function of the transverse coordinates `y`, and
each numerator monomial contributes `x^alpha = t^(alpha.d) y^alpha`. The
complete integrand is graded by the `t`-level:

    F_S * N = sum_L t^L * [ sum_{alpha, j : P + alpha.d + j = L}
                            c_alpha y^alpha c_j(y) ].

**Lemma 1 (graded lowest layer — no hidden cancellation).** Let
`v = min_{alpha in supp N} alpha.d`. The lowest populated level is
`L_min = P + v` and its coefficient is `c_0(y) * N_init(y)` where
`N_init = sum_{alpha.d = v} c_alpha y^alpha` is the initial form of `N`.
Since the layer coefficients live in an integral domain and `c_0 != 0`,
`N_init != 0`, this coefficient never vanishes. Under the STRICT RULE
(every level `<= 0` must vanish) it follows that

    F_S * N is LF on ray d
        <=>  score(S, d) + v > 0
        <=>  every monomial piece J(S + (alpha|0)), alpha in supp(N),
             is individually LF on ray d.

*Proof.* The level-`L_min` coefficient is the lowest graded piece of a product
of nonzero elements of a domain, hence nonzero; so LF requires `L_min >= 1`,
i.e. `P + alpha.d >= 1 - min(...)` for the minimizing alphas; monotonicity in
`alpha.d` extends this to all of `supp(N)`. Conversely if every piece has
positive score, every level of the sum sits at positive level. ∎

Consequences (all machine-checked by the implementation):

- Within the strict toric-ray gate there is **no** same-sector numerator
  cancellation: a polynomial numerator is admissible iff its support lies in
  the *clearing set* `A_LF(S, d) = { alpha : |alpha|_1 <= d, S + (alpha|0)
  is individually LF }`.
- The leading-cancellation linear system of Phase 2 step 5 (impose vanishing
  of all non-integrable leading coefficients of the complete integrand,
  solve for `c_alpha`) is still implemented honestly — its kernel provably
  equals `span{ e_alpha : alpha in A_LF }`, and the implementation asserts
  this equality (`lemma_consistent`) instead of assuming it.
- A "combined LF from separately divergent monomial pieces" object cannot
  exist under this gate; any test constructing one must (and does) assert its
  rejection. Genuine sub-toric numerator cancellation (numerators vanishing
  on the singular surface `A_l0(y) = 0` inside a boundary facet) is invisible
  to a toric scaling gate and belongs to boundary-aware analysis, which is
  **out of scope for this pass** by instruction.

## 4. Ray-sign obstruction (Lemma 2) and the Int2 sectors

**Lemma 2.** If every failing ray `d` of sector `S` satisfies `d <= 0`
componentwise, then for every nonzero polynomial `N` with `supp(N) >= 0` and
every degree, `v = min alpha.d <= 0`, so no polynomial numerator can make
`F_S * N` locally finite. (Scores are integers at ep=0; failing means
`score <= 0`; a cure needs `score + v >= 1` with `v <= 0`.)

Empirical status of the six certified normal-form sectors of External Int2
(probe, 12 candidate rays):

| Sector | LF | failing rays (direction, score) |
|---|---|---|
| `1/(x2*G0*G1)` | True | — |
| `1/(x2*G1*G3)` | True | — |
| `1/(G0*G3)` | True | — |
| `1/G1` | **False** | `(-1,0,0): -1`, `(-1,0,-1): 0`, `(-1,-1,0): 0`, `(-1,-1,-1): 0` |
| `1/G2` | **False** | `(-1,0,0): -1`, `(0,-1,0): 0`, `(-1,-1,0): -1`, `(-1,-1,-1): 0` |
| `x7/(G0*G3)` | True | — |

All failing rays of both non-LF sectors are componentwise `<= 0` (pure
`x -> oo` directions). **By Lemma 2, no polynomial numerator of any degree
makes `1/G1` or `1/G2` locally finite** — the divergence is lack of decay at
infinity, which polynomial numerators only worsen. The MVP still runs the
honest per-sector search (degrees 1, 2, all sectors near the normal form,
including deeper-m sectors where mixed-sign failing rays may occur and
numerators may genuinely help).

## 5. Phases

- **Phase 2 (module `finite_numerator.py`).** Per sector: enumerate ansatz
  monomials `|alpha|_1 <= d`; build the leading-cancellation linear system
  from the exact Laurent layers of the complete integrand on every failing
  ray; solve; cross-check the kernel against the clearing set (Lemma 1);
  verify every accepted candidate with the independent full-integrand gate;
  emit `FiniteNumeratorIntegral` objects (single-integrand export + defining
  expansion). Report per sector: failing rays, kernel dimension, clearing
  monomials, `numerator_cure_impossible_any_degree` (Lemma 2 flag).
- **Phase 3 (feasibility).** Because every admissible `M_a` expands into
  individually-LF monomial labels (Lemma 1), span feasibility through the
  new masters equals LF-constrained feasibility on the **numerator-extended
  label box** (n-range raised by the ansatz degree). Reuse the Method.1
  machinery (`lf_reduction_feasible_mod_p`) on box `[0,2]^3 x [-2,0]^4`
  against the Int2 target; defining rows are exactly the monomial expansions.
  Statuses: Feasible / Obstructed / BadSpecialization per (sample, prime).
- **Honest failure handling.** If no finite-numerator basis reaches the
  target: report `NoFiniteNumeratorWithinDegree` (or the Lemma 2 impossibility
  where it applies) and name the next candidates — raised propagator /
  quasi-finite basis (m-shifts *do* change the failing-ray scores via
  `f_l * val_d(G_l)`, e.g. `val_{(-1,0,0)}(G1) = 0` explains why G1-dots
  alone cannot fix `(-1,0,0)` either, so mixed dots on growing polynomials
  G0/G3 are the lever) or dimension shift. A partial non-LF result is never
  returned as Success.

## 6. Acceptance criteria

A candidate is accepted only if all of:

1. the complete integrand `F_S * N` independently passes
   `is_locally_finite = True` (not just the linear solver);
2. no cancellation after integration is used anywhere;
3. it is exported as ONE numerator-decorated integrand
   (`(N) * F_S` Wolfram/HyperInt text) plus its defining expansion rows.

The Int2 problem counts as solved by this route only if the target is
expressible through such masters alone with reconstructed rational
coefficients and a verified row-span certificate.
