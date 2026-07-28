# Toric surface validation (External Int2 Method.13)

Exact boundary criteria for the two IBP row families over the positive orthant
`R_+^N`, with the integrand

```
F(x) = prod_i x_i^{e_i} * prod_l G_l(x)^{f_l},
```

`G_l` polynomials with positive-orthant Newton data, `e_i = a_i + n_i`,
`f_l = b_l + m_l` (family base + integer label). Package ray convention:
a boundary ray is a primitive integer direction `d`, approached as
`x_i ~ lambda^{d_i}`, `lambda -> 0` (so `d_i > 0` means `x_i -> 0`,
`d_i < 0` means `x_i -> Infinity`).

Throughout, `score(H, d)` is the **full-measure scaling score** of an
integrand `H = x^p F`:

```
score(x^p F, d) = sum_i d_i (e_i + p_i + 1) + sum_l f_l * val_d(G_l),
val_d(G_l) = min_{a in supp G_l} (a . d)     (tropical valuation).
```

`score` is positively homogeneous and piecewise linear in `d`; it is linear on
every cone of the common refinement of the normal fans of the `G_l` (the only
non-linearity is the argmin switch inside `val_d`). Hence:

**Ray-sufficiency lemma.** If `score(H, d) > 0` for every ray `d` of a
complete fan on whose cones `score` is linear, then `score(H, d) > 0` for
every nonzero integer direction. The set produced by
`complete_polyhedral_rays` (all pairwise-independent intersections of `N-1`
wall normals `(a - a') . d = 0`, `a, a'` monomials of one `G_l`, plus the
coordinate rays) contains all extreme rays of such a fan — every extreme ray
of a complete fan in `R^N` lies on at least `N-1` independent walls — so
positivity on that finite set certifies positivity everywhere. A superset only
adds redundant (never wrong) conditions. If the enumeration exceeds its
budget, the set is not provably complete and every would-be `True` must
degrade to `Unknown`.

## 1. Setup: log coordinates

Substitute `x_k = e^{u_k}`. Then `R_+^N` becomes `R^N` and, for any polynomial
vector field `W = (W_1, ..., W_N)`,

```
Integral_{R_+^N} sum_k d/dx_k (W_k) dx  =  Integral_{R^N} sum_k d/du_k (V_k) du,
V_k(u) = W_k(x) / x_k * prod_j x_j .
```

*(Check: `d/dx_k = e^{-u_k} d/du_k`, `dx = e^{sum u} du`, and
`e^{-u_k} [d/du_k W_k] e^{sum u} = d/du_k (W_k e^{sum_{j != k} u_j})`.)*

By the divergence theorem in `u`-space the integral equals the total flux of
`V` through the sphere at infinity. The direction "`-d/|d|` at infinity in
`u`" is exactly the package ray `d`: `u = -R d + O(1)`, `R -> +Infinity`
means `x_k = lambda^{d_k}` with `lambda = e^{-R} -> 0`. The area element grows
only polynomially in `R = log(1/lambda)`, so a flux density `~ lambda^s`
vanishes iff `s > 0` **strictly**; `s = 0` is a finite or log-divergent flux
(Failed), `s < 0` diverges.

## 2. Coordinate primitive `d/dx_i (P F)`, `P = x^p`

Here `W = P F e_i` (only component `i` nonzero), so
`V = (P F / x_i) * prod_j x_j * e_i` and the flux density in direction `-d`
at infinity is

```
V . (-d)/|d| = -(d_i/|d|) * (P F / x_i) * prod_j x_j .
```

If `d_i = 0` the density vanishes identically — the field is tangential to
that face. For `d_i != 0`, under `x = lambda^d y`:

```
(P F / x_i) prod_j x_j  ~  lambda^{ val_d(P F) + sum_j d_j - d_i }
                        =  lambda^{ score(P F, d) - d_i } .
```

**Criterion (coordinate).** The primitive is surface-free iff for every ray
`d` of the complete set with `d_i != 0`:

```
surface_score(d, i) = score(P F, d) - d_i > 0     (strict).
```

Sanity: on `d = +e_i` this is `score - 1 = p_i + e_i + sum_l f_l
min_pow_i(G_l) > 0` — the historical `exp_zero > 0`; on `d = -e_i` it is
`-(p_i + e_i + sum_l f_l max_pow_i(G_l)) > 0` — the historical
`exp_inf < 0`. The mixed rays (`d_i != 0`, some other `d_k != 0`) are the new
content: they demand that the **transverse `(N-1)`-dimensional facet
integral** of `P F` decays, which pointwise decay in `x_i` alone does not
give. Method.12R counterexample (External Int2, label `[-1,0,0,-3,0,-3,0]`,
`i = x5`, `P = x5^2`): both component exponents pass, but on
`d = (1,-1,0)` the surface score is `-1` — the facet integral at
`x5 -> Infinity` blows up along `x2 -> 0` and the "vanishing" boundary term
is actually `0 * Infinity`.

## 3. Vector/tangent primitive `div(Q F)`

Now `W_k = Q_k F`, so `V_k = (Q_k F / x_k) prod_j x_j` and

```
V . (-d) = -(prod_j x_j) * F * N_d(x),      N_d(x) = sum_k d_k Q_k(x) / x_k .
```

`N_d` is the **normal component** of `Q` for the face of `d`; components with
`d_k = 0` are tangential and drop out. Crucially `N_d` must be assembled as a
single Laurent polynomial with exact coefficient arithmetic: monomials
contributed by different `Q_k` can cancel exactly (e.g. `Q = (x, -y)`,
`d = (-1,-1)`: `N_d = -1 + 1 = 0`, zero flux although each component in
isolation looks marginal). Only exact cancellation counts — a cancellation
that holds only modulo a prime or numerically is not a cancellation.

For each surviving monomial `x^m` of `N_d` (`m` may have a `-1` entry):

```
(prod_j x_j) F x^m  ~  lambda^{ score(x^m F, d) } .
```

**Criterion (flux).** `div(Q F)` is surface-free iff for every ray `d` of the
complete set, every monomial of the exactly-assembled `N_d` with nonzero
coefficient satisfies `score(x^m F, d) > 0` (strict). Since
`min_m score > 0` iff all `score > 0`, this is exactly the statement that the
initial form of the normal flux is suppressed.

Consistency: decomposing `div(Q F) = sum_i d/dx_i (Q_i F)` term by term and
applying §2 to each monomial `x^c` of `Q_i` gives
`score(x^c F, d) - d_i = score(x^{c - e_i} F, d)` — the per-term version of
the same criterion. §3 is strictly sharper only through the exact
cancellation in `N_d` and the dropping of tangential (`d_i = 0`) components;
it is never more permissive on a surviving term.

## 4. Strictness rules (Method.13 item 8)

Exact arithmetic only, no floating signs. For each tested ray:

- decided negative **or exactly zero** score -> `Failed` (a zero score is a
  finite/log boundary flux, not a vanishing one; STRICT RULE, same as the LF
  gate);
- unresolved symbolic sign (non-regulator parameter, or a regulator missing
  from the chamber point) -> `Unknown`;
- ray-enumeration budget exhausted -> the candidate set is not provably
  complete -> `Unknown`, never `True`;
- the random safety net (LF gate) may still find additional `False`
  witnesses, but can never certify `True`.

Both `SurfacePolicy` modes (`limit`: signs read as `epsilon -> 0^-/+` limits;
`chamber`: exact rational substitution) run over the SAME complete ray set —
the policy chooses only how the sign of the score expression is read.

## 5. What changed in Method.13

| | pre-13 (defective) | Method.13 |
|---|---|---|
| coordinate | component-local `exp_zero > 0`, `exp_inf < 0` only | all complete rays with `d_i != 0`, `score - d_i > 0` |
| tangent | every monomial of every `Q_i` on heuristic rays, no cancellation, `d_i = 0` included | exact `N_d` with cancellation, complete rays, tangential components dropped |
| rays | heuristic 12-ray set (Int2), joint faces missing | complete polyhedral set (18 for Int2), budget-guarded |
| Unknown | on undecidable signs | also on incomplete ray enumeration |

The pre-13 coordinate acceptance admitted invalid rows (Method.12R: 21% of
the production limit-policy Level-0 rows, 43% of the Method.10 chamber-only
rows have no valid justification); the pre-13 tangent test was additionally
over-strict on cancelling fields and under-scoped on mixed faces.
