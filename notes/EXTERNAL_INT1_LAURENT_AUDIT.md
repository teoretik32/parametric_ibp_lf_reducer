# External Int1 — Laurent structure audit

**Status: PASSED — no mismatch at any Laurent order through `ep^0`.**

- Script: `scripts/audit_external_int1_laurent.py` (standalone; does **not** import
  `parametric_ibp_lf_reducer` and does not touch the certified reduction or its
  coefficients — pure symbolic/numeric cross-examination).
- Machine-readable results: `validation/external_int1_laurent_audit.json`.
- Precision: `mp.dps = 45`.

## Setup

With `ep` the dimensional regulator and the two Feynman-parameter integrals taken
over the positive orthant:

```
J1(ep) = ∫∫ (1+x2)^ep (1+x6)^(ep-1) (1+x2+x6)^(ep-2) dx2 dx6
J2(ep) = ∫∫ (1+x2)^ep (1+x6)^ep     (1+x2+x6)^(ep-3) dx2 dx6

A(ep) = (4·ep − 1) / (3·(3·ep + 1))
B(ep) = (ep − 2)·(5·ep − 2) / (3·ep·(3·ep + 1))

C(ep)    = A·J1 + B·J2
Full(ep) = P(ep) · C(ep)        (units of 1/(s·t²))
```

where `P(ep)` is the external Γ/kinematic prefactor of External Int1 (applied
outside the reducer; its closed form is embedded in the script and cross-checked
numerically, see below).

Target Laurent series to audit against:

```
Full(ep) = 1/ep^4 − (π²/12)/ep² − (43·ζ3/6)/ep − π⁴/180 + O(ep)
```

## Method

1. **Inner-integral reduction to a 2F1 kernel.** Substituting `w = 1/(1+x)`:
   `∫₀^∞ (1+x)^a (1+x+y)^c dx = ∫₀¹ w^(−a−c−2) (1+y·w)^c dw`, which gives
   - J1 inner: `2F1(2−ep, 1−2ep; 2−2ep; −y)/(1−2ep)`
   - J2 inner: `2F1(3−ep, 2−2ep; 3−2ep; −y)/(2−2ep)`

   Validated against direct 2-D quadrature at `ep = 0.03`
   (|diff| = 1.85e-11 for J1, 5.71e-12 for J2 — at the 2-D quadrature's own accuracy).

2. **Exact prefactor Laurent series** via sympy (`series` + rational ζ-substitutions);
   the reduced closed form is additionally cross-checked numerically at
   `ep = 0.0173` and at complex `ep = 0.011 + 0.007i` (rel. diff ≤ 1.7e-45):

   ```
   P(ep) = (3/2)/ep³ + (9/2)/ep² − (π²/4)/ep + (−16·ζ3 − 3·π²/4)
           + (−48·ζ3 − 19·π⁴/80)·ep + O(ep²)
   ```

3. **Cauchy-circle Taylor extraction** for J1 (orders 0–3) and J2 (orders 0–4):
   circle `|ep| = 1/32`, 44 nodes, conjugate symmetry halves the evaluations;
   per-node 1-D quadrature error ≤ 1e-62; max imaginary residue of the extracted
   coefficients ≤ 1.04e-42.

4. **PSLQ identification** in the weight-graded basis
   `{1, π², ζ3, π⁴, ζ5, π²·ζ3}` with per-order weight cutoff `w ≤ k+1`
   (`tol = 1e-34`, `maxcoeff = 1e8`). Largest residual: 2.56e-41 (`J2[ep^4]`);
   acceptance gate 1e-32. The weight-5 slots (ζ5, π²·ζ3) were offered to
   `J2[ep^4]` and came back with zero coefficients.

5. **Assembly and comparison.** `C = A·J1 + B·J2` as exact series, multiplied by
   the exact `P(ep)`; per-order comparison against the target; numeric sanity
   check at `ep = 0.02` (direct minus truncated-Laurent = 0.103, consistent with
   the first omitted order `O(ep)`).

## Identified series (all PSLQ-certified, residuals ≤ 2.6e-41)

```
J1(ep) = 1 + (3 + π²/6)·ep + (7 + π²/6 + 7·ζ3)·ep² + (15 + π²/6 + 7·ζ3 + π⁴/3)·ep³ + O(ep⁴)

J2(ep) = 1/2 + (7/4)·ep + (35/8 + π²/12)·ep² + (155/16 + π²/8 + 7·ζ3/2)·ep³
         + (651/32 + 7·π²/48 + 21·ζ3/4 + π⁴/6)·ep⁴ + O(ep⁵)

C(ep)  = (2/3)/ep − 2 + (6 + π²/18)·ep + (−18 − π²/6 + 7·ζ3/3)·ep²
         + (54 + π²/2 − 7·ζ3 + π⁴/9)·ep³ + O(ep⁴)      (trustworthy through ep³)
```

## Per-order mismatch report

| Order   | Audit result   | Target         | Verdict |
|---------|----------------|----------------|---------|
| `ep^-4` | `1`            | `1`            | OK      |
| `ep^-3` | `0`            | `0`            | OK      |
| `ep^-2` | `−π²/12`       | `−π²/12`       | OK      |
| `ep^-1` | `−43·ζ3/6`     | `−43·ζ3/6`     | OK      |
| `ep^0`  | `−π⁴/180`      | `−π⁴/180`      | OK      |

**Conclusion:** the Laurent structure of External Int1 produced by the certified
reduction is independently confirmed at every order through the finite part.
No mismatch found; no changes to reducer code required (and none were made).

## Methodological pitfall found during the audit (audit-side, not reducer-side)

The first version of the audit gated the PSLQ basis by weight `w ≤ k` at Taylor
order `k`. That is wrong for these integrals: the two leftover parametric
integrations can raise the transcendental weight by one, e.g.
`J1'(0) = 3 + ζ2` carries weight 2 at order 1 (PSLQ then correctly "failed"
rather than misidentify). The cutoff is now `w ≤ k+1` and the numeric basis is
derived directly from the symbolic basis list, so the two cannot drift apart.
