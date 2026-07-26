# External Int2 Method.11c — offline reconstruction audit

## Verdict

The Stage-4 failure was **not** caused by too few sample points.  It was caused by
misclassifying the support-changing point

\[
(\epsilon,r)=\left(6,\frac{57}{11}\right)
\]

as a genuine coefficient zero and zero-filling it into the interpolation table.

All 152 records have the same RREF rank 26984.  However, at this one point the
normal-form basis specializes: one term disappears **and the other three coefficients
also differ from the generic rational coefficient functions**.  Thus it is not a
special zero of one coefficient; the whole sample is a basis/pivot specialization and
must be excluded from generic reconstruction.

## Data audit

- cache records: 152;
- distinct parameter samples: 38;
- primes per sample: 4;
- generic rank: 26984 for all records;
- modal support: four LF labels;
- full-support samples: 37;
- support-changing samples: one;
- stable modular comparisons: no mismatches;
- support-changing sample: all 16 coefficient/prime comparisons disagree with the
  generic functions (including the missing coefficient).

## Reconstructed generic coefficients

For label order `(x2,x5,x7,G0,G1,G2,G3)`:

\[
C_{(-1,0,0,-1,-1,0,0)}=
-\frac{2\epsilon r-\epsilon+1}{6r(3\epsilon+1)},
\]

\[
C_{(-1,0,0,0,-1,0,-1)}=
-\frac{(\epsilon-1)^2(r+1)}{6\epsilon r(3\epsilon+1)},
\]

\[
C_{(0,0,0,-1,0,0,-1)}=
\frac{(\epsilon-1)(4\epsilon^2r+\epsilon^2-2\epsilon r-1)}
{6\epsilon^2r(3\epsilon+1)},
\]

\[
C_{(0,0,1,-1,0,0,-1)}=
\frac{(\epsilon-1)(2\epsilon-1)(2\epsilon r+\epsilon+1)}
{6\epsilon^2r(3\epsilon+1)}.
\]

The existing dense interpolator reconstructs these functions in about four seconds
once the support-changing sample is removed.  The accepted degree pairs are `(2,2)`,
`(3,3)`, `(4,4)`, `(4,4)`, with 35 fit points and two holdouts.

## Required core correction

Pre-reconstruction support loss at generic rank cannot be classified as a true
`special_zero` from support agreement across primes alone.  Conservative semantics:

1. exclude every support-changing sample from fitting;
2. reconstruct generic coefficient functions from modal-support samples;
3. evaluate the reconstructed functions at the excluded sample;
4. classify it as a genuine special zero only if every present coefficient agrees and
   every missing coefficient evaluates to zero;
5. otherwise classify it as a basis/pivot specialization.

## Remaining acceptance step

The four-term formula is a candidate, not yet a final `Success`.  It requires an
independent off-sample row-span certificate using the chamber-valid row set.  Final
acceptance also requires all four labels to remain LF=True at epsilon=0 (already true in
all records) and, preferably, a direct numerical original-vs-LF-RHS comparison at
`epsilon=-3/5, r=1`.
