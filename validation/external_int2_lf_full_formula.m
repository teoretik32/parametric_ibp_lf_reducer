(* External Int2: certified four-term LF-basis reduction (Method.11/.12A). *)
(* Method.11 Phase C: independent row-span certificate for the reconstructed four-term LF reduction relation at NEW off-sample (ep, r) points. 'Certified' claims only that the claimed relation vector lies in the chamber-policy row span at the listed exact rational points modulo the listed primes, with the generic rank reproduced. It does NOT claim the analytic Laurent value; that cross-check is a separate later phase. *)
(* Variables: x2, x5, x7; parameters: ep, r. *)
(* The equality below holds AFTER integration over the positive orthant; *)
(* no pointwise integrand identity is implied. *)
G0 = x2 + 1;
G1 = x5 + 1;
G2 = x7 + 1;
G3 = r*x2*x5 + x2*x7 + x7 + 1;

J[f_] := Inactive[Integrate][f, {x2, 0, Infinity}, {x5, 0, Infinity}, {x7, 0, Infinity}];

C1 = -(2*ep*r - ep + 1)/(6*r*(3*ep + 1));  (* label [-1, 0, 0, -1, -1, 0, 0], LocallyFinite -> True *)
LFIntegrand1 = x2^(ep) * G0^(ep - 1) * G1^(ep - 1) * G2^(-ep - 1) * G3^(ep - 1);
C2 = -(ep - 1)^2*(r + 1)/(6*ep*r*(3*ep + 1));  (* label [-1, 0, 0, 0, -1, 0, -1], LocallyFinite -> True *)
LFIntegrand2 = x2^(ep) * G0^(ep) * G1^(ep - 1) * G2^(-ep - 1) * G3^(ep - 2);
C3 = (ep - 1)*(4*ep^2*r + ep^2 - 2*ep*r - 1)/(6*ep^2*r*(3*ep + 1));  (* label [0, 0, 0, -1, 0, 0, -1], LocallyFinite -> True *)
LFIntegrand3 = x2^(ep + 1) * G0^(ep - 1) * G1^(ep) * G2^(-ep - 1) * G3^(ep - 2);
C4 = (ep - 1)*(2*ep - 1)*(2*ep*r + ep + 1)/(6*ep^2*r*(3*ep + 1));  (* label [0, 0, 1, -1, 0, 0, -1], LocallyFinite -> True *)
LFIntegrand4 = x2^(ep + 1) * x7 * G0^(ep - 1) * G1^(ep) * G2^(-ep - 1) * G3^(ep - 2);

(* target label [0, 0, 0, 0, 0, 0, 0] *)
TargetIntegrand = x2^(ep + 1) * G0^(ep) * G1^(ep) * G2^(-ep - 1) * G3^(ep - 1);

ReductionIdentity =
  J[TargetIntegrand] == C1*J[LFIntegrand1] + C2*J[LFIntegrand2] + C3*J[LFIntegrand3] + C4*J[LFIntegrand4];

(* Method.12A physical wrapper (dimensionful External Int2 normalization). *)
r = s/t;

ExternalPrefactor2 =
  Exp[2*ep*EulerGamma]
  *t^(-3-ep)
  *Gamma[1-ep]*Gamma[-ep]^3*Gamma[ep]/
   (Gamma[-1-3*ep]*Gamma[-2*ep]);

FullPhysicalReduction =
  ExternalPrefactor2 * ReductionIdentity /. r -> s/t;
