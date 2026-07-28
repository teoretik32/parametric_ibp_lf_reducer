(* External Int2: four-term LF-basis reduction (Method.11/.12A) -- REVOKED by Method.12R. *)
(* The identity below does NOT hold: numerically at ep=-3/5, r=1 (all five integrals converge     *)
(* there) target vs RHS differ by ~107% relative, and symbolically the RHS has an ep^-3 pole      *)
(* -1/(6 r^2) while the target starts at ep^-2. LFIntegrand4 is not locally finite: its scaling   *)
(* score is exactly 0 on the joint infinity ray (0,-1,-1), so J[LFIntegrand4] = -1/(r ep) + O(1). *)
(* Kept verbatim for provenance -- see validation/external_int2_method12r.json.                   *)
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
C4 = (ep - 1)*(2*ep - 1)*(2*ep*r + ep + 1)/(6*ep^2*r*(3*ep + 1));  (* label [0, 0, 1, -1, 0, 0, -1], LocallyFinite -> False (Method.12R: score 0 on ray (0,-1,-1)) *)
LFIntegrand4 = x2^(ep + 1) * x7 * G0^(ep - 1) * G1^(ep) * G2^(-ep - 1) * G3^(ep - 2);

(* target label [0, 0, 0, 0, 0, 0, 0] *)
TargetIntegrand = x2^(ep + 1) * G0^(ep) * G1^(ep) * G2^(-ep - 1) * G3^(ep - 1);

(* REVOKED (Method.12R): this equality is false. Kept only as the audited claim. *)
ReductionIdentity =
  J[TargetIntegrand] == C1*J[LFIntegrand1] + C2*J[LFIntegrand2] + C3*J[LFIntegrand3] + C4*J[LFIntegrand4];

Method12RRefutation = <|
  "Status" -> "Refuted",
  "Point" -> <| "ep" -> -3/5, "r" -> 1 |>,
  "TargetNumeric" -> 3.9267,
  "RHSNumeric" -> -0.2891,
  "RelativeResidual" -> 1.07,
  "RHSLaurent" -> "-1/(6*r^2*ep^3) + O(1/ep^2)",
  "TargetLaurent" -> "-2/(3*r*ep^2) + O(1/ep)",
  "Artifact" -> "validation/external_int2_method12r.json"
|>;

(* Method.12A physical wrapper (dimensionful External Int2 normalization). *)
r = s/t;

ExternalPrefactor2 =
  Exp[2*ep*EulerGamma]
  *t^(-3-ep)
  *Gamma[1-ep]*Gamma[-ep]^3*Gamma[ep]/
   (Gamma[-1-3*ep]*Gamma[-2*ep]);

FullPhysicalReduction =
  ExternalPrefactor2 * ReductionIdentity /. r -> s/t;
