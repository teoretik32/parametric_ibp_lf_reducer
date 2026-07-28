(* REVOKED by External Int2 Method.12R (mixed-boundary contradiction audit).                     *)
(* The four-term relation below is NOT a valid integral identity: at ep=-3/5, r=1 (inside the     *)
(* convergence chamber, where all five integrals converge) the target and the right-hand side     *)
(* differ by ~107% relative, and symbolically the RHS carries an ep^-3 pole (-1/(6 r^2)) that the *)
(* target (-2/(3 r ep^2) + O(1/ep)) cannot have. Root cause: master 4 is NOT locally finite --    *)
(* base_score == 0 on the joint infinity ray (0,-1,-1) (x5, x7 -> Infinity at fixed x2), a ray    *)
(* the old candidate-ray generator never emitted. Labels and coefficients are kept verbatim for   *)
(* provenance; nothing below is re-derived. See validation/external_int2_method12r.json and       *)
(* notes/EXTERNAL_INT2_METHOD12R_CONTRADICTION_AUDIT.md.                                          *)
<|
  "Status" -> "Revoked(Method.12R)",
  "TargetLabel" -> {0,0,0,0,0,0,0},
  "AllLocallyFinite" -> False,
  "FormalRowSpanCertificate" -> "Passed",
  "IntegralIdentityStatus" -> "Revoked",
  "SurfaceValidationStatus" -> "Failed",
  "Terms" -> {
    <| "Integrand" -> 1/(x2*G0*G1), "Coefficient" -> -(2*ep*r - ep + 1)/(6*r*(3*ep + 1)), "LocallyFinite" -> True |>,
    <| "Integrand" -> 1/(x2*G1*G3), "Coefficient" -> -(ep - 1)^2*(r + 1)/(6*ep*r*(3*ep + 1)), "LocallyFinite" -> True |>,
    <| "Integrand" -> 1/(G0*G3), "Coefficient" -> (ep - 1)*(4*ep^2*r + ep^2 - 2*ep*r - 1)/(6*ep^2*r*(3*ep + 1)), "LocallyFinite" -> True |>,
    <| "Integrand" -> x7/(G0*G3), "Coefficient" -> (ep - 1)*(2*ep - 1)*(2*ep*r + ep + 1)/(6*ep^2*r*(3*ep + 1)), "LocallyFinite" -> False |>
  },
  "Diagnostics" -> <|
    "FormalSuccess" -> True,
    "ReconstructionVerified" -> True,
    "IndependentValidationPassed" -> False,
    "NumTerms" -> 4,
    "NonLFTerms" -> {{0,0,1,-1,0,0,-1}},
    "UnknownLFTerms" -> {}
  |>
|>

(* Method.12A: certificate provenance (recorded, no re-derivation). *)
CertificateProvenance = <|
  "CertificateStatus" -> "Passed",
  "CertificateArtifact" -> "validation/external_int2_lf_certificate.json",
  "CertificatePoints" -> 4,
  "CertificateChecks" -> 8,
  "GenericRank" -> 26984,
  "SurfacePolicy" -> "convergence_chamber",
  "SurfaceChamber" -> "ep=-3/5"
|>;

(* Method.12R: contradiction audit outcome (this artifact is no longer a Success claim). *)
Method12RInvalidation = <|
  "Method12RStatus" -> "Invalidated",
  "Artifact" -> "validation/external_int2_method12r.json",
  "Notes" -> "notes/EXTERNAL_INT2_METHOD12R_CONTRADICTION_AUDIT.md",
  "NonLFMaster" -> {0,0,1,-1,0,0,-1},
  "CounterexampleRay" -> {0,-1,-1},
  "BaseScoreOnRay" -> 0,
  "MasterPole" -> "J4(ep) = -1/(r*ep) + O(1)",
  "RHSLaurentOrder" -> -3,
  "TargetLaurentOrder" -> -2,
  "NumericResidualRelative" -> "~1.07 at ep=-3/5, r=1",
  "CertificateScope" -> "formal row span only (modular, at the recorded points); NOT an integral identity"
|>;
