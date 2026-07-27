<|
  "Status" -> "Success",
  "TargetLabel" -> {0,0,0,0,0,0,0},
  "AllLocallyFinite" -> True,
  "Terms" -> {
    <| "Integrand" -> 1/(x2*G0*G1), "Coefficient" -> -(2*ep*r - ep + 1)/(6*r*(3*ep + 1)), "LocallyFinite" -> True |>,
    <| "Integrand" -> 1/(x2*G1*G3), "Coefficient" -> -(ep - 1)^2*(r + 1)/(6*ep*r*(3*ep + 1)), "LocallyFinite" -> True |>,
    <| "Integrand" -> 1/(G0*G3), "Coefficient" -> (ep - 1)*(4*ep^2*r + ep^2 - 2*ep*r - 1)/(6*ep^2*r*(3*ep + 1)), "LocallyFinite" -> True |>,
    <| "Integrand" -> x7/(G0*G3), "Coefficient" -> (ep - 1)*(2*ep - 1)*(2*ep*r + ep + 1)/(6*ep^2*r*(3*ep + 1)), "LocallyFinite" -> True |>
  },
  "Diagnostics" -> <|
    "FormalSuccess" -> True,
    "ReconstructionVerified" -> True,
    "IndependentValidationPassed" -> True,
    "NumTerms" -> 4,
    "NonLFTerms" -> {},
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
