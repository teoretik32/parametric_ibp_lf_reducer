<|
  "Status" -> "Failure",
  "TargetLabel" -> {0,0,0,0,0,0,0},
  "AllLocallyFinite" -> False,
  "Terms" -> {
    <| "Integrand" -> 1/(x2*G0*G1), "Coefficient" -> (ep - 1)/(6*ep*r), "LocallyFinite" -> True |>,
    <| "Integrand" -> 1/(x2*G1*G3), "Coefficient" -> -(ep - 1)^2*(r + 1)/(6*ep^2*r), "LocallyFinite" -> True |>,
    <| "Integrand" -> 1/(G0*G3), "Coefficient" -> (ep - 1)*(2*ep^2*r + ep^2 - 1)/(6*ep^3*r), "LocallyFinite" -> True |>,
    <| "Integrand" -> 1/G1, "Coefficient" -> 1, "LocallyFinite" -> False |>,
    <| "Integrand" -> 1/G2, "Coefficient" -> -(ep + 1)/ep, "LocallyFinite" -> False |>,
    <| "Integrand" -> x7/(G0*G3), "Coefficient" -> (ep - 1)*(ep + 1)*(2*ep - 1)/(6*ep^3*r), "LocallyFinite" -> True |>
  },
  "Error" -> "NormalFormNotLocallyFinite",
  "ErrorDetail" -> "normal form contains non-locally-finite terms: ['{0,0,0,0,-1,0,0}', '{0,0,0,0,0,-1,0}']",
  "Diagnostics" -> <|
    "FormalSuccess" -> True,
    "ReconstructionVerified" -> True,
    "IndependentValidationPassed" -> True,
    "NumTerms" -> 6,
    "NonLFTerms" -> {{0,0,0,0,-1,0,0},{0,0,0,0,0,-1,0}},
    "UnknownLFTerms" -> {}
  |>
|>