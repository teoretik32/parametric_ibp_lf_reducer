<|
  "Example" -> "Example4StarCorrected",
  "Status" -> "Success",
  "Error" -> None,
  "Target" -> "15*ep*J[{0,0,0,0,0,0,0}] + 24*ep*J[{0,1,0,0,0,0,0}]",
  "LHSTerms" -> {
    <| "Label" -> {0,0,0,0,0,0,0}, "Coefficient" -> 15*ep |>,
    <| "Label" -> {0,1,0,0,0,0,0}, "Coefficient" -> 24*ep |>
  },
  "AllLocallyFinite" -> True,
  "Terms" -> {
    <| "Label" -> {1,1,0,-1,0,0,0}, "Coefficient" -> (47703*ep^3-521*ep^2-57*ep-1)/(3300*ep^2), "LocallyFinite" -> True |>,
    <| "Label" -> {1,1,0,0,0,-1,0}, "Coefficient" -> (816*ep^3+881*ep^2+66*ep+1)/(3300*ep^2), "LocallyFinite" -> True |>
  },
  "CertificateStatus" -> "Passed",
  "SelectedRank" -> 9924,
  "Note" -> "Known-value-only example (docs/05): certified row-span reduction of the corrected integrand via linearity over one shared row system; no reference LF decomposition exists in the notebook."
|>
