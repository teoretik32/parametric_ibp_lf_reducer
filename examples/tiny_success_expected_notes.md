# Expected behavior: `tiny_success_input.wl.txt`

A minimal 2-variable explicit family that runs the FULL pipeline (row generation,
modular reduction, reconstruction, independent validation, row-span certificate)
in ~1-2 seconds. This is the same document used by `tests/test_cli.py`.

## Command

```
python -m parametric_ibp_lf_reducer reduce examples/tiny_success_input.wl.txt \
    --out result.m --diagnostics-json diagnostics.json
```

## Expected outcome

- exit code `0`
- `result.m` (Wolfram-like association text):
  - `"Status" -> "Success"`
  - `"TargetLabel" -> {0,0,0,0}`
  - `"AllLocallyFinite" -> True`
  - exactly one term: `<| "Integrand" -> u/P0, "Coefficient" -> (ep - 1)/ep, "LocallyFinite" -> True |>`
- `diagnostics.json`:
  - `"status": "Success"`, `"success": true`, `"all_locally_finite": true`
  - `"terms"`: one entry with `"label": [1, 0, -1, 0]`, `"coefficient": "(ep - 1)/ep"`,
    `"integrand": "u/P0"`, `"locally_finite": true`
  - `"certificate_status": "Passed"` (row-span certificate gate is default-ON)
  - `"diagnostics"` block with `formal_success` / `reconstruction_verified` /
    `independent_validation_passed` all `true`

The mathematical content: the target `J[{0,0,0,0}]` (pure `x^a (1+u)^{-1+ep} (1+v)^{-2-ep}`
integrand) reduces to a single locally finite master with label `{1,0,-1,0}`
(extra power of `u`, one inverse power of `P0 = 1 + u`), with rational coefficient
`(ep - 1)/ep` in the regulator.
