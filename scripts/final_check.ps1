# Release sanity: fast suite + lint; -Heavy adds the ~25-30 min D4 acceptance runs.
param([switch]$Heavy)
$ErrorActionPreference = "Stop"
Set-Location (Join-Path $PSScriptRoot "..")
python -m pytest
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
ruff check .
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
if ($Heavy) {
    $env:RUN_D4_FULL = "1"
    python -m pytest tests/test_d4_vertical.py tests/test_d4_cli_e2e.py
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
}
Write-Host "final_check: OK"
