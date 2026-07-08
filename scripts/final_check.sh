#!/usr/bin/env bash
# Release sanity: fast suite + lint; --heavy adds the ~25-30 min D4 acceptance runs.
set -euo pipefail
cd "$(dirname "$0")/.."
python -m pytest
ruff check .
if [[ "${1:-}" == "--heavy" ]]; then
  RUN_D4_FULL=1 python -m pytest tests/test_d4_vertical.py tests/test_d4_cli_e2e.py
fi
echo "final_check: OK"
