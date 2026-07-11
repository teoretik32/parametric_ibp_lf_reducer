"""CLI.1 tests: thin ``python -m parametric_ibp_lf_reducer reduce`` wrapper.

The CLI is pure plumbing, so these tests only check wiring: exit codes, file
outputs, JSON shape. The tiny generic family runs the FULL pipeline (certificate
gate included) in ~1-2 s, so no math is mocked or bypassed.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from parametric_ibp_lf_reducer.cli import EXIT_FAILURE, EXIT_SUCCESS, EXIT_USAGE, main
from parametric_ibp_lf_reducer.result import (
    FAILURE_PARSER_NEEDS_EXPLICIT_FAMILY,
    STATUS_SUCCESS,
)

REPO_ROOT = Path(__file__).resolve().parents[1]

TINY_SUCCESS_TEXT = """
IBPInput = <|
  "Variables" -> {u, v},
  "Parameters" -> {ep},
  "Regulators" -> {ep},
  "Domain" -> "PositiveOrthant",
  "Polynomials" -> <| "P0" -> 1 + u, "P1" -> 1 + v |>,
  "MonomialExponents" -> <| u -> -1 - ep, v -> ep |>,
  "PolynomialExponents" -> <| "P0" -> -1 + ep, "P1" -> -2 - ep |>,
  "TargetMultiplier" -> 1
|>;
"""

# Same document without the explicit-family keys -> honest typed failure, not a crash.
NEEDS_FAMILY_TEXT = """
IBPInput = <|
  "Variables" -> {u, v},
  "Parameters" -> {ep},
  "Regulators" -> {ep},
  "Domain" -> "PositiveOrthant",
  "Integrand" -> 1
|>;
"""


@pytest.fixture()
def tiny_input(tmp_path):
    p = tmp_path / "input.m"
    p.write_text(TINY_SUCCESS_TEXT, encoding="utf-8")
    return p


# --- success path ----------------------------------------------------------------------------
def test_reduce_success_writes_out_and_json_and_exits_zero(tiny_input, tmp_path):
    out = tmp_path / "result.m"
    diag = tmp_path / "diagnostics.json"
    rc = main(["reduce", str(tiny_input), "--out", str(out), "--diagnostics-json", str(diag)])
    assert rc == EXIT_SUCCESS

    text = out.read_text(encoding="utf-8")
    assert '"Status" -> "Success"' in text
    assert '"Terms" ->' in text
    assert "**" not in text  # Wolfram-style powers only

    payload = json.loads(diag.read_text(encoding="utf-8"))
    assert payload["status"] == STATUS_SUCCESS
    assert payload["success"] is True
    assert payload["error"] is None
    assert payload["terms"], "success result must carry master terms"
    for term in payload["terms"]:
        assert set(term) == {"label", "coefficient", "integrand", "locally_finite"}
    inner = payload["diagnostics"]
    assert inner["formal_success"] is True
    assert inner["n_terms"] == len(payload["terms"])


def test_reduce_success_stdout_default(tiny_input, capsys):
    rc = main(["reduce", str(tiny_input)])
    assert rc == EXIT_SUCCESS
    captured = capsys.readouterr()
    assert '"Status" -> "Success"' in captured.out


def test_reduce_accepts_option_overrides(tiny_input, tmp_path):
    out = tmp_path / "result.m"
    rc = main(
        [
            "reduce",
            str(tiny_input),
            "--out",
            str(out),
            "--max-ibp-degree",
            "2",
            "--min-valid-records",
            "1",
        ]
    )
    assert rc == EXIT_SUCCESS
    assert '"Status" -> "Success"' in out.read_text(encoding="utf-8")


def test_reduce_accepts_rref_backend_flag(tiny_input, tmp_path):
    # Perf.11: backend selection only — identical results, so Success either way.
    out = tmp_path / "result.m"
    rc = main(["reduce", str(tiny_input), "--out", str(out), "--rref-backend", "dict"])
    assert rc == EXIT_SUCCESS
    assert '"Status" -> "Success"' in out.read_text(encoding="utf-8")


def test_reduce_accepts_rref_backend_auto(tiny_input, tmp_path):
    # Perf.12: "auto" is selection-only plumbing — the tiny input resolves to the
    # dict backend, so the Wolfram-style output is byte-identical to an explicit
    # dict run and the certificate gate still passes.
    out_dict = tmp_path / "result_dict.m"
    out_auto = tmp_path / "result_auto.m"
    diag = tmp_path / "diagnostics_auto.json"

    rc = main(["reduce", str(tiny_input), "--out", str(out_dict), "--rref-backend", "dict"])
    assert rc == EXIT_SUCCESS
    rc = main(
        [
            "reduce",
            str(tiny_input),
            "--out",
            str(out_auto),
            "--rref-backend",
            "auto",
            "--diagnostics-json",
            str(diag),
        ]
    )
    assert rc == EXIT_SUCCESS

    auto_text = out_auto.read_text(encoding="utf-8")
    assert '"Status" -> "Success"' in auto_text
    assert auto_text == out_dict.read_text(encoding="utf-8")

    payload = json.loads(diag.read_text(encoding="utf-8"))
    assert payload["status"] == STATUS_SUCCESS
    assert payload["certificate_status"] == "Passed"


def test_reduce_rejects_unknown_rref_backend(tiny_input, capsys):
    with pytest.raises(SystemExit) as exc:
        main(["reduce", str(tiny_input), "--rref-backend", "nonsense"])
    assert exc.value.code == EXIT_USAGE
    assert "invalid choice" in capsys.readouterr().err


# --- failure paths ---------------------------------------------------------------------------
def test_reduce_typed_failure_exits_one_with_reason(tmp_path, capsys):
    inp = tmp_path / "needs_family.m"
    inp.write_text(NEEDS_FAMILY_TEXT, encoding="utf-8")
    out = tmp_path / "result.m"
    diag = tmp_path / "diagnostics.json"
    rc = main(["reduce", str(inp), "--out", str(out), "--diagnostics-json", str(diag)])
    assert rc == EXIT_FAILURE

    # Result and diagnostics are still written: failures are honest outputs, not crashes.
    text = out.read_text(encoding="utf-8")
    assert '"Status" -> "Failure"' in text
    assert FAILURE_PARSER_NEEDS_EXPLICIT_FAMILY in text
    payload = json.loads(diag.read_text(encoding="utf-8"))
    assert payload["status"] == FAILURE_PARSER_NEEDS_EXPLICIT_FAMILY
    assert payload["success"] is False
    assert FAILURE_PARSER_NEEDS_EXPLICIT_FAMILY in capsys.readouterr().err


def test_reduce_malformed_document_exits_two(tmp_path, capsys):
    inp = tmp_path / "garbage.m"
    inp.write_text("this is not an association", encoding="utf-8")
    rc = main(["reduce", str(inp)])
    assert rc == EXIT_USAGE
    assert "malformed input document" in capsys.readouterr().err


def test_reduce_missing_input_file_exits_two(tmp_path, capsys):
    rc = main(["reduce", str(tmp_path / "does_not_exist.m")])
    assert rc == EXIT_USAGE
    assert "cannot read input file" in capsys.readouterr().err


def test_help_exits_zero():
    with pytest.raises(SystemExit) as exc:
        main(["--help"])
    assert exc.value.code == 0


# --- python -m entry point -------------------------------------------------------------------
def test_python_dash_m_entry(tiny_input, tmp_path):
    out = tmp_path / "result.m"
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "parametric_ibp_lf_reducer",
            "reduce",
            str(tiny_input),
            "--out",
            str(out),
        ],
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
        timeout=300,
    )
    assert proc.returncode == EXIT_SUCCESS, proc.stderr
    assert '"Status" -> "Success"' in out.read_text(encoding="utf-8")
