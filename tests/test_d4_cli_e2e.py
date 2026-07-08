"""D4.6 / Release.1 — opt-in heavy end-to-end CLI run on the D4 example document.

Runs ``python -m parametric_ibp_lf_reducer reduce`` (in-process ``cli.main``) on
``examples/d4_cli_example_input.wl.txt``, whose ``"Options"`` association carries the
FULL verified D4.4 configuration (scattered non-lattice samples, label box, preferred
masters M1..M5, off-sample rank-generic certificate points). Nothing is mocked and no
gate is bypassed: this is exactly what a user gets from the shipped example.

Heavy (~10-15 min): skips unless ``RUN_D4_FULL=1``.

The reference 5-term basis is NOT required: the recorded, certificate-verified outcome
is the 3-term LF basis {M1, M2, M3}, so we only require the terms to be a subset of it.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from parametric_ibp_lf_reducer.cli import EXIT_SUCCESS, main
from parametric_ibp_lf_reducer.result import STATUS_SUCCESS

REPO_ROOT = Path(__file__).resolve().parents[1]
D4_INPUT = REPO_ROOT / "examples" / "d4_cli_example_input.wl.txt"

# label = (x1, x2, x3, x4, G0, G1, G2); certified D4.4/D4.5 output support.
M1 = (0, 1, 1, 0, -2, -1, 0)
M2 = (1, 1, 0, 0, -2, -1, 0)
M3 = (0, 1, 1, 0, -3, -1, 0)
EXPECTED_SUPPORT = {M1, M2, M3}

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(
        not os.environ.get("RUN_D4_FULL"),
        reason="heavy ~10-15min D4 CLI end-to-end run; set RUN_D4_FULL=1 to run",
    ),
]


@pytest.fixture(scope="module")
def d4_cli_run(tmp_path_factory):
    """Run the heavy D4 CLI reduction ONCE and share (rc, result text, JSON payload)."""
    tmp = tmp_path_factory.mktemp("d4_cli_e2e")
    out = tmp / "result.m"
    diag = tmp / "diagnostics.json"
    rc = main(
        [
            "reduce",
            str(D4_INPUT),
            "--out",
            str(out),
            "--diagnostics-json",
            str(diag),
        ]
    )
    return rc, out, diag


def test_d4_cli_exits_zero(d4_cli_run):
    rc, _, _ = d4_cli_run
    assert rc == EXIT_SUCCESS


def test_d4_cli_result_text(d4_cli_run):
    _, out, _ = d4_cli_run
    assert out.is_file()
    text = out.read_text(encoding="utf-8")
    assert '"Status" -> "Success"' in text
    assert '"AllLocallyFinite" -> True' in text
    assert '"TargetLabel" -> {0,0,0,0,0,0,0}' in text
    assert '"LocallyFinite" -> False' not in text


def test_d4_cli_diagnostics_json(d4_cli_run):
    _, _, diag = d4_cli_run
    assert diag.is_file()
    payload = json.loads(diag.read_text(encoding="utf-8"))

    assert payload["status"] == STATUS_SUCCESS
    assert payload["success"] is True
    assert payload["all_locally_finite"] is True
    assert payload["target_label"] == [0, 0, 0, 0, 0, 0, 0]

    # Row-span certificate gate: default-ON and must have genuinely passed.
    assert payload["certificate_status"] == "Passed"
    assert payload["certificate"]["n_certificate_points_failed"] == 0
    assert payload["certificate"]["n_certificate_points_passed"] >= 1

    # Certified LF basis: subset of {M1,M2,M3}; the 5-term reference basis is NOT forced.
    labels = {tuple(t["label"]) for t in payload["terms"]}
    assert labels, "Success result must carry master terms"
    assert labels <= EXPECTED_SUPPORT
    assert all(t["locally_finite"] is True for t in payload["terms"])

    inner = payload["diagnostics"]
    assert inner["formal_success"] is True
    assert inner["reconstruction_verified"] is True
    assert inner["independent_validation_passed"] is True
    assert inner["non_lf_terms"] == []
    assert inner["unknown_lf_terms"] == []
