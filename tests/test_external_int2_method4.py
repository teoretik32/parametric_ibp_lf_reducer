"""Tests for ``scripts/run_external_int2_method4.py`` (Method.4 completeness audit).

Tiny inputs only: the script's default Int2 input is exercised solely in the
``RUN_EXTERNAL_INT2=1``-gated integration test — everything else runs on
``examples/tiny_success_input.wl.txt`` with hard caps, so the suite stays fast and never
launches anything resembling a production run.
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

SCRIPT = REPO_ROOT / "scripts" / "run_external_int2_method4.py"
TINY_INPUT = REPO_ROOT / "examples" / "tiny_success_input.wl.txt"
FINAL_JSON = REPO_ROOT / "validation" / "external_int2_method4.json"


def _load_script():
    spec = importlib.util.spec_from_file_location("run_external_int2_method4", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def script():
    return _load_script()


# --- usage errors -----------------------------------------------------------------------


def test_missing_input_is_usage_error(script, tmp_path):
    rc = script.main(
        ["--input", str(tmp_path / "nope.wl.txt"), "--out", str(tmp_path / "out.json")]
    )
    assert rc == 2


def test_bad_caps_are_usage_errors(script, tmp_path):
    rc = script.main(["--samples", "0", "--out", str(tmp_path / "out.json")])
    assert rc == 2
    rc = script.main(["--primes", "0", "--out", str(tmp_path / "out.json")])
    assert rc == 2


def test_bad_blocks_are_usage_errors(script, tmp_path):
    for bad in ("3", "3,3,3", "a,b", "-1,2", ""):
        # ``--blocks=<val>`` form: a plain ``-1,2`` argument would be eaten by argparse
        # itself as an unknown option (SystemExit) before ``_parse_blocks`` runs.
        rc = script.main([f"--blocks={bad}", "--out", str(tmp_path / "out.json")])
        assert rc == 2, bad


# --- _parse_blocks ----------------------------------------------------------------------


def test_parse_blocks(script):
    assert script._parse_blocks("3,3;4,4") == [(3, 3), (4, 4)]
    assert script._parse_blocks(" 1 , 2 ; ; 0,0 ") == [(1, 2), (0, 0)]
    with pytest.raises(ValueError):
        script._parse_blocks("3;4")
    with pytest.raises(ValueError):
        script._parse_blocks("-1,1")
    with pytest.raises(ValueError):
        script._parse_blocks(";")


# --- phase s (tiny) ---------------------------------------------------------------------


def test_phase_s_only_tiny(script, tmp_path):
    out = tmp_path / "m4_s.json"
    rc = script.main(
        [
            "--input",
            str(TINY_INPUT),
            "--out",
            str(out),
            "--phases",
            "s",
            "--blocks",
            "1,1",
        ]
    )
    assert rc == 0
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["phases"] == "s"
    assert payload["level_name"] == "deep"
    assert "never reported as an LF basis" in payload["scope_note"]
    (entry,) = payload["phase_s"]["blocks"]
    assert entry["block"] == [1, 1]
    assert entry["status"] == "ok"
    assert entry["n_fields"] >= 0 and entry["elapsed_s"] >= 0
    assert "phase_d" not in payload


def test_phase_s_budget_skip(script, tmp_path):
    out = tmp_path / "m4_budget.json"
    rc = script.main(
        [
            "--input",
            str(TINY_INPUT),
            "--out",
            str(out),
            "--phases",
            "s",
            "--blocks",
            "1,1;1,1;1,1",
            "--field-budget",
            "0",
        ]
    )
    assert rc == 0
    payload = json.loads(out.read_text(encoding="utf-8"))
    statuses = [b["status"] for b in payload["phase_s"]["blocks"]]
    # the first block always runs; once over the zero budget the rest are recorded skips
    assert statuses[0] == "ok"
    assert "skipped_budget" in statuses


# --- phases sd (tiny) -------------------------------------------------------------------


def test_phases_sd_tiny(script, tmp_path):
    out = tmp_path / "m4_sd.json"
    rc = script.main(
        [
            "--input",
            str(TINY_INPUT),
            "--out",
            str(out),
            "--phases",
            "sd",
            "--blocks",
            "1,1",
            "--samples",
            "1",
            "--primes",
            "1",
        ]
    )
    assert rc == 0
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["n_labels"] >= 1
    assert payload["n_lf_true"] >= 0
    d = payload["phase_d"]
    assert d["n_base_rows"] > 0
    assert (
        d["n_extra_rows_offered"] == d["n_extra_rows_new"] + d["n_extra_rows_duplicate"]
    )
    assert len(d["points"]) == 1
    point = d["points"][0]
    for key in ("baseline", "enriched"):
        assert point[key]["status"] in (
            "Feasible",
            "Obstructed",
            "BadSpecialization",
        )
    allowed = (
        script.VERDICT_FEASIBLE,
        script.VERDICT_OBSTRUCTED,
        script.VERDICT_MIXED,
        script.VERDICT_INCONCLUSIVE,
    )
    assert d["baseline_verdict"] in allowed
    assert d["enriched_verdict"] in allowed
    assert d["n_points_flipped_to_feasible"] >= 0


# --- real External Int2 family (heavy, gated) -------------------------------------------


@pytest.mark.skipif(
    os.environ.get("RUN_EXTERNAL_INT2") != "1",
    reason="heavy External Int2 Method.4 audit; set RUN_EXTERNAL_INT2=1",
)
def test_int2_method4_verdicts_reproduce(script, tmp_path):
    stored = json.loads(FINAL_JSON.read_text(encoding="utf-8"))
    out = tmp_path / "m4_full.json"
    rc = script.main(
        [
            "--out",
            str(out),
            "--phases",
            "sd",
            "--blocks",
            "3,3;4,4",
            "--samples",
            str(len(stored["samples"])),
            "--primes",
            str(len(stored["primes"])),
        ]
    )
    assert rc == 0
    fresh = json.loads(out.read_text(encoding="utf-8"))
    for key in ("baseline_verdict", "enriched_verdict", "n_points_flipped_to_feasible"):
        assert fresh["phase_d"][key] == stored["phase_d"][key]
    assert [p["baseline"]["status"] for p in fresh["phase_d"]["points"]] == [
        p["baseline"]["status"] for p in stored["phase_d"]["points"]
    ]
    assert [p["enriched"]["status"] for p in fresh["phase_d"]["points"]] == [
        p["enriched"]["status"] for p in stored["phase_d"]["points"]
    ]
