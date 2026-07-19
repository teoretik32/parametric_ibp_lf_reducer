"""Tests for ``scripts/run_external_int2_method1.py`` (Method.1 diagnostics).

Tiny inputs only: the script's default Int2 input is never exercised here — these tests
run the phases on ``examples/tiny_success_input.wl.txt`` with hard caps, so the suite
stays fast and never launches anything resembling a production run.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

SCRIPT = REPO_ROOT / "scripts" / "run_external_int2_method1.py"
TINY_INPUT = REPO_ROOT / "examples" / "tiny_success_input.wl.txt"


def _load_script():
    spec = importlib.util.spec_from_file_location("run_external_int2_method1", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def script():
    return _load_script()


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


def test_aggregate_verdicts(script):
    fea = script.STATUS_FEASIBLE
    obs = script.STATUS_OBSTRUCTED
    assert script._aggregate([fea, fea]) == script.VERDICT_FEASIBLE
    assert script._aggregate([obs]) == script.VERDICT_OBSTRUCTED
    assert script._aggregate([fea, obs]) == script.VERDICT_MIXED
    assert script._aggregate([]) == script.VERDICT_INCONCLUSIVE
    assert script._aggregate(["BadSpecialization"]) == script.VERDICT_INCONCLUSIVE


def test_phase_a_only_tiny(script, tmp_path):
    out = tmp_path / "m1.json"
    rc = script.main(
        [
            "--input", str(TINY_INPUT),
            "--out", str(out),
            "--phases", "a",
            "--audit-cap", "3",
            "--audit-trials", "8",
        ]
    )
    assert rc == 0
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["phases"] == "a"
    assert "phase_b" not in payload and "phase_c" not in payload
    pa = payload["phase_a"]
    assert payload["target_lf_verdict"] in (True, False, "Unknown")
    assert pa["target"]["verdict"] == payload["target_lf_verdict"]
    counts = pa["lf_counts"]
    assert sum(counts.values()) == payload["n_labels"]
    assert len(pa["detailed_non_lf_labels"]) <= 3
    assert pa["n_detailed_omitted"] == max(0, pa["n_non_lf_labels"] - 3)
    for report in [pa["target"], *pa["detailed_non_lf_labels"]]:
        assert report["n_failing_rays"] == len(report["failing_rays"])
        assert len(report["recommended_shifts"]) <= 8
    assert "never a global impossibility claim" in payload["scope_note"]


def test_phases_abc_tiny_consistent(script, tmp_path):
    out = tmp_path / "m1.json"
    rc = script.main(
        [
            "--input", str(TINY_INPUT),
            "--out", str(out),
            "--phases", "abc",
            "--samples", "2",
            "--primes", "1",
            "--audit-cap", "2",
            "--audit-trials", "8",
        ]
    )
    assert rc == 0
    payload = json.loads(out.read_text(encoding="utf-8"))
    pb, pc = payload["phase_b"], payload["phase_c"]
    assert pb["n_points"] == 2 == pc["n_points"]
    assert payload["n_rows"] >= 1
    assert pb["verdict"] in (
        script.VERDICT_FEASIBLE,
        script.VERDICT_OBSTRUCTED,
        script.VERDICT_MIXED,
        script.VERDICT_INCONCLUSIVE,
    )
    assert pc["consistent_with_phase_b"] is True
    for pt_b, pt_c in zip(pb["points"], pc["points"]):
        assert pt_b["status"] == pt_c["status"]
        if pt_c["status"] == script.STATUS_FEASIBLE:
            assert pt_c["n_support"] == len(pt_c["support"])
    if pb["verdict"] == script.VERDICT_FEASIBLE:
        assert pc["support_stable"] is True
        assert pc["common_support_size"] is not None


def test_level_b_deepens_box(script, tmp_path):
    n_labels = {}
    for level in ("A", "B"):
        out = tmp_path / f"{level}.json"
        rc = script.main(
            [
                "--input", str(TINY_INPUT),
                "--out", str(out),
                "--phases", "a",
                "--level", level,
                "--audit-cap", "0",
                "--audit-trials", "4",
            ]
        )
        assert rc == 0
        n_labels[level] = json.loads(out.read_text(encoding="utf-8"))["n_labels"]
    assert n_labels["B"] > n_labels["A"]
