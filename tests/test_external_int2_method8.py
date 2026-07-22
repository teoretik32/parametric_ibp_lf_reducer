"""Tests for ``scripts/run_external_int2_method8.py`` (Method.8: targeted Level-0 re-elim).

Tiny inputs only (``examples/tiny_success_input.wl.txt`` on a small level box) except the
``RUN_EXTERNAL_INT2_M8=1``-gated integration test, so the fast suite never launches the
HEAVY ~62k-row rerun. Module loading mirrors test_external_int2_method7.py: ``@dataclass``
under ``from __future__ import annotations`` requires ``sys.modules`` registration before
``exec_module``.
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

SCRIPT = REPO_ROOT / "scripts" / "run_external_int2_method8.py"
TINY_INPUT = REPO_ROOT / "examples" / "tiny_success_input.wl.txt"


def _load_script():
    name = "run_external_int2_method8"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod  # dataclass annotation resolution needs the module registered
    spec.loader.exec_module(mod)
    return mod


M8 = _load_script()
T2 = M8.T2

STATUS_FEASIBLE = M8.STATUS_FEASIBLE
STATUS_WITNESS = M8.STATUS_WITNESS


# --- shared tiny fixtures ----------------------------------------------------------------------


@pytest.fixture(scope="module")
def tiny_family():
    family = T2.parse_family_text(TINY_INPUT.read_text(encoding="utf-8"))
    _target, base_config = T2.build_reducer_config(family)
    spec = T2.LevelSpec(
        n_ranges=((-1, 1),) * family.nvars, m_range=(-1, 0), heavy=False, omit_special=False
    )
    return family, base_config, spec


@pytest.fixture(scope="module")
def tiny_rerun_payload(tiny_family):
    family, base_config, spec = tiny_family
    return M8.run_rerun_phase(
        family, base_config, primes=(10007,), spec=spec, out_path=None, recorded_dir=None
    )


# --- primes parsing ----------------------------------------------------------------------------


def test_parse_primes_good():
    assert M8._parse_primes("10007,10009") == (10007, 10009)
    assert M8._parse_primes(" 3 ") == (3,)


@pytest.mark.parametrize("bad", ["", ",", "0", "-7", "2,10007", "10007,1"])
def test_parse_primes_bad(bad):
    with pytest.raises(ValueError):
        M8._parse_primes(bad)


# --- recorded-baseline cross-checks ------------------------------------------------------------


def test_recorded_expectations_absent(tmp_path):
    out = M8._recorded_expectations(None)
    assert out == {"witness_artifact": None, "screen_artifact": None}
    out = M8._recorded_expectations(tmp_path)  # empty dir: files missing
    assert out == {"witness_artifact": None, "screen_artifact": None}


def test_recorded_expectations_extraction(tmp_path):
    wit = {
        "n_labels": 5,
        "n_rows_total": 42,
        "points": [
            {"sample": [["n1", -1]], "prime": 10007, "rank": 7, "nullity": 2,
             "support_size": 3, "extra": "ignored"}
        ],
    }
    scr = {
        "families": [
            {"family": "candidate_other", "n_candidates": 1, "n_breaks_total": 0},
            {"family": "candidate_ibp_deg3", "n_candidates": 6972, "n_breaks_total": 41832},
        ]
    }
    (tmp_path / M8.WITNESS_RECORDED_NAME).write_text(json.dumps(wit), encoding="utf-8")
    (tmp_path / M8.SCREEN_RECORDED_NAME).write_text(json.dumps(scr), encoding="utf-8")
    out = M8._recorded_expectations(tmp_path)
    assert out["witness_artifact"]["n_rows_total"] == 42
    pt = out["witness_artifact"]["points"][0]
    assert set(pt) == {"sample", "prime", "rank", "nullity", "support_size"}
    assert out["screen_artifact"] == {"n_candidates": 6972, "n_breaks_total": 41832}


def test_baseline_point_matching():
    expect = {
        "witness_artifact": {
            "points": [
                {"sample": [["n1", -1], ["m", 0]], "prime": 10007, "rank": 7, "nullity": 2},
                {"sample": [["n1", 1], ["m", 0]], "prime": 10007, "rank": 8, "nullity": 1},
            ]
        }
    }
    hit = M8._baseline_point(expect, [["n1", 1], ["m", 0]], 10007)
    assert hit is not None and hit["rank"] == 8
    assert M8._baseline_point(expect, [["n1", 1], ["m", 0]], 10009) is None
    assert M8._baseline_point(expect, [["n1", 2], ["m", 0]], 10007) is None
    assert M8._baseline_point({"witness_artifact": None}, [["n1", 1]], 10007) is None


# --- rerun phase on the tiny box ---------------------------------------------------------------


def test_rerun_schema_and_decision(tiny_rerun_payload):
    data = tiny_rerun_payload
    for key in (
        "script", "method", "scope_note", "level", "target", "n_labels", "n_rows_baseline",
        "ibp_deg3", "n_rows_enlarged", "samples", "excluded_special_sample", "primes",
        "recorded_baseline", "points", "decision", "elapsed_seconds",
    ):
        assert key in data, key
    assert data["level"] == M8.LEVEL
    g3 = data["ibp_deg3"]
    assert g3["degree"] == M8.IBP_DEGREE
    assert 0 <= g3["n_new_vs_baseline"] <= g3["n_generated"]
    assert g3["matches_screen_artifact"] is None  # no recorded_dir passed
    assert data["n_rows_enlarged"] == data["n_rows_baseline"] + g3["n_new_vs_baseline"]

    dec = data["decision"]
    n = dec["n_points"]
    assert n == len(data["points"]) == len(data["samples"]) > 0
    assert dec["n_feasible"] + dec["n_witness"] <= n
    expected = (
        "cured_all_points" if dec["n_feasible"] == n
        else "mixed" if dec["n_feasible"] > 0
        else "still_obstructed_all_points" if dec["n_witness"] == n
        else "inconclusive"
    )
    assert dec["outcome"] == expected
    assert len(dec["cured_points"]) == dec["n_feasible"]


def test_rerun_points_fields(tiny_rerun_payload):
    for pt in tiny_rerun_payload["points"]:
        assert pt["status"] in (STATUS_FEASIBLE, STATUS_WITNESS)
        assert pt["solve_seconds"] >= 0
        assert pt["support_size"] == len(pt["witness"])
        if pt["status"] == STATUS_WITNESS:
            # every included ibp_deg3 row must be annihilated by the NEW witness
            assert pt["included_ibp3_breaks_new_witness"] == 0
            assert pt["included_rows_annihilated"] is True
        else:
            assert "included_ibp3_breaks_new_witness" not in pt


def test_rerun_excluded_special_sample(tiny_rerun_payload):
    special = tiny_rerun_payload["excluded_special_sample"]
    assert special not in tiny_rerun_payload["samples"]


def test_rerun_max_points_cap_and_artifact(tiny_family, tmp_path):
    family, base_config, spec = tiny_family
    out_path = tmp_path / M8.OUT_NAME
    data = M8.run_rerun_phase(
        family, base_config, primes=(10007,), max_points=1, spec=spec, out_path=out_path
    )
    assert data["decision"]["n_points"] == 1
    on_disk = json.loads(out_path.read_text(encoding="utf-8"))
    assert on_disk["decision"] == data["decision"]


def test_rerun_recorded_baseline_deltas(tiny_family, tiny_rerun_payload, tmp_path):
    """Fake Method.7 artifacts drive matches_screen_artifact + per-point rank deltas."""
    family, base_config, spec = tiny_family
    first = tiny_rerun_payload["points"][0]
    n_new = tiny_rerun_payload["ibp_deg3"]["n_new_vs_baseline"]
    wit = {
        "n_labels": tiny_rerun_payload["n_labels"],
        "n_rows_total": tiny_rerun_payload["n_rows_baseline"],
        "points": [
            {
                "sample": first["sample"],
                "prime": first["prime"],
                "rank": max(first["rank"] - 1, 0),
                "nullity": first["nullity"] + 1,
                "support_size": 2,
            }
        ],
    }
    scr = {"families": [
        {"family": "candidate_ibp_deg3", "n_candidates": n_new, "n_breaks_total": 0}
    ]}
    (tmp_path / M8.WITNESS_RECORDED_NAME).write_text(json.dumps(wit), encoding="utf-8")
    (tmp_path / M8.SCREEN_RECORDED_NAME).write_text(json.dumps(scr), encoding="utf-8")
    data = M8.run_rerun_phase(
        family, base_config, primes=(10007,), max_points=1, spec=spec,
        out_path=None, recorded_dir=tmp_path,
    )
    assert data["ibp_deg3"]["matches_screen_artifact"] is True
    pt = data["points"][0]
    if pt["status"] == STATUS_WITNESS:
        assert pt["delta_rank"] == pt["rank"] - wit["points"][0]["rank"]
        assert pt["delta_nullity"] == pt["nullity"] - wit["points"][0]["nullity"]
    else:  # Feasible points never carry baseline deltas
        assert "delta_rank" not in pt


# --- CLI ---------------------------------------------------------------------------------------


def test_cli_nothing_runs_by_default(tmp_path, capsys):
    assert M8.main(["--out-dir", str(tmp_path)]) == 2
    assert "--phase is required" in capsys.readouterr().err


def test_cli_rerun_needs_allow_heavy(tmp_path, capsys):
    assert M8.main(["--phase", "rerun", "--out-dir", str(tmp_path)]) == 2
    assert "--allow-heavy" in capsys.readouterr().err


def test_cli_rejects_bad_primes(tmp_path, capsys):
    rc = M8.main(
        ["--phase", "rerun", "--allow-heavy", "--primes", "2,4", "--out-dir", str(tmp_path)]
    )
    assert rc == 2
    assert "bad primes" in capsys.readouterr().err


def test_cli_rejects_bad_max_points(tmp_path, capsys):
    rc = M8.main(
        ["--phase", "rerun", "--allow-heavy", "--max-points", "0", "--out-dir", str(tmp_path)]
    )
    assert rc == 2
    assert "--max-points" in capsys.readouterr().err


def test_cli_describe(tmp_path, capsys):
    artifact = {
        "decision": {"outcome": "still_obstructed_all_points", "n_feasible": 0, "n_points": 6},
        "n_rows_enlarged": 62047,
    }
    (tmp_path / M8.OUT_NAME).write_text(json.dumps(artifact), encoding="utf-8")
    assert M8.main(["--describe", "--out-dir", str(tmp_path)]) == 0
    out = capsys.readouterr().out
    assert "Method.8" in out
    assert "still_obstructed_all_points" in out


# --- gated integration -------------------------------------------------------------------------


@pytest.mark.skipif(
    os.environ.get("RUN_EXTERNAL_INT2_M8") != "1",
    reason="heavy integration run (hours); set RUN_EXTERNAL_INT2_M8=1",
)
def test_integration_full_method8(tmp_path):
    rc = M8.main(
        ["--phase", "rerun", "--allow-heavy", "--out-dir", str(tmp_path), "--max-points", "1"]
    )
    assert rc == 0
    data = json.loads((tmp_path / M8.OUT_NAME).read_text(encoding="utf-8"))
    assert data["decision"]["n_points"] == 1
