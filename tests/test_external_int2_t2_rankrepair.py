"""Tests for ``scripts/run_external_int2_t2_rankrepair.py`` (Method.6 reproduction + witness).

Tiny inputs only (``examples/tiny_success_input.wl.txt`` with tight caps) except the
``RUN_EXTERNAL_INT2=1``-gated integration test, so the fast suite never launches anything
resembling a heavy Level run. The script defines ``@dataclass`` types under
``from __future__ import annotations``; loading it via ``importlib`` therefore requires the module
to be registered in ``sys.modules`` before ``exec_module`` (dataclass annotation resolution).
"""

from __future__ import annotations

import importlib.util
import json
import os
import sys
from pathlib import Path

import pytest

from conftest import load_validation

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

SCRIPT = REPO_ROOT / "scripts" / "run_external_int2_t2_rankrepair.py"
TINY_INPUT = REPO_ROOT / "examples" / "tiny_success_input.wl.txt"

RECORDED = {
    0: "external_int2_t2_rankrepair_level0.json",
    1: "external_int2_t2_rankrepair_level1.json",
    2: "external_int2_t2_rankrepair_level2.json",
}


def _load_script():
    spec = importlib.util.spec_from_file_location("run_external_int2_t2_rankrepair", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod  # dataclass annotation resolution needs the module registered
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def script():
    return _load_script()


@pytest.fixture(scope="module")
def tiny_family(script):
    from parametric_ibp_lf_reducer import build_reducer_config, parse_family_text

    family = parse_family_text(TINY_INPUT.read_text(encoding="utf-8"))
    _target, base_config = build_reducer_config(family)
    return family, base_config


def _obstructed_spec(script):
    return script.LevelSpec(n_ranges=((0, 0), (0, 0)), m_range=(-1, 0), heavy=False, omit_special=False)


def _feasible_spec(script):
    return script.LevelSpec(n_ranges=((0, 1), (0, 1)), m_range=(-1, 0), heavy=False, omit_special=False)


# --- level table / config ---------------------------------------------------------------


def test_level_table_matches_recorded_configs(script):
    assert script.LEVELS[0].n_ranges == ((-1, 1), (0, 1), (0, 1))
    assert script.LEVELS[1].n_ranges == ((-1, 2), (0, 1), (0, 1))
    assert script.LEVELS[2].n_ranges == ((-1, 2), (0, 1), (0, 1))
    assert script.LEVELS[0].m_range == (-3, 0)
    assert script.LEVELS[1].m_range == (-3, 0)
    assert script.LEVELS[2].m_range == (-4, 0)
    assert script.LEVELS[0].heavy is False
    assert script.LEVELS[1].heavy is True and script.LEVELS[2].heavy is True
    assert script.LEVELS[2].omit_special is True
    # BASELINE must equal the recorded baseline dict (structurally).
    baseline = {
        "max_ibp_degree": script.BASELINE["max_ibp_degree"],
        "tangent_degree_blocks": [list(b) for b in script.BASELINE["tangent_degree_blocks"]],
        "extra_block": list(script.BASELINE["extra_block"]),
    }
    assert baseline == {
        "max_ibp_degree": 2,
        "tangent_degree_blocks": [[1, 1], [2, 2]],
        "extra_block": [3, 3],
    }


# --- CLI usage / gating -----------------------------------------------------------------


def test_describe_mode_no_side_effects(script, tmp_path):
    rc = script.main(["--describe", "--out-dir", str(tmp_path)])
    assert rc == 0
    assert list(tmp_path.iterdir()) == []  # describe writes nothing


def test_missing_levels_usage_error(script):
    assert script.main([]) == 2


def test_heavy_levels_gated(script, tmp_path):
    assert script.main(["--levels", "1", "--out-dir", str(tmp_path)]) == 2
    assert script.main(["--levels", "2", "--out-dir", str(tmp_path)]) == 2
    assert list(tmp_path.iterdir()) == []  # nothing written when gated out


# --- _analyze_points genericity ---------------------------------------------------------


def test_generic_special_classification(script):
    P = 2147483647
    points = [
        {"prime": P, "sample": [["ep", "1"]],
         "result": {"rank": 10, "status": "Obstructed", "residual_support": [[0, 0]]}},
        {"prime": P, "sample": [["ep", "2"]],
         "result": {"rank": 10, "status": "Obstructed", "residual_support": [[0, 0]]}},
        {"prime": P, "sample": [["ep", "3"]],
         "result": {"rank": 5, "status": "Feasible", "residual_support": []}},
    ]
    analysis = script._analyze_points(points, 2)
    assert points[2]["classification"] == "rank_deficient_special"
    assert points[0]["classification"] == "generic"
    assert analysis["generic_verdict"] == script.VERDICT_OBSTRUCTED
    assert analysis["special_only_feasible_rejected"] is True


# --- run_level (tiny) -------------------------------------------------------------------


def test_run_level_tiny_schema(script, tiny_family, tmp_path):
    family, base_config = tiny_family
    spec = _feasible_spec(script)
    samples = script._level_samples(family, spec, 3)
    out = tmp_path / "level_repro.json"
    payload = script.run_level(
        family, base_config, spec, samples, [30011], min_generic=1, out_path=out
    )
    for key in ("label_box", "baseline", "target", "n_labels", "n_lf_true", "n_rows_total",
                "points", "analysis", "verdict", "purpose"):
        assert key in payload
    assert json.loads(json.dumps(payload, sort_keys=True, default=str))  # JSON-safe roundtrip
    assert "residual_support=[target] does not imply quotient dimension one" in payload["purpose"]
    assert "codimension-1" not in payload["purpose"]
    assert out.exists()


def test_witness_mode_tiny(script, tiny_family):
    family, base_config = tiny_family
    spec = _obstructed_spec(script)
    samples = script._level_samples(family, spec, 3)
    witness_cfg = script.WitnessConfig(samples=2, primes=(30011, 30013))
    payload = script.run_level(
        family, base_config, spec, samples, [30011], min_generic=1, witness_cfg=witness_cfg
    )
    assert payload["verdict"] == script.VERDICT_OBSTRUCTED
    assert "witness" in payload
    pts = payload["witness"]["points"]
    assert pts and all(e["status"] == script.STATUS_WITNESS for e in pts)
    assert all(e["checks_pass"] for e in pts)
    assert payload["witness"]["witness_obstruction_consistent"] is True
    # Witness <=> Obstructed: obstructed verdict here, every witness point is a Witness.
    assert payload["witness"]["all_checks_pass"] is True


def test_probe_rows_tiny(script, tiny_family):
    from parametric_ibp_lf_reducer import ParamExpr, Row
    import parametric_ibp_lf_reducer.lf_obstruction_witness as low

    family, base_config = tiny_family

    def _row(terms):
        return Row("test", {}, {lab: ParamExpr.from_int(c, ("ep",)) for lab, c in terms.items()})

    labels = [(0, 0), (1, 0), (2, 0)]
    wit = low.lf_obstruction_witness_mod_p(
        [_row({(1, 0): 1, (2, 0): 1})], labels, (0, 0), {}, {"ep": 5}, 30011
    )
    payload_wit = low.witness_to_payload(wit)
    candidate_families = {
        "probe": [_row({(0, 0): 5}), _row({(1, 0): 1, (2, 0): -1})]  # breaking, annihilating
    }
    spec = _obstructed_spec(script)
    result = script.probe_rows(
        family, base_config, spec, payload_wit, [(3, 3)], candidate_families=candidate_families
    )
    (fam_entry,) = result["families"]
    assert fam_entry["n_candidates"] == 2
    assert fam_entry["n_breaks"] == 1
    assert fam_entry["n_annihilate"] == 1
    assert result["rerun_justified"] is True


# --- recorded artifact consistency (skips if untracked files absent) --------------------


def test_recorded_artifacts_consistency():
    try:
        level0 = load_validation(RECORDED[0])
        level1 = load_validation(RECORDED[1])
        level2 = load_validation(RECORDED[2])
    except FileNotFoundError:
        pytest.skip("recorded T2 rank-repair artifacts absent (untracked validation/)")

    assert level0["verdict"] == "Mixed"
    assert level1["verdict"] == "Mixed"
    assert level2["verdict"] == "Obstructed"
    assert (level0["n_labels"], level1["n_labels"], level2["n_labels"]) == (3072, 4096, 10000)
    assert (level0["n_rows_total"], level1["n_rows_total"], level2["n_rows_total"]) == (
        46737, 59605, 155298,
    )
    max_rank = {
        0: max(p["rank"] for p in level0["points"]),
        1: max(p["rank"] for p in level1["points"]),
        2: max(p["rank"] for p in level2["points"]),
    }
    assert max_rank == {0: 24617, 1: 30807, 2: 70827}
    assert level2["note"]  # ep=3 omission note present
    # generic (max-rank) points carry residual_support == [target] = [[0,0,0,0,0,0,0]]
    for data in (level0, level1, level2):
        mr = max(p["rank"] for p in data["points"])
        for p in data["points"]:
            if p["rank"] == mr and p["status"] == "Obstructed":
                assert p["residual_support"] == [[0, 0, 0, 0, 0, 0, 0]]


# --- real External Int2 Level 0 (heavy, gated) ------------------------------------------


@pytest.mark.integration
@pytest.mark.skipif(
    os.environ.get("RUN_EXTERNAL_INT2") != "1",
    reason="heavy External Int2 T2 Level 0 rerun; set RUN_EXTERNAL_INT2=1",
)
def test_heavy_rerun_level0(script):
    from parametric_ibp_lf_reducer import build_reducer_config, parse_family_text

    recorded = load_validation(RECORDED[0])
    family = parse_family_text(script.DEFAULT_INPUT.read_text(encoding="utf-8"))
    _target, base_config = build_reducer_config(family)
    spec = script.LEVELS[0]
    samples = script._level_samples(family, spec, 4)
    payload = script.run_level(family, base_config, spec, samples, [script.DEFAULT_PRIME])
    assert payload["n_labels"] == 3072 == recorded["n_labels"]
    assert payload["n_rows_total"] == 46737 == recorded["n_rows_total"]
    assert max(p["result"]["rank"] for p in payload["points"]) == 24617
    assert payload["verdict"] == "Mixed" == recorded["verdict"]
