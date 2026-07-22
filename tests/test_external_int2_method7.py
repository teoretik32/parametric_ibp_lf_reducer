"""Tests for ``scripts/run_external_int2_method7.py`` (Method.7: witness stability + screening).

Tiny inputs only (``examples/tiny_success_input.wl.txt`` on a small level box) except the
``RUN_EXTERNAL_INT2_M7=1``-gated integration test, so the fast suite never launches anything
resembling a heavy Level run. The script (and the Method.6 runner it reuses) define
``@dataclass`` types under ``from __future__ import annotations``; loading via ``importlib``
therefore requires registration in ``sys.modules`` before ``exec_module``.
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

SCRIPT = REPO_ROOT / "scripts" / "run_external_int2_method7.py"
TINY_INPUT = REPO_ROOT / "examples" / "tiny_success_input.wl.txt"

WITNESS_ARTIFACT = "external_int2_method7_witness.json"
SCREEN_ARTIFACT = "external_int2_method7_screen.json"


def _load_script():
    name = "run_external_int2_method7"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod  # dataclass annotation resolution needs the module registered
    spec.loader.exec_module(mod)
    return mod


M7 = _load_script()
T2 = M7.T2


# --- shared tiny fixtures (module scope: the screen phase costs ~15s once) ---------------------


@pytest.fixture(scope="module")
def tiny_family():
    family = T2.parse_family_text(TINY_INPUT.read_text(encoding="utf-8"))
    _target, base_config = T2.build_reducer_config(family)
    spec = T2.LevelSpec(
        n_ranges=((-1, 1),) * family.nvars, m_range=(-1, 0), heavy=False, omit_special=False
    )
    return family, base_config, spec


@pytest.fixture(scope="module")
def tiny_witness_payload(tiny_family):
    family, base_config, spec = tiny_family
    return M7.run_witness_phase(
        family, base_config, primes=(10007, 10009), spec=spec, out_path=None
    )


@pytest.fixture(scope="module")
def tiny_level(tiny_family):
    family, base_config, spec = tiny_family
    labels, target_label, config = T2.build_level_config(family, base_config, spec)
    rows = T2.assemble_level_rows(family, labels, config)
    return labels, target_label, rows


@pytest.fixture(scope="module")
def synthetic_witness(tiny_witness_payload, tiny_level):
    """A hand-built Witness payload supported on the target + one real row label."""
    labels, target_label, rows = tiny_level
    some_label = None
    for row in rows["merged"]:
        for lab in row.terms:
            if lab != target_label:
                some_label = lab
                break
        if some_label:
            break
    assert some_label is not None
    return {
        "status": "Witness",
        "prime": 10007,
        "sample": tiny_witness_payload["points"][0]["sample"],
        "target_label": list(target_label),
        "witness": [[list(target_label), 1], [list(some_label), 5]],
        "rank": 1,
        "n_projected_cols": 2,
        "nullity": 1,
        "nrows": 1,
        "n_projected_rows": 1,
        "n_allowed": 1,
        "n_forbidden": 1,
        "check_annihilation": True,
        "check_target_unit": True,
        "detail": "",
    }


@pytest.fixture(scope="module")
def tiny_screen_payload(tiny_family, synthetic_witness):
    family, base_config, spec = tiny_family
    return M7.run_screen_phase(
        family,
        base_config,
        [synthetic_witness],
        spec=spec,
        out_path=None,
        tangent55_budget=120.0,
        max_candidates=5000,
    )


# --- hashing / stability helpers ---------------------------------------------------------------


def test_support_hash_order_invariant():
    a = M7.support_hash([(1, 2, 0), (0, 0, 1)])
    b = M7.support_hash([(0, 0, 1), (1, 2, 0)])
    assert a == b
    assert M7.support_hash([(0, 0, 1)]) != a


def test_vector_hash_sensitive_to_coeffs():
    base = [[(0, 0, 1), 3], [(1, 2, 0), 7]]
    changed = [[(0, 0, 1), 3], [(1, 2, 0), 8]]
    assert M7.vector_hash(base) == M7.vector_hash(list(reversed(base)))
    assert M7.vector_hash(base) != M7.vector_hash(changed)


def test_jaccard_basics():
    assert M7.jaccard(set(), set()) == 1.0
    assert M7.jaccard({1}, {2}) == 0.0
    assert M7.jaccard({1, 2}, {2, 3}) == pytest.approx(1 / 3)


def _witness_point(support_coeffs, *, status="Witness", rank=5, nullity=2, checks=True):
    return {
        "status": status,
        "rank": rank,
        "nullity": nullity,
        "witness": [[list(lab), c] for lab, c in support_coeffs],
        "check_annihilation": checks,
        "check_target_unit": checks,
    }


def test_stability_block_stable():
    pts = [
        _witness_point([((0, 0), 1), ((1, 2), 9)]),
        _witness_point([((1, 2), 4), ((0, 0), 1)]),
    ]
    st = M7.stability_block(pts)
    assert st["n_witness_points"] == 2
    assert st["all_points_witness"] is True
    assert st["n_distinct_support_patterns"] == 1
    assert st["support_pattern_stable"] is True
    assert st["pairwise_support_jaccard"] == [[1.0, 1.0], [1.0, 1.0]]
    assert st["common_core_size"] == 2
    assert st["rank_constant"] is True and st["nullity_constant"] is True
    assert st["all_checks_pass"] is True


def test_stability_block_divergent():
    pts = [
        _witness_point([((0, 0), 1), ((1, 2), 9)]),
        _witness_point([((0, 0), 1), ((3, 3), 2)], rank=6),
    ]
    st = M7.stability_block(pts)
    assert st["n_distinct_support_patterns"] == 2
    assert st["support_pattern_stable"] is False
    assert st["pairwise_support_jaccard"][0][1] == pytest.approx(1 / 3)
    assert st["common_core_size"] == 1
    assert st["rank_constant"] is False


def test_stability_block_no_witness_points():
    st = M7.stability_block([_witness_point([], status="Feasible")])
    assert st["n_witness_points"] == 0
    assert st["all_points_witness"] is False
    assert st["support_pattern_stable"] is False
    assert st["all_checks_pass"] is False


# --- Phase A on the tiny family ----------------------------------------------------------------


def test_witness_phase_tiny_schema(tiny_witness_payload):
    wp = tiny_witness_payload
    assert wp["stability"]["statuses"] == ["Feasible"]  # tiny family reduces: no obstruction
    assert len(wp["samples"]) >= M7.MIN_GENERIC
    assert wp["excluded_special_sample"]
    assert len(wp["points"]) == len(wp["samples"]) * len(wp["primes"])
    for p in wp["points"]:
        assert {"support_size", "support_hash", "vector_hash", "solve_seconds"} <= set(p)
    assert "Method.7" in wp["scope_note"]
    json.dumps(wp, default=str)  # artifact must be JSON-safe


# --- candidate generators ----------------------------------------------------------------------


def test_ray_generator_parity_with_core(tiny_family, tiny_level):
    """Ray rows (single-variable multipliers, deg<=2) must be a subset of the core deg<=2 rows."""
    family, _config, _spec = tiny_family
    labels, _target, _rows = tiny_level
    from parametric_ibp_lf_reducer.row_generation import Row, generate_coordinate_ibp_rows

    core = generate_coordinate_ibp_rows(family, labels, 2)
    core_keys = {r.dedup_key() for r in core.rows}
    rows, rejected, truncated = M7.generate_ray_multiplier_rows(family, labels, degrees=(1, 2))
    assert rows and not truncated
    for row in rows:
        assert isinstance(row, Row)
        assert row.dedup_key() in core_keys  # same surface filter, same primitive rows
    assert set(rejected) <= {"surface_unknown", "surface_not_free", "trivial_row"}


def test_ray_generator_cap(tiny_family, tiny_level):
    family, _config, _spec = tiny_family
    labels, _target, _rows = tiny_level
    rows, _rejected, truncated = M7.generate_ray_multiplier_rows(
        family, labels, degrees=(1, 2), max_candidates=1
    )
    assert len(rows) == 1 and truncated is True


def test_tangent55_budget_abort_is_recorded(tiny_family, tiny_level):
    family, _config, _spec = tiny_family
    labels, _target, _rows = tiny_level
    rows, meta = M7.generate_tangent55_rows(
        family, labels, budget_seconds=1e-9, max_candidates=100
    )
    assert meta["status"] == "aborted_budget"  # never silent
    assert rows == [] and meta["n_labels_processed"] == 0


def test_tangent_cap_truncates(tiny_family, tiny_level):
    family, _config, _spec = tiny_family
    labels, _target, _rows = tiny_level
    rows, meta = M7.generate_tangent55_rows(
        family, labels, budget_seconds=120.0, max_candidates=1, block=(1, 1)
    )
    assert len(rows) <= 1
    assert meta["truncated_by_cap"] is (len(rows) == 1)


# --- Phase B/C screening -----------------------------------------------------------------------


def test_screen_phase_schema_and_decision(tiny_screen_payload, synthetic_witness):
    sc = tiny_screen_payload
    fams = {f["family"]: f for f in sc["families"]}
    assert {"existing_algebraic", "candidate_ibp_deg3", "candidate_ray_multipliers_deg4-6",
            "candidate_tangent_(5,5)"} <= set(fams)
    for f in sc["families"]:
        assert f["category"] in ("existing", "candidate")
        for pw in f["per_witness"]:
            assert pw["n_breaks"] + pw["n_annihilate"] == f["n_candidates"]
    dec = sc["decision"]
    n_candidate = sum(f["n_breaks_total"] for f in sc["families"] if f["category"] == "candidate")
    n_existing = sum(f["n_breaks_total"] for f in sc["families"] if f["category"] == "existing")
    assert dec["n_candidate_breaks"] == n_candidate
    assert dec["n_existing_breaks"] == n_existing
    assert dec["rerun_justified"] is (n_candidate > 0)
    assert dec["rerun_caveat"] == M7.RERUN_CAVEAT
    # the synthetic witness overlaps a real baseline row -> pairing machinery must detect it
    assert n_existing >= 1 and dec["existing_rows_all_annihilate"] is False
    wp = sc["witness_points"][0]
    assert wp["support_hash"] == M7.support_hash(
        [lab for lab, _c in synthetic_witness["witness"]]
    )
    json.dumps(sc, default=str)


def test_screen_breaking_examples_are_json_safe_provenance(tiny_screen_payload):
    for f in tiny_screen_payload["families"]:
        for pw in f["per_witness"]:
            for ex in pw["breaking_examples"]:
                assert set(ex) == {"kind", "provenance"}
                assert all(isinstance(v, str) for v in ex["provenance"].values())


def test_screen_requires_witness_status(tiny_family):
    family, base_config, spec = tiny_family
    feasible_only = [{"status": "Feasible"}]
    with pytest.raises(ValueError, match="no Witness-status"):
        M7.run_screen_phase(family, base_config, feasible_only, spec=spec)


def test_screen_rejects_wrong_target(tiny_family, synthetic_witness):
    family, base_config, spec = tiny_family
    bad = dict(synthetic_witness, target_label=[9] * len(synthetic_witness["target_label"]))
    with pytest.raises(ValueError, match="witness target"):
        M7.run_screen_phase(family, base_config, [bad], spec=spec)


def test_rejected_rows_cannot_be_paired():
    """RejectedRow carries no terms: families are built from ``.rows`` only, by construction."""
    from parametric_ibp_lf_reducer.row_generation import RejectedRow

    rej = RejectedRow("coordinate_ibp", {"seed": (0, 0)}, "surface_not_free")
    assert not hasattr(rej, "terms")
    assert not hasattr(rej, "dedup_key")


# --- CLI ---------------------------------------------------------------------------------------


def test_cli_nothing_runs_by_default(tmp_path, capsys):
    assert M7.main(["--out-dir", str(tmp_path)]) == 2
    assert "--phase is required" in capsys.readouterr().err


def test_cli_witness_needs_allow_heavy(tmp_path, capsys):
    assert M7.main(["--out-dir", str(tmp_path), "--phase", "witness"]) == 2
    assert "--allow-heavy" in capsys.readouterr().err


def test_cli_rejects_bad_primes(tmp_path, capsys):
    assert M7.main(["--out-dir", str(tmp_path), "--phase", "screen", "--witness-primes", "7"]) == 2
    assert "distinct primes" in capsys.readouterr().err
    assert (
        M7.main(["--out-dir", str(tmp_path), "--phase", "screen", "--witness-primes", "7,x"]) == 2
    )


def test_cli_rejects_bad_budget_and_cap(tmp_path):
    base = ["--out-dir", str(tmp_path), "--phase", "screen"]
    assert M7.main(base + ["--tangent55-budget", "0"]) == 2
    assert M7.main(base + ["--max-candidates", "0"]) == 2


def test_cli_screen_requires_witness_artifact(tmp_path, capsys):
    assert M7.main(["--out-dir", str(tmp_path), "--phase", "screen"]) == 2
    assert "witness artifact" in capsys.readouterr().err


def test_cli_describe(tmp_path, capsys):
    assert M7.main(["--describe", "--out-dir", str(tmp_path)]) == 0
    out = capsys.readouterr().out
    assert "Method.7" in out and "absent" in out


# --- recorded artifacts (present only after the real runs) -------------------------------------


@pytest.mark.skipif(
    not (REPO_ROOT / "validation" / WITNESS_ARTIFACT).exists(),
    reason="Method.7 witness artifact not recorded yet",
)
def test_recorded_witness_artifact_schema():
    data = load_validation(WITNESS_ARTIFACT)
    assert data["level"] == 0 and "Method.7" in data["scope_note"]
    st = data["stability"]
    assert st["n_points"] == len(data["points"]) >= 6
    assert len(data["samples"]) >= 3 and len(data["primes"]) >= 2
    for p in data["points"]:
        assert p["support_hash"] and p["vector_hash"]
        assert p["check_target_unit"] is True and p["check_annihilation"] is True


@pytest.mark.skipif(
    not (REPO_ROOT / "validation" / SCREEN_ARTIFACT).exists(),
    reason="Method.7 screen artifact not recorded yet",
)
def test_recorded_screen_artifact_schema():
    data = load_validation(SCREEN_ARTIFACT)
    dec = data["decision"]
    assert dec["rerun_caveat"] == M7.RERUN_CAVEAT
    assert isinstance(dec["rerun_justified"], bool)
    assert dec["existing_rows_all_annihilate"] is True  # witness exactness on known rows
    cats = {f["category"] for f in data["families"]}
    assert cats == {"existing", "candidate"}
    for f in data["families"]:
        for pw in f["per_witness"]:
            assert pw["n_breaks"] + pw["n_annihilate"] == f["n_candidates"]


# --- gated integration -------------------------------------------------------------------------


@pytest.mark.skipif(
    os.environ.get("RUN_EXTERNAL_INT2_M7") != "1",
    reason="heavy integration run (hours); set RUN_EXTERNAL_INT2_M7=1",
)
def test_integration_full_method7(tmp_path):
    rc = M7.main(["--phase", "all", "--allow-heavy", "--out-dir", str(tmp_path)])
    assert rc == 0
    assert (tmp_path / WITNESS_ARTIFACT).exists() and (tmp_path / SCREEN_ARTIFACT).exists()
