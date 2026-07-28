"""Tests for External Int2 Method.12N (numeric epsilon expansion).

Fast checks only: master parsing, certified-coefficient wiring, oracle
isolation, freeze/guard discipline, and (if the production artifact exists)
its schema, including the honest divergent-termwise verdict for L4.
"""

import ast
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "numeric_external_int2_lf_expansion.py"
ARTIFACT = REPO_ROOT / "validation" / "external_int2_lf_numeric_expansion.json"

CERTIFIED_LABELS = [
    [-1, 0, 0, -1, -1, 0, 0],
    [-1, 0, 0, 0, -1, 0, -1],
    [0, 0, 0, -1, 0, 0, -1],
    [0, 0, 1, -1, 0, 0, -1],
]


def load_mod():
    spec = importlib.util.spec_from_file_location("m12n", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def mod():
    return load_mod()


@pytest.fixture(scope="module")
def masters(mod):
    return mod.parse_masters()


def test_masters_labels_orders_and_common_log(masters):
    assert [m["label"] for m in masters] == CERTIFIED_LABELS
    assert [m["eps_order"] for m in masters] == [2, 3, 4, 4]
    a0 = masters[0]["a"]
    assert all(m["a"] == a0 for m in masters)
    # normalized tables must not carry zero entries
    for m in masters:
        assert 0 not in m["a"].values()
        assert 0 not in m["b"].values()
    # structural ep=0 exponents of L4 (the divergent-face master)
    assert masters[3]["b"] == {"x2": 1, "x7": 1, "G0": -1, "G2": -1, "G3": -2}


def test_certified_coefficients_pole_structure(mod, masters):
    import sympy as sp

    from fractions import Fraction

    coeffs = mod.load_certified_coefficients(masters)
    assert len(coeffs) == 4
    r = Fraction(3, 7)
    lead = []
    for c in coeffs:
        ser = mod.coefficient_series(c, r, lo=-2, hi=1)
        lead.append(ser[-2])
    # C1, C2 are finite at ep=0; C3, C4 start at ep^-2
    assert lead[0] == 0 and lead[1] == 0
    assert lead[2] != 0 and lead[3] != 0
    # exact C1(ep=0) = -1/(6r)
    ser1 = mod.coefficient_series(coeffs[0], r, lo=0, hi=0)
    assert sp.nsimplify(ser1[0]) == sp.nsimplify(-1 / (6 * sp.Rational(3, 7)))


def test_no_rref_or_modular_imports():
    text = SCRIPT.read_text(encoding="utf-8")
    tree = ast.parse(text)
    banned = ("parametric_ibp_lf_reducer", "rref", "reduction", "records")
    for node in ast.walk(tree):
        names = []
        if isinstance(node, ast.Import):
            names = [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            names = [node.module or ""]
        for name in names:
            assert not any(b in name for b in banned), name


def test_oracle_read_only_in_comparison_stage():
    text = SCRIPT.read_text(encoding="utf-8")
    tree = ast.parse(text)
    spans = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef):
            spans[node.name] = (node.lineno, node.end_lineno)
    lo, hi = spans["run_comparison_stage"]
    for i, line in enumerate(text.splitlines(), start=1):
        if "full_laurent_audit" in line:
            assert lo <= i <= hi, f"oracle path outside comparison stage: line {i}"


def test_compare_guard_refuses_divergent_artifact(mod, tmp_path):
    frozen = {
        "points": [],
        "all_masters_termwise_convergent": False,
        "verdict": "FAILURE_TERMWISE_EXPANSION_DIVERGENT",
    }
    blob = json.dumps(frozen, sort_keys=True).encode()
    doc = {
        "frozen": frozen,
        "frozen_sha256": hashlib.sha256(blob).hexdigest(),
        "oracle_read_during_compute": False,
    }
    p = tmp_path / "artifact.json"
    p.write_text(json.dumps(doc), encoding="utf-8")
    with pytest.raises(SystemExit, match="refusing to compare"):
        mod.run_comparison_stage(p)


def test_compare_guard_refuses_tampered_freeze(mod, tmp_path):
    frozen = {"points": [], "all_masters_termwise_convergent": True}
    doc = {"frozen": frozen, "frozen_sha256": "0" * 64}
    p = tmp_path / "artifact.json"
    p.write_text(json.dumps(doc), encoding="utf-8")
    with pytest.raises(SystemExit, match="freeze hash mismatch"):
        mod.run_comparison_stage(p)


def test_quadrature_anchors_pi_squared(mod, masters):
    # exact anchors at r=1: j_{2,0} = pi^2/8, j_{3,0} = pi^2/6
    acc = mod.sweep(masters, 1.0, 8, 12, n_max=0)
    assert abs(acc[1, 0] - math.pi**2 / 8) < 5e-4
    assert abs(acc[2, 0] - math.pi**2 / 6) < 1e-3


@pytest.mark.skipif(not ARTIFACT.exists(), reason="production artifact absent")
def test_artifact_schema_and_honest_verdict():
    doc = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    frozen = doc["frozen"]
    blob = json.dumps(frozen, sort_keys=True).encode()
    assert hashlib.sha256(blob).hexdigest() == doc["frozen_sha256"]
    assert doc["oracle_read_during_compute"] is False
    assert frozen["verdict"] == "FAILURE_TERMWISE_EXPANSION_DIVERGENT"
    assert frozen["all_masters_termwise_convergent"] is False
    assert doc["comparison"]["status"] == "not_run"
    for point in frozen["points"]:
        pm = point["per_master"]
        assert len(pm) == 4
        div = [lab for lab, m in pm.items() if m["status"] == "divergent_termwise"]
        assert div == [str(CERTIFIED_LABELS[3])]
        for lab, m in pm.items():
            if m["status"] == "converged":
                for ent in m["coefficients"].values():
                    for key in ("value", "err_est", "stab_precision", "stab_subdivision"):
                        assert key in ent
            else:
                g = m["growth_diagnostic"]
                # measured log-slope of j_{4,0} must match the analytic 1/r
                assert abs(g["log_slope_n0"] - g["predicted_slope_1_over_r"]) < 0.02
                incs = g["0"]["increments"]
                assert abs(incs[-1] / incs[-2] - 1.0) < 0.05  # constant => log-divergent
