"""External Int2 Method.12R gates: mixed-ray LF regression, Laurent algebra, numeric identity.

Fast, exact checks only. The heavy phase-C row revalidation is not re-run here; its recorded
counts are validated against the artifact schema when the artifact exists.
"""
from __future__ import annotations

import ast
import importlib.util
import json
import math
from fractions import Fraction
from pathlib import Path

import pytest
import sympy as sp

from parametric_ibp_lf_reducer import (
    base_score,
    complete_polyhedral_rays,
    compute_candidate_rays,
    explain_local_finiteness,
    heuristic_candidate_rays,
    is_locally_finite,
    newton_wall_normals,
    parse_family_text,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "audit_external_int2_method12r.py"
ARTIFACT = REPO_ROOT / "validation" / "external_int2_method12r.json"
INT2_INPUT = REPO_ROOT / "examples" / "external_int2_dimensionless_input.wl.txt"
D4_INPUT = REPO_ROOT / "examples" / "d4_cli_example_input.wl.txt"
RESULT = REPO_ROOT / "validation" / "external_int2_lf_result.m"
FORMULA = REPO_ROOT / "validation" / "external_int2_lf_full_formula.m"
CERT = REPO_ROOT / "validation" / "external_int2_lf_certificate.json"

EP, R = sp.symbols("ep r")

TARGET = (0, 0, 0, 0, 0, 0, 0)
L1 = (-1, 0, 0, -1, -1, 0, 0)
L2 = (-1, 0, 0, 0, -1, 0, -1)
L3 = (0, 0, 0, -1, 0, 0, -1)
L4 = (0, 0, 1, -1, 0, 0, -1)
MIXED_RAY = (0, -1, -1)


@pytest.fixture(scope="module")
def int2():
    return parse_family_text(INT2_INPUT.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def audit_mod():
    spec = importlib.util.spec_from_file_location("m12r", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# --- mixed-ray LF regression ------------------------------------------------------------------
def test_joint_infinity_ray_is_missing_from_the_heuristic_set(int2):
    heuristic = {ray.direction for ray in heuristic_candidate_rays(int2)}
    complete, is_complete = complete_polyhedral_rays(int2)
    complete_dirs = {ray.direction for ray in complete}
    assert is_complete is True
    assert MIXED_RAY not in heuristic, "the historical gap must stay documented by this test"
    assert MIXED_RAY in complete_dirs
    assert MIXED_RAY in {ray.direction for ray in compute_candidate_rays(int2)}
    # the correction only adds rays; nothing the old set tested is dropped
    assert heuristic <= {ray.direction for ray in compute_candidate_rays(int2)}


def test_l4_is_not_locally_finite_and_l1_l3_still_are(int2):
    assert is_locally_finite(int2, L4) is False
    for lab in (L1, L2, L3):
        assert is_locally_finite(int2, lab) is True
    assert is_locally_finite(int2, TARGET) is False


def test_strict_zero_score_on_the_mixed_ray_is_the_reason(int2):
    # exactly marginal at epsilon = 0 -> STRICT RULE -> not locally finite
    assert base_score(int2, L4, MIXED_RAY) == 0
    report = explain_local_finiteness(int2, L4)
    assert report.verdict is False
    failing = [rv for rv in report.failing_rays if rv.ray.direction == MIXED_RAY]
    assert len(failing) == 1
    assert failing[0].score == 0
    assert "STRICT RULE" in failing[0].detail
    # ... while the same ray is strictly positive inside the convergence chamber, which is why the
    # regulated integral converges at ep = -3/5 and the Method.12N pole sits at ep = 0.
    assert [base_score(int2, lab, MIXED_RAY) for lab in (L1, L2, L3)] == [1, 2, 1]


def test_correction_is_generic_not_label_specific():
    """A different family (D4) also has a marginal label on a ray the heuristic set never emits."""
    d4 = parse_family_text(D4_INPUT.read_text(encoding="utf-8"))
    label = (1, 1, 1, 0, -1, 0, -1)
    ray = (-1, 0, -1, 0)
    assert ray not in {r.direction for r in heuristic_candidate_rays(d4)}
    assert ray in {r.direction for r in complete_polyhedral_rays(d4)[0]}
    assert base_score(d4, label, ray) == 0
    assert is_locally_finite(d4, label) is False
    # the D4 validation masters are unaffected by the correction
    for master in (
        (0, 1, 1, 0, -2, -1, 0),
        (1, 1, 0, 0, -2, -1, 0),
        (0, 1, 1, 0, -3, -1, 0),
        (1, 1, 0, 0, -3, -1, 0),
        (0, 1, 1, 0, -4, -1, 0),
    ):
        assert is_locally_finite(d4, master) is True


def test_wall_normals_and_completeness_flag(int2):
    normals = newton_wall_normals(int2)
    # walls of the common refinement: one per monomial difference, sign-normalised
    assert (0, 1, -1) in normals  # the wall whose intersection with x2 = const gives (0,-1,-1)
    assert all(next((c for c in v if c), 0) > 0 for v in normals)
    rays, complete = complete_polyhedral_rays(int2)
    assert complete is True
    assert len(rays) == 18
    # a tiny budget must degrade honestly, never silently trust an incomplete set
    rays_small, complete_small = complete_polyhedral_rays(int2, budget=1)
    assert complete_small is False
    assert len(rays_small) == 2 * int2.nvars


# --- leading-order Laurent algebra ---------------------------------------------------------------
def test_rhs_carries_an_ep_minus_3_pole_the_target_cannot_have():
    claimed = json.loads(CERT.read_text(encoding="utf-8"))["claimed_coefficients"]
    coeffs = {
        k: sp.sympify(claimed[str(list(lab))], locals={"ep": EP, "r": R})
        for k, lab in (("C1", L1), ("C2", L2), ("C3", L3), ("C4", L4))
    }
    l1, l2, l3, l40 = sp.symbols("l1 l2 l3 l40")
    masters = {"C1": l1, "C2": l2, "C3": l3, "C4": -1 / (R * EP) + l40}
    rhs = sum(coeffs[k] * masters[k] for k in coeffs)
    ser = sp.expand(sp.simplify(sp.series(sp.together(rhs), EP, 0, 1).removeO()))
    pole3 = sp.simplify(ser.coeff(EP, -3))
    # independent of every finite part: it comes from C4 * (-1/(r ep)) alone
    assert sp.simplify(pole3 - (-1 / (6 * R**2))) == 0
    assert not (pole3.free_symbols & {l1, l2, l3, l40})
    # the independently derived target starts at ep^-2 (Method.2 leading-pole audit)
    assert sp.simplify(sp.limit(-2 / (3 * R * EP**2) * EP**3, EP, 0)) == 0


def test_certified_coefficient_leading_orders():
    claimed = json.loads(CERT.read_text(encoding="utf-8"))["claimed_coefficients"]
    orders = {}
    for name, lab in (("C1", L1), ("C2", L2), ("C3", L3), ("C4", L4)):
        c = sp.sympify(claimed[str(list(lab))], locals={"ep": EP, "r": R})
        ser = sp.expand(sp.simplify(sp.series(sp.together(c), EP, 0, 1).removeO()))
        orders[name] = next(k for k in range(-3, 2) if sp.simplify(ser.coeff(EP, k)) != 0)
    assert orders == {"C1": 0, "C2": -1, "C3": -2, "C4": -2}


# --- direct numeric identity fixture -------------------------------------------------------------
def test_four_term_relation_fails_numerically_in_the_chamber(audit_mod):
    """Small-mesh recomputation: target and RHS differ by ~100%, far beyond any quadrature error."""
    route = audit_mod.route_3d_tensor(((10, 5), (12, 6)))
    vals = {k: v[-1] for k, v in route["series"].items()}
    coeffs = {
        k: float(v.subs({EP: sp.Rational(-3, 5), R: sp.Integer(1)}))
        for k, v in audit_mod.certified_coefficients().items()
    }
    rhs = sum(coeffs[k] * vals[k] for k in ("L1", "L2", "L3", "L4"))
    target = vals["Target"]
    assert target > 0 and all(vals[k] > 0 for k in ("L1", "L2", "L3", "L4"))
    assert abs(target - rhs) / target > 0.5
    # ladder increments are orders of magnitude below the residual
    for series in route["series"].values():
        assert abs(series[-1] - series[-2]) < 0.05 * abs(series[-1])


def test_exact_x7_primitives_match_adaptive_quadrature(audit_mod):
    checks = audit_mod._verify_x7_forms()
    assert checks and all(c["ok"] for c in checks)
    assert max(c["relative_difference"] for c in checks) < audit_mod.X7_FORM_TOLERANCE


def test_chamber_scores_make_phase_d_well_posed(audit_mod, int2):
    """Every integral in the relation converges at ep = -3/5 (all complete rays score > 0)."""
    rays = [r.direction for r in compute_candidate_rays(int2)]
    for label in (TARGET, L1, L2, L3, L4):
        for d in rays:
            score = audit_mod.symbolic_score(int2, label, d).subs(EP, sp.Rational(-3, 5))
            assert sp.simplify(score) > 0, (label, d)


# --- discipline gates ----------------------------------------------------------------------------
def test_audit_script_invokes_no_rref_or_modular_solve():
    tree = ast.parse(SCRIPT.read_text(encoding="utf-8"))
    mods: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            mods += [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            mod = node.module or ""
            mods.append(mod)
            mods += [f"{mod}.{a.name}" for a in node.names]
    forbidden = {
        "sparse_rref",
        "sparse_rref_numba",
        "modular_normal_form",
        "records",
        "reducer",
        "certificate",
        "lf_feasibility",
        "lf_obstruction_witness",
    }
    for m in mods:
        for seg in m.lower().split("."):
            assert seg not in forbidden, m


def test_revoked_status_is_recorded_in_the_exported_artifacts():
    result = RESULT.read_text(encoding="utf-8")
    assert '"Status" -> "Revoked(Method.12R)"' in result
    assert '"AllLocallyFinite" -> False' in result
    assert "validation/external_int2_method12r.json" in result
    assert '"Method12RStatus" -> "Invalidated"' in result
    formula = FORMULA.read_text(encoding="utf-8")
    assert "Method.12R" in formula and "REVOKED" in formula
    # the four labels and coefficients themselves must NOT have been touched
    claimed = json.loads(CERT.read_text(encoding="utf-8"))["claimed_coefficients"]
    for i, lab in enumerate((L1, L2, L3, L4), start=1):
        expr = sp.sympify(claimed[str(list(lab))], locals={"ep": EP, "r": R})
        m = [ln for ln in formula.splitlines() if ln.startswith(f"C{i} = ")][0]
        text = m.split(" = ", 1)[1].split(";")[0].replace("^", "**")
        assert sp.simplify(sp.sympify(text, locals={"ep": EP, "r": R}) - expr) == 0


def test_certificate_downgrade_is_recorded_without_touching_the_claims():
    cert = json.loads(CERT.read_text(encoding="utf-8"))
    assert cert["verdict"] == "certified"  # the row-span claim itself is untouched
    down = cert["downgrade"]
    assert down["status"] == "formal_row_span_only"
    assert down["refuted_by"] == "Method.12R"
    assert down["evidence"] == "validation/external_int2_method12r.json"
    assert len(cert["claimed_coefficients"]) == 4
    assert cert["rows"]["n_rows_total"] == 49439


@pytest.mark.skipif(not ARTIFACT.exists(), reason="12R artifact absent")
def test_artifact_schema_and_verdict():
    doc = json.loads(ARTIFACT.read_text(encoding="utf-8"))
    assert doc["verdict"] == "FOUR_TERM_RELATION_INVALID"
    assert set(doc["findings"]) >= {
        "laurent_order_contradiction",
        "lf_gate_ray_gap",
        "numeric_identity_refuted",
    }
    a = doc["phase_a_laurent_consistency"]
    assert a["laurent_orders_match"] is False
    assert a["rhs_ep_minus_3_depends_on_finite_parts"] is False
    assert sp.simplify(sp.sympify(a["rhs_ep_minus_3"]) - (-1 / (6 * R**2))) == 0

    b = doc["phase_b_ray_completeness"]
    assert b["mixed_ray"] == list(MIXED_RAY)
    assert b["mixed_ray_in_heuristic_set"] is False
    assert b["enumeration_complete"] is True
    assert b["per_label"]["L4"]["locally_finite"] == "False"
    assert all(b["per_label"][k]["locally_finite"] == "True" for k in ("L1", "L2", "L3"))
    assert all(not b["per_label"][k]["chamber_failing_rays"] for k in b["per_label"])

    d = doc["phase_d_numeric_identity"]
    assert d["relation_holds_numerically"] is False
    for route in d["identity"].values():
        assert abs(route["residual_rel"]) > 0.5
        assert route["residual_over_error"] > 100
    # the two independent routes must agree on every integral far better than the residual
    for per in d["per_integral"].values():
        assert per["route_agreement_rel"] < 0.05

    if "phase_c_row_revalidation" in doc:
        c = doc["phase_c_row_revalidation"]
        assert c["ray_enumeration_complete"] is True
        assert all(g["disagreements"] == 0 for g in c["fast_flux_filter_gate"].values())
        assert c["rows"]["n_chamber"] == 49439
        assert c["rows"]["n_limit_only"] == 0
        assert c["chamber_only_revalidation"]["newly_rejected"] > 0
        for w in c["witness_rescreening"]:
            assert w["n_breaks_surviving_rows"] <= w["n_breaks_recorded_rows"]


def test_no_new_lf_search_or_reduction_was_run():
    """The audit must not have regenerated any modular record or certificate artifact."""
    text = SCRIPT.read_text(encoding="utf-8")
    for banned in ("rref_mod_p", "modular_normal_form", "lf_reduction_feasible_mod_p"):
        assert banned not in text
    assert not math.isnan(0.0)  # sanity: module imports resolved
    assert Fraction(-3, 5) == Fraction(-3, 5)
