"""Method.12A gates: integral-identity formula, provenance, prefactor, no-RREF audit."""
from __future__ import annotations

import ast
import json
import re
from pathlib import Path

import sympy as sp

REPO_ROOT = Path(__file__).resolve().parents[1]
FORMULA = REPO_ROOT / "validation" / "external_int2_lf_full_formula.m"
RESULT = REPO_ROOT / "validation" / "external_int2_lf_result.m"
CERT = REPO_ROOT / "validation" / "external_int2_lf_certificate.json"
LINRED = REPO_ROOT / "validation" / "external_int2_lf_linear_reducibility.json"
SCRIPTS = [
    REPO_ROOT / "scripts" / "export_external_int2_lf_masters.py",
    REPO_ROOT / "scripts" / "audit_external_int2_lf_linear_reducibility.py",
]

EP, R = sp.symbols("ep r")


def _wolfram_to_sympy(text: str) -> sp.Expr:
    return sp.sympify(text.replace("^", "**"), locals={"ep": EP, "r": R})


def test_formula_is_integral_identity_not_pointwise():
    text = FORMULA.read_text(encoding="utf-8")
    assert "Inactive[Integrate]" in text
    assert "J[TargetIntegrand] ==" in text
    for i in (1, 2, 3, 4):
        assert f"LFIntegrand{i} =" in text
        assert f"C{i}*J[LFIntegrand{i}]" in text
    assert "TargetIntegrand =" in text
    # the old pointwise forms must be gone
    assert re.search(r"\bJTarget\b", text) is None
    assert re.search(r"\bI[1-4]\s*=", text) is None
    assert "no pointwise integrand identity" in text


def test_labels_and_coefficients_unchanged():
    cert = json.loads(CERT.read_text(encoding="utf-8"))
    claimed = cert["claimed_coefficients"]
    text = FORMULA.read_text(encoding="utf-8")
    labels = sorted(claimed, key=lambda k: tuple(int(p) for p in k.strip("[] ").split(",")))
    assert len(labels) == 4
    for i, lab in enumerate(labels, start=1):
        assert f"label {lab}," in text, lab
        m = re.search(rf"^C{i} = (.+?);", text, flags=re.MULTILINE)
        assert m, f"C{i} missing"
        diff = sp.simplify(_wolfram_to_sympy(m.group(1)) - _wolfram_to_sympy(claimed[lab]))
        assert diff == 0, f"C{i} changed for label {lab}"


def test_certificate_provenance_present():
    text = RESULT.read_text(encoding="utf-8")
    for needle in (
        '"CertificateStatus" -> "Passed"',
        '"CertificateArtifact" -> "validation/external_int2_lf_certificate.json"',
        '"CertificatePoints" -> 4',
        '"CertificateChecks" -> 8',
        '"GenericRank" -> 26984',
        '"SurfacePolicy" -> "convergence_chamber"',
        '"SurfaceChamber" -> "ep=-3/5"',
    ):
        assert needle in text, needle


def test_external_prefactor_present():
    text = FORMULA.read_text(encoding="utf-8")
    assert "Exp[2*ep*EulerGamma]" in text
    assert "t^(-3-ep)" in text
    assert "Gamma[1-ep]*Gamma[-ep]^3*Gamma[ep]" in text
    assert "(Gamma[-1-3*ep]*Gamma[-2*ep])" in text
    assert "FullPhysicalReduction" in text
    assert "r = s/t;" in text


def test_no_rref_or_modular_solve_imported():
    for path in SCRIPTS:
        tree = ast.parse(path.read_text(encoding="utf-8"))
        mods: list[str] = []
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                mods += [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom):
                mod = node.module or ""
                mods.append(mod)
                mods += [f"{mod}.{a.name}" for a in node.names]
        forbidden_segments = {
            "sparse_rref",
            "sparse_rref_numba",
            "modular_normal_form",
            "records",
            "reducer",
        }
        for m in mods:
            for seg in m.lower().split("."):
                assert seg not in forbidden_segments, (path.name, m)


def test_linear_reducibility_report():
    rep = json.loads(LINRED.read_text(encoding="utf-8"))
    assert rep["all_linearly_reducible"] is True
    assert len(rep["masters"]) == 4
    for m in rep["masters"]:
        assert m["linearly_reducible"] is True
        assert m["chosen_order"] is not None
        assert m["n_valid_orders"] >= 1
        assert m["obstruction"] is None if "obstruction" in m else True
        hyperint = REPO_ROOT / m["hyperint_input"]
        assert hyperint.is_file()
        text = hyperint.read_text(encoding="utf-8")
        assert f"epsOrder := {m['required_ep_order']}" in text
        assert "# result := hyperInt(" in text  # integration NOT started
