"""Method.2 leading-pole audit of External Int2: fast guards + gated audit rerun.

Fast tests pin the audit contract: the audit script never touches the reducer core
(no ``parametric_ibp_lf_reducer`` import anywhere in its source; the ``reducer_core``
flags in the JSON artifact are all ``False``), the committed audit JSON is a
full-precision ``Success`` with every check passed and none skipped, the corrected
wrapper prefactor (``Exp[2*ep*EulerGamma]`` times the Gamma ratio) is byte-identical
between the wrapper script and the audit artifact, the prefactor Laurent data is
EulerGamma-free with leading ``6/ep^2`` (independent sympy re-derivation), and the
pole algebra ``P2_leading * J2_leading`` under ``r = s/t`` reproduces the supplied
AnsvInt2 leading pole ``-4/(s*t^2*ep^4)`` exactly.  The audit rerun is opt-in:
``RUN_INT2_POLE_AUDIT=1`` (fast mode, bounded).
"""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
import sympy as sp

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "audit_external_int2_leading_pole.py"
WRAPPER = REPO_ROOT / "scripts" / "run_external_int2.py"
AUDIT_JSON = REPO_ROOT / "validation" / "external_int2_leading_pole_audit.json"

EP, R, S, T = sp.symbols("ep r s t")

CHECK_NAMES = {
    "x7_identity_symbolic",
    "x7_identity_numeric",
    "hyp2f1_connection_formulas",
    "decomposition_consistency",
    "prefactor_series",
    "leading_pole_exact",
    "leading_pole_numeric",
    "j2_2d_cross_check",
}


@pytest.fixture(scope="module")
def payload():
    return json.loads(AUDIT_JSON.read_text(encoding="utf-8"))


def _sympify(text, names):
    return sp.sympify(text.replace("^", "**"), locals=names)


def test_audit_script_never_imports_reducer():
    """The docstring may NAME the reducer core; no code line may IMPORT it."""
    source = SCRIPT.read_text(encoding="utf-8")
    assert "import parametric_ibp_lf_reducer" not in source
    assert "from parametric_ibp_lf_reducer" not in source
    for line in source.splitlines():
        stripped = line.strip()
        if stripped.startswith(("import ", "from ")):
            assert "parametric_ibp_lf_reducer" not in stripped


def test_audit_json_reducer_core_untouched(payload):
    core = payload["reducer_core"]
    assert core["imported"] is False
    assert core["modified"] is False
    assert core["rref_run"] is False


def test_audit_json_full_precision_success(payload):
    assert payload["fast_mode"] is False
    assert payload["status"] == "Success"
    checks = payload["checks"]
    assert {c["name"] for c in checks} == CHECK_NAMES
    assert len(checks) == len(CHECK_NAMES)
    for check in checks:
        assert check["passed"] is True, check
        assert not check.get("skipped", False), check


def test_corrected_prefactor_matches_wrapper(payload):
    spec = importlib.util.spec_from_file_location("run_external_int2", WRAPPER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    corrected = payload["corrected_prefactor_text"]
    assert corrected == mod.EXTERNAL_PREFACTOR_TEXT
    assert corrected.startswith("Exp[2*ep*EulerGamma]*")
    assert payload["old_prefactor_text"] != corrected
    # The correction is exactly the Exp factor: old text is the corrected tail.
    assert corrected == "Exp[2*ep*EulerGamma]*" + payload["old_prefactor_text"]


def test_prefactor_leading_term_and_eulergamma_free():
    """Independent sympy re-derivation: leading 6/ep^2, EulerGamma-free data."""
    expr = (
        sp.exp(2 * EP * sp.EulerGamma)
        * sp.gamma(1 - EP)
        * sp.gamma(-EP) ** 3
        * sp.gamma(EP)
        / (sp.gamma(-1 - 3 * EP) * sp.gamma(-2 * EP))
    )
    lead = sp.limit(expr * EP**2, EP, 0)
    assert lead == 6
    ser = sp.series(expr, EP, 0, 1).removeO()
    # Structural check: EulerGamma cancels only after simplification (the old
    # ``free_symbols`` check was vacuous — EulerGamma is not a Symbol).
    ser = sp.simplify(sp.expand(ser))
    assert not ser.has(sp.EulerGamma)
    assert sp.simplify(ser - (6 / EP**2 + 18 / EP - 4 * sp.pi**2)) == 0


def test_pole_algebra_reproduces_supplied_ansv(payload):
    names = {"ep": EP, "r": R, "s": S, "t": T}
    p2 = _sympify(payload["p2_leading"], names)
    j2 = _sympify(payload["j2_leading"], names)
    prod = _sympify(payload["product_leading"], names)
    ansv = _sympify(payload["ansv_int2"]["supplied_leading_pole"], names)
    assert sp.simplify((p2 * j2).subs(R, S / T) - prod) == 0
    assert sp.simplify(prod - ansv) == 0
    assert payload["ansv_int2"]["leading_pole_match"] is True
    # The full AnsvInt2 value is metadata only (examples/) and is never invented.
    assert payload["ansv_int2"]["full_value_available"] is False


@pytest.mark.skipif(
    os.environ.get("RUN_INT2_POLE_AUDIT") != "1",
    reason="heavy audit rerun: set RUN_INT2_POLE_AUDIT=1",
)
def test_gated_audit_rerun_fast_mode(tmp_path):
    out = tmp_path / "audit.json"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--fast", "--json", str(out)],
        capture_output=True,
        text=True,
        timeout=3600,
        cwd=str(REPO_ROOT),
    )
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    data = json.loads(out.read_text(encoding="utf-8"))
    assert data["status"] == "Success"
    assert data["fast_mode"] is True
