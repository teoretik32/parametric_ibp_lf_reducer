"""Fast regression tests for the independent External Int2 analytic oracle."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "audit_external_int2_full_laurent.py"


def _load_module():
    spec = importlib.util.spec_from_file_location("external_int2_full_laurent", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_x7_preintegration_identity() -> None:
    module = _load_module()
    check = module.verify_x7_preintegration()
    assert check["passed"] is True


def test_ode_recurrence_s0_through_s3_is_exact() -> None:
    module = _load_module()
    check = module.verify_recurrence()
    assert check["passed"] is True
    assert set(check["residual_coefficients"].values()) == {"0"}


def test_full_laurent_coefficients_through_ep0_are_exact() -> None:
    module = _load_module()
    check = module.derive_full_coefficients()
    assert check["passed"] is True
    assert set(check["residuals"].values()) == {"0"}
    assert check["coefficients"]["ep^-4"] == "-4"


def test_subleading_pole_has_expected_log_structure() -> None:
    module = _load_module()
    r, t = sp.symbols("r t", positive=True)
    coefficients = module.compact_laurent_coefficients(r, t)
    assert sp.simplify(coefficients.f_m3 - (3 * sp.log(r) + 4 * sp.log(t))) == 0


def test_audit_json_is_machine_readable(tmp_path: Path) -> None:
    module = _load_module()
    output = tmp_path / "audit.json"
    result = module.run_audit(numeric=False, dps=30)
    output.write_text(json.dumps(result, ensure_ascii=False), encoding="utf-8")
    loaded = json.loads(output.read_text(encoding="utf-8"))
    assert loaded["passed"] is True
    assert loaded["schema"] == "external-int2-full-laurent-audit-v1"


def test_script_does_not_import_reducer_core() -> None:
    source = SCRIPT.read_text(encoding="utf-8")
    assert "import parametric_ibp_lf_reducer" not in source
    assert "from parametric_ibp_lf_reducer" not in source
