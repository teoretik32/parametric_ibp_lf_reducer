"""External Int1 (corrected): fast guards + one gated heavy run.

Fast tests pin the contract of the validation case: the input document is the pure family
``F1`` only (no MB/z1 interpretation, no Gamma/EulerGamma anywhere the reducer looks), the
external prefactor text matches the notebook's Gamma-ratio exactly, the wrapper artifact
keeps the prefactor strictly OUTSIDE the certified reduction, and the certificate/LF gates
are mandatory. The end-to-end adaptive reduction + numeric quadrature check is heavy:
opt in with ``RUN_EXTERNAL_INT1=1``.
"""

from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path

import pytest
import sympy as sp

from conftest import load_example
from parametric_ibp_lf_reducer import parse_family_text
from parametric_ibp_lf_reducer.api import build_reducer_config
from parametric_ibp_lf_reducer.reducer import CERTIFICATE_PASSED
from parametric_ibp_lf_reducer.result import (
    FAILURE_RESOURCE_LIMIT_REACHED,
    STATUS_SUCCESS,
    ReductionDiagnostics,
    ReductionResult,
    ReductionTerm,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "run_external_int1_corrected.py"
INPUT_NAME = "external_int1_corrected_input.wl.txt"

EP, S, T = sp.symbols("ep s t")


def _load_script():
    spec = importlib.util.spec_from_file_location("run_external_int1_corrected", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def script():
    return _load_script()


@pytest.fixture(scope="module")
def input_text() -> str:
    return load_example(INPUT_NAME)


@pytest.fixture(scope="module")
def family(input_text):
    return parse_family_text(input_text)


# --- fast guards -------------------------------------------------------------------------------


def test_input_has_no_mb_interpretation(input_text, family):
    """The document is the pure family F1; no Mellin-Barnes leftovers of any kind."""
    for token in ("z1", "Mellin", "MB"):
        assert token not in input_text, f"forbidden token {token!r} in the input document"

    assert family.variables == ("x2", "x6")
    assert family.poly_names == ("G0", "G1", "G2")
    # Polynomials exactly as specified (asserted at text level: the claim lives with the data).
    assert '"G0" -> 1 + x2' in input_text
    assert '"G1" -> 1 + x6' in input_text
    assert '"G2" -> 1 + x2 + x6' in input_text
    # Monomial exponents are all zero; polynomial exponents are (ep, ep, -1+ep).
    for exponent in family.monomial_exponents:
        assert sp.simplify(exponent.to_sympy()) == 0
    expected = (EP, EP, EP - 1)
    assert len(family.polynomial_exponents) == len(expected)
    for exponent, want in zip(family.polynomial_exponents, expected):
        assert sp.simplify(exponent.to_sympy() - want) == 0


def test_family_is_pure_no_prefactor_in_core(input_text, family):
    """TargetMultiplier == 1 and nothing Gamma-flavoured ever reaches the reducer."""
    assert sp.simplify(family.target_multiplier.to_sympy() - 1) == 0
    for token in ("Gamma", "EulerGamma", "Exp[", "Zeta"):
        assert token not in input_text, f"prefactor token {token!r} leaked into the family"

    target, config = build_reducer_config(family)
    assert tuple(target) == (0, 0, 0, 0, 0)
    assert config.require_certificate_for_success is True
    n_range, m_range = config.label_box
    assert tuple(tuple(r) for r in n_range) == ((0, 1), (0, 1))
    assert tuple(tuple(r) for r in m_range) == ((-2, 0), (-2, 0), (-2, 0))
    assert config.max_ibp_degree == 2
    assert tuple(tuple(b) for b in config.tangent_degree_blocks) == ((1, 1), (2, 2))


def test_prefactor_text_and_gamma_ratio(script):
    """The wrapper prefactor matches the notebook's Gamma-ratio (task/notebook notation)."""
    text = script.EXTERNAL_PREFACTOR_TEXT
    assert "Exp[2*ep*EulerGamma]" in text

    py = (
        text.replace("^", "**")
        .replace("Exp[", "exp(")
        .replace("Gamma[", "gamma(")
        .replace("]", ")")
    )
    expr = sp.sympify(
        py,
        locals={
            "gamma": sp.gamma,
            "exp": sp.exp,
            "EulerGamma": sp.EulerGamma,
            "ep": EP,
            "s": S,
            "t": T,
        },
    )
    expected = (
        sp.exp(2 * EP * sp.EulerGamma)
        * sp.gamma(1 - EP)
        * sp.gamma(-EP) ** 2
        * sp.gamma(EP)
        * sp.gamma(2 * EP)
        / (S * T**2 * sp.gamma(-1 - 3 * EP) * sp.gamma(1 + EP))
    )
    assert sp.simplify(expr - expected) == 0

    # The Laurent series is a reference text only (expansion around ep=0).
    assert script.REFERENCE_LAURENT_TEXT.startswith("1/(s*t^2)*")


def test_wrapper_preserves_prefactor(script):
    """build_full_formula_text keeps the prefactor OUTSIDE the pure reduction."""
    result = ReductionResult(
        status=STATUS_SUCCESS,
        target_label=(0, 0, 0, 0, 0),
        all_locally_finite=True,
        terms=(
            ReductionTerm(
                label=(0, 0, -1, 0, 0),
                coefficient_text="(3*ep + 1)/ep",
                integrand_text="1/G0",
                locally_finite=True,
            ),
            ReductionTerm(
                label=(1, 0, 0, 0, -1),
                coefficient_text="-1/2",
                integrand_text="x2/G2",
                locally_finite=True,
            ),
        ),
        formal_success=True,
    )
    text = script.build_full_formula_text(result)

    assert f"ExternalPrefactor1 = {script.EXTERNAL_PREFACTOR_TEXT};" in text
    assert "FullIntegralReduction = ExternalPrefactor1*PureReduction;" in text
    assert f"ReferenceLaurentSeries = {script.REFERENCE_LAURENT_TEXT};" in text

    pure = text.split("PureReduction = ", 1)[1].split(";", 1)[0]
    for token in ("Gamma", "EulerGamma", "s*t"):
        assert token not in pure, f"prefactor token {token!r} multiplied into PureReduction"
    assert pure.count("Int[") == 2
    assert "((3*ep + 1)/ep)*Int[" in pure
    assert "(-1/2)*Int[" in pure
    assert "{x2, 0, Infinity}, {x6, 0, Infinity}" in pure
    # label (0,0,-1,0,0) -> (1+x2)^(ep-1)*(1+x6)^ep*(1+x2+x6)^(ep-1)
    assert "(x2 + 1)^(ep - 1)" in pure
    # label (1,0,0,0,-1) -> x2*(1+x2)^ep*(1+x6)^ep*(1+x2+x6)^(ep-2)
    assert "(x2 + x6 + 1)^(ep - 2)" in pure


def test_certificate_and_lf_gates_mandatory(script, family):
    """Certificate gate is default-ON and the script's success gate never weakens it."""
    _, config = build_reducer_config(family)
    assert config.require_certificate_for_success is True

    def make(status=STATUS_SUCCESS, lf=True, cert=CERTIFICATE_PASSED):
        extra = {} if cert is None else {"certificate": {"certificate_status": cert}}
        return ReductionResult(
            status=status,
            target_label=(0, 0, 0, 0, 0),
            all_locally_finite=lf,
            terms=(),
            formal_success=status == STATUS_SUCCESS,
            diagnostics=ReductionDiagnostics(extra=extra),
        )

    assert script.certified_success(make()) is True
    assert script.certified_success(make(cert=None)) is False  # certificate never ran
    assert script.certified_success(make(cert="NotRun")) is False
    assert script.certified_success(make(cert="Failed")) is False
    assert script.certified_success(make(lf="Unknown")) is False
    assert script.certified_success(make(lf=False)) is False
    assert script.certified_success(make(status=FAILURE_RESOURCE_LIMIT_REACHED)) is False


# --- gated heavy end-to-end run ------------------------------------------------------------


@pytest.mark.integration
@pytest.mark.skipif(
    not os.environ.get("RUN_EXTERNAL_INT1"), reason="set RUN_EXTERNAL_INT1=1 to run"
)
def test_full_run_certified_and_numeric(script, tmp_path):
    out = tmp_path / "reduction.m"
    diag = tmp_path / "diagnostics.json"
    full = tmp_path / "full_formula.m"

    rc = script.main(["--out", str(out), "--json", str(diag), "--full-formula", str(full)])
    assert rc == 0

    payload = json.loads(diag.read_text(encoding="utf-8"))
    assert payload["status"] == "Success"
    assert payload["all_locally_finite"] is True
    assert payload["certificate_status"] == "Passed"
    assert payload["adaptive"], "adaptive per-level history missing from diagnostics"

    check = payload["numeric_check"]
    assert check["ran"] is True
    assert check["passed"] is True
    assert check["rel_diff"] < script.NUMERIC_REL_TOL

    reduction = out.read_text(encoding="utf-8")
    assert '"Status" -> "Success"' in reduction
    assert '"AllLocallyFinite" -> True' in reduction

    formula = full.read_text(encoding="utf-8")
    assert "FullIntegralReduction = ExternalPrefactor1*PureReduction;" in formula
    assert "ReferenceLaurentSeries" in formula
