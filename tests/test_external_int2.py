"""External Int2 (dimensionless): fast guards + one gated heavy run.

Fast tests pin the contract of the validation case: the input document is the pure family
``F2`` only (dimensionless rewrite ``r = s/t``; no Gamma/``t``-scaling anywhere the reducer
looks), the external prefactor text matches the task's Gamma-ratio ``P2`` exactly, the
wrapper artifact keeps the prefactor strictly OUTSIDE the certified reduction, the 2F1
kernel of the numeric check matches direct quadrature, no reference value is invented,
and the certificate/LF gates are mandatory. The end-to-end adaptive reduction + numeric
quadrature check is heavy: opt in with ``RUN_EXTERNAL_INT2=1``.
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
SCRIPT = REPO_ROOT / "scripts" / "run_external_int2.py"
INPUT_NAME = "external_int2_dimensionless_input.wl.txt"

EP, R, T = sp.symbols("ep r t")


def _load_script():
    spec = importlib.util.spec_from_file_location("run_external_int2", SCRIPT)
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


def test_input_is_dimensionless_pure_family(input_text, family):
    """The document is the dimensionless pure family F2 with parameters (ep, r) only."""
    for token in ("z1", "Mellin", "MB", "s/t"):
        assert token not in input_text.replace("r = s/t", ""), f"forbidden token {token!r}"

    assert family.variables == ("x2", "x5", "x7")
    assert family.parameters == ("ep", "r")
    assert family.regulators == ("ep",)
    assert family.poly_names == ("G0", "G1", "G2", "G3")
    assert any("r" in a and ">" in a and "0" in a for a in family.assumptions), (
        "assumption r > 0 missing"
    )
    # Polynomials exactly as specified (asserted at text level: the claim lives with the data).
    assert '"G0" -> 1 + x2' in input_text
    assert '"G1" -> 1 + x5' in input_text
    assert '"G2" -> 1 + x7' in input_text
    assert '"G3" -> 1 + x7 + x2*x7 + r*x2*x5' in input_text
    # Monomial exponents (1+ep, 0, 0); polynomial exponents (ep, ep, -1-ep, -1+ep).
    expected_mono = (1 + EP, sp.Integer(0), sp.Integer(0))
    assert len(family.monomial_exponents) == len(expected_mono)
    for exponent, want in zip(family.monomial_exponents, expected_mono):
        assert sp.simplify(exponent.to_sympy() - want) == 0
    expected_poly = (EP, EP, -1 - EP, -1 + EP)
    assert len(family.polynomial_exponents) == len(expected_poly)
    for exponent, want in zip(family.polynomial_exponents, expected_poly):
        assert sp.simplify(exponent.to_sympy() - want) == 0


def test_family_is_pure_no_prefactor_in_core(input_text, family):
    """TargetMultiplier == 1 and nothing Gamma/t-flavoured ever reaches the reducer."""
    assert sp.simplify(family.target_multiplier.to_sympy() - 1) == 0
    for token in ("Gamma", "EulerGamma", "Exp[", "Zeta", "t^"):
        assert token not in input_text, f"prefactor token {token!r} leaked into the family"

    target, config = build_reducer_config(family)
    assert tuple(target) == (0, 0, 0, 0, 0, 0, 0)
    assert config.require_certificate_for_success is True
    assert config.rref_backend == "auto"
    n_range, m_range = config.label_box
    assert tuple(tuple(r_) for r_ in n_range) == ((0, 1), (0, 1), (0, 1))
    assert tuple(tuple(r_) for r_ in m_range) == ((-2, 0), (-2, 0), (-2, 0), (-2, 0))
    assert config.max_ibp_degree == 2
    assert tuple(tuple(b) for b in config.tangent_degree_blocks) == ((1, 1), (2, 2))
    # Default samples are deterministic scattered points over BOTH parameters (ep, r).
    assert all(set(pt) == {"ep", "r"} for pt in config.samples)
    assert len({pt["ep"] for pt in config.samples}) == len(config.samples)


def test_prefactor_text_matches_p2(script):
    """The wrapper prefactor matches P2 = t^(-3-ep)*Gamma-ratio from the task metadata."""
    text = script.EXTERNAL_PREFACTOR_TEXT
    py = text.replace("^", "**").replace("Gamma[", "gamma(").replace("]", ")")
    expr = sp.sympify(py, locals={"gamma": sp.gamma, "ep": EP, "t": T})
    expected = (
        T ** (-3 - EP)
        * sp.gamma(1 - EP)
        * sp.gamma(-EP) ** 3
        * sp.gamma(EP)
        / (sp.gamma(-1 - 3 * EP) * sp.gamma(-2 * EP))
    )
    assert sp.simplify(expr - expected) == 0


def test_reference_value_not_invented(script):
    """AnsvInt2 stays absent metadata; no Laurent text, no numeric stand-in anywhere."""
    meta = script.REFERENCE_METADATA
    assert meta["name"] == "AnsvInt2"
    assert meta["available"] is False
    assert meta["compared_numerically"] is False
    assert "GPL" in meta["note"] and "not invented" in meta["note"]
    assert not hasattr(script, "REFERENCE_LAURENT_TEXT")


def test_i7_kernel_matches_direct_quadrature(script):
    """The exact 2F1 x7-kernel agrees with direct quadrature for k = 0 and k = 1."""
    import mpmath as mp

    old = mp.mp.dps
    mp.mp.dps = 25
    try:
        ep = mp.mpf(-3) / 5
        for a_exp, b_exp, k in (
            (-1 - ep, -1 + ep, 0),  # target-like label (m2 = m3 = 0, n7 = 0)
            (-3 - ep, -2 + ep, 0),  # deepened m-shifts
            (-2 - ep, -1 + ep, 1),  # n7 = 1 with m2 = -1
        ):
            c, b = mp.mpf(2.31), mp.mpf(1.47)
            direct = mp.quad(
                lambda x: x**k * (1 + x) ** a_exp * (c + b * x) ** b_exp,
                [0, 1, 10, 100, mp.inf],
            )
            val = script.i7_kernel(mp, c, b, a_exp, b_exp, k)
            assert abs(val - direct) / abs(direct) < mp.mpf(10) ** -18
        with pytest.raises(ValueError):  # divergent case is refused, never silently wrong
            script.i7_kernel(mp, mp.mpf(2), mp.mpf(1), -1 - ep, -1 + ep, 1)
    finally:
        mp.mp.dps = old


def test_wrapper_preserves_prefactor(script):
    """build_full_formula_text keeps the prefactor OUTSIDE the pure reduction."""
    result = ReductionResult(
        status=STATUS_SUCCESS,
        target_label=(0, 0, 0, 0, 0, 0, 0),
        all_locally_finite=True,
        terms=(
            ReductionTerm(
                label=(0, 0, 0, -1, 0, 0, 0),
                coefficient_text="(3*ep*r + 1)/ep",
                integrand_text="x2^(1+ep)/G0",
                locally_finite=True,
            ),
            ReductionTerm(
                label=(1, 0, 1, 0, 0, -1, -1),
                coefficient_text="-r/2",
                integrand_text="x2^(2+ep)*x7/(G2*G3)",
                locally_finite=True,
            ),
        ),
        formal_success=True,
    )
    text = script.build_full_formula_text(result)

    assert f"ExternalPrefactor2 = {script.EXTERNAL_PREFACTOR_TEXT};" in text
    assert "FullIntegralReduction = ExternalPrefactor2*PureReduction;" in text
    assert "AnsvInt2" in text and "NOT" in text  # absence of the reference is explicit

    pure = text.split("PureReduction = ", 1)[1].split(";", 1)[0]
    for token in ("Gamma", "EulerGamma", "t^"):
        assert token not in pure, f"prefactor token {token!r} multiplied into PureReduction"
    assert pure.count("Int[") == 2
    assert "((3*ep*r + 1)/ep)*Int[" in pure
    assert "(-r/2)*Int[" in pure
    assert "{x2, 0, Infinity}, {x5, 0, Infinity}, {x7, 0, Infinity}" in pure
    # label (0,0,0,-1,0,0,0): x2^(1+ep)*(1+x2)^(ep-1)*(1+x5)^ep*(1+x7)^(-1-ep)*G3^(-1+ep)
    assert "(x2 + 1)^(ep - 1)" in pure
    # label (1,0,1,0,0,-1,-1): x2^(2+ep)*x7*...*(1+x7)^(-2-ep)*G3^(-2+ep)
    assert "(x7 + 1)^(-ep - 2)" in pure
    assert "(r*x2*x5 + x2*x7 + x7 + 1)^(ep - 2)" in pure


def test_certificate_and_lf_gates_mandatory(script, family):
    """Certificate gate is default-ON and the script's success gate never weakens it."""
    _, config = build_reducer_config(family)
    assert config.require_certificate_for_success is True

    def make(status=STATUS_SUCCESS, lf=True, cert=CERTIFICATE_PASSED):
        extra = {} if cert is None else {"certificate": {"certificate_status": cert}}
        return ReductionResult(
            status=status,
            target_label=(0, 0, 0, 0, 0, 0, 0),
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
    not os.environ.get("RUN_EXTERNAL_INT2"), reason="set RUN_EXTERNAL_INT2=1 to run"
)
def test_full_run_certified_and_numeric(script, tmp_path):
    out = tmp_path / "result.m"
    diag = tmp_path / "diagnostics.json"
    full = tmp_path / "full_formula.m"

    rc = script.main(["--out", str(out), "--json", str(diag), "--full-formula", str(full)])
    assert rc == 0

    payload = json.loads(diag.read_text(encoding="utf-8"))
    assert payload["status"] == "Success"
    assert payload["all_locally_finite"] is True
    assert payload["certificate_status"] == "Passed"
    assert payload["adaptive"], "adaptive per-level history missing from diagnostics"
    assert payload["reference_value"]["available"] is False

    check = payload["numeric_check"]
    assert check["ran"] is True
    assert check["passed"] is True
    assert check["rel_diff"] < script.NUMERIC_REL_TOL
    assert check["kernel_check_rel"] < 1e-18

    reduction = out.read_text(encoding="utf-8")
    assert '"Status" -> "Success"' in reduction
    assert '"AllLocallyFinite" -> True' in reduction

    formula = full.read_text(encoding="utf-8")
    assert "FullIntegralReduction = ExternalPrefactor2*PureReduction;" in formula
    assert "AnsvInt2" in formula
