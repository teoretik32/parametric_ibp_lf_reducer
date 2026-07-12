"""Adaptive.2 tests: the adaptive schedule validated on a REAL explicit family.

Example 2 (``I3exampl2``, ``Examples_for_IBP_parametric.nb``) has a known 5-term reduction
(``validation/notebook_example2_n3_five_term_expected.json``). Starting from a deliberately
shallow label box, the **default MVP schedule** (no hand-crafted levels) must escalate
``base`` -> ``expand-1`` and certify exactly that basis with the notebook coefficients.

* ``test_default_schedule_certifies_real_five_term_family`` — API, normal suite, ~25 s: the
  single slowest normal test, kept because it is the only end-to-end proof that adaptive mode
  reproduces a known real-family result (everything else in Adaptive.1 runs tiny/stub math).
* ``test_cli_adaptive_certifies_real_family_medium`` — the same run through the CLI with the
  config carried in the document ``Options``; medium (~30-60 s), set ``RUN_ADAPTIVE_MEDIUM=1``
  to run (same convention as ``RUN_D4_FULL``).

No stubs, no mocks: every level runs the full pipeline including the certificate gate.
"""

from __future__ import annotations

import json
import os
from fractions import Fraction

import pytest
import sympy as sp

from conftest import load_example, load_validation
from parametric_ibp_lf_reducer import (
    FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE,
    STATUS_SUCCESS,
    AdaptiveSearchConfig,
    ReducerConfig,
    parse_family_text,
    reduce_family_adaptive,
)
from parametric_ibp_lf_reducer.api import _DEFAULT_PRIMES, default_scattered_samples
from parametric_ibp_lf_reducer.cli import EXIT_SUCCESS, main
from parametric_ibp_lf_reducer.reducer import CERTIFICATE_PASSED

EX2_FILE = "notebook_example2_n3_five_term_explicit_family.wl.txt"
EX2_FIXTURE = "notebook_example2_n3_five_term_expected.json"

TARGET = (0, 0, 0, 0, 0, 0, 0)

# Known 5-term basis, in the order of the fixture's ``relative_factors`` /
# ``coefficients`` lists: y/(G0*G1), y/(G0*G1^2), y/(G0*G1*G2), y/(G1^2*G2), y/(G1^2*G2^2).
EXPECTED_BASIS = (
    (0, 1, 0, -1, -1, 0, 0),
    (0, 1, 0, -1, -2, 0, 0),
    (0, 1, 0, -1, -1, -1, 0),
    (0, 1, 0, 0, -2, -1, 0),
    (0, 1, 0, 0, -2, -2, 0),
)

# Deliberately shallow base box: one m-depth short of a certifying configuration, so level 0
# must fail honestly and the default schedule's ``expand-1`` (m-ranges deepened by one,
# max_ibp_degree=2, tangent ((1,1),)) lands exactly on a box that certifies.
SHALLOW_BOX = (((0, 0), (0, 1), (0, 0)), ((0, 0), (-1, 0), (-1, 0), (0, 0)))

N_SAMPLES = 16  # 12 (default count) is not enough for interpolation here; 16 certifies.


def _base_config(family) -> ReducerConfig:
    return ReducerConfig(
        primes=list(_DEFAULT_PRIMES),
        samples=default_scattered_samples(family.parameters, N_SAMPLES),
        label_box=SHALLOW_BOX,
        preferred_masters=list(EXPECTED_BASIS),
        rref_backend="auto",
    )


def _assert_coefficients_match_fixture(terms) -> None:
    """Every reconstructed coefficient must equal the notebook value symbolically."""
    fixture = load_validation(EX2_FIXTURE)
    by_label = {t.label: t for t in terms}
    for label, expected_text in zip(EXPECTED_BASIS, fixture["coefficients"]):
        got = sp.sympify(by_label[label].coefficient_text.replace("^", "**"))
        want = sp.sympify(expected_text.replace("^", "**"))
        assert sp.simplify(got - want) == 0, (
            f"coefficient mismatch for {label}: got {got}, expected {want}"
        )


# --- fast real case (normal suite) -------------------------------------------------------------
def test_default_schedule_certifies_real_five_term_family():
    family = parse_family_text(load_example(EX2_FILE))
    result = reduce_family_adaptive(
        family, TARGET, _base_config(family), AdaptiveSearchConfig(max_levels=2)
    )

    assert result.status == STATUS_SUCCESS
    ad = result.diagnostics.extra["adaptive"]
    assert ad["stop_reason"] == "success"
    levels = ad["levels"]
    assert [lvl["name"] for lvl in levels] == ["base", "expand-1"]

    # Level 0: ran for real and failed honestly, with the *correct* recommendation.
    base = levels[0]
    assert base["ran"] is True
    assert base["status"] == FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE
    assert "expand the label box" in base["recommendation"]

    # Level 1: certified Success on the expanded box.
    exp1 = levels[1]
    assert exp1["ran"] is True
    assert exp1["status"] == STATUS_SUCCESS
    assert exp1["certificate_status"] == CERTIFICATE_PASSED
    assert exp1["max_ibp_degree"] == 2

    # Exactly the known basis, with the notebook coefficients.
    assert sorted(t.label for t in result.terms) == sorted(EXPECTED_BASIS)
    _assert_coefficients_match_fixture(result.terms)


# --- gated medium case: same run through the CLI ------------------------------------------------
def _wolfram_document_with_adaptive_options(family) -> str:
    """The example document with the shallow-box config carried in ``Options``."""
    samples = default_scattered_samples(family.parameters, N_SAMPLES)

    def _frac(v: Fraction) -> str:
        return str(v)

    sample_text = ",\n      ".join(
        "<| " + ", ".join(f'"{name}" -> {_frac(val)}' for name, val in pt.items()) + " |>"
        for pt in samples
    )
    masters_text = ", ".join("{" + ",".join(map(str, lab)) + "}" for lab in EXPECTED_BASIS)
    n_text = ", ".join("{%d,%d}" % r for r in SHALLOW_BOX[0])
    m_text = ", ".join("{%d,%d}" % r for r in SHALLOW_BOX[1])
    injected = (
        '"Options" -> <|\n'
        f'    "LabelBox" -> {{{{{n_text}}}, {{{m_text}}}}},\n'
        f'    "PreferredMasters" -> {{ {masters_text} }},\n'
        '    "RREFBackend" -> "auto",\n'
        f'    "Samples" -> {{\n      {sample_text}\n    }},\n'
    )
    text = load_example(EX2_FILE)
    assert '"Options" -> <|' in text
    return text.replace('"Options" -> <|', injected, 1)


@pytest.mark.skipif(
    not os.environ.get("RUN_ADAPTIVE_MEDIUM"),
    reason="medium ~30-60s real-family adaptive CLI e2e; set RUN_ADAPTIVE_MEDIUM=1 to run",
)
def test_cli_adaptive_certifies_real_family_medium(tmp_path):
    family = parse_family_text(load_example(EX2_FILE))
    doc = tmp_path / "ex2_adaptive.m"
    doc.write_text(_wolfram_document_with_adaptive_options(family), encoding="utf-8")
    out = tmp_path / "result.m"
    diag = tmp_path / "diag.json"

    rc = main(
        [
            "reduce",
            str(doc),
            "--out",
            str(out),
            "--adaptive",
            "--adaptive-max-levels",
            "2",
            "--diagnostics-json",
            str(diag),
        ]
    )
    assert rc == EXIT_SUCCESS

    text = out.read_text(encoding="utf-8")
    assert '"Status" -> "Success"' in text
    assert "**" not in text  # Wolfram-style powers only

    payload = json.loads(diag.read_text(encoding="utf-8"))
    assert payload["status"] == STATUS_SUCCESS
    assert sorted(tuple(t["label"]) for t in payload["terms"]) == sorted(EXPECTED_BASIS)

    ad = payload["adaptive"]
    assert ad["stop_reason"] == "success"
    assert len(ad["levels"]) == 2
    base, exp1 = ad["levels"]
    assert base["ran"] is True
    assert base["status"] != STATUS_SUCCESS
    assert base["recommendation"], "honest level-0 failure must carry a recommendation"
    assert exp1["status"] == STATUS_SUCCESS
    assert exp1["certificate_status"] == CERTIFICATE_PASSED
