"""Corrected Example 4* (linearity over one shared row system) — fast unit tests + opt-in run.

Fast tests pin: the corrected input document parses and its ``LHSTerms`` claim matches the
``15*ep + 24*ep*x7`` multiplier; symbolic combination is exact (zeros drop); and the generic
``lhs_terms`` extension of the row-span certificate certifies a combined LHS honestly (InSpan
for the true relation, NotInSpan for a wrong one, classic single-target default unchanged).

The end-to-end corrected reduction is heavy (two reductions over a 972-label box; the
single-target 648-label run took ~45-50 min): opt in with ``RUN_EXAMPLE4_STAR=1``.
"""

from __future__ import annotations

import importlib.util
import json
import os
from fractions import Fraction
from pathlib import Path

import pytest
import sympy as sp

from parametric_ibp_lf_reducer import (
    STATUS_IN_SPAN,
    STATUS_NOT_IN_SPAN,
    algebraic_row,
    parse_family_text,
    verify_reduction_relation_mod_p,
    zero_label,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "run_example4_star_corrected.py"
INPUT = REPO_ROOT / "examples" / "example4_star_corrected_input.wl.txt"

PRIME = 2_147_483_647

ONE_VAR = """
IBPInput = <|
  "Variables" -> {x}, "Parameters" -> {ep}, "Regulators" -> {ep},
  "Polynomials" -> <| "G0" -> 1 + x |>,
  "MonomialExponents" -> <| x -> 0 |>,
  "PolynomialExponents" -> <| "G0" -> -2 |>
|>
"""


def _load_script():
    spec = importlib.util.spec_from_file_location("run_example4_star_corrected", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def script():
    return _load_script()


@pytest.fixture(scope="module")
def corrected_family():
    return parse_family_text(INPUT.read_text(encoding="utf-8"))


# --- corrected input document -----------------------------------------------------------------
def test_corrected_document_parses_and_matches_family(corrected_family):
    fam = corrected_family
    assert fam.nvars == 3 and fam.npolys == 4
    assert [str(p) for p in fam.parameters] == ["ep"]


def test_lhs_terms_read_from_document(script, corrected_family):
    lhs = script.lhs_terms_from_document(corrected_family)
    ep = sp.Symbol("ep")
    assert set(lhs) == {(0, 0, 0, 0, 0, 0, 0), (0, 1, 0, 0, 0, 0, 0)}
    assert sp.simplify(lhs[(0, 0, 0, 0, 0, 0, 0)] - 15 * ep) == 0
    assert sp.simplify(lhs[(0, 1, 0, 0, 0, 0, 0)] - 24 * ep) == 0


def test_label_box_gives_x7_headroom(corrected_family):
    # n7 target shift is +1; the box must extend to n7 = 2 so IBP shifts have room.
    from parametric_ibp_lf_reducer.api import build_reducer_config

    _, config = build_reducer_config(corrected_family)
    n_range, m_range = config.label_box
    assert tuple(n_range[1]) == (0, 2)
    assert len(m_range) == 4


# --- symbolic combination ---------------------------------------------------------------------
def test_combine_coefficients_is_linear_and_exact(script):
    ep = sp.Symbol("ep")
    parts = [
        (15 * ep, {(0, 0): sp.Integer(1) / ep, (1, 0): ep}),
        (24 * ep, {(1, 0): 1, (2, 0): sp.Integer(2)}),
    ]
    out = script.combine_coefficients(parts)
    assert sp.simplify(out[(0, 0)] - 15) == 0
    assert sp.simplify(out[(1, 0)] - (15 * ep**2 + 24 * ep)) == 0
    assert sp.simplify(out[(2, 0)] - 48 * ep) == 0


def test_combine_coefficients_drops_exact_zeros(script):
    out = script.combine_coefficients(
        [(sp.Integer(2), {(0, 0): sp.Integer(3)}), (sp.Integer(3), {(0, 0): sp.Integer(-2)})]
    )
    assert out == {}


# --- generic lhs_terms certificate ------------------------------------------------------------
def test_lhs_terms_true_combined_relation_in_span():
    fam = parse_family_text(ONE_VAR)
    row = algebraic_row(fam, zero_label(1, 1), 0)  # J(0,0) - J(0,-1) - J(1,-1) = 0
    ep = sp.Symbol("ep")
    for s in (Fraction(2), Fraction(7, 3)):
        cert = verify_reduction_relation_mod_p(
            fam,
            [row],
            (0, 0),
            {(0, -1): 15 * ep, (1, -1): 15 * ep},
            {"ep": s},
            PRIME,
            lhs_terms={(0, 0): 15 * ep},
        )
        assert cert.status == STATUS_IN_SPAN and cert.in_span


def test_lhs_terms_wrong_relation_not_in_span():
    fam = parse_family_text(ONE_VAR)
    row = algebraic_row(fam, zero_label(1, 1), 0)
    cert = verify_reduction_relation_mod_p(
        fam,
        [row],
        (0, 0),
        {(0, -1): 2, (1, -1): 3},
        {"ep": Fraction(2)},
        PRIME,
        lhs_terms={(0, 0): 2},
    )
    assert cert.status == STATUS_NOT_IN_SPAN and not cert.in_span
    assert cert.residual


def test_default_lhs_terms_unchanged_classic_behavior():
    fam = parse_family_text(ONE_VAR)
    row = algebraic_row(fam, zero_label(1, 1), 0)
    cert = verify_reduction_relation_mod_p(
        fam, [row], (0, 0), {(0, -1): 1, (1, -1): 1}, {"ep": Fraction(2)}, PRIME
    )
    assert cert.status == STATUS_IN_SPAN and cert.in_span and cert.residual == {}


# --- heavy end-to-end run (opt-in) ------------------------------------------------------------
RUN_FULL = os.environ.get("RUN_EXAMPLE4_STAR") == "1"


@pytest.mark.skipif(not RUN_FULL, reason="heavy (~1-2 h): set RUN_EXAMPLE4_STAR=1 to run")
def test_corrected_example4_star_end_to_end(tmp_path, script):
    out_m = tmp_path / "result.m"
    out_json = tmp_path / "diagnostics.json"
    rc = script.main(["--out", str(out_m), "--json", str(out_json)])
    assert rc == 0
    doc = json.loads(out_json.read_text(encoding="utf-8"))
    assert doc["status"] == "Success" and doc["success"] is True
    assert doc["combined"]["all_locally_finite"] is True
    assert doc["combined"]["certificate"]["certificate_status"] == "Passed"
    assert all(sr["status"] == "Success" for sr in doc["subruns"].values())
    assert out_m.read_text(encoding="utf-8").startswith("<|")
