"""Tests for matrix assembly and single-sample modular normal-form extraction (Pass 2F)."""

from __future__ import annotations

import sympy as sp

import pytest

from parametric_ibp_lf_reducer import (
    BadSpecialization,
    ParamExpr,
    algebraic_row,
    assemble_rows_mod_p,
    is_locally_finite,
    modular_normal_form,
    parse_family_text,
    zero_label,
)
from parametric_ibp_lf_reducer.row_generation import Row

PRIME = 2_147_483_647

ONE_VAR = """
IBPInput = <|
  "Variables" -> {x}, "Parameters" -> {ep}, "Regulators" -> {ep},
  "Polynomials" -> <| "G0" -> 1 + x |>,
  "MonomialExponents" -> <| x -> 0 |>,
  "PolynomialExponents" -> <| "G0" -> -2 |>
|>
"""


def _fam():
    return parse_family_text(ONE_VAR)


def _row(terms_int, params=("ep",)):
    return Row("test", {}, {lab: ParamExpr.from_int(v, params) for lab, v in terms_int.items()})


def test_assemble_specializes_to_integers():
    fam = _fam()
    row = algebraic_row(fam, zero_label(1, 1), 0)  # J(0,0) - J(0,-1) - J(1,-1) = 0
    matrix = assemble_rows_mod_p(fam, [row], {"ep": 5}, PRIME)
    assert matrix == [{(0, 0): 1, (0, -1): PRIME - 1, (1, -1): PRIME - 1}]


def test_target_reduces_to_master_labels():
    fam = _fam()
    row = algebraic_row(fam, zero_label(1, 1), 0)
    res = modular_normal_form(fam, [row], (0, 0), {"ep": 5}, PRIME)
    assert res.status == "Reduced"
    assert res.formal_success is True
    assert res.pivot_label == (0, 0)
    # J(0,0) = J(0,-1) + J(1,-1)
    assert res.terms == {(0, -1): 1, (1, -1): 1}
    assert res.all_terms_lf is True and res.non_lf_terms == []


def test_bad_specialization_is_rejected_not_patched():
    fam = _fam()
    # coefficient 1/ep is singular at ep = 0.
    bad = _row({(0, 0): 1})
    bad.terms[(1, 0)] = ParamExpr.from_sympy(sp.sympify("1/ep"), ("ep",))
    res = modular_normal_form(fam, [bad], (0, 0), {"ep": 0}, PRIME)
    assert res.status == "BadSpecialization"
    assert res.formal_success is False
    assert res.terms == {}
    # the low-level assembler raises rather than skipping the row
    with pytest.raises(BadSpecialization):
        assemble_rows_mod_p(fam, [bad], {"ep": 0}, PRIME)


def test_empty_system_and_absent_target():
    fam = _fam()
    assert modular_normal_form(fam, [], (0, 0), {"ep": 3}, PRIME).status == "EmptySystem"
    row = algebraic_row(fam, zero_label(1, 1), 0)  # labels: (0,0),(0,-1),(1,-1)
    res = modular_normal_form(fam, [row], (9, 9), {"ep": 3}, PRIME)  # target absent
    assert res.status == "TargetNotReducible"
    assert res.formal_success is False


def test_lf_diagnostic_flags_non_lf_master_without_success():
    fam = _fam()
    # (0,1) is a non-LF label; force a reduction target=(0,0) -> (0,1).
    assert is_locally_finite(fam, (0, 1)) is False
    row = _row({(0, 0): 1, (0, 1): -1})  # J(0,0) - J(0,1) = 0
    res = modular_normal_form(fam, [row], (0, 0), {"ep": 5}, PRIME)
    assert res.status == "Reduced" and res.formal_success is True
    assert res.terms == {(0, 1): 1}
    # formal reduction succeeded, but a term is not locally finite -> NOT physical success.
    assert res.all_terms_lf is False
    assert res.non_lf_terms == [(0, 1)]


def test_output_is_deterministic():
    fam = _fam()
    row = algebraic_row(fam, zero_label(1, 1), 0)
    r1 = modular_normal_form(fam, [row], (0, 0), {"ep": 7}, PRIME)
    r2 = modular_normal_form(fam, [row], (0, 0), {"ep": 7}, PRIME)
    assert r1.terms == r2.terms and r1.status == r2.status
