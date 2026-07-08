"""Fast unit tests for the modular row-span certificate helper (Pass D4.4).

``verify_reduction_relation_mod_p`` certifies a *claimed* relation against a generated row
system at one exact ``(sample, prime)`` point. These tests pin: a true relation is InSpan, a
wrong one is NotInSpan with an honest residual, bad specializations (row or claimed coefficient)
reject the point, coefficient types (SymPy / ParamExpr / int / Fraction) evaluate exactly, and
no ``Success`` is involved anywhere.
"""

from __future__ import annotations

from fractions import Fraction

import pytest
import sympy as sp

from parametric_ibp_lf_reducer import (
    STATUS_IN_SPAN,
    STATUS_NOT_IN_SPAN,
    ParamExpr,
    algebraic_row,
    parse_family_text,
    verify_reduction_relation_mod_p,
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


def _coeff_row(sympy_coeff: str) -> Row:
    """Row  J(0,0) - C(ep) J(1,0) = 0  with C(ep) given by a SymPy string."""
    ep = ("ep",)
    return Row("test", {}, {
        (0, 0): ParamExpr.from_int(1, ep),
        (1, 0): ParamExpr.from_sympy(sp.sympify(f"-({sympy_coeff})"), ep),
    })


def test_true_relation_is_in_span():
    fam = _fam()
    row = algebraic_row(fam, zero_label(1, 1), 0)  # J(0,0) - J(0,-1) - J(1,-1) = 0
    cert = verify_reduction_relation_mod_p(
        fam, [row], (0, 0), {(0, -1): 1, (1, -1): 1}, {"ep": Fraction(2)}, PRIME
    )
    assert cert.status == STATUS_IN_SPAN and cert.in_span
    assert cert.residual == {}
    assert cert.nrows == 1 and cert.rank == 1


def test_wrong_relation_is_not_in_span_with_honest_residual():
    fam = _fam()
    row = algebraic_row(fam, zero_label(1, 1), 0)
    cert = verify_reduction_relation_mod_p(
        fam, [row], (0, 0), {(0, -1): 1, (1, -1): 2}, {"ep": Fraction(2)}, PRIME
    )
    assert cert.status == STATUS_NOT_IN_SPAN and not cert.in_span
    assert cert.residual  # surviving columns reported, never silently zeroed


def test_parametric_coefficient_certified_across_samples():
    fam = _fam()
    f = "(2*ep - 1)/(ep + 3)"
    row = _coeff_row(f)
    for s in (Fraction(2), Fraction(5), Fraction(7, 3)):
        cert = verify_reduction_relation_mod_p(
            fam, [row], (0, 0), {(1, 0): sp.sympify(f)}, {"ep": s}, PRIME
        )
        assert cert.in_span, f"true parametric relation rejected at ep={s}"
    wrong = verify_reduction_relation_mod_p(
        fam, [row], (0, 0), {(1, 0): sp.sympify("(2*ep - 1)/(ep + 4)")}, {"ep": Fraction(2)}, PRIME
    )
    assert not wrong.in_span


def test_mixed_coefficient_types_are_exact():
    fam = _fam()
    row = algebraic_row(fam, zero_label(1, 1), 0)
    terms = {(0, -1): ParamExpr.from_int(1, ("ep",)), (1, -1): Fraction(1)}
    cert = verify_reduction_relation_mod_p(fam, [row], (0, 0), terms, {"ep": Fraction(3)}, PRIME)
    assert cert.in_span


def test_bad_row_specialization_rejects_point():
    fam = _fam()
    row = Row("test", {}, {
        (0, 0): ParamExpr.from_int(1, ("ep",)),
        (1, 0): ParamExpr.from_sympy(sp.sympify("1/ep"), ("ep",)),  # singular at ep=0
    })
    cert = verify_reduction_relation_mod_p(
        fam, [row], (0, 0), {(1, 0): 1}, {"ep": Fraction(0)}, PRIME
    )
    assert cert.status == "BadSpecialization" and not cert.in_span


def test_bad_claimed_coefficient_rejects_point():
    fam = _fam()
    row = algebraic_row(fam, zero_label(1, 1), 0)
    terms = {(0, -1): sp.sympify("1/(ep - 2)"), (1, -1): 1}  # pole at the sample
    cert = verify_reduction_relation_mod_p(fam, [row], (0, 0), terms, {"ep": Fraction(2)}, PRIME)
    assert cert.status == "BadSpecialization" and not cert.in_span


def test_empty_system_reported():
    fam = _fam()
    cert = verify_reduction_relation_mod_p(
        fam, [], (0, 0), {(1, 0): 1}, {"ep": Fraction(2)}, PRIME
    )
    assert cert.status == "EmptySystem" and not cert.in_span


def test_non_rational_coefficient_is_a_caller_error():
    fam = _fam()
    row = algebraic_row(fam, zero_label(1, 1), 0)
    with pytest.raises(ValueError):
        verify_reduction_relation_mod_p(
            fam, [row], (0, 0), {(0, -1): sp.sympify("q + 1")}, {"ep": Fraction(2)}, PRIME
        )
