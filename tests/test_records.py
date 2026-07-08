"""Tests for the multi-sample modular normal-form record collector (Pass 2G.1).

The collector runs the real single-sample ``modular_normal_form`` over a grid of primes x samples
and records every point honestly. These tests also confirm the records feed reconstruction
end-to-end (still univariate; no Success, no multivariate here).
"""

from __future__ import annotations

from fractions import Fraction

import sympy as sp

from parametric_ibp_lf_reducer import (
    NormalFormRecord,
    ParamExpr,
    algebraic_row,
    collect_normal_form_records,
    parse_family_text,
    reconstruct_coefficients,
    record_from_result,
    summarize_records,
    zero_label,
)
from parametric_ibp_lf_reducer.modular_normal_form import NormalFormResult
from parametric_ibp_lf_reducer.row_generation import Row

PRIMES = [2_147_483_647, 2_147_483_629, 2_147_483_587]

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


def _coeff_row(sympy_master_coeff: str):
    """Row  J(0,0) - C(ep) J(1,0) = 0  with C(ep) given by a SymPy string."""
    ep = ("ep",)
    return Row(
        "test",
        {},
        {
            (0, 0): ParamExpr.from_int(1, ep),
            (1, 0): ParamExpr.from_sympy(sp.sympify(f"-({sympy_master_coeff})"), ep),
        },
    )


def _samples(values):
    return [{"ep": Fraction(v)} for v in values]


def test_record_from_result_maps_fields():
    res = NormalFormResult(
        status="Reduced", target_label=(0, 0), prime=101, sample={"ep": Fraction(2)},
        formal_success=True, terms={(1, 0): 5, (2, 0): 7}, pivot_label=(0, 0),
        all_terms_lf=True, non_lf_terms=[], unknown_lf_terms=[], nrows=3, rank=1,
    )
    rec = record_from_result(res)
    assert isinstance(rec, NormalFormRecord)
    assert rec.coeffs == {(1, 0): 5, (2, 0): 7}
    assert rec.support == ((1, 0), (2, 0))
    assert rec.status == "Reduced" and rec.formal_success is True
    assert rec.all_terms_lf is True and rec.rank == 1
    assert rec.diagnostics == {"nrows": 3, "pivot_label": (0, 0)}


def test_collect_records_over_grid_all_reduced():
    fam = _fam()
    row = algebraic_row(fam, zero_label(1, 1), 0)  # J(0,0)-J(0,-1)-J(1,-1)=0
    samples = _samples([2, 3, 4])
    records = collect_normal_form_records(fam, [row], (0, 0), PRIMES, samples)
    assert len(records) == len(PRIMES) * len(samples)
    assert all(r.status == "Reduced" and r.formal_success for r in records)
    assert all(r.coeffs == {(0, -1): 1, (1, -1): 1} for r in records)
    assert summarize_records(records)["reduced"] == len(records)


def test_bad_specialization_recorded_not_dropped():
    fam = _fam()
    row = Row("test", {}, {
        (0, 0): ParamExpr.from_int(1, ("ep",)),
        (1, 0): ParamExpr.from_sympy(sp.sympify("1/ep"), ("ep",)),  # singular at ep=0
    })
    records = collect_normal_form_records(fam, [row], (0, 0), PRIMES, _samples([0, 5]))
    assert len(records) == 2 * len(PRIMES)  # nothing dropped
    for r in records:
        if r.sample["ep"] == Fraction(0):
            assert r.status == "BadSpecialization" and r.coeffs == {} and not r.formal_success
        else:
            assert r.status == "Reduced" and r.formal_success


def test_target_not_reducible_recorded():
    fam = _fam()
    row = algebraic_row(fam, zero_label(1, 1), 0)  # labels (0,0),(0,-1),(1,-1)
    records = collect_normal_form_records(fam, [row], (9, 9), PRIMES[:1], _samples([3]))
    assert len(records) == 1
    assert records[0].status == "TargetNotReducible" and not records[0].formal_success


def test_records_feed_univariate_reconstruction():
    fam = _fam()
    ep = sp.Symbol("ep")
    row = _coeff_row("(2*ep-1)/(ep+3)")  # J(0,0) = C(ep) J(1,0)
    records = collect_normal_form_records(fam, [row], (0, 0), PRIMES, _samples(range(2, 12)))
    coeffs = reconstruct_coefficients(records, ["ep"])
    assert sp.simplify(coeffs[(1, 0)] - (2 * ep - 1) / (ep + 3)) == 0


def test_skipped_bad_sample_does_not_corrupt_reconstruction():
    fam = _fam()
    ep = sp.Symbol("ep")
    row = _coeff_row("(2*ep-1)/(ep+3)")  # denominator vanishes at ep=-3
    samples = _samples(list(range(2, 12)) + [-3])
    records = collect_normal_form_records(fam, [row], (0, 0), PRIMES, samples)
    bad = [r for r in records if r.sample["ep"] == Fraction(-3)]
    assert bad and all(r.status == "BadSpecialization" for r in bad)
    # bad points are skipped (not dropped from the list, not patched); good points still recover
    coeffs = reconstruct_coefficients(records, ["ep"])
    assert sp.simplify(coeffs[(1, 0)] - (2 * ep - 1) / (ep + 3)) == 0


def test_union_support_zero_fill_through_collector():
    fam = _fam()
    ep = sp.Symbol("ep")
    row = _coeff_row("ep - 5")  # coefficient is exactly 0 at ep=5 -> term absent there
    samples = _samples(range(1, 10))
    records = collect_normal_form_records(fam, [row], (0, 0), PRIMES, samples)
    recs5 = [r for r in records if r.sample["ep"] == Fraction(5)]
    assert recs5 and all(r.status == "Reduced" and r.support == () for r in recs5)
    # union support across the other samples fills 0 at ep=5, recovering ep-5 exactly
    coeffs = reconstruct_coefficients(records, ["ep"])
    assert sp.simplify(coeffs[(1, 0)] - (ep - 5)) == 0
