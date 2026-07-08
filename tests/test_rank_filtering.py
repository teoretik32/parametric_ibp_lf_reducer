"""Tests for rank-consistency record selection before reconstruction (Pass D4.3).

D4.2 showed that rank-deficient sample points (whose RREF rank drops below the generic rank)
produce normal-form records with a shrunken/shifted support; union-support 0-fill then poisons
interpolation. These tests pin the fix: only ``Reduced``/formal-success records at the maximal
observed rank feed reconstruction by default (``rank_policy="max_rank"``), rank-deficient records
are skipped + counted, and a *special zero* inside a max-rank record is still an honest zero.
``rank_policy="all"`` keeps the old behaviour for tests/debugging. Nothing here weakens the LF
gate or hardcodes validation cases.
"""

from __future__ import annotations

from fractions import Fraction

import pytest
import sympy as sp

from parametric_ibp_lf_reducer import (
    RANK_POLICY_ALL,
    InterpolationFailed,
    ParamExpr,
    parse_family_text,
    reconstruct_coefficients,
    reduce_rows_once,
    select_records_for_reconstruction,
)
from parametric_ibp_lf_reducer.modular_normal_form import NormalFormResult
from parametric_ibp_lf_reducer.result import STATUS_SUCCESS
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


def _to_res(value: Fraction, prime: int) -> int:
    return value.numerator % prime * pow(value.denominator % prime, -1, prime) % prime


def _reduced(sample: dict, prime: int, coeffs: dict, rank: int) -> NormalFormResult:
    return NormalFormResult(
        status="Reduced", target_label=(0, 0), prime=prime, sample=dict(sample),
        formal_success=True, terms=dict(coeffs), rank=rank,
    )


def _records_for_function(f, samples, primes, label=(1, 0), rank=2, param="ep"):
    """Reduced records for ``target = f(param) * J[label]`` with an explicit RREF rank."""
    ep = sp.Symbol(param)
    records = []
    for s in samples:
        val = sp.Rational(f.subs(ep, s))
        frac = Fraction(int(val.p), int(val.q))
        for p in primes:
            res = _to_res(frac, p)
            coeffs = {label: res} if res else {}
            records.append(_reduced({param: Fraction(s)}, p, coeffs, rank))
    return records


# --- selection helper unit tests --------------------------------------------------------------
def test_select_records_max_rank_policy_and_diagnostics():
    good = [_reduced({"ep": Fraction(s)}, PRIMES[0], {(1, 0): 5}, rank=7) for s in (2, 3, 4)]
    low = [_reduced({"ep": Fraction(9)}, PRIMES[0], {(9, 9): 1}, rank=4)]
    invalid = [
        NormalFormResult(status="BadSpecialization", target_label=(0, 0), prime=PRIMES[0],
                         sample={"ep": Fraction(0)}, formal_success=False),
        NormalFormResult(status="TargetNotReducible", target_label=(0, 0), prime=PRIMES[0],
                         sample={"ep": Fraction(1)}, formal_success=False, rank=9),
    ]
    selected, diag = select_records_for_reconstruction(good + low + invalid)
    assert selected == good  # only max-rank valid records; invalid never coefficient records
    assert diag["rank_policy"] == "max_rank"
    assert diag["n_records"] == 6
    assert diag["n_valid_records_before_rank_filter"] == 4
    assert diag["selected_rank"] == 7
    assert diag["n_selected_records"] == 3
    assert diag["n_rank_filtered_records"] == 1
    assert diag["rank_histogram"] == {4: 1, 7: 3}  # invalid records never enter the histogram
    assert diag["support_after_rank_filter"] == ((1, 0),)  # low-rank support (9,9) excluded


def test_select_records_all_policy_keeps_every_valid_record():
    good = [_reduced({"ep": Fraction(2)}, PRIMES[0], {(1, 0): 5}, rank=7)]
    low = [_reduced({"ep": Fraction(9)}, PRIMES[0], {(9, 9): 1}, rank=4)]
    selected, diag = select_records_for_reconstruction(good + low, rank_policy=RANK_POLICY_ALL)
    assert selected == good + low
    assert diag["selected_rank"] is None
    assert diag["n_rank_filtered_records"] == 0
    assert diag["rank_histogram"] == {4: 1, 7: 1}
    assert diag["support_after_rank_filter"] == ((1, 0), (9, 9))


def test_select_records_unknown_policy_rejected():
    with pytest.raises(ValueError):
        select_records_for_reconstruction([], rank_policy="best_effort")


def test_select_records_no_valid_records():
    invalid = [
        NormalFormResult(status="BadSpecialization", target_label=(0, 0), prime=PRIMES[0],
                         sample={"ep": Fraction(0)}, formal_success=False),
    ]
    selected, diag = select_records_for_reconstruction(invalid)
    assert selected == []
    assert diag["selected_rank"] is None
    assert diag["rank_histogram"] == {}
    assert diag["support_after_rank_filter"] == ()


# --- reconstruction behaviour -----------------------------------------------------------------
def test_low_rank_records_are_filtered_but_poison_the_all_policy():
    """High-rank records reconstruct exactly; injected low-rank records with shrunk support are
    skipped by the default policy, while the old ``"all"`` behaviour is poisoned by them."""
    ep = sp.Symbol("ep")
    f = (2 * ep - 1) / (ep + 3)
    records = _records_for_function(f, range(2, 12), PRIMES, label=(1, 0), rank=2)
    # a rank-deficient point: different support {(9,9)} and a wrong "coefficient" value 7
    for p in PRIMES:
        records.append(_reduced({"ep": Fraction(13)}, p, {(9, 9): 7}, rank=1))

    coeffs = reconstruct_coefficients(records, ["ep"])  # default: max_rank
    assert set(coeffs) == {(1, 0)}  # the low-rank support never enters the union support
    assert sp.simplify(coeffs[(1, 0)] - f) == 0

    # old behaviour: the shrunken-support point 0-fills (1,0) at ep=13 and injects (9,9)=7,
    # so the value sequences are incoherent as rational functions -> honest refusal
    with pytest.raises(InterpolationFailed):
        reconstruct_coefficients(records, ["ep"], rank_policy=RANK_POLICY_ALL)


def test_special_zero_in_max_rank_record_is_still_zero_not_skipped():
    """A coefficient that vanishes at a *max-rank* point stays an honest 0-fill (union support)."""
    ep = sp.Symbol("ep")
    f = ep - 5  # exact zero at ep=5 -> term absent from that record's support
    records = _records_for_function(f, range(1, 10), PRIMES, label=(1, 0), rank=2)
    selected, diag = select_records_for_reconstruction(records)
    assert diag["n_rank_filtered_records"] == 0  # the special-zero record is NOT rank-filtered
    assert len(selected) == len(records)
    coeffs = reconstruct_coefficients(records, ["ep"])
    assert sp.simplify(coeffs[(1, 0)] - f) == 0


def test_uniform_rank_behavior_unchanged():
    """If all valid records share one rank, max_rank selects everything = old behaviour."""
    ep = sp.Symbol("ep")
    f = (3 * ep + 1) / ep
    records = _records_for_function(f, range(2, 10), PRIMES, label=(1, 0), rank=5)
    selected, diag = select_records_for_reconstruction(records)
    assert selected == records and diag["n_rank_filtered_records"] == 0
    new = reconstruct_coefficients(records, ["ep"])
    old = reconstruct_coefficients(records, ["ep"], rank_policy=RANK_POLICY_ALL)
    assert set(new) == set(old) == {(1, 0)}
    assert sp.simplify(new[(1, 0)] - old[(1, 0)]) == 0


# --- reducer end-to-end: a real rank-deficient sample point -----------------------------------
def test_reducer_filters_rank_deficient_sample_end_to_end():
    """A genuine rank drop (a row vanishing at ep=2) shrinks that sample's system; the reducer
    must skip it via the max-rank selection and report the filter in its diagnostics."""
    fam = parse_family_text(ONE_VAR)
    ep_syms = ("ep",)
    T, A, B = (0, 0), (1, 0), (2, 0)
    row1 = Row("test", {}, {
        T: ParamExpr.from_int(1, ep_syms),
        A: ParamExpr.from_int(-1, ep_syms),
        B: ParamExpr.from_int(-1, ep_syms),
    })  # T = A + B
    row2 = Row("test", {}, {
        B: ParamExpr.from_sympy(sp.sympify("ep - 2"), ep_syms),
        A: ParamExpr.from_sympy(sp.sympify("-(ep - 2)"), ep_syms),
    })  # (ep-2) * (B - A) = 0: vanishes entirely at ep=2 -> rank drops 2 -> 1
    samples = [{"ep": Fraction(v)} for v in (3, 4, 5, 6, 7, 2)]
    res = reduce_rows_once(
        fam, T, [T, A, B], [row1, row2], PRIMES[:2], samples,
        lf_flags={T: False, A: True, B: True},
    )
    # generic points give T = 2*A (rank 2); ep=2 gives T = A + B (rank 1) and must be filtered
    assert res.status == STATUS_SUCCESS
    assert res.all_locally_finite is True
    assert [t.label for t in res.terms] == [A]
    assert sp.simplify(sp.sympify(res.terms[0].coefficient_text) - 2) == 0

    ex = res.diagnostics.extra
    assert ex["n_reduced_records"] == 12  # every point reduced, incl. the deficient one
    assert ex["n_selected_records"] == 10
    sel = ex["record_selection"]
    assert sel["rank_histogram"] == {1: 2, 2: 10}
    assert sel["selected_rank"] == 2
    assert sel["n_rank_filtered_records"] == 2
    assert sel["n_valid_records_before_rank_filter"] == 12
    assert sel["support_after_rank_filter"] == (A,)
    assert any("rank filter" in m for m in res.diagnostics.messages)
