"""Tests for the dual LF-obstruction certificate (Method.6, External Int2).

Hand-built row systems with known answers exercise: no witness for Feasible systems; a valid
right-nullspace witness with exact ``A_projected @ w == 0`` and ``w[target] == 1`` for Obstructed
systems (both the target-free-column and target-pivot cases); the codimension-one refutation
(``residual_support == (target,)`` with nullity 2); determinism / JSON-safety; the row-pairing
helper (breaking vs annihilating candidate rows); BadSpecialization reporting; agreement with
``lf_reduction_feasible_mod_p``; and a golden guard that the reused feasibility semantics are
unchanged.

The witness helper ``test_rows_against_obstruction_witness`` matches pytest's ``test_*`` pattern,
so this module imports the witness MODULE as ``low`` and never ``from ... import`` that helper.
"""

from __future__ import annotations

import json

import sympy as sp
import parametric_ibp_lf_reducer.lf_obstruction_witness as low
from parametric_ibp_lf_reducer import (
    ParamExpr,
    Row,
    lf_reduction_feasible_mod_p,
)

EP = ("ep",)
PRIME = 30011
SAMPLE = {"ep": 5}

T = (0, 0)
A = (1, 0)
F1 = (1, 0)
F2 = (2, 0)


def _row(terms: dict) -> Row:
    return Row("test", {}, {lab: ParamExpr.from_int(c, EP) for lab, c in terms.items()})


def test_feasible_system_no_witness():
    # T = A with A allowed (LF-True): projected system is {T:1} alone -> feasible, no witness.
    rows = [_row({T: 1, A: 1})]
    res = low.lf_obstruction_witness_mod_p(rows, [T, A], T, {A: True}, SAMPLE, PRIME)
    assert res.status == low.STATUS_FEASIBLE
    assert res.witness == ()
    feas = lf_reduction_feasible_mod_p(rows, [T, A], T, {A: True}, SAMPLE, PRIME)
    assert feas.status == "Feasible"


def test_obstructed_target_free_column():
    rows = [_row({F1: 1, F2: 1})]
    res = low.lf_obstruction_witness_mod_p(rows, [T, F1, F2], T, {}, SAMPLE, PRIME)
    assert res.status == low.STATUS_WITNESS
    assert res.witness == ((T, 1),)
    assert res.check_annihilation is True
    assert res.check_target_unit is True


def test_obstructed_target_pivot_column():
    rows = [_row({T: 1, F1: 2})]
    res = low.lf_obstruction_witness_mod_p(rows, [T, F1], T, {}, SAMPLE, PRIME)
    assert res.status == low.STATUS_WITNESS
    inv_neg2 = pow(-2, PRIME - 2, PRIME)
    assert dict(res.witness) == {T: 1, F1: inv_neg2}
    # exact row . w == 0 : 1*1 + 2*inv(-2) == 0
    assert res.check_annihilation is True
    assert res.check_target_unit is True


def test_nullity_gt_one_target_only_residual():
    # codimension-one refutation in miniature.
    labels = [(0, 0), (1, 0), (2, 0)]
    rows = [_row({(1, 0): 1, (2, 0): 1})]
    feas = lf_reduction_feasible_mod_p(rows, labels, (0, 0), {}, SAMPLE, PRIME)
    assert feas.status == "Obstructed"
    assert feas.residual_support == ((0, 0),)  # target only
    res = low.lf_obstruction_witness_mod_p(rows, labels, (0, 0), {}, SAMPLE, PRIME)
    assert res.status == low.STATUS_WITNESS
    assert res.witness == (((0, 0), 1),)
    assert res.rank == 1
    assert res.n_projected_cols == 3
    assert res.nullity == 2  # NOT one — residual_support=[target] does not bound the quotient
    assert res.check_annihilation and res.check_target_unit


def test_witness_deterministic_json_safe():
    labels = [(0, 0), (1, 0), (2, 0)]
    r1 = low.lf_obstruction_witness_mod_p(
        [_row({(1, 0): 1, (2, 0): 1})], labels, (0, 0), {}, SAMPLE, PRIME
    )
    r2 = low.lf_obstruction_witness_mod_p(
        [_row({(1, 0): 1, (2, 0): 1})], labels, (0, 0), {}, SAMPLE, PRIME
    )
    assert r1 == r2
    payload = low.witness_to_payload(r1)
    assert json.loads(json.dumps(payload, sort_keys=True)) == payload
    assert low.witness_from_payload(payload) == r1
    for _lab, coeff in r1.witness:
        assert 0 <= coeff < PRIME


def test_target_pivot_free_tiebreak_stable():
    rows = [_row({T: 1, F1: 1, F2: 1})]
    base = low.lf_obstruction_witness_mod_p(rows, [T, F1, F2], T, {}, SAMPLE, PRIME)
    assert base.status == low.STATUS_WITNESS
    # chosen free column is F1 (first in sorted order); F2 stays a pivot -> not in witness support
    assert set(dict(base.witness)) == {T, F1}
    # permuting labels / duplicating rows leaves the witness unchanged (determinism rules 2 & 4)
    permuted = low.lf_obstruction_witness_mod_p(rows, [F2, F1, T], T, {}, SAMPLE, PRIME)
    duplicated = low.lf_obstruction_witness_mod_p(
        [_row({T: 1, F1: 1, F2: 1}), _row({T: 1, F1: 1, F2: 1})], [T, F1, F2], T, {}, SAMPLE, PRIME
    )
    assert permuted.witness == base.witness
    assert duplicated.witness == base.witness


def test_pairing_detects_breaking_row():
    labels = [(0, 0), (1, 0), (2, 0)]
    wit = low.lf_obstruction_witness_mod_p(
        [_row({(1, 0): 1, (2, 0): 1})], labels, (0, 0), {}, SAMPLE, PRIME
    )
    (pairing,) = low.test_rows_against_obstruction_witness(
        [_row({(0, 0): 5})], wit, SAMPLE, PRIME
    )
    assert pairing.pairing == 5
    assert pairing.breaks is True


def test_pairing_flags_annihilating_rows():
    labels = [(0, 0), (1, 0), (2, 0)]
    wit = low.lf_obstruction_witness_mod_p(
        [_row({(1, 0): 1, (2, 0): 1})], labels, (0, 0), {}, SAMPLE, PRIME
    )
    candidates = [_row({(1, 0): 1, (2, 0): -1}), _row({(1, 0): 7})]
    pairings = low.test_rows_against_obstruction_witness(candidates, wit, SAMPLE, PRIME)
    assert [p.pairing for p in pairings] == [0, 0]
    assert all(p.breaks is False for p in pairings)


def test_bad_specialization_status():
    singular = Row(
        "test",
        {},
        {T: ParamExpr.from_sympy(sp.sympify("1/ep"), EP), A: ParamExpr.from_int(1, EP)},
    )
    res = low.lf_obstruction_witness_mod_p([singular], [T, A], T, {A: True}, {"ep": 0}, PRIME)
    assert res.status == low.STATUS_BAD_SPECIALIZATION
    assert res.witness == ()
    assert res.detail


def test_witness_feasibility_agreement():
    cases = [
        ([_row({T: 1, A: 1})], [T, A], {A: True}),  # feasible
        ([_row({F1: 1, F2: 1})], [T, F1, F2], {}),  # obstructed (target free)
        ([_row({T: 1, F1: 2})], [T, F1], {}),  # obstructed (target pivot)
        ([_row({(1, 0): 1, (2, 0): 1})], [(0, 0), (1, 0), (2, 0)], {}),  # obstructed nullity 2
    ]
    for rows, labels, flags in cases:
        target = labels[0]
        wit = low.lf_obstruction_witness_mod_p(rows, labels, target, flags, SAMPLE, PRIME)
        feas = lf_reduction_feasible_mod_p(rows, labels, target, flags, SAMPLE, PRIME)
        assert (wit.status == low.STATUS_WITNESS) == (feas.status == "Obstructed")
        assert (wit.status == low.STATUS_FEASIBLE) == (feas.status == "Feasible")
        assert wit.rank == feas.rank


def test_existing_reducer_results_unchanged():
    # golden mini-check: the reused feasibility path keeps its pre-Method.6 verdicts/residuals.
    feasible = lf_reduction_feasible_mod_p([_row({T: 1, A: 1})], [T, A], T, {A: True}, SAMPLE, PRIME)
    assert feasible.status == "Feasible"
    assert feasible.residual_support == ()
    obstructed = lf_reduction_feasible_mod_p([_row({T: 1, A: 1})], [T, A], T, {A: False}, SAMPLE, PRIME)
    assert obstructed.status == "Obstructed"
    assert obstructed.residual_support == (T,)
