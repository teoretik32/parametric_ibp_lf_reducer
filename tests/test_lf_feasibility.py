"""Tests for the LF-constrained feasibility span test (Method.1, External Int2).

Hand-built row systems with known answers: the test asserts the projection semantics
(allowed = LF-True only, Unknown/missing = forbidden), ranking-independence, determinism,
BadSpecialization reporting, and that extracted coefficients certify via the standard
row-span certificate.
"""

from __future__ import annotations

import sympy as sp
from parametric_ibp_lf_reducer import (
    STATUS_IN_SPAN,
    ParamExpr,
    Row,
    feasibility_to_payload,
    lf_reduction_coefficients_mod_p,
    lf_reduction_feasible_mod_p,
    verify_reduction_relation_mod_p,
)

EP = ("ep",)
PRIME = 30011
SAMPLE = {"ep": 5}

T = (0, 0)
A = (-1, 0)
B = (0, -1)
C = (-2, -1)  # negative n-shift label (generator parity)


def _row(terms: dict) -> Row:
    return Row("test", {}, {lab: ParamExpr.from_int(c, EP) for lab, c in terms.items()})


def test_feasible_through_lf_label_only():
    # T = A and A = B; A is forbidden (non-LF), B is allowed -> the LF-only path T = B exists
    # even though a ranking keeping A as master would report a non-LF normal form.
    rows = [_row({T: 1, A: -1}), _row({A: 1, B: -1})]
    flags = {T: False, A: False, B: True}
    res = lf_reduction_feasible_mod_p(rows, [T, A, B], T, flags, SAMPLE, PRIME)
    assert res.status == "Feasible"
    assert res.residual_support == ()
    assert res.n_allowed == 1
    assert res.n_forbidden >= 1


def test_obstructed_when_target_needs_forbidden_label():
    rows = [_row({T: 1, A: -1})]  # the only relation goes through forbidden A
    flags = {T: False, A: False, B: True}
    res = lf_reduction_feasible_mod_p(rows, [T, A, B], T, flags, SAMPLE, PRIME)
    assert res.status == "Obstructed"
    assert res.residual_support  # a nonzero canonical residual is the witness


def test_unknown_and_missing_flags_are_forbidden():
    rows = [_row({T: 1, A: -1}), _row({A: 1, B: -1})]
    # B Unknown -> forbidden -> no allowed labels at all -> Obstructed.
    res_unknown = lf_reduction_feasible_mod_p(
        rows, [T, A, B], T, {T: False, A: False, B: "Unknown"}, SAMPLE, PRIME
    )
    assert res_unknown.status == "Obstructed"
    assert res_unknown.n_allowed == 0
    # B missing from lf_flags entirely -> same conservative answer.
    res_missing = lf_reduction_feasible_mod_p(
        rows, [T, A, B], T, {T: False, A: False}, SAMPLE, PRIME
    )
    assert res_missing.status == "Obstructed"


def test_ranking_independence_via_column_order_override():
    rows = [_row({T: 1, A: -1}), _row({A: 1, B: -1})]
    flags = {T: False, A: False, B: True}
    res_default = lf_reduction_feasible_mod_p(rows, [T, A, B], T, flags, SAMPLE, PRIME)
    res_override = lf_reduction_feasible_mod_p(
        rows, [T, A, B], T, flags, SAMPLE, PRIME, column_order=[T, A]
    )
    assert res_default.status == res_override.status == "Feasible"


def test_negative_shift_label_feasible():
    rows = [_row({T: 1, C: -1})]
    res = lf_reduction_feasible_mod_p(rows, [T, C], T, {T: False, C: True}, SAMPLE, PRIME)
    assert res.status == "Feasible"


def test_bad_specialization_is_reported_not_raised():
    singular = Row("test", {}, {T: ParamExpr.from_sympy(sp.sympify("1/ep"), EP), A: ParamExpr.from_int(1, EP)})
    res = lf_reduction_feasible_mod_p(
        [singular], [T, A], T, {A: True}, {"ep": 0}, PRIME
    )
    assert res.status == "BadSpecialization"
    assert res.detail


def test_determinism_and_payload_roundtrip():
    import json

    rows = [_row({T: 1, A: -1}), _row({A: 1, B: -1})]
    flags = {T: False, A: False, B: True}
    r1 = lf_reduction_feasible_mod_p(rows, [T, A, B], T, flags, SAMPLE, PRIME)
    r2 = lf_reduction_feasible_mod_p(rows, [T, A, B], T, flags, SAMPLE, PRIME)
    assert r1 == r2
    p1 = feasibility_to_payload(r1)
    assert json.loads(json.dumps(p1, sort_keys=True)) == p1


def test_coefficients_certify_via_row_span_certificate():
    # T = A, A = B, plus a parametric spectator relation to make the combination nontrivial.
    rows = [
        _row({T: 1, A: -1}),
        _row({A: 1, B: -1}),
        Row("test", {}, {A: ParamExpr.from_sympy(sp.sympify("ep"), EP), C: ParamExpr.from_int(1, EP)}),
    ]
    flags = {T: False, A: False, B: True, C: True}
    res, coeffs = lf_reduction_coefficients_mod_p(rows, [T, A, B, C], T, flags, SAMPLE, PRIME)
    assert res.status == "Feasible"
    assert coeffs  # nonempty explicit reduction
    assert set(coeffs) <= {B, C}  # support only on allowed labels
    cert = verify_reduction_relation_mod_p(None, rows, T, coeffs, SAMPLE, PRIME)
    assert cert.status == STATUS_IN_SPAN


def test_coefficients_match_feasibility_statuses():
    rows = [_row({T: 1, A: -1})]
    flags = {T: False, A: False, B: True}
    res, coeffs = lf_reduction_coefficients_mod_p(rows, [T, A, B], T, flags, SAMPLE, PRIME)
    assert res.status == "Obstructed"
    assert coeffs == {}
