"""Tests for the all-row-support LF feasibility mode (Method Audit.1 Phase A).

Regression target: `lf_reduction_feasible_mod_p` builds its allowed set from the seed
box only, so an LF-True row-support column OUTSIDE the seed box is silently forbidden
and a genuinely feasible reduction is reported `Obstructed`. The all-support mode must
recover `Feasible` on such a fixture WITHOUT generating any new rows.
"""

from fractions import Fraction

from parametric_ibp_lf_reducer.coefficients import ParamExpr
from parametric_ibp_lf_reducer.labels import Label
from parametric_ibp_lf_reducer.lf_feasibility import (
    STATUS_FEASIBLE,
    STATUS_OBSTRUCTED,
    classify_row_support_lf,
    collect_row_support_labels,
    lf_reduction_feasible_all_support_mod_p,
    lf_reduction_feasible_mod_p,
)
from parametric_ibp_lf_reducer.row_generation import Row

PRIME = 2147483647
SAMPLE = {"ep": Fraction(15, 7), "r": Fraction(32, 11)}

T = Label((2, 1, 1, 0))  # target
A = Label((1, 1, 1, 0))  # in seed box, LF-True
B = Label((1, 1, 0, 0))  # in seed box, LF-False
X = Label((0, 1, 1, 0))  # OUTSIDE seed box, LF-True — the audited defect column
U = Label((0, 0, 1, 0))  # outside seed box, LF-Unknown

PARAMS = ("ep", "r")


def _c(k: int) -> ParamExpr:
    return ParamExpr.from_int(k, PARAMS)


def _row(terms) -> Row:
    return Row(kind="test", provenance={"origin": "test_lf_feasibility_all_support"}, terms=terms)


# Single row: T - A - X = 0 (constant coefficients keep the fixture exact).
ROW_TX = _row({T: _c(1), A: _c(-1), X: _c(-1)})
# Row touching only forbidden support: T - B = 0.
ROW_TB = _row({T: _c(1), B: _c(-1)})

SEED_BOX = [T, A, B]
SEED_LF = {T: True, A: True, B: False, X: True, U: "Unknown"}


def test_collect_row_support_labels_sorted_union():
    assert collect_row_support_labels([ROW_TX, ROW_TB]) == tuple(sorted({T, A, B, X}))
    assert collect_row_support_labels([]) == ()


def test_classify_row_support_lf_uses_known_flags_without_family():
    verdicts = classify_row_support_lf(None, [ROW_TX], known_lf_flags=SEED_LF)
    assert verdicts == {T: True, A: True, X: True}


def test_classify_row_support_lf_missing_flag_without_family_raises():
    try:
        classify_row_support_lf(None, [ROW_TX], known_lf_flags={T: True, A: True})
    except ValueError as exc:
        assert "support label" in str(exc)
    else:
        raise AssertionError("expected ValueError for uncovered support label")


def test_seed_only_mode_obstructed_by_out_of_box_lf_true_column():
    # Production behaviour (unchanged): X is not in the seed box => forbidden => Obstructed.
    res = lf_reduction_feasible_mod_p(
        [ROW_TX], SEED_BOX, T, {lab: SEED_LF[lab] for lab in SEED_BOX}, SAMPLE, PRIME
    )
    assert res.status == STATUS_OBSTRUCTED


def test_all_support_mode_recovers_feasible_same_rows():
    res, verdicts = lf_reduction_feasible_all_support_mod_p(
        None, [ROW_TX], T, SEED_LF, SAMPLE, PRIME
    )
    assert res.status == STATUS_FEASIBLE
    assert verdicts[X] is True  # the recovered column, LF-True but out-of-box


def test_all_support_mode_keeps_false_and_unknown_forbidden():
    # T = B with B LF-False must stay Obstructed even in all-support mode.
    res, _ = lf_reduction_feasible_all_support_mod_p(None, [ROW_TB], T, SEED_LF, SAMPLE, PRIME)
    assert res.status == STATUS_OBSTRUCTED
    # Same with an Unknown support column: T = U.
    row_tu = _row({T: _c(1), U: _c(-1)})
    res_u, verdicts_u = lf_reduction_feasible_all_support_mod_p(
        None, [row_tu], T, SEED_LF, SAMPLE, PRIME
    )
    assert res_u.status == STATUS_OBSTRUCTED
    assert verdicts_u[U] == "Unknown"


def test_all_support_mode_target_absent_from_rows_gets_flag_not_crash():
    # Rows that never touch the target: reduction is vacuously infeasible, no KeyError.
    row_ax = _row({A: _c(1), X: _c(-1)})
    res, verdicts = lf_reduction_feasible_all_support_mod_p(
        None, [row_ax], T, SEED_LF, SAMPLE, PRIME
    )
    assert T in verdicts
    assert res.status == STATUS_OBSTRUCTED
