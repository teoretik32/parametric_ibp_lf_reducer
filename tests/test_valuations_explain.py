"""Tests for the Method.1 explainable local-finiteness audit (explain_local_finiteness).

The report must agree with ``is_locally_finite`` verbatim (it delegates and cross-checks), and
its recommendations must be honest: every recommended unit shift strictly increases
``base_score`` along every failing ray (recomputed here directly, not via the report).
"""

from __future__ import annotations

from conftest import load_example
from parametric_ibp_lf_reducer import (
    base_score,
    explain_local_finiteness,
    is_locally_finite,
    make_label,
    parse_family_text,
    report_to_payload,
    zero_label,
)
from test_valuations import (
    BULK_SINGULAR,
    ONE_VAR_CONVERGENT,
    ONE_VAR_EPS_REGULATED,
    ONE_VAR_LOG_DIVERGENT,
    ONE_VAR_SYMBOLIC_EXPONENT,
)


def _apply_shift(label, shift):
    return tuple(a + b for a, b in zip(label, shift))


def test_verdict_matches_is_locally_finite_matrix():
    cases = [
        (ONE_VAR_CONVERGENT, zero_label(1, 1)),
        (ONE_VAR_LOG_DIVERGENT, zero_label(1, 1)),
        (ONE_VAR_EPS_REGULATED, zero_label(1, 1)),
        (ONE_VAR_SYMBOLIC_EXPONENT, zero_label(1, 1)),
        (BULK_SINGULAR, zero_label(2, 2)),
    ]
    for text, lab in cases:
        fam = parse_family_text(text)
        report = explain_local_finiteness(fam, lab)
        assert report.verdict == is_locally_finite(fam, lab)
        assert report.label == lab


def test_d4_masters_and_target_verdicts():
    fam = parse_family_text(load_example("d4_explicit_family.wl.txt"))
    m1 = make_label([0, 1, 1, 0], [-2, -1, 0])
    rep_m1 = explain_local_finiteness(fam, m1)
    assert rep_m1.verdict is True
    assert rep_m1.failing_rays == ()
    assert rep_m1.recommended_shifts == ()

    rep_t = explain_local_finiteness(fam, zero_label(4, 3))
    assert rep_t.verdict is False
    assert rep_t.failing_rays  # at least one nonpos witness ray
    # Every recommended shift strictly improves base_score on every failing ray.
    for rec in rep_t.recommended_shifts:
        assert rec.improves_all
        shifted = _apply_shift(zero_label(4, 3), rec.shift)
        for rv in rep_t.failing_rays:
            assert base_score(fam, shifted, rv.ray) - base_score(fam, zero_label(4, 3), rv.ray) > 0


def test_log_divergent_failing_ray_is_marginal_at_infinity():
    fam = parse_family_text(ONE_VAR_LOG_DIVERGENT)
    report = explain_local_finiteness(fam, zero_label(1, 1))
    assert report.verdict is False
    dirs = {rv.ray.direction for rv in report.failing_rays}
    assert (-1,) in dirs
    marginal = [rv for rv in report.failing_rays if rv.ray.direction == (-1,)]
    assert marginal[0].score == 0
    assert "STRICT RULE" in marginal[0].detail
    # Improving shifts: n_x -= 1 (delta = 1) and m_G0 -= 1 (delta = 1 since val_{-1}(G0) = -1).
    rec_shifts = {rec.shift for rec in report.recommended_shifts}
    assert (-1, 0) in rec_shifts
    assert (0, -1) in rec_shifts
    for rec in report.recommended_shifts:
        shifted = _apply_shift(zero_label(1, 1), rec.shift)
        for rv in report.failing_rays:
            assert base_score(fam, shifted, rv.ray) - base_score(fam, zero_label(1, 1), rv.ray) > 0


def test_eps_regulated_zero_score_is_failing():
    fam = parse_family_text(ONE_VAR_EPS_REGULATED)
    report = explain_local_finiteness(fam, zero_label(1, 1))
    assert report.verdict is False
    assert any(rv.score == 0 for rv in report.failing_rays)


def test_symbolic_exponent_reports_unknown_rays():
    fam = parse_family_text(ONE_VAR_SYMBOLIC_EXPONENT)
    report = explain_local_finiteness(fam, zero_label(1, 1))
    assert report.verdict == "Unknown"
    assert report.failing_rays == ()
    assert report.unknown_rays  # the undecidable rays are reported, not swallowed
    assert any("undecidable" in note for note in report.notes)


def test_bulk_singular_reports_bulk_note():
    fam = parse_family_text(BULK_SINGULAR)
    report = explain_local_finiteness(fam, zero_label(2, 2))
    assert report.verdict == "Unknown"
    assert report.failing_rays == ()
    assert report.unknown_rays == ()
    assert report.bulk_safe is False
    assert any("bulk" in note for note in report.notes)


def test_convergent_true_report_shape():
    fam = parse_family_text(ONE_VAR_CONVERGENT)
    report = explain_local_finiteness(fam, zero_label(1, 1))
    assert report.verdict is True
    assert report.failing_rays == ()
    assert report.unknown_rays == ()
    assert report.bulk_safe is True
    # All 2*(nvars+npolys) unit shifts are tabulated even when nothing fails.
    assert len(report.shift_deltas) == 2 * (1 + 1)
    assert report.recommended_shifts == ()


def test_payload_is_deterministic_and_json_safe():
    import json

    fam = parse_family_text(ONE_VAR_LOG_DIVERGENT)
    p1 = report_to_payload(explain_local_finiteness(fam, zero_label(1, 1)))
    p2 = report_to_payload(explain_local_finiteness(fam, zero_label(1, 1)))
    assert p1 == p2
    s = json.dumps(p1, sort_keys=True)
    assert json.loads(s) == p1
    assert p1["verdict"] is False
    assert all(isinstance(r["score"], (str, type(None))) for r in p1["rays"])
