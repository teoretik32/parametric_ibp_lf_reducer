"""Tests for label ranking (Pass 2D): elimination order, no label removal, non-LF before LF."""

from __future__ import annotations

from parametric_ibp_lf_reducer import (
    is_locally_finite,
    parse_family_text,
    rank_labels,
)

# G0 = 1 + x, integrand x^0 G0^-2.
#   (0,0) LF ; (0,-1) LF ; (2,-2) LF (complex) ; (0,1) non-LF (simple) ; (1,1) non-LF.
ONE_VAR = """
IBPInput = <|
  "Variables" -> {x}, "Parameters" -> {ep}, "Regulators" -> {ep},
  "Polynomials" -> <| "G0" -> 1 + x |>,
  "MonomialExponents" -> <| x -> 0 |>,
  "PolynomialExponents" -> <| "G0" -> -2 |>
|>
"""

LABELS = [(0, 0), (0, -1), (2, -2), (0, 1), (1, 1)]


def _fam():
    return parse_family_text(ONE_VAR)


def test_ranking_never_removes_labels():
    fam = _fam()
    ranked = rank_labels(fam, LABELS)
    assert sorted(ranked.ordered) == sorted(set(LABELS))  # a permutation, nothing dropped
    assert len(ranked.ordered) == len(LABELS)


def test_target_label_is_eliminated_first():
    fam = _fam()
    ranked = rank_labels(fam, LABELS, target=(0, -1))
    assert ranked.ordered[0] == (0, -1)
    assert ranked.tiers[(0, -1)] == 0


def test_all_non_lf_ranked_before_all_lf():
    fam = _fam()
    ranked = rank_labels(fam, LABELS)
    order = {lab: i for i, lab in enumerate(ranked.ordered)}
    non_lf = [lab for lab in LABELS if ranked.lf[lab] is not True]
    lf = [lab for lab in LABELS if ranked.lf[lab] is True]
    assert non_lf and lf  # the fixture actually contains both
    assert max(order[lab] for lab in non_lf) < min(order[lab] for lab in lf)


def test_simple_non_lf_does_not_stay_free():
    fam = _fam()
    # (0,1) is a *simple* non-LF label; (0,0) is a simple LF label.
    assert is_locally_finite(fam, (0, 1)) is False
    assert is_locally_finite(fam, (0, 0)) is True
    ranked = rank_labels(fam, LABELS)
    order = {lab: i for i, lab in enumerate(ranked.ordered)}
    # The simple non-LF is eliminated earlier than the simple LF -> it is NOT left free.
    assert order[(0, 1)] < order[(0, 0)]
    # The free-most (last) label is a locally finite one, never a non-LF one.
    assert ranked.lf[ranked.ordered[-1]] is True


def test_complex_lf_eliminated_before_simple_lf():
    fam = _fam()
    ranked = rank_labels(fam, LABELS)
    order = {lab: i for i, lab in enumerate(ranked.ordered)}
    # (2,-2) is a complex LF label, (0,0) is the simplest LF label.
    assert ranked.lf[(2, -2)] is True and ranked.lf[(0, 0)] is True
    assert ranked.complexity[(2, -2)] > ranked.complexity[(0, 0)]
    assert order[(2, -2)] < order[(0, 0)]  # complex eliminated first, simple kept freer


def test_preferred_masters_are_freest():
    fam = _fam()
    ranked = rank_labels(fam, LABELS, preferred_masters=[(0, -1)])
    order = {lab: i for i, lab in enumerate(ranked.ordered)}
    generic_lf = [lab for lab in LABELS if ranked.lf[lab] is True and lab != (0, -1)]
    assert ranked.tiers[(0, -1)] == 3
    assert order[(0, -1)] > max(order[lab] for lab in generic_lf)  # ranked after generic LF


def test_ranking_no_hardcode_custom_names():
    fam = parse_family_text("""
    IBPInput = <|
      "Variables" -> {u, v}, "Parameters" -> {}, "Regulators" -> {},
      "Polynomials" -> <| "A" -> 1 + u + v |>,
      "MonomialExponents" -> <| u -> 0, v -> 0 |>,
      "PolynomialExponents" -> <| "A" -> -3 |>
    |>
    """)
    labels = [(0, 0, 0), (1, 0, 0), (0, 0, 2)]
    ranked = rank_labels(fam, labels)
    assert sorted(ranked.ordered) == sorted(labels)
