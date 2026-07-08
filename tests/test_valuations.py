"""Tests for rays, base_score, and the strict local-finiteness decision at epsilon=0."""

from __future__ import annotations

from conftest import load_example
from parametric_ibp_lf_reducer import (
    base_score,
    compute_candidate_rays,
    is_locally_finite,
    make_label,
    parse_family_text,
    valuation_poly,
    zero_label,
)

# --- small synthetic families (no hardcoded validation polynomials in the core) ----------
ONE_VAR_CONVERGENT = """
IBPInput = <|
  "Variables" -> {x}, "Parameters" -> {ep}, "Regulators" -> {ep},
  "Polynomials" -> <| "G0" -> 1 + x |>,
  "MonomialExponents" -> <| x -> 0 |>,
  "PolynomialExponents" -> <| "G0" -> -2 |>
|>
"""

ONE_VAR_LOG_DIVERGENT = """
IBPInput = <|
  "Variables" -> {x}, "Parameters" -> {ep}, "Regulators" -> {ep},
  "Polynomials" -> <| "G0" -> 1 + x |>,
  "MonomialExponents" -> <| x -> 0 |>,
  "PolynomialExponents" -> <| "G0" -> -1 |>
|>
"""

# base_score == 0 at eps=0 but the eps-coefficient is negative (would regulate for eps<0).
ONE_VAR_EPS_REGULATED = """
IBPInput = <|
  "Variables" -> {x}, "Parameters" -> {ep}, "Regulators" -> {ep},
  "Polynomials" -> <| "G0" -> 1 + x |>,
  "MonomialExponents" -> <| x -> -1 - ep |>,
  "PolynomialExponents" -> <| "G0" -> 0 |>
|>
"""

# exponent depends on a non-regulator parameter r with no assumption -> undecidable.
ONE_VAR_SYMBOLIC_EXPONENT = """
IBPInput = <|
  "Variables" -> {x}, "Parameters" -> {ep, r}, "Regulators" -> {ep},
  "Polynomials" -> <| "G0" -> 1 + x |>,
  "MonomialExponents" -> <| x -> r |>,
  "PolynomialExponents" -> <| "G0" -> -2 |>
|>
"""

# boundary-convergent but a denominator can vanish in the interior (coeff -1, no positivity) ->
# the boundary-ray test is insufficient -> Unknown.
BULK_SINGULAR = """
IBPInput = <|
  "Variables" -> {x, y}, "Parameters" -> {ep}, "Regulators" -> {ep},
  "Polynomials" -> <| "G0" -> 1 + x + y, "G1" -> 2 - x |>,
  "MonomialExponents" -> <| x -> 0, y -> 0 |>,
  "PolynomialExponents" -> <| "G0" -> -2, "G1" -> -2 |>
|>
"""


def test_candidate_rays_include_coordinate_and_mixed():
    fam = parse_family_text("""
    IBPInput = <|
      "Variables" -> {x, y}, "Parameters" -> {ep}, "Regulators" -> {ep},
      "Polynomials" -> <| "G0" -> 1 + x + y |>,
      "MonomialExponents" -> <| x -> 0, y -> 0 |>,
      "PolynomialExponents" -> <| "G0" -> -2 |>
    |>
    """)
    dirs = {ray.direction for ray in compute_candidate_rays(fam)}
    assert (1, 0) in dirs and (0, 1) in dirs  # coordinate zero rays
    assert (-1, 0) in dirs and (0, -1) in dirs  # infinity rays
    assert (1, 1) in dirs or (-1, -1) in dirs  # a mixed/diagonal toric ray
    # valuation_poly delegates to SparsePoly.valuation (min convention)
    g0 = fam.polynomial("G0")
    assert valuation_poly(g0, (1, 1)) == 0
    assert valuation_poly(g0, (-1, -1)) == -1


def test_base_score_positive_is_locally_finite():
    fam = parse_family_text(ONE_VAR_CONVERGENT)
    lab = zero_label(1, 1)
    assert base_score(fam, lab, (1,)) == 1  # x->0
    assert base_score(fam, lab, (-1,)) == 1  # x->inf
    assert is_locally_finite(fam, lab) is True


def test_base_score_zero_is_not_locally_finite():
    fam = parse_family_text(ONE_VAR_LOG_DIVERGENT)
    lab = zero_label(1, 1)
    assert base_score(fam, lab, (-1,)) == 0  # marginal at infinity
    assert is_locally_finite(fam, lab) is False  # strict rule: 0 is not LF


def test_base_score_zero_with_negative_eps_coeff_is_not_locally_finite():
    fam = parse_family_text(ONE_VAR_EPS_REGULATED)
    lab = zero_label(1, 1)
    # At eps=0 the score is exactly 0; the eps-term (coefficient < 0) would regulate for eps<0,
    # but that must NOT count as local finiteness.
    assert base_score(fam, lab, (1,)) == 0
    assert is_locally_finite(fam, lab) is False


def test_symbolic_exponent_gives_unknown():
    fam = parse_family_text(ONE_VAR_SYMBOLIC_EXPONENT)
    assert is_locally_finite(fam, zero_label(1, 1)) == "Unknown"


def test_bulk_singularity_without_assumptions_gives_unknown():
    fam = parse_family_text(BULK_SINGULAR)
    assert is_locally_finite(fam, zero_label(2, 2)) == "Unknown"


def test_d4_master_is_locally_finite():
    fam = parse_family_text(load_example("d4_explicit_family.wl.txt"))
    # Master M1 relative factor x2*x3/(G0^2*G1) -> n=(0,1,1,0), m=(-2,-1,0).
    m1 = make_label([0, 1, 1, 0], [-2, -1, 0])
    assert is_locally_finite(fam, m1) is True


def test_d4_base_target_is_not_locally_finite():
    # The base integrand (target) has base_score == 0 along x4 -> infinity at eps=0
    # (score = -ep), i.e. it is only epsilon-regulated -> NOT locally finite.
    fam = parse_family_text(load_example("d4_explicit_family.wl.txt"))
    assert is_locally_finite(fam, zero_label(4, 3)) is False
