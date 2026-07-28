"""Tests for the explicit SurfacePolicy API (External Int2 Method.11, Phase A).

Covers: immutable constructors (limit / chamber), exact-rational chamber values (floats
rejected), strict chamber sign semantics (zero / unknown reject), default equivalence of the
threaded filters with the historical ``epsilon -> 0^-`` limit reading, and the reducer-level
threading (``ReducerConfig.surface_policy`` + diagnostics record of mode and values).
"""

from __future__ import annotations

import dataclasses
from fractions import Fraction

import pytest
import sympy as sp

from parametric_ibp_lf_reducer import (
    ParamExpr,
    SparsePoly,
    SurfacePolicy,
    coordinate_primitive_surface_free,
    parse_family_text,
    regulated_sign,
    vector_field_surface_free,
    zero_label,
)
from parametric_ibp_lf_reducer.reducer import ReducerConfig, _generate_rows

# G0 = 1+x+y ; x boundaries are fine, but y (and the diagonal) diverge.
ASYMMETRIC = """
IBPInput = <|
  "Variables" -> {x, y}, "Parameters" -> {ep}, "Regulators" -> {ep},
  "Polynomials" -> <| "G0" -> 1 + x + y |>,
  "MonomialExponents" -> <| x -> 0, y -> -3 |>,
  "PolynomialExponents" -> <| "G0" -> -2 |>
|>
"""

# x-exponent depends on the regulator: the x-check with P = x is decided by the x -> 0 facet
# score 1 + ep, so chamber points decide it differently. (Method.13: y-exponent 0, not -3, so
# no mixed transverse ray interferes with the decision — the deciding facet stays +e_x.)
EP_FAMILY = """
IBPInput = <|
  "Variables" -> {x, y}, "Parameters" -> {ep}, "Regulators" -> {ep},
  "Polynomials" -> <| "G0" -> 1 + x + y |>,
  "MonomialExponents" -> <| x -> ep, y -> 0 |>,
  "PolynomialExponents" -> <| "G0" -> -2 |>
|>
"""


# ---- constructors / immutability -------------------------------------------
def test_limit_constructor_and_validation():
    pol = SurfacePolicy.limit()
    assert pol.mode == "limit" and pol.direction == "minus" and pol.point == {}
    assert SurfacePolicy.limit("plus").direction == "plus"
    with pytest.raises(ValueError):
        SurfacePolicy.limit("sideways")
    with pytest.raises(ValueError):
        SurfacePolicy(mode="both")
    with pytest.raises(ValueError):  # limit mode takes no chamber point
        SurfacePolicy(mode="limit", point_items=(("ep", Fraction(1)),))


def test_chamber_constructor_exact_values_only():
    pol = SurfacePolicy.chamber({"ep": Fraction(-3, 5)})
    assert pol.mode == "chamber"
    assert pol.point == {"ep": Fraction(-3, 5)}
    # int and string forms are exact and accepted
    assert SurfacePolicy.chamber({"ep": -2}).point == {"ep": Fraction(-2)}
    assert SurfacePolicy.chamber({"ep": "-3/5"}).point == {"ep": Fraction(-3, 5)}
    # floats are rejected: no floating-point sign decisions
    with pytest.raises(TypeError):
        SurfacePolicy.chamber({"ep": -0.6})
    # a chamber must actually pin at least one regulator
    with pytest.raises(ValueError):
        SurfacePolicy.chamber({})


def test_policy_is_immutable():
    pol = SurfacePolicy.chamber({"ep": Fraction(-3, 5)})
    with pytest.raises(dataclasses.FrozenInstanceError):
        pol.mode = "limit"
    # the exposed point is a fresh dict; mutating it does not touch the policy
    d = pol.point
    d["ep"] = Fraction(1)
    assert pol.point == {"ep": Fraction(-3, 5)}


def test_describe_records_mode_and_values():
    assert SurfacePolicy.limit().describe() == {"mode": "limit", "direction": "minus"}
    assert SurfacePolicy.chamber({"ep": Fraction(-3, 5)}).describe() == {
        "mode": "chamber",
        "point": {"ep": "-3/5"},
    }


# ---- sign semantics --------------------------------------------------------
def test_sign_at_limit_matches_regulated_sign():
    ep = sp.Symbol("ep")
    for expr in (ep, -ep, 1 + 2 * ep, sp.Integer(3), sp.Integer(0)):
        for direction in ("minus", "plus"):
            pol = SurfacePolicy.limit(direction)
            assert pol.sign_at(expr, ["ep"]) == regulated_sign(expr, ["ep"], direction)


def test_sign_at_chamber_is_strict_exact_substitution():
    ep = sp.Symbol("ep")
    expr = 1 + 2 * ep
    assert SurfacePolicy.chamber({"ep": Fraction(-1, 4)}).sign_at(expr, ["ep"]) == "pos"
    assert SurfacePolicy.chamber({"ep": Fraction(-3, 5)}).sign_at(expr, ["ep"]) == "neg"
    # exact zero is "zero" (marginal), never coerced to a sign
    assert SurfacePolicy.chamber({"ep": Fraction(-1, 2)}).sign_at(expr, ["ep"]) == "zero"
    # leftover symbols (non-regulator parameter or missing regulator) are "unknown"
    delta = sp.Symbol("delta")
    pol = SurfacePolicy.chamber({"ep": Fraction(-3, 5)})
    assert pol.sign_at(1 + delta, ["ep"]) == "unknown"
    assert pol.sign_at(1 + 2 * sp.Symbol("ep2"), ["ep", "ep2"]) == "unknown"


# ---- filter threading: default equivalence --------------------------------
def test_filters_default_equals_explicit_limit_policy():
    fam = parse_family_text(ASYMMETRIC)
    lab = zero_label(2, 1)
    limit_pol = SurfacePolicy.limit("minus")
    for var_index in range(2):
        for p in ((0, 0), (1, 0), (0, 1), (1, 1)):
            assert coordinate_primitive_surface_free(
                fam, lab, var_index, multiplier_exps=p
            ) == coordinate_primitive_surface_free(
                fam, lab, var_index, multiplier_exps=p, policy=limit_pol
            )
    one = ParamExpr.one(())
    field = [SparsePoly.monomial(2, (1, 0), one), SparsePoly.monomial(2, (0, 1), one)]
    assert vector_field_surface_free(fam, lab, field) == vector_field_surface_free(
        fam, lab, field, policy=limit_pol
    )


def test_chamber_point_decides_coordinate_check():
    fam = parse_family_text(EP_FAMILY)
    lab = zero_label(2, 1)
    xi = fam.variables.index("x")
    # limit minus: exponents 1 + ep (pos) and ep - 1 (neg) -> surface-free
    assert coordinate_primitive_surface_free(fam, lab, xi, multiplier_exps=(1, 0)) is True
    # a chamber point near 0 on the same side agrees
    near = SurfacePolicy.chamber({"ep": Fraction(-1, 4)})
    assert (
        coordinate_primitive_surface_free(fam, lab, xi, multiplier_exps=(1, 0), policy=near) is True
    )
    # a far chamber flips the x -> 0 exponent negative -> rejected
    far = SurfacePolicy.chamber({"ep": Fraction(-3, 2)})
    assert (
        coordinate_primitive_surface_free(fam, lab, xi, multiplier_exps=(1, 0), policy=far) is False
    )
    # on the chamber wall the exponent is exactly zero -> strict reject (marginal)
    wall = SurfacePolicy.chamber({"ep": Fraction(-1)})
    assert (
        coordinate_primitive_surface_free(fam, lab, xi, multiplier_exps=(1, 0), policy=wall)
        is False
    )


# ---- reducer threading + diagnostics ---------------------------------------
def test_generate_rows_records_policy_and_default_is_limit():
    fam = parse_family_text(ASYMMETRIC)
    lab = zero_label(2, 1)
    default_cfg = ReducerConfig(labels=[lab], max_ibp_degree=1)
    rows_default, diag_default = _generate_rows(fam, [lab], default_cfg)
    assert diag_default["surface_policy"] == {"mode": "limit", "direction": "minus"}

    explicit_cfg = ReducerConfig(
        labels=[lab], max_ibp_degree=1, surface_policy=SurfacePolicy.limit("minus")
    )
    rows_explicit, diag_explicit = _generate_rows(fam, [lab], explicit_cfg)
    assert diag_explicit["surface_policy"] == diag_default["surface_policy"]
    assert [r.dedup_key() for r in rows_explicit] == [r.dedup_key() for r in rows_default]

    chamber_cfg = ReducerConfig(
        labels=[lab],
        max_ibp_degree=1,
        surface_policy=SurfacePolicy.chamber({"ep": Fraction(-3, 5)}),
    )
    _, diag_chamber = _generate_rows(fam, [lab], chamber_cfg)
    assert diag_chamber["surface_policy"] == {"mode": "chamber", "point": {"ep": "-3/5"}}
