"""Tests for tangent (logarithmic) vector fields via the syzygy ansatz (Pass 2B)."""

from __future__ import annotations

import sympy as sp

from parametric_ibp_lf_reducer import (
    ParamExpr,
    SparsePoly,
    generate_tangent_fields,
    parse_family_text,
    verify_tangent,
)

# G = 1 + x + y  (the spec 11.1 sanity family)
G_XY = """
IBPInput = <|
  "Variables" -> {x, y}, "Parameters" -> {}, "Regulators" -> {},
  "Polynomials" -> <| "G" -> 1 + x + y |>,
  "MonomialExponents" -> <| x -> 0, y -> 0 |>,
  "PolynomialExponents" -> <| "G" -> -1 |>
|>
"""


def _mono(exps, coeff_int=1, params=()):
    return SparsePoly.monomial(len(exps), exps, ParamExpr.from_int(coeff_int, params))


def test_verify_tangent_field_11_1():
    fam = parse_family_text(G_XY)
    # Q = (x*y, -x*y) is tangent to G = 1+x+y with H = 0.
    qx = _mono((1, 1), 1)
    qy = _mono((1, 1), -1)
    ok, hs = verify_tangent(fam, [qx, qy])
    assert ok is True
    assert len(hs) == 1 and hs[0].is_zero  # Q.grad G = 0


def test_verify_rejects_non_tangent_field():
    fam = parse_family_text(G_XY)
    # Q = (x, 0): Q.grad G = x, not divisible by 1+x+y.
    ok, hs = verify_tangent(fam, [_mono((1, 0), 1), SparsePoly.zero(2, ())])
    assert ok is False and hs is None


def test_generate_fields_all_verified_zero_dropped_and_deduped():
    fam = parse_family_text(G_XY)
    fields = generate_tangent_fields(fam, [(1, 0), (2, 0)])
    assert len(fields) >= 1
    # Every returned field is genuinely tangent to all G_l.
    assert all(f.is_tangent(fam) for f in fields)
    # The zero field is never returned.
    assert all(not f.is_zero for f in fields)
    # No two returned fields are parameter-scalar proportional (dedup worked).
    from parametric_ibp_lf_reducer.tangent_fields import _proportional

    xsyms = [sp.Symbol("x"), sp.Symbol("y")]
    sym = [tuple(q.to_sympy(fam.variables) for q in f.components) for f in fields]
    for i in range(len(sym)):
        for j in range(i + 1, len(sym)):
            assert not _proportional(sym[i], sym[j], xsyms)


def test_generated_field_reproduces_xy_field_in_span():
    # The (x*y, -x*y) field must be expressible from the degree-2 tangent basis.
    fam = parse_family_text(G_XY)
    fields = generate_tangent_fields(fam, [(2, 0)])
    xsyms = [sp.Symbol("x"), sp.Symbol("y")]
    target = [sp.sympify("x*y"), sp.sympify("-x*y")]
    # Solve target = sum_k c_k * field_k over the (constant) scalars c_k.
    cs = sp.symbols(f"c0:{len(fields)}")
    combo = [sp.Integer(0), sp.Integer(0)]
    for c, f in zip(cs, fields):
        comp = [q.to_sympy(fam.variables) for q in f.components]
        combo = [combo[0] + c * comp[0], combo[1] + c * comp[1]]
    eqs = []
    for comp_expr, tgt in zip(combo, target):
        for _, coeff in sp.Poly(sp.expand(comp_expr - tgt), *xsyms).terms():
            eqs.append(coeff)
    sol = sp.linsolve(eqs, cs)
    assert sol != sp.EmptySet  # a solution exists -> (xy,-xy) is in the span


def test_parametric_tangent_field_has_param_coefficients():
    fam = parse_family_text("""
    IBPInput = <|
      "Variables" -> {x, y}, "Parameters" -> {s}, "Regulators" -> {},
      "Polynomials" -> <| "G" -> 1 + s*x + y |>,
      "MonomialExponents" -> <| x -> 0, y -> 0 |>,
      "PolynomialExponents" -> <| "G" -> -1 |>
    |>
    """)
    fields = generate_tangent_fields(fam, [(1, 0)])
    assert fields and all(f.is_tangent(fam) for f in fields)
    s = sp.Symbol("s")
    has_param = any(
        s in q.to_sympy(fam.variables).free_symbols for f in fields for q in f.components
    )
    assert has_param  # parameter-dependent (rational) coefficients are supported


def test_custom_names_no_hardcode():
    fam = parse_family_text("""
    IBPInput = <|
      "Variables" -> {u, v}, "Parameters" -> {}, "Regulators" -> {},
      "Polynomials" -> <| "A" -> 1 + u + v |>,
      "MonomialExponents" -> <| u -> 0, v -> 0 |>,
      "PolynomialExponents" -> <| "A" -> -1 |>
    |>
    """)
    ok, _ = verify_tangent(fam, [_mono((1, 1), 1), _mono((1, 1), -1)])  # (u*v, -u*v)
    assert ok is True
    fields = generate_tangent_fields(fam, [(1, 0)])
    assert fields and all(f.is_tangent(fam) for f in fields)
