"""Tests for the finite-numerator LF basis search (single-integrand semantics).

Hand-verified synthetic model ``G0 = x + x^2 + x^3`` with sector ``1/G0``
(label ``(0, -1)``, measure weight ``x^0``):

* ray ``+1`` (x -> 0): ``score = 1 - val_{+1}(G0) = 1 - 1 = 0`` -> divergent
  (``G0 ~ x`` at the origin, ``1/G0`` is a log divergence);
* ray ``-1`` (x -> oo): ``score = -1 + 3 = 2`` -> fine.

Numerator monomials: ``x`` clears both rays (``x/G0`` is LF everywhere), while
``1`` fails at the origin and ``x^2`` fails at infinity (``score = 0``). So the
degree-2 ansatz ``{1, x, x^2}`` has exactly one admissible support ``{x}`` and —
per Lemma 1 (graded lowest layer) — NO linear combination containing ``1`` or
``x^2`` can be locally finite as a complete integrand: separately divergent
pieces never combine into an accepted candidate. This is asserted below (the
corrected direction: no cancellation after integration is ever used).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from parametric_ibp_lf_reducer import (
    Row,
    lf_reduction_feasible_mod_p,
    parse_family_text,
)
from parametric_ibp_lf_reducer.coefficients import ParamExpr
from parametric_ibp_lf_reducer.finite_numerator import (
    STATUS_ALREADY_LF,
    STATUS_FOUND,
    STATUS_IMPOSSIBLE,
    decorated_label,
    finite_numerator_search,
    full_integrand_lf,
    impossible_any_degree,
    leading_cancellation_kernel,
    numerator_ansatz,
    search_sector,
)
from parametric_ibp_lf_reducer.sparse_poly import SparsePoly
from parametric_ibp_lf_reducer.valuations import is_locally_finite

CURABLE = """
IBPInput = <|
  "Variables" -> {x}, "Parameters" -> {ep}, "Regulators" -> {ep},
  "Polynomials" -> <| "G0" -> x + x^2 + x^3 |>,
  "MonomialExponents" -> <| x -> 0 |>,
  "PolynomialExponents" -> <| "G0" -> 0 |>
|>
"""

INCURABLE = """
IBPInput = <|
  "Variables" -> {x}, "Parameters" -> {ep}, "Regulators" -> {ep},
  "Polynomials" -> <| "G0" -> 1 + x |>,
  "MonomialExponents" -> <| x -> 0 |>,
  "PolynomialExponents" -> <| "G0" -> 0 |>
|>
"""

SECTOR = (0, -1)  # 1/G0


@pytest.fixture(scope="module")
def curable():
    return parse_family_text(CURABLE)


@pytest.fixture(scope="module")
def incurable():
    return parse_family_text(INCURABLE)


@pytest.fixture(scope="module")
def int2_family():
    path = (
        Path(__file__).resolve().parents[1]
        / "examples"
        / "external_int2_dimensionless_input.wl.txt"
    )
    return parse_family_text(path.read_text(encoding="utf-8"))


INT2_SECTORS = {
    "1/(x2*G0*G1)": (-1, 0, 0, -1, -1, 0, 0),
    "1/(x2*G1*G3)": (-1, 0, 0, 0, -1, 0, -1),
    "1/(G0*G3)": (0, 0, 0, -1, 0, 0, -1),
    "1/G1": (0, 0, 0, 0, -1, 0, 0),
    "1/G2": (0, 0, 0, 0, 0, -1, 0),
    "x7/(G0*G3)": (0, 0, 1, -1, 0, 0, -1),
    "1/(G1*G3)": (0, 0, 0, 0, -1, 0, -1),
}


def _mono(family, alpha):
    return SparsePoly.monomial(
        family.nvars, alpha, ParamExpr.from_int(1, tuple(family.parameters))
    )


# --------------------------------------------------------------------------- #
# Basic structure


def test_ansatz_and_decorated_label():
    assert numerator_ansatz(1, 2) == [(0,), (1,), (2,)]
    assert decorated_label((0, -1), (2,)) == (2, -1)
    assert decorated_label((-1, 0, 0, 0, -1, 0, 0), (1, 0, 2)) == (0, 0, 2, 0, -1, 0, 0)


# --------------------------------------------------------------------------- #
# Synthetic curable sector: candidate found, exported as ONE integrand


def test_synthetic_clearing_monomial_accepted(curable):
    rep = search_sector(curable, SECTOR, degree=2)
    assert rep.bare_lf is False
    assert rep.status == STATUS_FOUND
    assert rep.impossible_any_degree is False
    assert ((1,), "1") not in [(a, "x") for a in rep.clearing]  # sanity of shape
    assert rep.clearing == [(1,)]
    mono_cands = [c for c in rep.candidates if len(c.numerator.support()) == 1]
    assert len(mono_cands) == 1
    cand = mono_cands[0]
    assert cand.lf_verdict is True
    assert cand.degree == 1
    # ONE numerator-decorated integrand, HyperInt-ready
    assert cand.to_wolfram_text() == "(x)*(1/G0)"
    # exact defining expansion into monomial labels
    assert cand.defining_expansion() == [((1, -1), "1")]
    # every expansion label is individually LF (never a divergent member)
    for lab, _c in cand.defining_expansion():
        assert is_locally_finite(curable, lab) is True


def test_lemma1_kernel_matches_clearing(curable):
    alphas = numerator_ansatz(1, 2)
    kernel, n_conditions = leading_cancellation_kernel(curable, SECTOR, alphas)
    # honest solver: conditions exist and kernel == span{e_x} exactly (Lemma 1)
    assert n_conditions > 0
    assert len(kernel) == 1
    vec = kernel[0]
    assert vec[alphas.index((1,))] != 0
    assert vec[alphas.index((0,))] == 0
    assert vec[alphas.index((2,))] == 0
    rep = search_sector(curable, SECTOR, degree=2)
    assert rep.lemma_consistent is True
    assert rep.kernel_dim == 1


def test_separately_divergent_pieces_rejected(curable):
    # pieces: 1/G0 divergent at 0, x^2/G0 divergent at oo, x/G0 finite
    assert is_locally_finite(curable, (0, -1)) is False
    assert is_locally_finite(curable, (2, -1)) is False
    assert is_locally_finite(curable, (1, -1)) is True
    # No polynomial mixing a divergent piece is ever accepted as one integrand:
    for terms in [{(0,): 1, (1,): 1}, {(1,): 1, (2,): 1}, {(0,): 1, (2,): -1}]:
        num = SparsePoly(
            1,
            tuple(curable.parameters),
            {a: ParamExpr.from_int(c, tuple(curable.parameters)) for a, c in terms.items()},
        )
        assert full_integrand_lf(curable, SECTOR, num) is False


# --------------------------------------------------------------------------- #
# Lemma 2: incurable sector (failing ray componentwise <= 0)


def test_lemma2_impossible_flag(incurable):
    assert is_locally_finite(incurable, SECTOR) is False  # 1/(1+x) at x -> oo
    assert impossible_any_degree(incurable, SECTOR) is True
    rep = search_sector(incurable, SECTOR, degree=2)
    assert rep.status == STATUS_IMPOSSIBLE
    assert rep.clearing == []
    assert rep.kernel_dim == 0
    assert rep.candidates == []
    assert rep.lemma_consistent is True


# --------------------------------------------------------------------------- #
# External Int2: honest verdicts for the six normal-form sectors


def test_int2_sector_search(int2_family):
    payload = finite_numerator_search(int2_family, INT2_SECTORS, degrees=(0, 1, 2))
    by_name = {}
    for rep in payload["reports"]:
        by_name.setdefault(rep["sector_name"], []).append(rep)
    for name in ("1/(x2*G0*G1)", "1/(x2*G1*G3)", "1/(G0*G3)", "x7/(G0*G3)"):
        assert all(r["status"] == STATUS_ALREADY_LF for r in by_name[name])
    for name in ("1/G1", "1/G2", "1/(G1*G3)"):
        for r in by_name[name]:
            assert r["status"] == STATUS_IMPOSSIBLE
            assert r["numerator_cure_impossible_any_degree"] is True
            assert r["candidates"] == []
            # all failing rays componentwise <= 0: x -> oo divergences
            assert all(
                all(c <= 0 for c in f["direction"]) for f in r["failing_rays"]
            )
    assert payload["new_lf_masters_found"] is False
    assert payload["lemma_consistent_everywhere"] is True


def test_int2_reducer_facing_verdicts_unchanged(int2_family):
    # The search is read-only: the gate verdicts it builds on stay as certified.
    expected = {
        "1/(x2*G0*G1)": True,
        "1/(x2*G1*G3)": True,
        "1/(G0*G3)": True,
        "1/G1": False,
        "1/G2": False,
        "x7/(G0*G3)": True,
        "1/(G1*G3)": False,
    }
    for name, sector in INT2_SECTORS.items():
        assert is_locally_finite(int2_family, sector) is expected[name]


def test_module_never_imports_reducer_core():
    import parametric_ibp_lf_reducer.finite_numerator as fn

    src = Path(fn.__file__).read_text(encoding="utf-8")
    for banned in ("from .reducer", "from .certificate", "import reducer"):
        assert banned not in src


# --------------------------------------------------------------------------- #
# Offset convention: total exponent = base exponent + label shift

BASED = """
IBPInput = <|
  "Variables" -> {x}, "Parameters" -> {ep}, "Regulators" -> {ep},
  "Polynomials" -> <| "G0" -> x + x^2 + x^3 |>,
  "MonomialExponents" -> <| x -> 1 |>,
  "PolynomialExponents" -> <| "G0" -> -1 |>
|>
"""


def test_offset_convention_base_plus_shift(curable):
    based = parse_family_text(BASED)
    # Label (0, 0) on the based family carries the same TOTAL exponents as
    # label (1, -1) on the zero-base family: x^1 / G0 (LF everywhere).
    assert is_locally_finite(based, (0, 0)) is True
    assert is_locally_finite(based, (0, 0)) == is_locally_finite(curable, (1, -1))
    # Shifting back down reproduces the divergent bare sector 1/G0.
    assert is_locally_finite(based, (-1, 0)) is False
    assert is_locally_finite(based, (-1, 0)) == is_locally_finite(curable, SECTOR)
    rep_based = search_sector(based, (-1, 0), 2)
    rep_zero = search_sector(curable, SECTOR, 2)
    assert rep_based.status == rep_zero.status == STATUS_FOUND
    assert rep_based.clearing == rep_zero.clearing


# --------------------------------------------------------------------------- #
# Defining rows bridge into the Method.1 modular span test


def test_defining_rows_feed_modular_feasibility(curable):
    rep = search_sector(curable, SECTOR, 2)
    cand = next(c for c in rep.candidates if c.degree == 1)
    exp = cand.defining_expansion()
    assert exp == [((1, -1), "1")]
    # Lemma 1: every defining-expansion label is individually LF ...
    for lab, _ in exp:
        assert is_locally_finite(curable, lab) is True
    # ... so the span test may mark it allowed: the target reduces through the
    # master's expansion label although the bare sector itself is forbidden.
    master = exp[0][0]
    row = Row(
        "embed",
        {},
        {
            SECTOR: ParamExpr.from_int(1, ("ep",)),
            master: ParamExpr.from_int(-1, ("ep",)),
        },
    )
    res = lf_reduction_feasible_mod_p(
        [row],
        [SECTOR, master],
        SECTOR,
        {SECTOR: False, master: True},
        {"ep": 5},
        30011,
    )
    assert res.status == "Feasible"
    assert res.residual_support == ()
