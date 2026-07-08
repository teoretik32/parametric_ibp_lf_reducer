"""Tests for coefficient reconstruction (Pass 2G): CRT + rational recon + validated interp."""

from __future__ import annotations

from fractions import Fraction

import sympy as sp

import pytest

from parametric_ibp_lf_reducer import (
    InterpolationFailed,
    interpolate_multivariate,
    interpolate_univariate,
    rational_reconstruction,
    reconstruct_coefficients,
    reconstruct_rational,
)
from parametric_ibp_lf_reducer.modular_normal_form import NormalFormResult
from parametric_ibp_lf_reducer.reconstruction import collect_value_table

PRIMES = [2_147_483_647, 2_147_483_629, 2_147_483_587]


def _to_res(value: Fraction, prime: int) -> int:
    return value.numerator % prime * pow(value.denominator % prime, -1, prime) % prime


def test_rational_reconstruction_roundtrip():
    m = 2_147_483_647  # bound sqrt(m/2) ~ 32767 comfortably covers the fractions below
    for v in (Fraction(-3, 7), Fraction(5, 2), Fraction(0), Fraction(11), Fraction(-1, 999)):
        a = _to_res(v, m)
        assert rational_reconstruction(a, m) == v


def test_rational_reconstruction_none_when_modulus_too_small():
    # -1/999 cannot be recovered when sqrt(m/2) < 999.
    m = 1_000_003
    assert rational_reconstruction(_to_res(Fraction(-1, 999), m), m) is None


def test_reconstruct_rational_needs_multiple_primes():
    v = Fraction(-34, 135)
    residues = {p: _to_res(v, p) for p in PRIMES}
    assert reconstruct_rational(residues) == v


def test_interpolate_univariate_recovers_known_function():
    ep = sp.Symbol("ep")
    f = (2 * ep - 1) / ep**2
    values = {
        Fraction(k): Fraction(int(f.subs(ep, k).p), int(f.subs(ep, k).q)) for k in range(2, 10)
    }
    got = interpolate_univariate(values, "ep")
    assert sp.simplify(got - f) == 0


def test_interpolate_from_one_point_fails():
    with pytest.raises(InterpolationFailed):
        interpolate_univariate({Fraction(2): Fraction(1)}, "ep")


def _synthetic_records(functions, samples, primes, param="ep"):
    """Build Reduced NormalFormResults for label->C(ep) across (prime, sample)."""
    ep = sp.Symbol(param)
    records = []
    for s in samples:
        for p in primes:
            terms = {}
            for label, f in functions.items():
                val = sp.Rational(f.subs(ep, s))
                res = _to_res(Fraction(int(val.p), int(val.q)), p)
                if res:
                    terms[label] = res
            records.append(
                NormalFormResult(
                    status="Reduced", target_label=(0,), prime=p, sample={param: Fraction(s)},
                    formal_success=True, terms=terms,
                )
            )
    return records


def test_reconstruct_coefficients_end_to_end_univariate():
    ep = sp.Symbol("ep")
    functions = {(1, -2): (2 * ep - 1) / ep**2, (0, -1): sp.Rational(3) * (ep + 1) / ep}
    records = _synthetic_records(functions, samples=range(2, 12), primes=PRIMES)
    coeffs = reconstruct_coefficients(records, ["ep"])
    for label, f in functions.items():
        assert sp.simplify(coeffs[label] - f) == 0


def test_bad_specialization_records_are_skipped_not_dropped():
    ep = sp.Symbol("ep")
    functions = {(2, 0): (ep - 1) / (ep + 2)}
    records = _synthetic_records(functions, samples=range(2, 12), primes=PRIMES)
    # inject bad-specialization records that must be ignored (not corrupt the result)
    records.append(
        NormalFormResult(status="BadSpecialization", target_label=(0,), prime=PRIMES[0],
                         sample={"ep": Fraction(3)}, formal_success=False)
    )
    coeffs = reconstruct_coefficients(records, ["ep"])
    assert sp.simplify(coeffs[(2, 0)] - functions[(2, 0)]) == 0


def _synthetic_records_multi(functions, points, primes, params=("ep", "r")):
    """Build Reduced NormalFormResults for label->C(ep,r) across (prime, point)."""
    syms = [sp.Symbol(p) for p in params]
    records = []
    for pt in points:
        subs = {s: sp.Rational(v) for s, v in zip(syms, pt)}
        for p in primes:
            terms = {}
            for label, f in functions.items():
                val = sp.Rational(f.subs(subs))
                res = _to_res(Fraction(int(val.p), int(val.q)), p)
                if res:
                    terms[label] = res
            records.append(
                NormalFormResult(
                    status="Reduced", target_label=(0,), prime=p,
                    sample={params[i]: Fraction(pt[i]) for i in range(len(params))},
                    formal_success=True, terms=terms,
                )
            )
    return records


def _grid(ep_vals, r_vals):
    return [(a, b) for a in ep_vals for b in r_vals]


def test_interpolate_multivariate_recovers_known_function():
    ep, r = sp.symbols("ep r")
    f = (ep + 2 * r) / (1 + ep + r)
    pts = _grid([1, 2, 3, 4], [1, 2, 3])
    values = {
        (Fraction(a), Fraction(b)): Fraction(int(f.subs({ep: a, r: b}).p), int(f.subs({ep: a, r: b}).q))
        for a, b in pts
    }
    got = interpolate_multivariate(values, ["ep", "r"])
    assert sp.simplify(got - f) == 0


def test_reconstruct_coefficients_multivariate_end_to_end():
    ep, r = sp.symbols("ep r")
    functions = {(1, 0): (ep + 2 * r) / (1 + ep + r), (0, -1): (2 * ep - r) / (ep + r + 3)}
    records = _synthetic_records_multi(functions, _grid([1, 2, 3, 4], [1, 2, 3, 4]), PRIMES)
    coeffs = reconstruct_coefficients(records, ["ep", "r"])
    for label, f in functions.items():
        assert sp.simplify(coeffs[label] - f) == 0


def test_multivariate_insufficient_points_refused():
    ep, r = sp.symbols("ep r")
    # only 3 distinct points -> cannot pin down/validate a 2-var rational function -> refuse
    records = _synthetic_records_multi({(0,): (ep + r) / (ep + 1)}, _grid([1, 2, 3], [1]), PRIMES)
    with pytest.raises(InterpolationFailed):
        reconstruct_coefficients(records, ["ep", "r"])


def test_union_support_zero_fill_for_special_zero():
    # A coefficient that is 0 at some samples (special zero) must still reconstruct as that fn.
    ep = sp.Symbol("ep")
    functions = {(5, 0): ep - 3}  # zero at ep=3 -> term absent there, union support fills 0
    records = _synthetic_records(functions, samples=range(1, 10), primes=PRIMES)
    coeffs = reconstruct_coefficients(records, ["ep"])
    assert sp.simplify(coeffs[(5, 0)] - (ep - 3)) == 0


# --- Pass 2H.1: multivariate reconstruction hardening ---------------------------------------


def _mv_values(expr, pts, params=("ep", "r")):
    """Exact-rational value table ``{(Fraction, ...): Fraction}`` for ``expr`` over ``pts``."""
    syms = [sp.Symbol(p) for p in params]
    values = {}
    for pt in pts:
        subs = {s: sp.Rational(v) for s, v in zip(syms, pt)}
        val = sp.Rational(expr.subs(subs))
        values[tuple(Fraction(v) for v in pt)] = Fraction(int(val.p), int(val.q))
    return values


def test_multivariate_recovers_bivariate_polynomial():
    ep, r = sp.symbols("ep r")
    f = 2 * ep + 3 * r - 5
    values = _mv_values(f, _grid([1, 2, 3, 4], [1, 2, 3]))
    got = interpolate_multivariate(values, ["ep", "r"], max_deg=3)
    assert sp.simplify(got - f) == 0


def test_multivariate_recovers_bivariate_rational_squared_denominator():
    ep, r = sp.symbols("ep r")
    f = (2 * ep * r - ep - 4 * r) / (ep**2 * (r + 1))
    # >= 4 distinct values per variable so the total-degree-3 denominator (needs r^3) is resolved
    values = _mv_values(f, _grid([1, 2, 3, 4, 5], [1, 2, 3, 4]))
    got = interpolate_multivariate(values, ["ep", "r"], max_deg=3)
    assert sp.simplify(got - f) == 0


def test_multivariate_denominator_without_constant_term():
    ep, r = sp.symbols("ep r")
    f = (ep + r) / (ep * r)  # denominator has no constant term
    values = _mv_values(f, _grid([1, 2, 3, 4], [1, 2, 3]))
    got = interpolate_multivariate(values, ["ep", "r"], max_deg=3)
    assert sp.simplify(got - f) == 0


def test_multivariate_rational_numeric_prefactor():
    ep, r = sp.symbols("ep r")
    f = (3 * ep + r / 2) / (ep - r)  # rational numeric coefficient; must stay exact (no float)
    values = _mv_values(f, _grid([2, 4, 6, 8], [1, 3, 5]))  # ep != r everywhere -> denom != 0
    got = interpolate_multivariate(values, ["ep", "r"], max_deg=2)
    assert sp.simplify(got - f) == 0
    assert not got.atoms(sp.Float)  # exact rationals only, never floats


def test_multivariate_corrupted_holdout_is_refused():
    ep, r = sp.symbols("ep r")
    f = (ep + 2 * r) / (1 + ep + r)
    values = _mv_values(f, _grid([1, 2, 3, 4], [1, 2, 3]))
    worst = max(values)  # largest tuple -> lands in the sorted holdout slice
    values[worst] = values[worst] + Fraction(1)  # corrupt one validation point
    with pytest.raises(InterpolationFailed):
        interpolate_multivariate(values, ["ep", "r"], max_deg=2)


def test_multivariate_insufficient_points_direct_refused():
    values = {
        (Fraction(1), Fraction(1)): Fraction(1),
        (Fraction(2), Fraction(1)): Fraction(2),
        (Fraction(3), Fraction(1)): Fraction(3),
    }  # only 3 distinct points -> below min_validation + 2
    with pytest.raises(InterpolationFailed):
        interpolate_multivariate(values, ["ep", "r"])


def test_multivariate_union_support_special_zero():
    ep, r = sp.symbols("ep r")
    functions = {(1, 0): ep - r, (0, -1): ep + r + 1}  # (1,0) vanishes on the diagonal ep == r
    records = _synthetic_records_multi(functions, _grid([1, 2, 3, 4], [1, 2, 3, 4]), PRIMES)
    coeffs = reconstruct_coefficients(records, ["ep", "r"])
    assert sp.simplify(coeffs[(1, 0)] - (ep - r)) == 0  # special-zero term reconstructed
    assert sp.simplify(coeffs[(0, -1)] - (ep + r + 1)) == 0


def test_multivariate_bad_and_nonpivot_records_skipped_not_patched():
    ep, r = sp.symbols("ep r")
    functions = {(1, 0): (ep + 2 * r) / (1 + ep + r)}
    records = _synthetic_records_multi(functions, _grid([1, 2, 3, 4], [1, 2, 3, 4]), PRIMES)
    # inject records that must be ignored: a non-pivot target and a bad specialization
    records.append(
        NormalFormResult(status="TargetNotReducible", target_label=(0,), prime=PRIMES[0],
                         sample={"ep": Fraction(2), "r": Fraction(2)}, formal_success=False)
    )
    records.append(
        NormalFormResult(status="BadSpecialization", target_label=(0,), prime=PRIMES[0],
                         sample={"ep": Fraction(3), "r": Fraction(3)}, formal_success=False)
    )
    _, _, _, n_skipped = collect_value_table(records)
    assert n_skipped == 2  # both counted as skipped, never consumed
    coeffs = reconstruct_coefficients(records, ["ep", "r"])
    assert sp.simplify(coeffs[(1, 0)] - functions[(1, 0)]) == 0  # unaffected by injected records
