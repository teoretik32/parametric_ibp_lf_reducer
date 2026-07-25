"""Tests for Method.11b: special-zero vs unstable support classification (Phase A).

A ``special_zero`` sample (all primes at generic rank agree on the same missing label set)
must be retained and zero-filled; an ``unstable`` sample (rank drop or per-prime support
disagreement) must be rejected from the value table entirely.
"""

from __future__ import annotations

import sys
from fractions import Fraction as F
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

import sympy as sp  # noqa: E402

from parametric_ibp_lf_reducer import (  # noqa: E402
    SUPPORT_SPECIAL_ZERO,
    SUPPORT_STABLE,
    SUPPORT_UNSTABLE,
    InterpolationFailed,
    classify_sample_supports,
    collect_value_table,
    interpolate_multivariate,
    reconstruct_coefficients,
)
from parametric_ibp_lf_reducer.modular_normal_form import STATUS_REDUCED  # noqa: E402
from parametric_ibp_lf_reducer.records import NormalFormRecord  # noqa: E402

PRIMES = (2147483647, 2147483629, 2147483587)
TARGET = (1, 1, 1, 1, 1, 1, -2)
LAB_A = (0, 0, 1, -1, 0, 0, -1)
LAB_B = (0, 1, 0, 0, 0, 0, -1)
GENERIC_RANK = 100


# Exact coefficient functions used to synthesize records.
def _coeff_a(ep: F, r: F) -> F:
    return (ep - 6) / (r + 1)  # vanishes exactly at ep=6


def _coeff_b(ep: F, r: F) -> F:
    return 1 / (ep + r)  # never vanishes on our positive samples


def _r_of(k: int) -> F:
    return F((3 * k) % 7 + 1, 3)  # scattered: keeps the samples off any single line


GENERIC_POINTS = [(F(k, 5), _r_of(k)) for k in range(1, 12)]
SPECIAL_POINT = (F(6), F(57, 11))  # sorts last (max ep) -> lands in the holdout set


def _mod(value: F, p: int) -> int:
    return value.numerator * pow(value.denominator, -1, p) % p


def _record(pt, prime, *, drop=(), rank=GENERIC_RANK, extra=()):  # noqa: ANN001
    ep, r = pt
    coeffs = {}
    for lab, fn in ((LAB_A, _coeff_a), (LAB_B, _coeff_b)):
        val = fn(ep, r)
        if lab in drop or val == 0:
            continue  # a missing label == exact zero / dropped support
        coeffs[lab] = _mod(val, prime)
    for lab in extra:
        coeffs[lab] = 12345 % prime
    return NormalFormRecord(
        prime=prime,
        sample={"ep": ep, "r": r},
        target_label=TARGET,
        status=STATUS_REDUCED,
        formal_success=True,
        coeffs=coeffs,
        support=tuple(sorted(coeffs)),
        rank=rank,
    )


def _base_records(points):
    return [_record(pt, p) for pt in points for p in PRIMES]


def _key_of(samples: dict, ep: F) -> object:
    (key,) = [k for k, s in samples.items() if F(s["ep"]) == ep]
    return key


def test_special_zero_sample_classified_and_zero_filled():
    records = _base_records(GENERIC_POINTS + [SPECIAL_POINT])
    cls = classify_sample_supports(records)
    summary = cls["summary"]
    assert summary["generic_rank"] == GENERIC_RANK
    assert summary["reference_support_size"] == 2
    assert summary["n_samples"] == 12
    assert summary["n_stable"] == 11
    assert summary["n_special_zero"] == 1
    assert summary["n_unstable"] == 0
    assert summary["special_zero_samples"] == ["ep=6,r=57/11"]
    entry = cls["by_sample"]["ep=6,r=57/11"]
    assert entry["classification"] == SUPPORT_SPECIAL_ZERO
    assert entry["missing_labels"] == [list(LAB_A)]
    assert entry["ranks"] == [GENERIC_RANK]

    labels, table, samples, n_skipped = collect_value_table(records)
    assert n_skipped == 0  # special-zero sample is retained
    assert len(samples) == 12
    special_key = _key_of(samples, F(6))
    assert table[LAB_A][special_key] == 0  # zero-filled exact zero
    assert table[LAB_B][special_key] == _coeff_b(*SPECIAL_POINT)


def test_rank_drop_sample_rejected_as_unstable():
    records = _base_records(GENERIC_POINTS)
    bad_pt = (F(7), F(2, 3))
    # One prime collapses to a lower rank (support loss), two stay at generic rank:
    # pre-11b the two generic-rank records would have been kept and A zero-filled.
    records.append(_record(bad_pt, PRIMES[0], rank=GENERIC_RANK - 1, drop=(LAB_A,)))
    records.append(_record(bad_pt, PRIMES[1]))
    records.append(_record(bad_pt, PRIMES[2]))
    cls = classify_sample_supports(records)
    assert cls["summary"]["n_unstable"] == 1
    assert cls["summary"]["unstable_samples"] == ["ep=7,r=2/3"]
    assert cls["by_sample"]["ep=7,r=2/3"]["classification"] == SUPPORT_UNSTABLE

    labels, table, samples, n_skipped = collect_value_table(records)
    assert len(samples) == 11  # the unstable sample is rejected entirely
    assert n_skipped == 3
    assert all(F(s["ep"]) != 7 for s in samples.values())


def test_prime_disagreement_rejected_as_unstable():
    records = _base_records(GENERIC_POINTS)
    bad_pt = (F(8), F(3, 5))
    records.append(_record(bad_pt, PRIMES[0], drop=(LAB_A,)))  # same rank, missing A
    records.append(_record(bad_pt, PRIMES[1]))
    records.append(_record(bad_pt, PRIMES[2]))
    cls = classify_sample_supports(records)
    entry = cls["by_sample"]["ep=8,r=3/5"]
    assert entry["classification"] == SUPPORT_UNSTABLE  # primes disagree -> not a zero

    labels, table, samples, n_skipped = collect_value_table(records)
    assert len(samples) == 11
    assert n_skipped == 3


def test_all_stable_supports_unchanged():
    records = _base_records(GENERIC_POINTS)
    cls = classify_sample_supports(records)
    summary = cls["summary"]
    assert summary["n_stable"] == summary["n_samples"] == 11
    assert summary["n_special_zero"] == 0 and summary["n_unstable"] == 0
    assert all(e["classification"] == SUPPORT_STABLE for e in cls["by_sample"].values())

    labels, table, samples, n_skipped = collect_value_table(records)
    assert n_skipped == 0
    assert len(samples) == 11
    key = _key_of(samples, F(1, 5))
    assert table[LAB_A][key] == _coeff_a(F(1, 5), _r_of(1))


def test_end_to_end_reconstruction_through_special_zero():
    records = _base_records(GENERIC_POINTS + [SPECIAL_POINT])
    report: dict = {}
    coeffs = reconstruct_coefficients(records, ["ep", "r"], report=report)
    ep, r = sp.symbols("ep r")
    assert sp.simplify(coeffs[LAB_A] - (ep - 6) / (r + 1)) == 0
    assert sp.simplify(coeffs[LAB_B] - 1 / (ep + r)) == 0
    # The reconstructed coefficient honours the special zero exactly.
    assert coeffs[LAB_A].subs({ep: 6, r: sp.Rational(57, 11)}) == 0
    info = report[str(LAB_A)]
    assert info["status"] == "validated"
    assert (info["num_deg"], info["den_deg"]) == (1, 1)
    assert info["n_fit"] + info["n_hold"] == 12


def test_underdetermined_failure_reports_point_deficit():
    # Degree-(3,1) truth with far too few points: every viable degree pair is either
    # underdetermined or fails holdout; the report must expose the exact deficits (#12).
    points = [(F(k, 3), F(k + 1, 4)) for k in range(1, 9)]
    values = {pt: (pt[0] ** 3 + 1) / (pt[1] + 1) for pt in points}
    info: dict = {}
    with pytest.raises(InterpolationFailed):
        interpolate_multivariate(values, ["ep", "r"], info=info)
    assert info["status"] == "failed"
    assert info["n_fit"] == 6 and info["n_hold"] == 2
    under = [a for a in info["attempts"] if a["status"].startswith("underdetermined")]
    assert under, "expected underdetermined attempts in the report"
    # The true degree pair (3,1) needs 10+3-1=12 fit points; deficit must be stated.
    (need_31,) = [a for a in under if (a["num_deg"], a["den_deg"]) == (3, 1)]
    assert "needs 12 fit points, have 6" in need_31["status"]
