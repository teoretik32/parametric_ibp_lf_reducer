"""Method.11c regression tests: post-reconstruction validation of support deviations.

The pre-reconstruction assumption ``same-rank + same missing support across primes =>
special_zero`` is gone.  A support-changing sample is excluded from fitting and marked
``support_deviation_pending_validation``; only ``validate_support_deviations`` may upgrade
it to ``GenuineSpecialZero`` (all present coefficients match AND all missing labels
evaluate to zero) -- otherwise it is a ``BasisSpecialization`` and stays excluded.

The real 152-record External Int2 cache must classify ep=6, r=57/11 as
``BasisSpecialization`` (generic rank 26984, one label disappears, all three surviving
coefficients disagree with the generic functions).
"""

from __future__ import annotations

import importlib.util
import json
import sys
from fractions import Fraction as F
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

import sympy as sp  # noqa: E402

from parametric_ibp_lf_reducer import (  # noqa: E402
    DEVIATION_BASIS_SPECIALIZATION,
    DEVIATION_GENUINE_SPECIAL_ZERO,
    SUPPORT_DEVIATION_PENDING,
    SUPPORT_UNSTABLE,
    classify_sample_supports,
    reconstruct_coefficients,
    validate_support_deviations,
)
from parametric_ibp_lf_reducer.modular_normal_form import STATUS_REDUCED  # noqa: E402
from parametric_ibp_lf_reducer.records import NormalFormRecord  # noqa: E402

PRIMES = (2147483647, 2147483629, 2147483587)
TARGET = (1, 1, 1, 1, 1, 1, -2)
LAB_A = (0, 0, 1, -1, 0, 0, -1)
LAB_B = (0, 1, 0, 0, 0, 0, -1)
GENERIC_RANK = 100

CACHE_PATH = REPO_ROOT / "validation" / "external_int2_method11_records.jsonl"
CACHE_SAMPLE = "ep=6,r=57/11"


def _coeff_a(ep: F, r: F) -> F:
    return (ep - 6) / (r + 1)  # vanishes exactly at ep=6


def _coeff_b(ep: F, r: F) -> F:
    return 1 / (ep + r)


def _r_of(k: int) -> F:
    return F((3 * k) % 7 + 1, 3)


GENERIC_POINTS = [(F(k, 5), _r_of(k)) for k in range(1, 12)]
SPECIAL_POINT = (F(6), F(57, 11))


def _mod(value: F, p: int) -> int:
    return value.numerator * pow(value.denominator, -1, p) % p


def _record(pt, prime, *, drop=(), rank=GENERIC_RANK, shift=0):  # noqa: ANN001
    ep, r = pt
    coeffs = {}
    for lab, fn in ((LAB_A, _coeff_a), (LAB_B, _coeff_b)):
        val = fn(ep, r)
        if lab in drop or val == 0:
            continue
        coeffs[lab] = (_mod(val, prime) + shift) % prime
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


def _fit_and_validate(records):
    coeffs = reconstruct_coefficients(records, ["ep", "r"])
    return coeffs, validate_support_deviations(coeffs, records)


def test_true_special_zero_upgraded_to_genuine():
    """Missing coefficient is exactly zero and every other coefficient agrees."""
    records = _base_records(GENERIC_POINTS + [SPECIAL_POINT])
    cls = classify_sample_supports(records)
    assert cls["by_sample"][CACHE_SAMPLE]["classification"] == SUPPORT_DEVIATION_PENDING
    coeffs, val = _fit_and_validate(records)
    ep, r = sp.symbols("ep r")
    assert sp.simplify(coeffs[LAB_A] - (ep - 6) / (r + 1)) == 0
    entry = val["by_sample"][CACHE_SAMPLE]
    assert entry["verdict"] == DEVIATION_GENUINE_SPECIAL_ZERO
    assert entry["n_present_mismatches"] == 0
    assert entry["missing_nonzero"] == []
    assert val["summary"]["n_genuine_special_zero"] == 1
    assert val["summary"]["n_basis_specialization"] == 0


def test_same_rank_basis_specialization_stays_excluded():
    """One label missing at generic rank but the surviving coefficient changed."""
    records = _base_records(GENERIC_POINTS)
    # Same rank at every prime, LAB_A consistently absent, LAB_B consistently shifted:
    # pre-11c this was classified special_zero and zero-filled (the poison).
    for p in PRIMES:
        records.append(_record(SPECIAL_POINT, p, drop=(LAB_A,), shift=1))
    cls = classify_sample_supports(records)
    assert cls["by_sample"][CACHE_SAMPLE]["classification"] == SUPPORT_DEVIATION_PENDING
    coeffs, val = _fit_and_validate(records)
    entry = val["by_sample"][CACHE_SAMPLE]
    assert entry["verdict"] == DEVIATION_BASIS_SPECIALIZATION
    assert entry["n_present_mismatches"] == len(PRIMES)  # LAB_B disagrees at every prime
    # LAB_A genuinely vanishes at ep=6, so the missing label itself is not the witness:
    assert entry["missing_nonzero"] == []
    assert val["summary"]["basis_specialization_samples"] == [CACHE_SAMPLE]


def test_rank_drop_remains_unstable_never_validated():
    records = _base_records(GENERIC_POINTS)
    bad_pt = (F(7), F(2, 3))
    records.append(_record(bad_pt, PRIMES[0], rank=GENERIC_RANK - 1, drop=(LAB_A,)))
    records.append(_record(bad_pt, PRIMES[1]))
    records.append(_record(bad_pt, PRIMES[2]))
    cls = classify_sample_supports(records)
    assert cls["by_sample"]["ep=7,r=2/3"]["classification"] == SUPPORT_UNSTABLE
    coeffs, val = _fit_and_validate(records)
    # unstable samples are never upgraded (no verdict entry at all):
    assert "ep=7,r=2/3" not in val["by_sample"]
    assert val["summary"]["n_pending"] == 0


def test_prime_disagreement_remains_unstable_never_validated():
    records = _base_records(GENERIC_POINTS)
    bad_pt = (F(8), F(3, 5))
    records.append(_record(bad_pt, PRIMES[0], drop=(LAB_A,)))
    records.append(_record(bad_pt, PRIMES[1]))
    records.append(_record(bad_pt, PRIMES[2]))
    cls = classify_sample_supports(records)
    assert cls["by_sample"]["ep=8,r=3/5"]["classification"] == SUPPORT_UNSTABLE
    coeffs, val = _fit_and_validate(records)
    assert "ep=8,r=3/5" not in val["by_sample"]
    assert val["summary"]["n_pending"] == 0


# --- the real 152-record External Int2 cache --------------------------------------------------
def _load_cache_records():
    path = REPO_ROOT / "scripts" / "run_external_int2_method11.py"
    spec = importlib.util.spec_from_file_location("m11c_cache_loader", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    lines = CACHE_PATH.read_text(encoding="utf-8").splitlines()
    return [mod._record_from_json(json.loads(ln)) for ln in lines[1:] if ln.strip()]


@pytest.fixture(scope="module")
def cache_fit():
    records = _load_cache_records()
    coeffs = reconstruct_coefficients(records, ["ep", "r"])
    return records, coeffs


def test_cache_sample_is_basis_specialization(cache_fit):
    """ep=6, r=57/11 in the real cache: generic rank kept, but a basis specialization."""
    records, coeffs = cache_fit
    assert len(records) == 152
    cls = classify_sample_supports(records)
    assert cls["summary"]["generic_rank"] == 26984
    assert cls["summary"]["support_deviation_samples"] == [CACHE_SAMPLE]
    assert cls["summary"]["n_stable"] == 37
    val = validate_support_deviations(coeffs, records, classification=cls)
    entry = val["by_sample"][CACHE_SAMPLE]
    assert entry["verdict"] == DEVIATION_BASIS_SPECIALIZATION
    # three surviving labels disagree at every one of the 4 primes...
    assert entry["n_present_mismatches"] == 12
    # ...and the missing label does NOT evaluate to zero there:
    assert len(entry["missing_nonzero"]) == 1
    assert entry["missing_nonzero"][0]["label"] == list(LAB_A)


def test_cache_reconstruction_recovers_generic_functions(cache_fit):
    """Excluding the deviating sample, the 37 stable samples give the exact functions."""
    _, coeffs = cache_fit
    ep, r = sp.symbols("ep r")
    expected = {
        (-1, 0, 0, -1, -1, 0, 0): -(2 * ep * r - ep + 1) / (6 * r * (3 * ep + 1)),
        (-1, 0, 0, 0, -1, 0, -1): -((ep - 1) ** 2) * (r + 1) / (6 * ep * r * (3 * ep + 1)),
        (0, 0, 0, -1, 0, 0, -1): (ep - 1)
        * (4 * ep**2 * r + ep**2 - 2 * ep * r - 1)
        / (6 * ep**2 * r * (3 * ep + 1)),
        (0, 0, 1, -1, 0, 0, -1): (ep - 1)
        * (2 * ep - 1)
        * (2 * ep * r + ep + 1)
        / (6 * ep**2 * r * (3 * ep + 1)),
    }
    assert set(coeffs) == set(expected)
    for lab, truth in expected.items():
        assert sp.simplify(coeffs[lab] - truth) == 0, lab
