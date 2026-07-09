"""Perf.3: parallel (prime, sample) record collection must be bit-identical to serial.

``jobs`` selects how many worker processes compute the independent ``(prime, sample)``
points. ``jobs=1`` is the exact serial path; for ``jobs>1`` the results come back via
``ProcessPoolExecutor.map`` (order-preserving), so records, statuses and coefficients
must not change at all. These tests pin: record-level equality, full-pipeline equality,
input validation, the timing caveat (per-point stage keys read 0.0 with ``jobs>1``),
and the config coercion for the ``jobs`` override.
"""

from __future__ import annotations

from fractions import Fraction

import pytest
import sympy as sp

from parametric_ibp_lf_reducer import (
    STATUS_SUCCESS,
    ParamExpr,
    Row,
    parse_family_text,
    reduce_rows_once,
)
from parametric_ibp_lf_reducer.api import build_reducer_config
from parametric_ibp_lf_reducer.records import collect_normal_form_records

PRIMES = [2_147_483_647, 2_147_483_629, 2_147_483_587]

GENERIC_FAMILY_TEXT = """
IBPInput = <|
  "Variables" -> {u, v},
  "Parameters" -> {ep},
  "Regulators" -> {ep},
  "Domain" -> "PositiveOrthant",
  "Polynomials" -> <| "P0" -> 1 + u, "P1" -> 1 + v |>,
  "MonomialExponents" -> <| u -> -1 - ep, v -> ep |>,
  "PolynomialExponents" -> <| "P0" -> -1 + ep, "P1" -> -2 - ep |>,
  "TargetMultiplier" -> 1
|>;
"""

T = (0, 0, 0, 0)
M = (0, 0, -1, 0)
LF_MAP = {T: True, M: True}


def _family():
    return parse_family_text(GENERIC_FAMILY_TEXT)


def _samples(vals):
    return [{"ep": Fraction(v)} for v in vals]


def _row(params, terms):
    row = Row(kind="synthetic", provenance={})
    for label, expr in terms.items():
        row.add_term(label, ParamExpr.from_sympy(sp.sympify(expr), params))
    return row


def _tiny_rows(fam):
    return [_row(fam.parameters, {T: 1, M: "ep + 3"})]


# --- record-level equality ---------------------------------------------------------------------
def test_jobs2_records_bit_identical_to_serial():
    """jobs=2 must return the same records in the same order as the serial path."""
    fam = _family()
    rows = _tiny_rows(fam)
    samples = _samples([1, 2, 3, 4, 5])

    serial = collect_normal_form_records(fam, rows, T, PRIMES, samples, lf_map=LF_MAP, jobs=1)
    parallel = collect_normal_form_records(fam, rows, T, PRIMES, samples, lf_map=LF_MAP, jobs=2)

    assert len(serial) == 15
    assert parallel == serial  # bit-for-bit: statuses, coeffs, support, LF flags, ranks


def test_jobs_exceeding_task_count_still_identical():
    """max_workers is clamped to len(tasks); oversized jobs must not change anything."""
    fam = _family()
    rows = _tiny_rows(fam)
    samples = _samples([1, 2])

    serial = collect_normal_form_records(fam, rows, T, PRIMES, samples, lf_map=LF_MAP)
    parallel = collect_normal_form_records(fam, rows, T, PRIMES, samples, lf_map=LF_MAP, jobs=64)
    assert parallel == serial


def test_single_task_stays_on_serial_path():
    """One (prime, sample) point never spawns a pool (jobs>1 requires len(tasks)>1)."""
    fam = _family()
    rows = _tiny_rows(fam)
    serial = collect_normal_form_records(fam, rows, T, PRIMES[:1], _samples([1]), lf_map=LF_MAP)
    parallel = collect_normal_form_records(
        fam, rows, T, PRIMES[:1], _samples([1]), lf_map=LF_MAP, jobs=8
    )
    assert len(parallel) == 1
    assert parallel == serial


# --- validation --------------------------------------------------------------------------------
@pytest.mark.parametrize("bad", [0, -1, True, 2.0, "2", None])
def test_jobs_rejects_non_positive_and_non_int(bad):
    fam = _family()
    with pytest.raises(ValueError):
        collect_normal_form_records(
            fam, _tiny_rows(fam), T, PRIMES[:1], _samples([1]), lf_map=LF_MAP, jobs=bad
        )


# --- full pipeline -----------------------------------------------------------------------------
def test_full_pipeline_jobs2_matches_serial():
    """reduce_rows_once(jobs=2) must match jobs=1 everywhere except wall-clock timings."""
    fam = _family()
    samples = _samples([1, 2, 3, 4, 5])

    r1 = reduce_rows_once(fam, T, [T, M], _tiny_rows(fam), PRIMES, samples, lf_flags=LF_MAP, jobs=1)
    r2 = reduce_rows_once(fam, T, [T, M], _tiny_rows(fam), PRIMES, samples, lf_flags=LF_MAP, jobs=2)

    assert r1.status == STATUS_SUCCESS
    assert r2.status == r1.status
    assert r2.target_label == r1.target_label
    assert r2.all_locally_finite == r1.all_locally_finite
    assert r2.formal_success == r1.formal_success
    assert r2.terms == r1.terms  # reconstructed coefficients are bit-identical

    d1, d2 = r1.diagnostics, r2.diagnostics
    assert d2.n_records == d1.n_records
    assert d2.n_skipped_records == d1.n_skipped_records
    for key in (
        "n_records",
        "n_reduced_records",
        "n_selected_records",
        "n_bad_specializations",
        "n_skipped_records",
        "record_selection",
        "reconstruction_diagnostics",
    ):
        assert d2.extra[key] == d1.extra[key]


def test_jobs2_timing_caveat_per_point_stages_zero():
    """With jobs>1, per-point stage keys accumulate in workers and read 0.0 at the caller;
    the caller-side records_total stage still measures the whole collection."""
    fam = _family()
    samples = _samples([1, 2, 3, 4, 5])
    res = reduce_rows_once(
        fam, T, [T, M], _tiny_rows(fam), PRIMES, samples, lf_flags=LF_MAP, jobs=2
    )
    assert res.status == STATUS_SUCCESS
    timings = res.diagnostics.extra["timings"]
    assert timings["records_total"] > 0.0
    for key in ("assemble_rows_mod_p", "rref_mod_p", "extract_normal_form"):
        assert timings[key] == 0.0


# --- config plumbing ---------------------------------------------------------------------------
def test_config_jobs_override_coerced_and_defaulted():
    fam = _family()
    _, cfg_default = build_reducer_config(fam, {})
    assert cfg_default.jobs == 1  # serial by default; parallelism is strictly opt-in

    _, cfg = build_reducer_config(fam, {"jobs": 4})
    assert cfg.jobs == 4

    _, cfg_str = build_reducer_config(fam, {"jobs": "4"})
    assert cfg_str.jobs == 4  # _as_int coercion (CLI/document overrides arrive as text)
