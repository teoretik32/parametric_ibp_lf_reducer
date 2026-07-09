"""Perf.1: ranking is hoisted out of the per-(prime, sample) record loop.

The ranking is a pure function of the label set (order is total per-label:
``(tier, -complexity, label)``), so building it once per run and reusing it at every
normal-form point must be bit-for-bit identical to the old per-record ranking. These tests
pin that: call count == 1, record-level equality with the un-hoisted path, and the timing
schema (``ranking_once`` populated, per-record ``ranking`` now 0.0).
"""

from __future__ import annotations

from fractions import Fraction

import sympy as sp

from parametric_ibp_lf_reducer import (
    STATUS_SUCCESS,
    ParamExpr,
    Row,
    parse_family_text,
    reduce_rows_once,
)
from parametric_ibp_lf_reducer.modular_normal_form import modular_normal_form
from parametric_ibp_lf_reducer.ranking import rank_labels
from parametric_ibp_lf_reducer.records import (
    collect_normal_form_records,
    record_from_result,
)
from parametric_ibp_lf_reducer.timing import STAGE_KEYS

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


def test_rank_labels_called_exactly_once_for_many_records(monkeypatch):
    """15 (prime, sample) records must trigger exactly one rank_labels call (the hoist)."""
    import importlib

    # NB: the package __init__ re-exports the *function* modular_normal_form, which shadows
    # the submodule under ``import pkg.mod as m`` — go through importlib for the real module.
    mnf_mod = importlib.import_module("parametric_ibp_lf_reducer.modular_normal_form")
    records_mod = importlib.import_module("parametric_ibp_lf_reducer.records")

    calls = {"n": 0}

    def counting_rank_labels(*args, **kwargs):
        calls["n"] += 1
        return rank_labels(*args, **kwargs)

    monkeypatch.setattr(records_mod, "rank_labels", counting_rank_labels)
    monkeypatch.setattr(mnf_mod, "rank_labels", counting_rank_labels)

    fam = _family()
    records = collect_normal_form_records(
        fam, _tiny_rows(fam), T, PRIMES, _samples([1, 2, 3, 4, 5]), lf_map=LF_MAP
    )
    assert len(records) == 15
    assert all(rec.status == "Reduced" for rec in records)
    assert calls["n"] == 1


def test_hoisted_records_identical_to_per_point_ranking():
    """collect (ranking once) == manual per-point modular_normal_form (ranking per record)."""
    fam = _family()
    rows = _tiny_rows(fam)
    samples = _samples([1, 2, 3, 4, 5])

    hoisted = collect_normal_form_records(fam, rows, T, PRIMES, samples, lf_map=LF_MAP)
    per_point = [
        record_from_result(
            modular_normal_form(fam, rows, T, dict(sample), prime, lf_map=LF_MAP)
        )
        for sample in samples
        for prime in PRIMES
    ]
    assert hoisted == per_point  # bit-for-bit: statuses, coeffs, support, LF flags, ranks


def test_ranking_timing_moved_out_of_records_loop():
    fam = _family()
    res = reduce_rows_once(
        fam, T, [T, M], _tiny_rows(fam), PRIMES, _samples([1, 2, 3, 4, 5]), lf_flags=LF_MAP
    )
    assert res.status == STATUS_SUCCESS  # gates and math untouched by the hoist
    timings = res.diagnostics.extra["timings"]
    assert set(timings) == set(STAGE_KEYS)  # schema intact, ranking_once included
    assert timings["ranking"] == 0.0  # per-record ranking never runs anymore
    assert timings["ranking_once"] > 0.0  # the hoisted ranking really ran, once
    assert timings["ranking_once"] <= timings["records_total"] + 1e-6
