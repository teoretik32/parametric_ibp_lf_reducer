"""Perf.0: lightweight stage timing diagnostics.

Timings are pure observability: every stage key is always present as a non-negative float,
the snapshot is JSON-safe (exported by the CLI payload), and adding them must not change
statuses, terms, or the certificate/LF gates. No math is asserted here beyond the existing
Success outcomes of the tiny fixtures.
"""

from __future__ import annotations

import json
from fractions import Fraction

import sympy as sp

from parametric_ibp_lf_reducer import (
    STATUS_SUCCESS,
    ParamExpr,
    ReducerConfig,
    Row,
    parse_family_text,
    reduce_family_once,
    reduce_rows_once,
)
from parametric_ibp_lf_reducer.cli import _diagnostics_payload
from parametric_ibp_lf_reducer.timing import STAGE_KEYS, StageTimings, new_stage_timings

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


def _family():
    return parse_family_text(GENERIC_FAMILY_TEXT)


def _samples(vals):
    return [{"ep": Fraction(v)} for v in vals]


def _row(params, terms):
    row = Row(kind="synthetic", provenance={})
    for label, expr in terms.items():
        row.add_term(label, ParamExpr.from_sympy(sp.sympify(expr), params))
    return row


def _tiny_success():
    fam = _family()
    row = _row(fam.parameters, {T: 1, M: "ep + 3"})
    return reduce_rows_once(
        fam, T, [T, M], [row], PRIMES, _samples([1, 2, 3, 4, 5]), lf_flags={T: True, M: True}
    )


# --- accumulator unit behaviour --------------------------------------------------------------
def test_new_stage_timings_preseeds_all_keys_to_zero():
    t = new_stage_timings()
    assert set(t) == set(STAGE_KEYS)
    assert all(v == 0.0 for v in t.values())


def test_stage_accumulates_and_survives_exceptions():
    t = StageTimings()
    with t.stage("x"):
        pass
    first = t["x"]
    assert first >= 0.0
    try:
        with t.stage("x"):
            raise RuntimeError("boom")
    except RuntimeError:
        pass
    assert t["x"] >= first  # accumulated (not overwritten), even on the exception path


# --- orchestration wiring --------------------------------------------------------------------
def test_reduce_rows_once_emits_full_timing_schema():
    res = _tiny_success()
    assert res.status == STATUS_SUCCESS  # gate untouched by timing instrumentation
    timings = res.diagnostics.extra["timings"]
    assert set(timings) == set(STAGE_KEYS)
    assert all(isinstance(v, float) and v >= 0.0 for v in timings.values())
    # ready-made rows: row-generation stages exist in the schema but never ran
    assert timings["row_generation_total"] == 0.0
    assert timings["tangent_fields"] == 0.0
    assert timings["tangent_rows"] == 0.0


def test_reduce_family_once_times_row_generation():
    fam = _family()
    cfg = ReducerConfig(
        primes=PRIMES[:2],
        samples=_samples([1, 2, 3, 4]),
        label_box=((0, 0), (-1, 0)),
        max_ibp_degree=1,
    )
    res = reduce_family_once(fam, T, cfg)
    timings = res.diagnostics.extra["timings"]
    assert set(timings) == set(STAGE_KEYS)
    assert timings["row_generation_total"] > 0.0
    assert timings["algebraic_rows"] >= 0.0
    assert timings["coordinate_rows"] >= 0.0
    # inner stages never exceed their enclosing totals (loose sanity, no flaky margins)
    assert (
        timings["algebraic_rows"]
        + timings["coordinate_rows"]
        + timings["tangent_fields"]
        + timings["tangent_rows"]
        <= timings["row_generation_total"] + 1e-6
    )
    assert (
        timings["assemble_rows_mod_p"]
        + timings["ranking"]
        + timings["rref_mod_p"]
        + timings["extract_normal_form"]
        <= timings["records_total"] + 1e-6
    )


def test_timings_are_json_safe_in_cli_payload():
    res = _tiny_success()
    payload = _diagnostics_payload(res)
    json.dumps(payload)  # must not raise
    assert set(payload["diagnostics"]["timings"]) == set(STAGE_KEYS)
    assert "certificate_total" in payload["diagnostics"]["timings"]
    # existing diagnostics contract unchanged
    for key in ("formal_success", "n_records", "n_skipped_records", "messages"):
        assert key in payload["diagnostics"]


# --- Perf.0.1: certificate stage timing -------------------------------------------------------
def test_certificate_timing_keys_on_tiny_success():
    assert "certificate_total" in STAGE_KEYS
    assert "certificate_points_total" in STAGE_KEYS
    res = _tiny_success()
    assert res.status == STATUS_SUCCESS  # gates untouched by certificate timing
    timings = res.diagnostics.extra["timings"]
    assert timings["certificate_total"] >= 0.0
    assert timings["certificate_points_total"] >= 0.0
    # per-point verification is an inner stage of the certificate step
    assert timings["certificate_points_total"] <= timings["certificate_total"] + 1e-6
    # tiny Success actually certifies (auto points), so the stage really ran
    assert res.diagnostics.extra.get("certificate_status") != "NotRun"
