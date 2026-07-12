"""Adaptive.1 tests: opt-in deterministic schedule over the fixed-pass reducer.

The tiny 2-var generic family runs the FULL pipeline (certificate gate included) in well under
a second per level, so real-math tests are not mocked. Loop-policy tests (best-partial order,
non-LF/Unknown honesty) use stubbed ``reduce_family_once`` results — they test the adaptive
orchestration only, never the math. Heavy D4/Example4 runs stay opt-in elsewhere.
"""

from __future__ import annotations

import json
from dataclasses import replace

import pytest

from parametric_ibp_lf_reducer import (
    FAILURE_INTERPOLATION_FAILED,
    FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE,
    FAILURE_RESOURCE_LIMIT_REACHED,
    FAILURE_TARGET_NOT_REDUCIBLE,
    FAILURE_VERIFICATION_FAILED,
    STATUS_SUCCESS,
    AdaptiveSearchConfig,
    SearchLevel,
    build_reducer_config,
    default_search_levels,
    parse_family_text,
    reduce_family_adaptive,
    reduce_family_once,
    reduce_wolfram_style_input_adaptive,
)
from parametric_ibp_lf_reducer import adaptive as adaptive_mod
from parametric_ibp_lf_reducer.adaptive import _extended_primes, _extended_samples
from parametric_ibp_lf_reducer.cli import EXIT_SUCCESS, EXIT_USAGE, main
from parametric_ibp_lf_reducer.finite_field import is_probable_prime
from parametric_ibp_lf_reducer.result import ReductionDiagnostics, ReductionResult

TINY_TEXT = """
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

# Same mathematics under arbitrary variable/polynomial/parameter names (no D4/Example4 or
# name hardcode anywhere in the adaptive layer).
RENAMED_TEXT = """
IBPInput = <|
  "Variables" -> {xx, yy},
  "Parameters" -> {mu},
  "Regulators" -> {mu},
  "Domain" -> "PositiveOrthant",
  "Polynomials" -> <| "QA" -> 1 + xx, "RB" -> 1 + yy |>,
  "MonomialExponents" -> <| xx -> -1 - mu, yy -> mu |>,
  "PolynomialExponents" -> <| "QA" -> -1 + mu, "RB" -> -2 - mu |>,
  "TargetMultiplier" -> 1
|>;
"""

GOOD_BOX = ((0, 0), (-1, 0))  # the package-default box; tiny family certifies here
OFF_BOX = ((1, 1), (0, 0))  # anchored away from the target -> honest TargetNotReducible
OFF_BOX_2 = ((2, 2), (0, 0))  # a second off-target box (for exhaustion tests)


def _level(name: str, box) -> SearchLevel:
    return SearchLevel(name=name, label_box=box, max_ibp_degree=1, tangent_degree_blocks=())


@pytest.fixture(scope="module")
def tiny():
    family = parse_family_text(TINY_TEXT)
    target, config = build_reducer_config(family)
    return family, target, config


def _strip_wall(adaptive_dict: dict) -> dict:
    """Deep copy of the adaptive payload with the (wall-clock) observability fields removed."""
    out = json.loads(json.dumps(adaptive_dict))
    for level in out["levels"]:
        level.pop("wall_seconds")
    return out


# --- escalation and stop-on-success -----------------------------------------------------------
def test_level0_fails_level1_succeeds(tiny):
    family, target, config = tiny
    search = AdaptiveSearchConfig(levels=(_level("off", OFF_BOX), _level("good", GOOD_BOX)))
    res = reduce_family_adaptive(family, target, config, search)

    assert res.success and res.status == STATUS_SUCCESS
    ad = res.diagnostics.extra["adaptive"]
    assert ad["stop_reason"] == "success"
    assert ad["n_levels_run"] == 2
    assert ad["best_level"] == 1
    assert ad["levels"][0]["status"] == FAILURE_TARGET_NOT_REDUCIBLE
    assert "expand the label box" in ad["levels"][0]["recommendation"]
    assert ad["levels"][1]["status"] == STATUS_SUCCESS
    assert ad["levels"][1]["recommendation"] is None
    # Success came from the ordinary strict gate: certified, reconstruction verified.
    assert res.diagnostics.extra["certificate"]["certificate_status"] == "Passed"
    assert res.diagnostics.reconstruction_verified


def test_stops_after_first_certified_success(tiny):
    family, target, config = tiny
    search = AdaptiveSearchConfig(levels=(_level("good", GOOD_BOX), _level("off", OFF_BOX)))
    res = reduce_family_adaptive(family, target, config, search)

    ad = res.diagnostics.extra["adaptive"]
    assert res.success
    assert ad["n_levels_planned"] == 2
    assert ad["n_levels_run"] == 1  # the second level never ran
    assert len(ad["levels"]) == 1
    assert ad["best_level"] == 0


def test_all_levels_fail_returns_deterministic_best_partial_with_history(tiny):
    family, target, config = tiny
    search = AdaptiveSearchConfig(levels=(_level("off1", OFF_BOX), _level("off2", OFF_BOX_2)))
    res = reduce_family_adaptive(family, target, config, search)

    assert not res.success
    assert res.status == FAILURE_TARGET_NOT_REDUCIBLE
    ad = res.diagnostics.extra["adaptive"]
    assert ad["stop_reason"] == "levels_exhausted"
    assert [lv["ran"] for lv in ad["levels"]] == [True, True]
    assert ad["best_level"] == 0  # equal quality -> deterministic earlier-level tie-break

    res2 = reduce_family_adaptive(family, target, config, search)
    assert _strip_wall(res2.diagnostics.extra["adaptive"]) == _strip_wall(ad)


# --- best-partial ordering (orchestration only; math stubbed) ---------------------------------
def _stub_result(status, *, non_lf=(), unknown=(), verified=False, cert=None):
    diag = ReductionDiagnostics(
        reconstruction_verified=verified,
        non_lf_terms=tuple(non_lf),
        unknown_lf_terms=tuple(unknown),
    )
    if cert is not None:
        diag.extra["certificate"] = {"certificate_status": cert}
    return ReductionResult(
        status=status,
        target_label=(0, 0, 0, 0),
        all_locally_finite="Unknown",
        terms=(),
        diagnostics=diag,
    )


def _patch_sequence(monkeypatch, results):
    remaining = list(results)

    def fake(family, target_label, config):
        return remaining.pop(0)

    monkeypatch.setattr(adaptive_mod, "reduce_family_once", fake)


def test_best_partial_prefers_reduced_target(monkeypatch, tiny):
    family, target, config = tiny
    stub0 = _stub_result(FAILURE_TARGET_NOT_REDUCIBLE)
    stub1 = _stub_result(FAILURE_INTERPOLATION_FAILED)  # target reduced -> strictly better
    _patch_sequence(monkeypatch, [stub0, stub1])

    search = AdaptiveSearchConfig(levels=(_level("a", OFF_BOX), _level("b", OFF_BOX_2)))
    res = reduce_family_adaptive(family, target, config, search)

    assert res is not stub0 and res.status == FAILURE_INTERPOLATION_FAILED
    ad = res.diagnostics.extra["adaptive"]
    assert ad["stop_reason"] == "levels_exhausted"
    assert ad["best_level"] == 1
    assert "samples and/or primes" in ad["levels"][1]["recommendation"]


def test_non_lf_and_unknown_terms_never_become_success(monkeypatch, tiny):
    family, target, config = tiny
    worse = _stub_result(
        FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE, non_lf=[(1, 0, 0, 0), (0, 1, 0, 0)], cert="Passed"
    )
    better = _stub_result(
        FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE, non_lf=[(1, 0, 0, 0)], cert="Passed"
    )
    _patch_sequence(monkeypatch, [worse, better])

    search = AdaptiveSearchConfig(levels=(_level("a", OFF_BOX), _level("b", OFF_BOX_2)))
    res = reduce_family_adaptive(family, target, config, search)

    # A passed certificate on a non-LF normal form must NOT be upgraded to Success.
    assert not res.success
    assert res.status == FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE
    ad = res.diagnostics.extra["adaptive"]
    assert ad["best_level"] == 1  # fewer non-LF terms wins
    assert "non-LF/Unknown labels" in ad["levels"][0]["recommendation"]


def test_verification_failed_recommendation_never_accepts_coefficients(monkeypatch, tiny):
    family, target, config = tiny
    stub0 = _stub_result(FAILURE_VERIFICATION_FAILED, cert="Failed")
    stub1 = _stub_result(FAILURE_VERIFICATION_FAILED, cert="Failed")
    _patch_sequence(monkeypatch, [stub0, stub1])

    search = AdaptiveSearchConfig(levels=(_level("a", OFF_BOX), _level("b", OFF_BOX_2)))
    res = reduce_family_adaptive(family, target, config, search)

    assert not res.success
    assert res.status == FAILURE_VERIFICATION_FAILED
    ad = res.diagnostics.extra["adaptive"]
    assert ad["best_level"] == 0  # identical quality -> earlier level
    assert "never accept the current coefficients" in ad["levels"][0]["recommendation"]


# --- resource limits ---------------------------------------------------------------------------
def test_max_labels_preflight_gives_honest_typed_failure(tiny):
    family, target, config = tiny
    search = AdaptiveSearchConfig(levels=(_level("good", GOOD_BOX),), max_labels=1)
    res = reduce_family_adaptive(family, target, config, search)

    assert not res.success
    assert res.status == FAILURE_RESOURCE_LIMIT_REACHED
    ad = res.diagnostics.extra["adaptive"]
    assert ad["stop_reason"] == "resource_limit"
    assert ad["n_levels_run"] == 0
    assert ad["resource_limit"] == {"kind": "max_labels", "limit": 1, "observed": 4, "level": 0}
    assert ad["levels"][0]["ran"] is False
    assert "max_labels" in ad["levels"][0]["skipped_reason"]


def test_max_rows_stops_escalation_after_level_completes(tiny):
    family, target, config = tiny
    search = AdaptiveSearchConfig(
        levels=(_level("off", OFF_BOX), _level("good", GOOD_BOX)), max_rows=1
    )
    res = reduce_family_adaptive(family, target, config, search)

    # Level 0 ran (and honestly failed); the good level was never attempted.
    assert not res.success
    assert res.status == FAILURE_TARGET_NOT_REDUCIBLE
    ad = res.diagnostics.extra["adaptive"]
    assert ad["stop_reason"] == "resource_limit"
    assert ad["n_levels_run"] == 1
    assert ad["resource_limit"]["kind"] == "max_rows"
    assert ad["resource_limit"]["limit"] == 1
    assert ad["resource_limit"]["observed"] > 1


def test_timeout_sec_checked_between_levels_only(tiny):
    family, target, config = tiny
    search = AdaptiveSearchConfig(
        levels=(_level("off1", OFF_BOX), _level("off2", OFF_BOX_2)), timeout_sec=0.0
    )
    res = reduce_family_adaptive(family, target, config, search)

    ad = res.diagnostics.extra["adaptive"]
    assert ad["n_levels_run"] == 1  # level 0 always runs; the budget gates level 1
    assert ad["stop_reason"] == "resource_limit"
    assert ad["resource_limit"]["kind"] == "timeout_sec"
    assert ad["resource_limit"]["level"] == 1
    assert not res.success  # honest failure from the completed level, not a fabricated one


# --- equivalence with the fixed pass -----------------------------------------------------------
def test_adaptive_success_equals_fixed_pass_on_same_level(tiny):
    family, target, config = tiny
    fixed_cfg = replace(
        config, label_box=GOOD_BOX, labels=None, max_ibp_degree=1, tangent_degree_blocks=None
    )
    fixed = reduce_family_once(family, target, fixed_cfg)
    res = reduce_family_adaptive(
        family, target, config, AdaptiveSearchConfig(levels=(_level("good", GOOD_BOX),))
    )

    assert fixed.success and res.success
    assert res.status == fixed.status
    assert [(t.label, t.coefficient_text) for t in res.terms] == [
        (t.label, t.coefficient_text) for t in fixed.terms
    ]
    assert (
        res.diagnostics.extra["certificate"]["certificate_status"]
        == fixed.diagnostics.extra["certificate"]["certificate_status"]
        == "Passed"
    )


# --- default schedule ---------------------------------------------------------------------------
def test_default_schedule_shape_and_determinism(tiny):
    family, _target, config = tiny
    levels = default_search_levels(family, config)
    assert levels == default_search_levels(family, config)  # deterministic
    assert [lv.name for lv in levels] == ["base", "expand-1", "deep"]

    base, expand, deep = levels
    assert base.max_ibp_degree == 1 and base.tangent_degree_blocks == ()
    assert expand.max_ibp_degree == 2 and expand.tangent_degree_blocks == ((1, 1),)
    assert deep.tangent_degree_blocks == ((1, 1), (2, 2))
    assert deep.extra_samples > 0 and deep.extra_primes > 0

    # m-ranges deepen by one per level; n-ranges never change.
    assert base.label_box == (((0, 0), (0, 0)), ((-1, 0), (-1, 0)))
    assert expand.label_box == (((0, 0), (0, 0)), ((-2, 0), (-2, 0)))
    assert deep.label_box == (((0, 0), (0, 0)), ((-3, 0), (-3, 0)))

    # An explicit label list cannot be grown deterministically.
    frozen = replace(config, labels=((0, 0, 0, 0),), label_box=None)
    with pytest.raises(ValueError, match="explicit"):
        default_search_levels(family, frozen)


def test_extended_primes_and_samples_are_deterministic_and_valid(tiny):
    family, _target, config = tiny
    primes = _extended_primes(config.primes, 2)
    assert primes == _extended_primes(config.primes, 2)
    assert primes[: len(list(config.primes))] == [int(p) for p in config.primes]
    assert len(primes) == len(list(config.primes)) + 2
    assert len(set(primes)) == len(primes)
    assert all(is_probable_prime(p) for p in primes)

    samples = _extended_samples(family, config.samples, 4)
    assert samples == _extended_samples(family, config.samples, 4)
    assert len(samples) == len(list(config.samples)) + 4
    keys = {tuple(sorted(pt.items())) for pt in samples}
    assert len(keys) == len(samples)  # no duplicate points
    assert all(set(pt) == set(family.parameters) for pt in samples)


# --- arbitrary names + API smoke ----------------------------------------------------------------
def test_arbitrary_names_through_adaptive_api():
    res = reduce_wolfram_style_input_adaptive(RENAMED_TEXT)
    assert res.success
    ad = res.diagnostics.extra["adaptive"]
    assert ad["stop_reason"] == "success"
    assert res.diagnostics.extra["certificate"]["certificate_status"] == "Passed"


def test_api_adaptive_escalation_with_explicit_schedule():
    search = AdaptiveSearchConfig(levels=(_level("off", OFF_BOX), _level("good", GOOD_BOX)))
    res = reduce_wolfram_style_input_adaptive(TINY_TEXT, search=search)
    assert res.success
    assert res.diagnostics.extra["adaptive"]["best_level"] == 1


# --- CLI ----------------------------------------------------------------------------------------
@pytest.fixture()
def tiny_input(tmp_path):
    p = tmp_path / "input.m"
    p.write_text(TINY_TEXT, encoding="utf-8")
    return p


def test_cli_adaptive_success_smoke(tiny_input, tmp_path):
    out = tmp_path / "result.m"
    diag = tmp_path / "diagnostics.json"
    rc = main(
        [
            "reduce",
            str(tiny_input),
            "--adaptive",
            "--out",
            str(out),
            "--diagnostics-json",
            str(diag),
        ]
    )
    assert rc == EXIT_SUCCESS
    assert '"Status" -> "Success"' in out.read_text(encoding="utf-8")
    payload = json.loads(diag.read_text(encoding="utf-8"))
    assert payload["status"] == STATUS_SUCCESS
    assert payload["adaptive"]["stop_reason"] == "success"
    assert payload["adaptive"]["n_levels_run"] >= 1


def test_cli_adaptive_max_levels_flag(tiny_input, tmp_path):
    diag = tmp_path / "d.json"
    rc = main(
        [
            "reduce",
            str(tiny_input),
            "--adaptive",
            "--adaptive-max-levels",
            "1",
            "--diagnostics-json",
            str(diag),
        ]
    )
    assert rc == EXIT_SUCCESS
    payload = json.loads(diag.read_text(encoding="utf-8"))
    assert payload["adaptive"]["n_levels_planned"] == 1


def test_cli_adaptive_max_levels_requires_adaptive(tiny_input):
    assert main(["reduce", str(tiny_input), "--adaptive-max-levels", "2"]) == EXIT_USAGE


def test_cli_adaptive_max_levels_must_be_positive(tiny_input):
    rc = main(["reduce", str(tiny_input), "--adaptive", "--adaptive-max-levels", "0"])
    assert rc == EXIT_USAGE


def test_cli_fixed_path_unchanged_without_adaptive_flag(tiny_input, tmp_path):
    diag = tmp_path / "d.json"
    rc = main(["reduce", str(tiny_input), "--diagnostics-json", str(diag)])
    assert rc == EXIT_SUCCESS
    payload = json.loads(diag.read_text(encoding="utf-8"))
    assert "adaptive" not in payload
