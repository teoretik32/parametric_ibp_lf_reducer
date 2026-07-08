"""Public text API tests (Pass 2I.3): reduce_wolfram_style_input / _to_text.

Covers: callability, a tiny explicit family end-to-end, the honest
ParserNeedsExplicitFamily failure for integrand-only input, Wolfram-style output
formatting (caret powers, association syntax, never Python ``**``), and the CLI
staying a non-success path.
"""

from __future__ import annotations

import pytest

import parametric_ibp_lf_reducer as pkg
from parametric_ibp_lf_reducer import (
    FAILURE_PARSER_NEEDS_EXPLICIT_FAMILY,
    reduce_wolfram_style_input,
    reduce_wolfram_style_input_to_text,
)
from parametric_ibp_lf_reducer.result import ReductionResult

TINY_EXPLICIT = """
IBPInput = <|
  "Variables" -> {x1, x2},
  "Parameters" -> {ep},
  "Regulators" -> {ep},
  "Polynomials" -> <| "G0" -> 1 + x1 + x2 |>,
  "MonomialExponents" -> <| x1 -> 0, x2 -> 0 |>,
  "PolynomialExponents" -> <| "G0" -> -1 - ep |>
|>
"""

INTEGRAND_ONLY = """
IBPInput = <|
  "Variables" -> {x1, x2},
  "Parameters" -> {ep},
  "Regulators" -> {ep},
  "Integrand" -> x1^(-ep)*(1 + x1 + x2)^(2*ep)
|>
"""


def test_api_functions_are_callable():
    assert callable(pkg.reduce_wolfram_style_input)
    assert callable(pkg.reduce_wolfram_style_input_to_text)
    assert callable(pkg.build_reducer_config)
    assert callable(pkg.default_scattered_samples)


def test_tiny_explicit_input_runs_through_pipeline():
    res = reduce_wolfram_style_input(TINY_EXPLICIT)
    assert isinstance(res, ReductionResult)
    assert res.target_label == (0, 0, 0)
    # The formal reduction happens; the strict LF gate then reports honestly.
    assert res.status == "NormalFormNotLocallyFinite"
    assert res.formal_success is True
    assert len(res.terms) == 3


def test_integrand_only_input_returns_typed_parser_failure():
    # No factorization is guessed (spec §3.2): honest typed Failure, no exception.
    res = reduce_wolfram_style_input(INTEGRAND_ONLY)
    assert res.status == FAILURE_PARSER_NEEDS_EXPLICIT_FAMILY == "ParserNeedsExplicitFamily"
    assert res.formal_success is False
    assert res.terms == ()
    text = reduce_wolfram_style_input_to_text(INTEGRAND_ONLY)
    assert '"Status" -> "Failure"' in text
    assert "ParserNeedsExplicitFamily" in text


def test_output_text_is_wolfram_style():
    text = reduce_wolfram_style_input_to_text(TINY_EXPLICIT)
    assert "<|" in text and "|>" in text
    assert "->" in text
    assert "^" in text  # powers rendered Wolfram-style ...
    assert "**" not in text  # ... never Python-style


def test_cli_is_still_not_a_success_path():
    from parametric_ibp_lf_reducer.__main__ import main

    with pytest.raises(SystemExit) as exc:
        main()
    assert exc.value.code != 0
