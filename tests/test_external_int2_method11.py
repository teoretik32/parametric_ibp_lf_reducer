"""Tests for External Int2 Method.11 Phase B (certified reduction under a chamber policy)."""

from __future__ import annotations

import importlib.util
import json
import sys
from fractions import Fraction
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from parametric_ibp_lf_reducer import SurfacePolicy, parse_family_text, zero_label  # noqa: E402
from parametric_ibp_lf_reducer.reducer import ReducerConfig  # noqa: E402
from parametric_ibp_lf_reducer.result import ReductionResult, ReductionTerm  # noqa: E402


def _load_script(name: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def m11():
    return _load_script("run_external_int2_method11")


@pytest.fixture(scope="module")
def t2(m11):
    return sys.modules["run_external_int2_t2_rankrepair"]


# x-exponent depends on the regulator: chamber points on different sides of ep = -1 decide
# the x-boundary check differently (same family as the Phase A tests).
EP_FAMILY = """
IBPInput = <|
  "Variables" -> {x, y}, "Parameters" -> {ep}, "Regulators" -> {ep},
  "Polynomials" -> <| "G0" -> 1 + x + y |>,
  "MonomialExponents" -> <| x -> ep, y -> -3 |>,
  "PolynomialExponents" -> <| "G0" -> -2 |>
|>
"""


def _small_setup(t2, surface_policy=None):
    fam = parse_family_text(EP_FAMILY)
    labels = [zero_label(2, 1), (1, 0, 0), (0, 1, 0), (0, 0, -1)]
    cfg = ReducerConfig(
        labels=tuple(labels),
        max_ibp_degree=t2.BASELINE["max_ibp_degree"],
        tangent_degree_blocks=t2.BASELINE["tangent_degree_blocks"],
        surface_policy=surface_policy,
    )
    return fam, labels, cfg


def _keys(asm) -> list:
    return sorted(r.dedup_key() for r in asm["merged"])


class TestChamberPolicy:
    def test_exact_point(self, m11):
        pol = m11.chamber_policy(Fraction(-3, 5))
        assert pol.mode == "chamber"
        assert pol.point == {"ep": Fraction(-3, 5)}

    def test_float_rejected(self, m11):
        with pytest.raises(TypeError):
            m11.chamber_policy(-0.6)

    def test_default_ep_matches_method10_chamber(self, m11):
        assert m11.DEFAULT_EP == Fraction(-3, 5)


class TestAssembleChamberRows:
    def test_limit_mode_equals_t2_assembly(self, m11, t2):
        """No policy set -> identical row set to the pre-policy T2 assembly (Phase A contract)."""
        fam, labels, cfg = _small_setup(t2)
        assert _keys(m11.assemble_chamber_rows(fam, labels, cfg)) == _keys(
            t2.assemble_level_rows(fam, labels, cfg)
        )

    def test_explicit_limit_equals_default(self, m11, t2):
        fam, labels, cfg = _small_setup(t2)
        _, _, cfg_limit = _small_setup(t2, surface_policy=SurfacePolicy.limit("minus"))
        assert _keys(m11.assemble_chamber_rows(fam, labels, cfg_limit)) == _keys(
            m11.assemble_chamber_rows(fam, labels, cfg)
        )

    def test_chamber_policy_recorded_and_counts_consistent(self, m11, t2):
        fam, labels, cfg = _small_setup(t2, surface_policy=m11.chamber_policy(Fraction(-3, 5)))
        asm = m11.assemble_chamber_rows(fam, labels, cfg)
        assert asm["row_diagnostics"]["surface_policy"] == {
            "mode": "chamber",
            "point": {"ep": "-3/5"},
        }
        assert asm["n_rows_total"] == asm["n_base_rows"] + asm["n_extra_new"]
        assert sum(asm["by_kind"].values()) == asm["n_rows_total"]

    def test_far_chamber_changes_row_set(self, m11, t2):
        """ep = -3/2 flips the x-boundary exponent sign -> different surviving rows."""
        fam, labels, cfg = _small_setup(t2)
        _, _, cfg_far = _small_setup(
            t2, surface_policy=SurfacePolicy.chamber({"ep": Fraction(-3, 2)})
        )
        k_limit = _keys(m11.assemble_chamber_rows(fam, labels, cfg))
        k_far = _keys(m11.assemble_chamber_rows(fam, labels, cfg_far))
        assert k_far != k_limit


class TestCrosscheckMethod10:
    def _artifact(self, tmp_path, **over):
        data = {"chamber_ep": "-3/5", "level": 0, "rows": {"n_chamber": 100}}
        data.update(over)
        p = tmp_path / "m10.json"
        p.write_text(json.dumps(data), encoding="utf-8")
        return p

    def test_missing_artifact(self, m11, tmp_path):
        cross = m11.crosscheck_method10(1, Fraction(-3, 5), tmp_path / "absent.json")
        assert cross["available"] is False and cross["matches"] is None

    def test_match_and_mismatch(self, m11, tmp_path):
        art = self._artifact(tmp_path)
        ok = m11.crosscheck_method10(100, Fraction(-3, 5), art)
        assert ok["comparable"] is True and ok["matches"] is True
        bad = m11.crosscheck_method10(101, Fraction(-3, 5), art)
        assert bad["comparable"] is True and bad["matches"] is False

    def test_other_chamber_point_not_comparable(self, m11, tmp_path):
        art = self._artifact(tmp_path)
        cross = m11.crosscheck_method10(100, Fraction(-1, 2), art)
        assert cross["available"] is True
        assert cross["comparable"] is False and cross["matches"] is None

    def test_recorded_artifact_pins_row_total(self, m11):
        """The committed Method.10 artifact fixes the Level-0 chamber row total at 49439."""
        if not m11.METHOD10_ARTIFACT.exists():
            pytest.skip("recorded Method.10 artifact not present")
        cross = m11.crosscheck_method10(49439, Fraction(-3, 5))
        assert cross["comparable"] is True and cross["matches"] is True
        assert cross["recorded_n_chamber"] == 49439


class TestResultPayload:
    def test_json_ready(self, m11):
        result = ReductionResult(
            status="Success",
            target_label=(0, 0, 0),
            all_locally_finite=True,
            terms=(
                ReductionTerm(
                    label=(0, 0, -1),
                    coefficient_text="1/(1 - ep)",
                    integrand_text="1/G0",
                    locally_finite=True,
                ),
            ),
        )
        payload = m11.result_payload(result)
        json.dumps(payload, default=str)  # must not raise
        assert payload["status"] == "Success"
        assert payload["n_terms"] == 1
        assert payload["terms"][0]["label"] == [0, 0, -1]
        assert payload["diagnostics"]["formal_success"] is False


class TestCLIGates:
    @pytest.mark.parametrize(
        "argv",
        [
            ["--primes", "2,5"],
            ["--primes", ""],
            ["--n-samples", "0"],
            ["--n-samples", "38"],
            ["--min-valid-records", "0"],
            ["--min-certificate-points", "0"],
            ["--jobs", "0"],
        ],
    )
    def test_bad_args_exit_before_any_work(self, m11, argv):
        with pytest.raises(SystemExit) as exc:
            m11.main(argv + ["--input", "does_not_exist.wl.txt"])
        assert exc.value.code == 2


class TestScopeNote:
    def test_no_analytic_value_claim(self, m11):
        assert "does NOT claim the analytic Laurent value" in m11.SCOPE_NOTE
        assert "default surface policy" in m11.SCOPE_NOTE
