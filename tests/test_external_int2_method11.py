"""Tests for External Int2 Method.11 Phase B (certified reduction under a chamber policy)."""

from __future__ import annotations

import importlib.util
import json
import sys
from dataclasses import replace
from fractions import Fraction
from pathlib import Path
from types import SimpleNamespace

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


class TestRecordCacheJson:
    """Perf.12: the JSONL record cache round-trips exactly and rejects corruption."""

    def _record(self, m11):
        return m11.NormalFormRecord(
            prime=101,
            sample={"ep": Fraction(-3, 5)},
            target_label=(0, 0, 0),
            status="Reduced",
            formal_success=True,
            coeffs={(0, 0, -1): 7, (1, 0, 0): 100},
            support=((0, 0, -1), (1, 0, 0)),
            all_terms_lf=True,
            non_lf_terms=(),
            unknown_lf_terms=((1, 0, 0),),
            rank=5,
            diagnostics={"nrows": 9, "pivot_label": (0, 0, 0)},
        )

    def test_roundtrip_is_faithful(self, m11):
        rec = self._record(m11)
        key = m11.sample_prime_key(rec.sample, rec.prime)
        back = m11._record_from_json(json.loads(json.dumps(m11._record_to_json(key, rec))))
        assert back == rec

    def test_key_is_canonical(self, m11):
        k = m11.sample_prime_key
        assert k({"ep": Fraction(2, 1)}, 7) == k({"ep": 2}, 7)
        assert k({"a": 1, "b": Fraction(1, 2)}, 7) == k({"b": Fraction(1, 2), "a": 1}, 7)
        assert k({"ep": 1}, 7) != k({"ep": 1}, 11)

    def test_load_missing_file_is_empty(self, m11, tmp_path):
        assert m11._load_record_cache(tmp_path / "none.jsonl", (0, 0, 0)) == {}

    def test_load_roundtrip_and_duplicate_last_wins(self, m11, tmp_path):
        rec = self._record(m11)
        key = m11.sample_prime_key(rec.sample, rec.prime)
        path = tmp_path / "cache.jsonl"
        lines = [
            json.dumps(m11._record_to_json(key, rec)),
            json.dumps(m11._record_to_json(key, replace(rec, rank=6))),
        ]
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        cache = m11._load_record_cache(path, (0, 0, 0))
        assert set(cache) == {key}
        assert cache[key].rank == 6

    def test_load_rejects_key_mismatch(self, m11, tmp_path):
        rec = self._record(m11)
        path = tmp_path / "cache.jsonl"
        path.write_text(
            json.dumps(m11._record_to_json("p=999;ep=-3/5", rec)) + "\n", encoding="utf-8"
        )
        with pytest.raises(SystemExit, match="key mismatch"):
            m11._load_record_cache(path, (0, 0, 0))

    def test_load_rejects_foreign_target(self, m11, tmp_path):
        rec = self._record(m11)
        key = m11.sample_prime_key(rec.sample, rec.prime)
        path = tmp_path / "cache.jsonl"
        path.write_text(json.dumps(m11._record_to_json(key, rec)) + "\n", encoding="utf-8")
        with pytest.raises(SystemExit, match="foreign target"):
            m11._load_record_cache(path, (9, 9, 9))


class TestCollectRecordCacheSeam:
    """Perf.12: cache hits skip the per-point solve; misses stream through ``on_record``."""

    @staticmethod
    def _install_stub(monkeypatch, records_mod, calls):
        def fake_mnf(family, rows, target_label, sample, prime, **kwargs):
            calls.append((tuple(sorted(sample.items())), prime))
            return SimpleNamespace(
                terms={(0, 0, -1): prime - 1},
                prime=prime,
                sample=dict(sample),
                target_label=target_label,
                status="Reduced",
                formal_success=True,
                all_terms_lf=True,
                non_lf_terms=(),
                unknown_lf_terms=(),
                rank=3,
                nrows=4,
                pivot_label=target_label,
            )

        monkeypatch.setattr(records_mod, "modular_normal_form", fake_mnf)

    def test_cache_and_stream_behavior(self, monkeypatch):
        import parametric_ibp_lf_reducer.records as records_mod

        target = (0, 0, 0)
        samples = [{"ep": Fraction(1, 3)}, {"ep": Fraction(2, 3)}]
        primes = [101, 103]
        calls: list = []
        self._install_stub(monkeypatch, records_mod, calls)
        kw = dict(preferred_masters=(), lf_map={}, ranking=object())

        baseline = records_mod.collect_normal_form_records(None, [], target, primes, samples, **kw)
        assert len(calls) == 4

        calls.clear()
        cache: dict = {}
        seen: list = []
        first = records_mod.collect_normal_form_records(
            None,
            [],
            target,
            primes,
            samples,
            record_cache=cache,
            on_record=lambda key, rec: seen.append(key),
            **kw,
        )
        assert first == baseline
        assert len(calls) == 4 and len(cache) == 4
        assert seen == [records_mod.sample_prime_key(s, p) for s in samples for p in primes]

        calls.clear()
        seen.clear()
        warm = records_mod.collect_normal_form_records(
            None,
            [],
            target,
            primes,
            samples,
            record_cache=cache,
            on_record=lambda key, rec: seen.append(key),
            **kw,
        )
        assert warm == baseline
        assert calls == [] and seen == []

    def test_partial_cache_computes_only_missing_points(self, monkeypatch):
        import parametric_ibp_lf_reducer.records as records_mod

        target = (0, 0, 0)
        samples = [{"ep": Fraction(1, 3)}, {"ep": Fraction(2, 3)}]
        calls: list = []
        self._install_stub(monkeypatch, records_mod, calls)
        kw = dict(preferred_masters=(), lf_map={}, ranking=object())

        cache: dict = {}
        records_mod.collect_normal_form_records(
            None, [], target, [101, 103], samples, record_cache=cache, **kw
        )
        calls.clear()
        out = records_mod.collect_normal_form_records(
            None, [], target, [101, 103, 107], samples, record_cache=cache, **kw
        )
        assert [p for _, p in calls] == [107, 107]  # only the new prime, per sample
        # deterministic order preserved: samples outer, primes inner
        assert [(r.sample["ep"], r.prime) for r in out] == [
            (s["ep"], p) for s in samples for p in [101, 103, 107]
        ]


class TestInterpolationFailureMessage:
    """The original reconstruction failure reason must surface in the diagnostics messages."""

    def test_reconstruction_error_text_is_propagated(self, m11, monkeypatch):
        import parametric_ibp_lf_reducer.reducer as reducer_mod

        records = [
            m11.NormalFormRecord(
                prime=p,
                sample={"ep": Fraction(1, 3)},
                target_label=(0, 0, 0),
                status="Reduced",
                formal_success=True,
                coeffs={(0, 0, -1): 2},
                support=((0, 0, -1),),
                rank=3,
            )
            for p in (101, 103)
        ]

        def boom(selected, parameters):
            raise reducer_mod.InterpolationFailed("boom-424242")

        monkeypatch.setattr(reducer_mod, "reconstruct_coefficients", boom)
        fam = parse_family_text(EP_FAMILY)
        run = reducer_mod.ReducerRunDiagnostics(n_labels=1, n_rows=1, row_diagnostics={})
        st = reducer_mod._select_and_reconstruct(
            fam, records, 1, run, reducer_mod.new_stage_timings()
        )
        assert st.interpolation_failed is True
        assert st.coeffs is None
        assert any("boom-424242" in msg for msg in st.messages)
