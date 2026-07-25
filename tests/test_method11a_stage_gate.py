"""Tests for scripts/method11a_stage_gate.py (Phase C stage-gate audit)."""

import importlib.util
import json
from fractions import Fraction
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]


def _load(name, relpath):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relpath)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def gate():
    return _load("gate_under_test", "scripts/method11a_stage_gate.py")


@pytest.fixture(scope="module")
def m11(gate):
    return gate.load_m11()


def _record(m11, ep, prime, rank=9, support=((0, 0, -1), (1, 0, 0)), formal_success=True):
    return m11.NormalFormRecord(
        prime=prime,
        sample={"ep": Fraction(ep, 7)},
        target_label=(0, 0, 0),
        status="Reduced" if formal_success else "Unknown",
        formal_success=formal_success,
        coeffs={t: 1 for t in support},
        support=tuple(support),
        all_terms_lf=True,
        non_lf_terms=(),
        unknown_lf_terms=(),
        rank=rank,
        diagnostics={},
    )


def _grid(m11, points=(1, 2), primes=(101, 103), **kw):
    return [_record(m11, ep, p, **kw) for ep in points for p in primes]


CAPACITY_FAIL = {
    "reconstruction": {"formal_success": False, "status": "Failure(InterpolationFailed)"}
}
OTHER_FAIL = {"reconstruction": {"formal_success": False, "status": "Failure(rank collapse)"}}
SUCCESS = {"reconstruction": {"formal_success": True, "status": "Success"}}


class TestAuditRecords:
    def test_uniform_grid_is_clean(self, gate, m11):
        rep = gate.audit_records(_grid(m11), m11, expected_points=2, expected_primes=2)
        assert rep["n_records"] == 4
        assert rep["n_formal_success"] == 4
        assert rep["max_rank"] == 9
        assert rep["n_rank_drops"] == 0
        assert rep["support_deviations"] == []
        assert rep["n_points_seen"] == 2
        assert rep["points_incomplete_primes"] == []

    def test_support_deviation_is_reported_with_diff(self, gate, m11):
        recs = _grid(m11)
        recs[3] = _record(m11, 2, 103, support=((0, 0, -1), (2, 0, 0)))
        rep = gate.audit_records(recs, m11, 2, 2)
        (dev,) = rep["support_deviations"]
        assert dev["key"] == m11.sample_prime_key({"ep": Fraction(2, 7)}, 103)
        assert dev["missing"] == [[1, 0, 0]] and dev["extra"] == [[2, 0, 0]]

    def test_rank_drop_excludes_record_from_top(self, gate, m11):
        recs = _grid(m11) + [_record(m11, 3, 101, rank=5)]
        rep = gate.audit_records(recs, m11, 3, 2)
        assert rep["n_rank_drops"] == 1
        assert rep["rank_histogram"] == {"5": 1, "9": 4}
        # the dropped point never reaches max rank -> counted as unseen
        assert rep["n_points_seen"] == 2

    def test_failed_records_are_ignored(self, gate, m11):
        recs = _grid(m11) + [_record(m11, 9, 101, formal_success=False)]
        rep = gate.audit_records(recs, m11, 2, 2)
        assert rep["n_records"] == 5 and rep["n_formal_success"] == 4

    def test_incomplete_prime_coverage_detected(self, gate, m11):
        rep = gate.audit_records(_grid(m11)[:-1], m11, 2, 2)
        assert len(rep["points_incomplete_primes"]) == 1


class TestGateVerdict:
    def _clean(self, gate, m11, **kw):
        return gate.audit_records(_grid(m11, **kw), m11, 2, 2)

    def test_capacity_failure_with_clean_records_proceeds(self, gate, m11):
        verdict, reasons = gate.gate_verdict(self._clean(gate, m11), CAPACITY_FAIL)
        assert verdict == "PROCEED"
        assert reasons == ["capacity-bound failure"]

    def test_reconstruction_success_is_done_even_with_gaps(self, gate, m11):
        rep = gate.audit_records(_grid(m11)[:-1], m11, 2, 2)
        verdict, _ = gate.gate_verdict(rep, SUCCESS)
        assert verdict == "DONE"

    def test_non_capacity_failure_stops(self, gate, m11):
        verdict, reasons = gate.gate_verdict(self._clean(gate, m11), OTHER_FAIL)
        assert verdict == "STOP"
        assert any("non-capacity" in r for r in reasons)

    def test_support_instability_stops_despite_capacity_failure(self, gate, m11):
        recs = _grid(m11)
        recs[0] = _record(m11, 1, 101, support=((0, 0, -1),))
        rep = gate.audit_records(recs, m11, 2, 2)
        verdict, reasons = gate.gate_verdict(rep, CAPACITY_FAIL)
        assert verdict == "STOP"
        assert any("unstable support" in r for r in reasons)

    def test_missing_artifact_section_stops(self, gate, m11):
        verdict, reasons = gate.gate_verdict(self._clean(gate, m11), {"no": "reconstruction"})
        assert verdict == "STOP"
        assert any("lacks a reconstruction section" in r for r in reasons)

    def test_no_artifact_records_only_proceeds(self, gate, m11):
        verdict, _ = gate.gate_verdict(self._clean(gate, m11), None)
        assert verdict == "PROCEED"


class TestEndToEnd:
    def _write_cache(self, m11, path, records):
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(json.dumps({"fingerprint": {"chamber_ep": "-3/5"}}) + "\n")
            for r in records:
                key = m11.sample_prime_key(r.sample, r.prime)
                fh.write(json.dumps(m11._record_to_json(key, r)) + "\n")

    def test_main_proceed_exit_zero_and_header_flag(self, gate, m11, tmp_path, capsys):
        cache = tmp_path / "cache.jsonl"
        self._write_cache(m11, cache, _grid(m11))
        art = tmp_path / "art.json"
        art.write_text(json.dumps(CAPACITY_FAIL), encoding="utf-8")
        out = tmp_path / "gate.json"
        code = gate.main(
            [
                "--cache",
                str(cache),
                "--artifact",
                str(art),
                "--expected-points",
                "2",
                "--expected-primes",
                "2",
                "--json-out",
                str(out),
            ]
        )
        assert code == 0
        rep = json.loads(out.read_text(encoding="utf-8"))
        assert rep["verdict"] == "PROCEED"
        assert rep["cache_has_v2_header"] is True
        assert rep["n_records"] == 4

    def test_main_stop_exit_two(self, gate, m11, tmp_path, capsys):
        cache = tmp_path / "cache.jsonl"
        self._write_cache(m11, cache, _grid(m11)[:-1])
        code = gate.main(
            ["--cache", str(cache), "--expected-points", "2", "--expected-primes", "2"]
        )
        assert code == 2
