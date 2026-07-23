"""Tests for External Int2 Method.10 (chamber-policy controlled LF-feasibility solve)."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

import sympy as sp  # noqa: E402

from parametric_ibp_lf_reducer import surface  # noqa: E402


def _load_script(name: str):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def m10():
    return _load_script("run_external_int2_method10")


class FakeRow:
    """Duck-typed stand-in: Method.10 helpers only use ``dedup_key()`` and ``kind``."""

    def __init__(self, key, kind="tangent_ibp"):
        self._key = key
        self.kind = kind

    def dedup_key(self):
        return self._key


class TestDiffByDedupKey:
    def test_disjoint_and_overlap(self, m10):
        base = [FakeRow(("a",)), FakeRow(("b",))]
        other = [FakeRow(("b",)), FakeRow(("c",)), FakeRow(("d",))]
        only_other, only_base = m10.diff_by_dedup_key(base, other)
        assert [r.dedup_key() for r in only_other] == [("c",), ("d",)]
        assert [r.dedup_key() for r in only_base] == [("a",)]

    def test_identical_sets_give_empty_diffs(self, m10):
        rows = [FakeRow((i,)) for i in range(4)]
        only_other, only_base = m10.diff_by_dedup_key(rows, list(reversed(rows)))
        assert only_other == [] and only_base == []

    def test_duplicate_keys_in_other_reported_once(self, m10):
        base = [FakeRow(("a",))]
        other = [FakeRow(("x",)), FakeRow(("x",))]
        only_other, _ = m10.diff_by_dedup_key(base, other)
        assert len(only_other) == 1


class TestOverallVerdict:
    def test_empty_is_screening_only(self, m10):
        assert m10.overall_verdict([]) == m10.VERDICT_SCREENING_ONLY

    def test_all_feasible(self, m10):
        st = [m10.T2.STATUS_FEASIBLE] * 3
        assert m10.overall_verdict(st) == m10.VERDICT_FEASIBLE

    def test_all_obstructed(self, m10):
        st = [m10.T2.STATUS_OBSTRUCTED]
        assert m10.overall_verdict(st) == m10.VERDICT_OBSTRUCTED

    def test_mixed(self, m10):
        st = [m10.T2.STATUS_FEASIBLE, m10.T2.STATUS_OBSTRUCTED]
        assert m10.overall_verdict(st) == m10.VERDICT_MIXED

    def test_bad_specialization_is_not_feasible(self, m10):
        st = [m10.T2.STATUS_FEASIBLE, "BadSpecialization"]
        assert m10.overall_verdict(st) == m10.VERDICT_MIXED


class TestCLIGates:
    def test_solves_require_allow_heavy(self, m10):
        with pytest.raises(SystemExit) as exc:
            m10.main(["--max-solves", "1"])
        assert exc.value.code == 2

    def test_negative_max_solves_rejected(self, m10):
        with pytest.raises(SystemExit) as exc:
            m10.main(["--allow-heavy", "--max-solves", "-1"])
        assert exc.value.code == 2

    def test_bad_primes_rejected(self, m10):
        with pytest.raises(SystemExit) as exc:
            m10.main(["--allow-heavy", "--max-solves", "1", "--primes", "2"])
        assert exc.value.code == 2


class TestDiagnosticIsolation:
    def test_module_load_leaves_production_sign_policy(self, m10):
        """Importing the runner must not swap ``surface.regulated_sign``."""
        assert surface.regulated_sign.__name__ == "regulated_sign"

    def test_chamber_swap_is_scoped(self, m10):
        chamber = m10.M9.chamber_sign_factory(m10.DEFAULT_EP)
        orig = surface.regulated_sign
        with m10.M9.sign_policy(chamber):
            assert surface.regulated_sign is chamber
        assert surface.regulated_sign is orig


class TestFamilyConventionAudit:
    """Phase A: the parsed family at target label zero IS the source integrand (symbolically)."""

    def test_target_label_zero_reconstructs_source_integrand(self, m10):
        import re

        import sympy as sp

        family = m10.T2.parse_family_text(
            (REPO_ROOT / "examples" / "external_int2_dimensionless_input.wl.txt").read_text(
                encoding="utf-8"
            )
        )
        # Reference from the checked-in transcription of the source notebook (NOT re-typed here):
        ref_text = (REPO_ROOT / "examples" / "external_int2_source_reference.wl.txt").read_text(
            encoding="utf-8"
        )
        match = re.search(r"ExternalInt2PureIntegrand\s*=\s*(.*?);", ref_text, re.S)
        assert match, "reference file lost the ExternalInt2PureIntegrand block"
        reference = sp.sympify(match.group(1).replace("^", "**"))

        target = tuple([0] * (family.nvars + family.npolys))
        var_exps, poly_exps = surface._label_exps_symbolic(family, target)
        xs = [sp.Symbol(str(v)) for v in family.variables]
        integrand = sp.Integer(1)
        for i in range(family.nvars):
            integrand *= xs[i] ** var_exps[i]
        for name, pexp in zip(family.polynomials, poly_exps):
            integrand *= family.polynomials[name].to_sympy(family.variables) ** pexp
        assert sp.simplify(integrand / reference) == 1


ARTIFACT = REPO_ROOT / "validation" / "external_int2_method10.json"


@pytest.mark.skipif(not ARTIFACT.exists(), reason="recorded Method.10 artifact absent")
class TestRecordedArtifact:
    def test_schema_and_gates(self):
        data = json.loads(ARTIFACT.read_text(encoding="utf-8"))
        assert data["script"] == "scripts/run_external_int2_method10.py"
        assert data["level"] == 0
        rows = data["rows"]
        # Chamber acceptance must strictly extend production acceptance (limit_only == 0):
        assert rows["n_limit_only"] == 0
        assert rows["n_chamber"] == rows["n_limit"] + rows["n_chamber_only"]
        scr = data["screening"]
        assert scr["rerun_justified"] == (scr["n_breaks_total"] > 0)
        verdict = data["solve"]["verdict"]
        statuses = [p["status"] for p in data["solve"]["points"]]
        if verdict == "Feasible(modular)":
            assert statuses and all(s == "Feasible" for s in statuses)
        # scope note must keep the not-an-LF-basis disclaimer:
        assert "NOT an LF basis" in data["scope_note"]


class TestChamberSignExact:
    """Exact sign-policy contract: ``score = -3*ep - 1`` separates the two policies."""

    @pytest.fixture
    def chamber(self):
        m9 = _load_script("run_external_int2_method9")
        return m9.chamber_sign_factory(sp.Rational(-3, 5))

    def test_limit_policy_score_is_neg(self):
        ep = sp.Symbol("ep")
        assert surface.regulated_sign(-3 * ep - 1, ("ep",), "minus") == "neg"

    def test_chamber_policy_score_is_pos(self, chamber):
        ep = sp.Symbol("ep")
        # -3*(-3/5) - 1 = 4/5 > 0: the chamber accepts what the limit rejects.
        assert chamber(-3 * ep - 1, ("ep",), "minus") == "pos"

    def test_leftover_external_symbol_is_lowercase_unknown(self, chamber):
        ep = sp.Symbol("ep")
        assert chamber(sp.Symbol("r") * ep + 1, ("ep",), "minus") == "unknown"

    def test_exact_zero(self, chamber):
        ep = sp.Symbol("ep")
        assert chamber(5 * ep + 3, ("ep",), "minus") == "zero"
