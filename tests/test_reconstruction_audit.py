"""Tests for External Int2 Method.11c offline reconstruction audit (scripts/audit_...py)."""

import importlib.util
import json
import sys
from fractions import Fraction
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "audit_external_int2_reconstruction.py"
VALUE_TABLE = REPO_ROOT / "validation" / "external_int2_method11_value_table.json"
AUDIT = REPO_ROOT / "validation" / "external_int2_method11_reconstruction_audit.json"


def _load_script(stem: str):
    if stem in sys.modules:
        return sys.modules[stem]
    spec = importlib.util.spec_from_file_location(stem, REPO_ROOT / "scripts" / f"{stem}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[stem] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def audit_mod():
    return _load_script("audit_external_int2_reconstruction")


@pytest.fixture(scope="module")
def loaded(audit_mod):
    header, records, raw_rows = audit_mod.load_records(audit_mod.RECORDS_PATH)
    from parametric_ibp_lf_reducer.reconstruction import collect_value_table

    labels, table, samples, n_skipped = collect_value_table(records)
    assert n_skipped == 0
    return header, records, raw_rows, labels, table, samples


@pytest.fixture(scope="module")
def value_doc():
    return json.loads(VALUE_TABLE.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def audit_doc():
    return json.loads(AUDIT.read_text(encoding="utf-8"))


class TestValueTable:
    def test_verification_ok(self, value_doc):
        v = value_doc["verification"]
        assert v["ok"], v["problems"]
        assert v["n_samples"] == 38
        assert v["n_records"] == 152
        assert v["selected_rank"] == 26984
        assert v["primes"] == [2147483579, 2147483587, 2147483629, 2147483647]

    def test_cache_table_roundtrip(self, audit_mod, loaded, value_doc):
        """On-disk table values match a fresh exact rebuild from the cache."""
        _, _, _, labels, table, samples = loaded
        by_canon = {e["canonical"]: e for e in value_doc["samples"]}
        assert len(by_canon) == len(samples)
        for key, smp in samples.items():
            canon = audit_mod._canon(smp)
            entry = by_canon[canon]
            for lab in labels:
                got = Fraction(entry["coefficients"][audit_mod._lab_str(lab)]["value"])
                assert got == table[lab][key]

    def test_crt_residue_consistency(self, value_doc):
        """CRT value reduces to the recorded residue at every one of the 4 primes."""
        for entry in value_doc["samples"]:
            for lab, cf in entry["coefficients"].items():
                v = Fraction(cf["value"])
                for p_str, res in cf["residues"].items():
                    p = int(p_str)
                    assert (v.numerator * pow(v.denominator, -1, p)) % p == res, (
                        entry["canonical"],
                        lab,
                    )

    def test_special_zero_classification(self, value_doc):
        v = value_doc["verification"]
        assert list(v["special_zero"]) == ["ep=6,r=57/11"]
        assert v["special_zero"]["ep=6,r=57/11"]["missing_labels"] == ["(0, 0, 1, -1, 0, 0, -1)"]
        entry = next(e for e in value_doc["samples"] if e["canonical"] == "ep=6,r=57/11")
        cf = entry["coefficients"]["(0, 0, 1, -1, 0, 0, -1)"]
        assert cf["zero_filled"] is True
        assert Fraction(cf["value"]) == 0
        assert all(res == 0 for res in cf["residues"].values())


class TestDegreeDiagnostics:
    def test_monomial_and_unknown_counts(self, audit_mod):
        assert len(audit_mod._monoms_total(6)) == 28
        assert len(audit_mod._monoms_total(7)) == 36
        assert len(audit_mod._monoms_tensor(2, 3)) == 12

    def test_dense_cell_statuses_on_real_data(self, audit_mod, loaded):
        _, _, _, labels, table, samples = loaded
        fit, hold = audit_mod.split_fit_hold(samples)
        assert (len(fit), len(hold)) == (36, 2)
        vals = audit_mod._values_by_point(table[labels[0]], samples)
        primes = list(audit_mod.EXPECTED_PRIMES)
        e00, _ = audit_mod.cell_diagnostics(vals, fit, (0, 0), primes)
        assert e00["status"].startswith("inconsistent (dim 0)")
        e66, _ = audit_mod.cell_diagnostics(vals, fit, (6, 6), primes)
        assert e66["unknowns"] == 56
        assert "needs 55 fit points" in e66["status"]

    def test_artifact_crosscheck_matches(self, audit_doc):
        xc = audit_doc["phase_b_degree_diagnostics"]["artifact_crosscheck"]
        assert xc["available"] and xc["matches"], xc

    def test_row_reduce_mod_smoke(self, audit_mod):
        p = 101
        rows = [[1, 2, 3], [2, 4, 6], [0, 1, 1]]
        rank, piv, _ = audit_mod._row_reduce_mod(rows, p)
        assert rank == 2 and piv == [0, 1]
        dim, vec = audit_mod._nullspace_dim_mod(rows, p)
        assert dim == 1
        for row in rows:
            assert sum(a * b for a, b in zip(row, vec)) % p == 0


class TestSchedule:
    def test_deterministic_schedule(self, audit_mod, loaded, audit_doc):
        """Phase D is reproducible and respects the exclusion rules."""
        _, _, _, labels, table, samples = loaded
        fit, _ = audit_mod.split_fit_hold(samples)
        b = audit_doc["phase_b_degree_diagnostics"]
        c = audit_doc["phase_c_structured_experiments"]
        d1 = audit_mod.phase_d(b, c, samples, False, fit)
        d2 = audit_mod.phase_d(b, c, samples, False, fit)
        assert d1 == d2
        sched = d1["recommended_new_samples"]
        assert len(sched) == d1["schedule_target"]["additional_samples"]
        pts = [(Fraction(s["ep"]), Fraction(s["r"])) for s in sched]
        assert len(set(pts)) == len(pts)
        existing = {(s["ep"], s["r"]) for s in samples.values()}
        for ep, r in pts:
            assert ep not in (3, 6) and r != Fraction(57, 11)
            assert (ep, r) not in existing

    def test_costs_are_positive_and_sorted(self, audit_doc):
        costs = audit_doc["phase_d_schedule"]["viable_ansatz_costs"]
        adds = [e["additional_samples"] for e in costs]
        assert adds == sorted(adds)
        assert all(a >= 1 for a in adds)


class TestNoRrefOrReducer:
    def test_static_no_rref_or_reducer_imports(self):
        imports = [
            line.strip()
            for line in SCRIPT.read_text(encoding="utf-8").splitlines()
            if line.strip().startswith(("import ", "from "))
        ]
        for line in imports:
            assert "sparse_rref" not in line, line
            assert ".reducer" not in line and "import reducer" not in line, line

    def test_audit_flags_no_rref(self, audit_doc):
        assert audit_doc["no_rref_or_reducer_solve"] is True

    def test_negative_results_recorded_honestly(self, audit_doc):
        rec = audit_doc["reconstructed_coefficients"]
        c = audit_doc["phase_c_structured_experiments"]
        for lab, got in rec.items():
            if got is None:
                assert lab in c["C_sparse_tensor_supports"]
        assert "statement" in c["A_univariate_ep_on_fixed_r"]
        assert "statement" in c["B_univariate_r_on_fixed_ep"]
