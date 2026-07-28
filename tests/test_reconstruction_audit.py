"""Tests for the External Int2 Method.11 reconstruction audits (Method.11c semantics).

Three artifacts / tools are covered:

* ``validation/external_int2_method11_value_table.json`` — the historical Phase-A value
  table written in the zero-fill era (schema ``m11c-value-table/v1``, 38 samples).  Its
  exact CRT values for the 37 full-support samples remain authoritative; its treatment of
  (ep, r) = (6, 57/11) as a zero-filled "special zero" is kept on disk as recorded
  history and asserted here *as history*.
* ``scripts/audit_external_int2_reconstruction.py`` — the Phase-A audit script, whose
  ``build_value_table`` now follows Method.11c semantics (schema ``m11c-value-table/v2``):
  the support-deviating sample is excluded from the value table (never zero-filled) and
  documented under ``verification["excluded_from_table"]``.
* ``validation/external_int2_method11_reconstruction_audit.json`` — the Method.11c audit
  artifact (schema ``external-int2-method11-reconstruction-audit/v1``) written by
  ``scripts/audit_external_int2_method11_reconstruction.py``: support classification,
  the four reconstructed generic rational coefficients (validated on hold-out points),
  and an exact modular cross-check of every cached record.  Verdict:
  ``basis_specialization_not_special_zero``.

The old audit's Phase B–D outputs (degree-diagnostics cross-check, ansatz experiments,
Stage-5 schedule) were superseded by the Method.11c fix and live in git history only;
their pure helpers (dense-cell diagnostics, local row reduction) are still tested here.
"""

from __future__ import annotations

import ast
import importlib.util
import json
import sys
from fractions import Fraction
from pathlib import Path

import pytest
import sympy as sp

REPO_ROOT = Path(__file__).resolve().parents[1]
OLD_SCRIPT = REPO_ROOT / "scripts" / "audit_external_int2_reconstruction.py"
NEW_SCRIPT = REPO_ROOT / "scripts" / "audit_external_int2_method11_reconstruction.py"
VALUE_TABLE = REPO_ROOT / "validation" / "external_int2_method11_value_table.json"
AUDIT_DOC = REPO_ROOT / "validation" / "external_int2_method11_reconstruction_audit.json"

DEVIATING = "ep=6,r=57/11"
MISSING_LABEL = "(0, 0, 1, -1, 0, 0, -1)"
EXPECTED_PRIMES = (2147483579, 2147483587, 2147483629, 2147483647)
GENERIC_RANK = 26984
N_RECORDS = 152
N_CACHED_SAMPLES = 38
N_TABLE_SAMPLES = 37

EXPECTED_COEFFS = {
    "(-1, 0, 0, -1, -1, 0, 0)": "-(2*ep*r - ep + 1)/(6*r*(3*ep + 1))",
    "(-1, 0, 0, 0, -1, 0, -1)": "-(ep - 1)**2*(r + 1)/(6*ep*r*(3*ep + 1))",
    "(0, 0, 0, -1, 0, 0, -1)": "(ep - 1)*(4*ep**2*r + ep**2 - 2*ep*r - 1)/(6*ep**2*r*(3*ep + 1))",
    "(0, 0, 1, -1, 0, 0, -1)": "(ep - 1)*(2*ep - 1)*(2*ep*r + ep + 1)/(6*ep**2*r*(3*ep + 1))",
}

EP, R = sp.symbols("ep r")


#: The 152-record modular cache these audits replay. It is gitignored
#: (``.gitignore``: ``validation/*_records.jsonl``) because it is regenerated
#: evidence, not source, so a clean checkout does not have it. Without it these
#: tests cannot run at all -- and they must not be silently weakened either:
#: they assert the exact record count, the generic rank and the per-prime CRT
#: residues, none of which a trimmed fixture could reproduce honestly. So they
#: skip with a precise reason instead of erroring at collection.
RECORDS_CACHE = REPO_ROOT / "validation" / "external_int2_method11_records.jsonl"

pytestmark = pytest.mark.skipif(
    not RECORDS_CACHE.is_file(),
    reason=(
        f"requires the regenerated modular cache {RECORDS_CACHE.name} "
        f"({N_RECORDS} records, gitignored via 'validation/*_records.jsonl', so "
        "absent in a clean checkout). Regenerate with "
        "'python scripts/run_external_int2_method11.py'. These audits replay the "
        "cache exactly (record count, generic rank, per-prime CRT residues) and "
        "cannot be run against a trimmed fixture."
    ),
)


@pytest.fixture(scope="module")
def audit_mod():
    src = str(REPO_ROOT / "src")
    if src not in sys.path:
        sys.path.insert(0, src)
    spec = importlib.util.spec_from_file_location("_m11_recon_audit", OLD_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def loaded(audit_mod):
    header, records, raw_rows = audit_mod.load_records(audit_mod.RECORDS_PATH)
    return header, records, raw_rows


@pytest.fixture(scope="module")
def fresh(audit_mod, loaded):
    _header, records, raw_rows = loaded
    return audit_mod.build_value_table(records, raw_rows)


@pytest.fixture(scope="module")
def value_doc():
    return json.loads(VALUE_TABLE.read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def audit_doc():
    return json.loads(AUDIT_DOC.read_text(encoding="utf-8"))


def _key_for(audit_mod, samples, canon):
    for key, sample in samples.items():
        if audit_mod._canon(sample) == canon:
            return key
    raise AssertionError(f"sample {canon} not found in value table")


# --------------------------------------------------------------------------------------
# Historical Phase-A artifact (zero-fill era) — asserted as recorded history.
# --------------------------------------------------------------------------------------
class TestHistoricalValueTable:
    def test_doc_shape_and_verification(self, value_doc):
        assert value_doc["schema"] == "m11c-value-table/v1"
        assert len(value_doc["labels"]) == 4
        assert set(value_doc["labels"]) == set(EXPECTED_COEFFS)
        assert len(value_doc["samples"]) == N_CACHED_SAMPLES
        ver = value_doc["verification"]
        assert ver["n_samples"] == N_CACHED_SAMPLES
        assert ver["n_records"] == N_RECORDS
        assert tuple(ver["primes"]) == EXPECTED_PRIMES
        assert ver["selected_rank"] == GENERIC_RANK
        assert ver["special_zero"] == {DEVIATING: {"missing_labels": [MISSING_LABEL]}}

    def test_historical_zero_fill_entry(self, value_doc):
        # The zero-fill era recorded (ep, r) = (6, 57/11) as a "special zero" with the
        # missing label zero-filled.  Method.11c reclassified this point as a basis
        # specialization (see the audit artifact); the disk doc is kept as history.
        entry = next(s for s in value_doc["samples"] if s["canonical"] == DEVIATING)
        assert entry["classification"] == "special_zero"
        coeff = entry["coefficients"][MISSING_LABEL]
        assert coeff["zero_filled"] is True
        assert coeff["value"] == "0"
        assert all(v == 0 for v in coeff["residues"].values())


# --------------------------------------------------------------------------------------
# Fresh rebuild under Method.11c semantics.
# --------------------------------------------------------------------------------------
class TestValueTableRebuild:
    def test_skip_count(self, audit_mod, loaded):
        _header, records, _raw = loaded
        labels, _table, samples, n_skipped = audit_mod.collect_value_table(records)
        # One support-deviating sample x four primes is excluded, never zero-filled.
        assert n_skipped == 4
        assert len(samples) == N_TABLE_SAMPLES
        assert [audit_mod._lab_str(lab) for lab in labels] == sorted(EXPECTED_COEFFS)

    def test_verification_ok(self, fresh):
        doc, _labels, _table, samples, _classification = fresh
        ver = doc["verification"]
        assert ver["problems"] == []
        assert ver["ok"] is True
        assert doc["schema"] == "m11c-value-table/v2"
        assert ver["n_samples"] == N_TABLE_SAMPLES == len(samples)
        assert ver["n_cached_samples"] == N_CACHED_SAMPLES
        assert ver["n_records"] == N_RECORDS

    def test_excluded_from_table(self, fresh):
        doc, *_ = fresh
        excluded = doc["verification"]["excluded_from_table"]
        assert len(excluded) == 1
        (entry,) = excluded
        assert entry["canonical"] == DEVIATING
        assert entry["classification"] == "support_deviation_pending_validation"
        assert entry["missing_labels"] == [MISSING_LABEL]
        assert entry["cache_keys"] == [f"p={p};{DEVIATING}" for p in EXPECTED_PRIMES]

    def test_roundtrip_against_disk(self, fresh, value_doc):
        doc, *_ = fresh
        assert doc["labels"] == value_doc["labels"]
        fresh_by = {s["canonical"]: s for s in doc["samples"]}
        disk_by = {s["canonical"]: s for s in value_doc["samples"]}
        # The only doc-level divergence is the deviating sample, dropped by Method.11c.
        assert sorted(set(disk_by) - set(fresh_by)) == [DEVIATING]
        assert sorted(set(fresh_by) - set(disk_by)) == []
        norm = lambda e: json.loads(json.dumps(e, default=str))  # noqa: E731
        for canon, entry in fresh_by.items():
            assert norm(entry["coefficients"]) == norm(disk_by[canon]["coefficients"]), canon


# --------------------------------------------------------------------------------------
# Exact arithmetic cross-checks between the artifacts.
# --------------------------------------------------------------------------------------
class TestCrtConsistency:
    def test_table_values_match_disk_residues(self, audit_mod, fresh, value_doc):
        _doc, labels, table, samples, _cls = fresh
        canon = "ep=15/7,r=32/11"
        key = _key_for(audit_mod, samples, canon)
        disk = next(s for s in value_doc["samples"] if s["canonical"] == canon)
        for lab in labels:
            value = table[lab][key]
            assert isinstance(value, Fraction)
            residues = disk["coefficients"][audit_mod._lab_str(lab)]["residues"]
            for p_str, res in residues.items():
                p = int(p_str)
                got = (value.numerator * pow(value.denominator, -1, p)) % p
                assert got == res, (canon, lab, p)

    def test_reconstructed_coefficients_match_table(self, audit_mod, fresh, audit_doc):
        # The Method.11c generic coefficients must reproduce the exact CRT value table
        # at a generic sample: ties the reconstruction audit to the Phase-A table.
        _doc, labels, table, samples, _cls = fresh
        key = _key_for(audit_mod, samples, "ep=15/7,r=32/11")
        point = {EP: sp.Rational(15, 7), R: sp.Rational(32, 11)}
        for lab in labels:
            expr = sp.sympify(
                audit_doc["coefficients"][audit_mod._lab_str(lab)]["expression"],
                locals={"ep": EP, "r": R},
            )
            got = expr.subs(point)
            want = table[lab][key]
            assert got == sp.Rational(want.numerator, want.denominator), lab


# --------------------------------------------------------------------------------------
# Pure helpers of the Phase-A script (kept after the Phase B-D retirement).
# --------------------------------------------------------------------------------------
class TestDegreeDiagnostics:
    def test_dense_cell_statuses_on_real_data(self, audit_mod, fresh):
        _doc, labels, table, samples, _cls = fresh
        fit, hold = audit_mod.split_fit_hold(samples)
        assert (len(fit), len(hold)) == (35, 2)
        vals = audit_mod._values_by_point(table[labels[0]], samples)
        primes = list(EXPECTED_PRIMES)
        e00, _ = audit_mod.cell_diagnostics(vals, fit, (0, 0), primes)
        assert e00["status"].startswith("inconsistent (dim 0)")
        # With the poisoned point excluded, the true (2,2) rational cell now fits:
        # the reconstructed coefficient for labels[0] has bidegree (2,2).
        e22, _ = audit_mod.cell_diagnostics(vals, fit, (2, 2), primes)
        assert e22["status"] == "candidate (dim 1)"
        e66, _ = audit_mod.cell_diagnostics(vals, fit, (6, 6), primes)
        assert e66["unknowns"] == 56
        assert e66["status"] == "underdetermined (needs 55 fit points, have 35)"

    def test_row_reduce_mod_smoke(self, audit_mod):
        rows = [[1, 2, 3], [2, 4, 6], [0, 1, 1]]
        rank, piv_cols, echelon = audit_mod._row_reduce_mod(rows, 7)
        assert rank == 2
        assert piv_cols == [0, 1]
        assert echelon[2] == [0, 0, 0]
        assert rows == [[1, 2, 3], [2, 4, 6], [0, 1, 1]]  # operates on a copy


# --------------------------------------------------------------------------------------
# Method.11c reconstruction audit artifact.
# --------------------------------------------------------------------------------------
class TestMethod11cAudit:
    def test_schema_and_counts(self, audit_doc):
        assert audit_doc["schema"] == "external-int2-method11-reconstruction-audit/v1"
        assert audit_doc["cache_schema"] == "m11-record-cache/v2"
        assert audit_doc["n_records"] == N_RECORDS
        assert audit_doc["n_samples"] == N_CACHED_SAMPLES
        assert audit_doc["n_primes"] == len(EXPECTED_PRIMES)
        assert audit_doc["generic_rank"] == GENERIC_RANK
        assert audit_doc["rank_histogram"] == {str(GENERIC_RANK): N_RECORDS}

    def test_support_classification(self, audit_doc):
        assert audit_doc["n_stable_samples"] == N_TABLE_SAMPLES
        assert audit_doc["n_changed_samples"] == 1
        modal = [str(tuple(lab)) for lab in audit_doc["modal_support"]]
        assert modal == sorted(EXPECTED_COEFFS)
        assert audit_doc["modal_support_record_count"] == 148
        (changed,) = audit_doc["changed_samples"]
        assert changed["sample"] == DEVIATING
        assert changed["ranks"] == [GENERIC_RANK]
        assert changed["n_records"] == 4
        (support,) = changed["supports"]
        assert len(support) == 3
        assert MISSING_LABEL not in {str(tuple(lab)) for lab in support}

    def test_coefficients_validated(self, audit_doc):
        coeffs = audit_doc["coefficients"]
        assert sorted(coeffs) == sorted(EXPECTED_COEFFS)
        for label, info in coeffs.items():
            assert info["status"] == "validated"
            assert (info["n_points"], info["n_fit"], info["n_hold"]) == (37, 35, 2)
            got = sp.sympify(info["expression"], locals={"ep": EP, "r": R})
            want = sp.sympify(EXPECTED_COEFFS[label], locals={"ep": EP, "r": R})
            assert sp.simplify(got - want) == 0, label
            num, den = sp.fraction(sp.together(got))
            assert sp.Poly(sp.expand(num), EP, R).total_degree() == info["num_deg"]
            assert sp.Poly(sp.expand(den), EP, R).total_degree() == info["den_deg"]

    def test_exact_crosscheck_bookkeeping(self, audit_doc):
        assert audit_doc["stable_record_mismatches"] == []
        assert audit_doc["stable_records_exactly_matched"] is True
        mism = audit_doc["changed_record_mismatches"]
        # The one basis-specialized sample deviates at all 4 labels x 4 primes.
        assert len(mism) == 16
        assert {m["sample"] for m in mism} == {DEVIATING}
        assert {m["prime"] for m in mism} == set(EXPECTED_PRIMES)
        missing = [m for m in mism if str(tuple(m["label"])) == MISSING_LABEL]
        assert len(missing) == 4
        for m in missing:
            # Recorded 0 (label absent from that basis), generically nonzero.
            assert m["recorded"] == 0
            assert m["predicted"] != 0
        assert audit_doc["changed_sample_is_not_generic_zero"] is True

    def test_verdict(self, audit_doc):
        assert audit_doc["classification"] == "basis_specialization_not_special_zero"
        assert "Exclude support-changing samples from fitting" in audit_doc["recommendation"]


# --------------------------------------------------------------------------------------
# Independence guard: neither audit script may import the reducer/RREF machinery.
# --------------------------------------------------------------------------------------
class TestNoRrefOrReducer:
    FORBIDDEN = {"reducer", "rref", "sparse_rref", "sparse_rref_numba"}

    @pytest.mark.parametrize("script", [OLD_SCRIPT, NEW_SCRIPT], ids=lambda p: p.name)
    def test_static_no_rref_or_reducer_imports(self, script):
        tree = ast.parse(script.read_text(encoding="utf-8"))
        imported: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported.update(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                imported.add(node.module)
                imported.update(f"{node.module}.{alias.name}" for alias in node.names)
        bad = {
            mod
            for mod in imported
            if mod.split(".")[-1] in self.FORBIDDEN
            or "rref" in mod.split(".")[-1]
        }
        assert not bad, f"{script.name} imports forbidden modules: {sorted(bad)}"
