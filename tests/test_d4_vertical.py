"""Pass D4.1 — targeted D4 vertical validation (canonical 11.3).

This is NOT full end-to-end reconstruction. Its job is to *separate two questions* that the
vertical audit surfaced:

1. Is the expected relation ``J[T] = sum_i C_i J[M_i]`` actually **in the row span** of the
   generated row system? (a modular certificate, ``test_d4_expected_relation_is_in_row_span_mod_p``)
2. Does ``reduce_family_once`` currently *select* that LF basis? (diagnostic,
   ``test_d4_reduce_family_once_current_config`` — Success not required yet)

If (1) holds but (2) does not pick M1..M5, that is a basis-selection issue, not a row-generation
failure. Nothing here fabricates ``Success``; a failed certificate is reported as the real blocker.
"""

from __future__ import annotations

import os
from fractions import Fraction

import pytest
import sympy as sp

from parametric_ibp_lf_reducer import (
    ALL_FAILURE_REASONS,
    FAILURE_INTERPOLATION_FAILED,
    FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE,
    STATUS_SUCCESS,
    ReducerConfig,
    enumerate_box,
    generate_algebraic_rows,
    generate_coordinate_ibp_rows,
    generate_tangent_fields,
    generate_tangent_ibp_rows,
    is_locally_finite,
    modular_normal_form,
    parse_family_text,
    reduce_family_once,
    verify_reduction_relation_mod_p,
)

from conftest import load_example, load_validation

# --- D4 fixtures -----------------------------------------------------------------------------
# label = (x1, x2, x3, x4, G0, G1, G2)
T = (0, 0, 0, 0, 0, 0, 0)
M1 = (0, 1, 1, 0, -2, -1, 0)  # x2*x3/(G0^2*G1)
M2 = (1, 1, 0, 0, -2, -1, 0)  # x1*x2/(G0^2*G1)
M3 = (0, 1, 1, 0, -3, -1, 0)  # x2*x3/(G0^3*G1)
M4 = (1, 1, 0, 0, -3, -1, 0)  # x1*x2/(G0^3*G1)
M5 = (0, 1, 1, 0, -4, -1, 0)  # x2*x3/(G0^4*G1)
MASTERS = [M1, M2, M3, M4, M5]

PRIME = 2_147_483_647
SAMPLE = {"ep": Fraction(2, 3), "r": Fraction(3)}

# Safe (fractional) certificate points, one per large prime — each is verified in-span for the
# REFERENCE relation. NOTE: (3/2, 4) is rank-DEFICIENT (rank 2011 vs generic 2041); the reference
# relation still reduces there, but reducer-output certificates need rank-GENERIC points.
CERT_POINTS = [
    (Fraction(2, 3), Fraction(3), 2_147_483_647),
    (Fraction(3, 2), Fraction(4), 2_147_483_629),
    (Fraction(5, 4), Fraction(6), 2_147_483_587),
]

# Certificate points for the reducer's OWN output (D4.4): probe-verified rank-generic
# (rank 2041, support {M1,M2,M3}) and NOT part of the reduction sample grid below, so they are
# genuinely independent off-sample checks of the reconstructed coefficient functions.
REDUCER_CERT_POINTS = [
    (Fraction(2, 3), Fraction(3), 2_147_483_647),
    (Fraction(5, 4), Fraction(6), 2_147_483_587),
    (Fraction(7, 3), Fraction(9, 2), 2_147_483_629),
]

# D4 label box: x1,x2,x3 in 0..1, x4 fixed 0; G0 in -4..0, G1 in -1..0, G2 fixed 0.
N_RANGE = [(0, 1), (0, 1), (0, 1), (0, 0)]
M_RANGE = [(-4, 0), (-1, 0), (0, 0)]


def d4_family():
    return parse_family_text(load_example("d4_explicit_family.wl.txt"))


def _expected_coefficients():
    """C1..C5 as SymPy expressions from the validation fixture."""
    data = load_validation("expected_d4_coefficients.json")["coefficients"]
    return [sp.sympify(data[f"C{i}"]) for i in range(1, 6)]


def _build_d4_rows(family):
    """Generate the full D4 row system (algebraic + coordinate-IBP deg2 + tangent) for the box."""
    seeds = list(enumerate_box(family.nvars, family.npolys, N_RANGE, M_RANGE))
    alg = generate_algebraic_rows(family, seeds)
    coord = generate_coordinate_ibp_rows(family, seeds, 2)
    fields = generate_tangent_fields(family, [(1, 1), (2, 2)])
    tang = generate_tangent_ibp_rows(family, seeds, fields)
    rows = [*alg.rows, *coord.rows, *tang.rows]
    counts = {
        "seeds": len(seeds),
        "algebraic": len(alg.rows),
        "coordinate_ibp": len(coord.rows),
        "coordinate_rejected": len(coord.rejected),
        "tangent_fields": len(fields),
        "tangent_ibp": len(tang.rows),
        "tangent_rejected": len(tang.rejected),
        "total_rows": len(rows),
    }
    return rows, counts


@pytest.fixture(scope="module")
def d4_rows():
    """Build the (expensive) D4 row system once and share it across the heavy tests."""
    family = d4_family()
    rows, counts = _build_d4_rows(family)
    return family, rows, counts


# --- Test 1: LF verdicts (fast) --------------------------------------------------------------
def test_d4_expected_labels_are_lf_and_target_is_not_lf():
    family = d4_family()
    assert is_locally_finite(family, T) is False
    for i, m in enumerate(MASTERS, start=1):
        assert is_locally_finite(family, m) is True, f"M{i} {m} expected locally finite"


# --- Test 2: modular row-span certificate (core), over multiple points ------------------------
@pytest.mark.integration
@pytest.mark.parametrize("ep, r, prime", CERT_POINTS)
def test_d4_expected_relation_is_in_row_span_mod_p(d4_rows, ep, r, prime):
    """The reference relation J[T] = sum C_i J[M_i] is certified in the row span at several
    points (guards against a one-point accident). Uses the generic certificate helper (D4.4)."""
    family, rows, counts = d4_rows
    assert counts["total_rows"] > 0

    ref_terms = dict(zip(MASTERS, _expected_coefficients()))
    cert = verify_reduction_relation_mod_p(family, rows, T, ref_terms, {"ep": ep, "r": r}, prime)
    assert cert.in_span, (
        f"expected D4 relation is NOT in the row span at ep={ep}, r={r}, p={prime}; "
        f"residual columns = {sorted(cert.residual)} (row counts: {counts})"
    )
    assert cert.rank > 0 and 0 < cert.nrows <= counts["total_rows"]


# --- Test 2b: reducer's own generic normal form is a locally-finite subset of M1..M5 ----------
@pytest.mark.integration
@pytest.mark.parametrize(
    "ep, r, prime",
    [(Fraction(3), Fraction(2), 2_147_483_647), (Fraction(4), Fraction(5), 2_147_483_629)],
)
def test_d4_target_reduces_to_lf_subset_of_masters(d4_rows, ep, r, prime):
    """At generic points the target reduces (with preferred_masters=M1..M5) to a subset of those
    masters, all locally finite — i.e. ``preferred_masters`` already frees only LF masters."""
    family, rows, _ = d4_rows
    nf = modular_normal_form(
        family, rows, T, {"ep": ep, "r": r}, prime, preferred_masters=MASTERS
    )
    assert nf.status == "Reduced"
    labels = set(nf.terms)
    assert labels, "target produced an empty normal form"
    assert labels <= set(MASTERS), f"free terms outside M1..M5: {sorted(labels - set(MASTERS))}"
    assert all(is_locally_finite(family, lab) is True for lab in labels)


# --- Test 3: row counts (informational) ------------------------------------------------------
@pytest.mark.integration
def test_d4_row_counts_are_reported(d4_rows, capsys):
    _, _, counts = d4_rows
    with capsys.disabled():
        print(f"\nD4 row counts @ box {N_RANGE}|{M_RANGE}, deg2, tangent[(1,1),(2,2)]: {counts}")
    assert counts["algebraic"] > 0
    assert counts["coordinate_ibp"] > 0
    assert counts["total_rows"] == (
        counts["algebraic"] + counts["coordinate_ibp"] + counts["tangent_ibp"]
    )


# --- Test 4: current reducer outcome (diagnostic; Success not required) -----------------------
def test_d4_reduce_family_once_current_config(capsys):
    family = d4_family()
    primes = [PRIME, 2_147_483_629]
    samples = [{"ep": Fraction(a), "r": Fraction(b)} for a in (2, 3, 5) for b in (2, 3, 5)]
    cfg = ReducerConfig(
        primes=primes, samples=samples, labels=[T, *MASTERS], max_ibp_degree=1
    )
    res = reduce_family_once(family, T, cfg)
    with capsys.disabled():
        print(
            f"\nD4 reduce_family_once (labels=T+M1..M5, deg1): status={res.status} "
            f"error={res.error} extra={res.diagnostics.extra.get('row_diagnostics')}"
        )
    # Honest typed outcome; we do NOT require Success yet (basis-selection is a later pass).
    assert res.status in ({STATUS_SUCCESS} | ALL_FAILURE_REASONS)


# --- Test 5: full-config reducer attempt with preferred masters (heavy; opt-in) ---------------
# History:
# - D4.2 (2026-07-07, 16 samples x 3 primes, ~350 s): InterpolationFailed — rank-deficient
#   samples (e.g. ep=2,r=3 -> rank 1995 vs generic 2041) shrank the normal-form support and the
#   union-support 0-fill poisoned interpolation.
# - D4.3 adds max-rank record selection before reconstruction, so that poisoning path is gone:
#   rank-deficient records are skipped + counted in diagnostics.extra["record_selection"].
#   The grid here is enlarged to 6x6 so the surviving max-rank samples can pin down the
#   coefficient degrees. Runs only when RUN_D4_FULL is set (else ~10-15 min).
# - Recorded D4.3 outcome (2026-07-07, 795 s, 6x6 INTEGER grid): formal Success — terms
#   {M1,M2,M3} all LF, rank_histogram {1995: 18, 2041: 90}, 90/108 records selected. BUT the
#   D4.4 row-span certificate then showed those interpolated coefficients are WRONG off-grid:
#   a 6x6 product lattice is degenerate for the dense degree search (prod(ep-k), k=2..7 has
#   degree 6 = max_deg and vanishes on the whole lattice, holdout included), so a wrong
#   candidate can pass on-lattice validation. True values at ep=2/3, r=3 are (-5, 5, -2) for
#   (M1, M2, M3) — the interpolants disagreed there.
# - D4.4 therefore samples a SCATTERED non-lattice grid (no product structure, no low-degree
#   curve through the points), shares the run via a module fixture, and adds row-span
#   certificates for the reducer's own output + equivalence with the reference relation at
#   independent off-sample rank-generic points.
@pytest.fixture(scope="module")
def d4_full_result():
    """Run the heavy full-config D4 reduction ONCE (opt-in) and share it across D4.4 tests."""
    if not os.environ.get("RUN_D4_FULL"):
        pytest.skip("heavy ~10-15min full-config D4 reduce; set RUN_D4_FULL=1 to run")
    family = d4_family()
    primes = [PRIME, 2_147_483_629, 2_147_483_587]
    # 35 scattered rational points (ep distinct with denominator 7, r jumps with denominator 6 —
    # no product-lattice degeneracy, no low-degree curve) + one known rank-deficient point (2,3)
    # to keep the D4.3 rank filter exercised in vivo.
    samples = [
        {"ep": Fraction(14 + k, 7), "r": Fraction(12 + ((11 * k + 5) % 36), 6)}
        for k in range(35)
    ]
    samples.append({"ep": Fraction(2), "r": Fraction(3)})
    cfg = ReducerConfig(
        primes=primes,
        samples=samples,
        label_box=(N_RANGE, M_RANGE),
        max_ibp_degree=2,
        tangent_degree_blocks=[(1, 1), (2, 2)],
        min_valid_records=16,
        preferred_masters=MASTERS,
        # D4.5: the reducer certifies its own relation at these probe-verified rank-generic
        # off-sample points; Success now REQUIRES the certificate to pass (default gate).
        certificate_points=[{"ep": ep, "r": r} for ep, r, _ in REDUCER_CERT_POINTS],
    )
    return family, cfg, reduce_family_once(family, T, cfg)


@pytest.mark.integration
def test_d4_reduce_family_once_full_config_preferred_masters(d4_full_result, capsys):
    family, cfg, res = d4_full_result
    ex = res.diagnostics.extra
    sel = ex.get("record_selection", {})
    with capsys.disabled():
        print(
            f"\nD4 full-config: status={res.status} all_lf={res.all_locally_finite} "
            f"terms={[t.label for t in res.terms]} "
            f"n_rows={ex.get('n_rows')} n_records={ex.get('n_records')} "
            f"n_reduced={ex.get('n_reduced_records')} "
            f"n_selected={ex.get('n_selected_records')} "
            f"n_bad_spec={ex.get('n_bad_specializations')} "
            f"n_target_not_pivot={ex.get('n_target_not_pivot')}\n"
            f"  record_selection: rank_histogram={sel.get('rank_histogram')} "
            f"selected_rank={sel.get('selected_rank')} "
            f"n_rank_filtered={sel.get('n_rank_filtered_records')} "
            f"support={sel.get('support_after_rank_filter')}\n"
            f"  messages={res.diagnostics.messages}"
        )

    assert res.status in ({STATUS_SUCCESS} | ALL_FAILURE_REASONS)  # honest typed outcome

    # D4.3: the rank filter must be wired in and engaged — the grid contains rank-deficient
    # points (e.g. ep=2,r=3), so mixed ranks appear and only the max-rank records may feed
    # reconstruction. This kills the old rank-poisoning interpolation symptom structurally.
    assert sel, "record_selection diagnostics missing from the reducer run"
    assert sel["n_valid_records_before_rank_filter"] == ex["n_reduced_records"]
    assert len(sel["rank_histogram"]) >= 2, "expected rank-deficient samples in this grid"
    assert sel["n_rank_filtered_records"] > 0
    assert sel["selected_rank"] == max(sel["rank_histogram"])

    if res.status == STATUS_SUCCESS:
        # If it succeeds, every returned term must be locally finite...
        assert res.all_locally_finite is True
        assert all(t.locally_finite is True for t in res.terms)
        term_labels = {t.label for t in res.terms}
        assert term_labels <= set(MASTERS)
        # ...and (D4.5) Success now implies the reducer's own row-span certificate passed:
        assert ex["certificate"]["certificate_status"] == "Passed"
        assert ex["certificate"]["n_certificate_points_failed"] == 0
        # ...and if the basis is exactly M1..M5, coefficients must match the reference symbolically.
        if term_labels == set(MASTERS):
            expected = dict(zip(MASTERS, _expected_coefficients()))
            for t in res.terms:
                got = sp.sympify(t.coefficient_text.replace("^", "**"))
                assert sp.simplify(got - expected[t.label]) == 0, f"coeff mismatch at {t.label}"
    elif res.status == FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE:
        # Would prove a protected/free-basis policy is needed: which non-LF labels remained.
        assert res.diagnostics.non_lf_terms or res.diagnostics.unknown_lf_terms
    elif res.status == FAILURE_INTERPOLATION_FAILED:
        # If reconstruction still fails, it must be for a NEW honest reason (e.g. coefficient
        # degree exceeding the sample budget), not the old rank-poisoning: the value table now
        # only sees max-rank records with a consistent support.
        assert ex.get("n_reduced_records", 0) > 0
        assert ex.get("n_target_not_pivot", 0) == 0
        assert ex.get("n_selected_records", 0) >= cfg.min_valid_records, (
            "rank filter starved reconstruction below min_valid_records — enlarge the grid"
        )


# --- Test 6 (D4.4): the reducer's OWN output is row-span certified (heavy; opt-in) -------------
@pytest.mark.integration
def test_d4_success_result_is_row_span_certified(d4_full_result, d4_rows):
    """D4 acceptance: the full-config result must be a strict-gate Success, all terms locally
    finite within {M1..M5}, and the reconstructed relation T = sum C_i * L_i must be certified
    in the generated row span at several independent rank-generic points that are NOT sample
    points of the reduction grid. (This is the check that exposed the 6x6-lattice interpolation
    failure — see the history note above.)"""
    _, _, res = d4_full_result
    family, rows, _ = d4_rows

    assert res.status == STATUS_SUCCESS
    assert res.all_locally_finite is True
    assert res.terms and {t.label for t in res.terms} <= set(MASTERS)
    assert all(t.locally_finite is True for t in res.terms)

    reducer_terms = {t.label: t.coefficient for t in res.terms}
    for ep, r, prime in REDUCER_CERT_POINTS:
        cert = verify_reduction_relation_mod_p(
            family, rows, T, reducer_terms, {"ep": ep, "r": r}, prime
        )
        assert cert.in_span, (
            f"reducer relation NOT in row span at ep={ep}, r={r}, p={prime}; "
            f"residual columns = {sorted(cert.residual)}"
        )


# --- Test 7 (D4.4): equivalence of the reducer result and the reference relation --------------
@pytest.mark.integration
def test_d4_reducer_relation_equivalent_to_reference(d4_full_result, d4_rows):
    """Both relations — reference ``T - sum C_ref_i M_i`` and reducer ``T - sum C_red_i L_i`` —
    vanish modulo the SAME generated row span at the SAME sampled points. Hence they are
    equivalent modulo the generated IBP/algebraic rows: their difference
    ``sum C_ref_i M_i - sum C_red_i L_i`` is itself in the row span. No coefficient-by-coefficient
    comparison is required (the bases differ: M4, M5 are reducible in this row system)."""
    _, _, res = d4_full_result
    family, rows, _ = d4_rows
    assert res.status == STATUS_SUCCESS

    reducer_terms = {t.label: t.coefficient for t in res.terms}
    reference_terms = dict(zip(MASTERS, _expected_coefficients()))
    for ep, r, prime in REDUCER_CERT_POINTS:
        sample = {"ep": ep, "r": r}
        cert_ref = verify_reduction_relation_mod_p(family, rows, T, reference_terms, sample, prime)
        cert_red = verify_reduction_relation_mod_p(family, rows, T, reducer_terms, sample, prime)
        assert cert_ref.in_span, f"reference relation not in span at {sample}, p={prime}"
        assert cert_red.in_span, f"reducer relation not in span at {sample}, p={prime}"


# --- Test 8 (D4.4, diagnostic): WHY the basis is 3-term — M4, M5 are themselves reducible -----
@pytest.mark.integration
@pytest.mark.parametrize("m_label", [M4, M5], ids=["M4", "M5"])
def test_d4_m4_m5_reduce_to_smaller_basis(d4_rows, m_label):
    """Diagnostic (no Success involved): with M1..M3 preferred free, M4 and M5 each reduce to a
    combination of {M1,M2,M3} inside the same row system — which is exactly why the reducer's
    LF basis is 3-term rather than the reference 5-term one."""
    family, rows, _ = d4_rows
    ep, r, prime = CERT_POINTS[0]
    nf = modular_normal_form(
        family, rows, m_label, {"ep": ep, "r": r}, prime, preferred_masters=[M1, M2, M3]
    )
    assert nf.status == "Reduced"
    assert set(nf.terms) <= {M1, M2, M3}, (
        f"{m_label} reduced to labels outside M1..M3: {sorted(set(nf.terms) - {M1, M2, M3})}"
    )
