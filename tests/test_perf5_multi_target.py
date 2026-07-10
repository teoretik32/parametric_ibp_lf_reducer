"""Perf.5: multi-target normal-form/certificate reuse must match the single-target paths.

``collect_normal_form_records_multi`` shares ONE assemble + RREF per ``(prime, sample)``
point across several targets; ``verify_reduction_relations_mod_p`` does the same for
row-span certificates. For a single target both must be bit-identical to the singular
APIs. For several targets we pin: per-target record lists equal to a multi run with the
*same shared ranking* (the honesty caveat — masters may differ from per-target rankings,
which is why equality is stated against the shared-ranking serial path, not against
independent per-target runs), point order, jobs>1 equality, and per-relation
BadSpecialization isolation in the plural certifier.
"""

from __future__ import annotations

from fractions import Fraction

import pytest
import sympy as sp

from parametric_ibp_lf_reducer import ParamExpr, Row, parse_family_text
from parametric_ibp_lf_reducer.certificate import (
    STATUS_BAD_SPECIALIZATION,
    verify_reduction_relation_mod_p,
    verify_reduction_relations_mod_p,
)
from parametric_ibp_lf_reducer.ranking import rank_labels
from parametric_ibp_lf_reducer.records import (
    collect_normal_form_records,
    collect_normal_form_records_multi,
)

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

T0 = (0, 0, 0, 0)
T1 = (0, 1, 0, 0)
M = (0, 0, -1, 0)
LF_MAP = {T0: True, T1: True, M: True}


def _family():
    return parse_family_text(GENERIC_FAMILY_TEXT)


def _samples(vals):
    return [{"ep": Fraction(v)} for v in vals]


def _row(params, terms):
    row = Row(kind="synthetic", provenance={})
    for label, expr in terms.items():
        row.add_term(label, ParamExpr.from_sympy(sp.sympify(expr), params))
    return row


def _two_target_rows(fam):
    return [
        _row(fam.parameters, {T0: 1, M: "ep + 3"}),
        _row(fam.parameters, {T1: 1, M: "2*ep - 5", T0: "ep"}),
    ]


def _shared_ranking(fam, rows, targets):
    labels = sorted({c for row in rows for c in row.terms} | set(targets))
    return rank_labels(fam, labels, targets=list(targets), lf_map=LF_MAP)


# --- records: single target is bit-identical to the singular collector --------------------------
def test_multi_single_target_bit_identical():
    fam = _family()
    rows = _two_target_rows(fam)
    samples = _samples([1, 2, 3])
    ranking = _shared_ranking(fam, rows, [T0])

    singular = collect_normal_form_records(
        fam, rows, T0, PRIMES, samples, lf_map=LF_MAP, ranking=ranking
    )
    multi = collect_normal_form_records_multi(
        fam, rows, [T0], PRIMES, samples, lf_map=LF_MAP, ranking=ranking
    )
    assert list(multi) == [T0]
    assert multi[T0] == singular  # bit-for-bit


# --- records: several targets match the shared-ranking serial baseline --------------------------
def test_multi_two_targets_match_shared_ranking_serial():
    """Multi run == per-target singular runs *with the same shared ranking* (honesty caveat)."""
    fam = _family()
    rows = _two_target_rows(fam)
    samples = _samples([1, 2, 3, 4])
    ranking = _shared_ranking(fam, rows, [T0, T1])

    multi = collect_normal_form_records_multi(
        fam, rows, [T0, T1], PRIMES, samples, lf_map=LF_MAP, ranking=ranking
    )
    for tgt in (T0, T1):
        singular = collect_normal_form_records(
            fam, rows, tgt, PRIMES, samples, lf_map=LF_MAP, ranking=ranking
        )
        assert multi[tgt] == singular
        assert len(multi[tgt]) == len(PRIMES) * len(samples)  # every point recorded


def test_multi_point_order_samples_outer_primes_inner():
    fam = _family()
    rows = _two_target_rows(fam)
    samples = _samples([1, 2])
    multi = collect_normal_form_records_multi(fam, rows, [T0, T1], PRIMES, samples, lf_map=LF_MAP)
    expected = [(s["ep"], p) for s in samples for p in PRIMES]
    for tgt in (T0, T1):
        got = [(Fraction(r.sample["ep"]), r.prime) for r in multi[tgt]]
        assert got == expected


# --- records: jobs>1 equality and validation ----------------------------------------------------
def test_multi_jobs2_bit_identical_to_serial():
    fam = _family()
    rows = _two_target_rows(fam)
    samples = _samples([1, 2, 3, 4, 5])

    serial = collect_normal_form_records_multi(
        fam, rows, [T0, T1], PRIMES, samples, lf_map=LF_MAP, jobs=1
    )
    parallel = collect_normal_form_records_multi(
        fam, rows, [T0, T1], PRIMES, samples, lf_map=LF_MAP, jobs=2
    )
    assert parallel == serial


@pytest.mark.parametrize("bad", [0, -1, True, 2.0, "2", None])
def test_multi_jobs_rejects_non_positive_and_non_int(bad):
    fam = _family()
    with pytest.raises(ValueError):
        collect_normal_form_records_multi(
            fam, _two_target_rows(fam), [T0], PRIMES[:1], _samples([1]), lf_map=LF_MAP, jobs=bad
        )


def test_multi_requires_at_least_one_target():
    fam = _family()
    with pytest.raises(ValueError):
        collect_normal_form_records_multi(
            fam, _two_target_rows(fam), [], PRIMES[:1], _samples([1]), lf_map=LF_MAP
        )


def test_multi_deduplicates_targets_preserving_order():
    fam = _family()
    rows = _two_target_rows(fam)
    multi = collect_normal_form_records_multi(
        fam, rows, [T1, T0, T1], PRIMES[:1], _samples([1]), lf_map=LF_MAP
    )
    assert list(multi) == [T1, T0]


# --- certificates: plural verifier matches the singular one -------------------------------------
def _cert_setup():
    fam = _family()
    rows = _two_target_rows(fam)
    relations = {
        T0: {M: sp.sympify("-(ep + 3)")},
        T1: {M: sp.sympify("-(2*ep - 5)"), T0: sp.sympify("-ep")},
    }
    sample = {"ep": Fraction(1, 3)}
    return fam, rows, relations, sample


def test_plural_certificates_match_singular():
    fam, rows, relations, sample = _cert_setup()
    for prime in PRIMES:
        plural = verify_reduction_relations_mod_p(fam, rows, relations, sample, prime)
        assert list(plural) == [T0, T1]
        for tgt, terms in relations.items():
            single = verify_reduction_relation_mod_p(fam, rows, tgt, terms, sample, prime)
            assert plural[tgt] == single
            assert plural[tgt].in_span is True


def test_plural_certificate_bad_relation_isolated_per_target():
    """A pole in ONE relation's claimed coefficients must not poison the other target."""
    fam, rows, relations, sample = _cert_setup()
    relations = dict(relations)
    relations[T1] = {M: sp.sympify("1/(3*ep - 1)")}  # pole exactly at ep = 1/3
    plural = verify_reduction_relations_mod_p(fam, rows, relations, sample, PRIMES[0])
    assert plural[T1].status == STATUS_BAD_SPECIALIZATION
    assert plural[T0].in_span is True


def test_plural_certificate_empty_rows_marks_all_targets():
    fam, _, relations, sample = _cert_setup()
    plural = verify_reduction_relations_mod_p(fam, [], relations, sample, PRIMES[0])
    assert {r.status for r in plural.values()} == {"EmptySystem"}
