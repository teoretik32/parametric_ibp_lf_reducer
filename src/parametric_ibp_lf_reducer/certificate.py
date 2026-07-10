"""Modular row-span certificates for *given* reduction relations (Pass D4.4).

Verifies that a claimed relation ``J[target] = sum_i C_i(params) * J[label_i]`` lies in the row
span of a generated row system at one exact ``(sample, prime)`` point:

1. assemble the parametric rows to integer rows modulo ``prime``;
2. evaluate every claimed coefficient ``C_i`` at the sample, exactly, modulo ``prime``;
3. reduce the relation vector ``J[target] - sum_i C_i * J[label_i]`` by the RREF pivot rows.

A zero residual certifies the relation is a linear combination of the generated rows at that
point; a nonzero residual is reported honestly with the surviving columns. This is *evidence*,
not a verdict: nothing here stamps ``Success`` (the strict gate in :mod:`result` remains the only
place), no local-finiteness is asserted, and a bad specialization (a vanishing denominator in a
row or claimed coefficient) rejects the point instead of being patched.

Coefficients may be SymPy expressions/numbers, :class:`ParamExpr` (anything with ``eval_mod_p``),
or exact ``int``/``Fraction`` values. Evaluation is exact — no floats. SymPy is used only for
this one-off coefficient evaluation (validation utility, not the row-generation/RREF hot loop).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from fractions import Fraction

import sympy as sp

from .family import ParametricFamily
from .labels import Label
from .modular_normal_form import (
    STATUS_BAD_SPECIALIZATION,
    STATUS_EMPTY_SYSTEM,
    BadSpecialization,
    assemble_rows_mod_p,
)
from .row_generation import Row
from .sparse_rref import rref_mod_p

STATUS_IN_SPAN = "InSpan"
STATUS_NOT_IN_SPAN = "NotInSpan"


@dataclass
class CertificateResult:
    """Outcome of one modular row-span certificate check (never a ``Success`` stamp)."""

    status: str  # InSpan | NotInSpan | BadSpecialization | EmptySystem
    in_span: bool
    target_label: Label
    prime: int
    sample: dict
    relation: dict = field(default_factory=dict)  # assembled relation vector mod prime
    residual: dict = field(default_factory=dict)  # nonzero leftover columns (empty iff in span)
    nrows: int = 0
    rank: int = 0


def _coeff_mod_p(coeff, sample: Mapping, prime: int) -> int:
    """Evaluate one claimed coefficient at ``sample`` modulo ``prime``, exactly.

    Raises :class:`BadSpecialization` if its denominator vanishes modulo ``prime`` (or at the
    sample), and ``ValueError`` if it does not evaluate to an exact rational (leftover symbols).
    """
    if hasattr(coeff, "eval_mod_p"):  # ParamExpr path
        v = coeff.eval_mod_p(dict(sample), prime)
        if v is None:
            raise BadSpecialization(
                f"claimed coefficient denominator vanishes mod {prime} at {dict(sample)}"
            )
        return v % prime
    expr = sp.sympify(coeff)
    subs = {sp.Symbol(str(k)): sp.Rational(Fraction(v)) for k, v in sample.items()}
    value = expr.subs(subs)
    if value.is_finite is False:  # pole at the sample (zoo/oo) -> reject the point, honestly
        raise BadSpecialization(
            f"claimed coefficient {coeff!r} has a pole at {dict(sample)}"
        )
    try:
        value = sp.Rational(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"claimed coefficient {coeff!r} did not evaluate to an exact rational "
            f"at {dict(sample)} (got {value!r})"
        ) from exc
    num, den = int(value.p), int(value.q)
    if den % prime == 0:
        raise BadSpecialization(
            f"claimed coefficient denominator vanishes mod {prime} at {dict(sample)}"
        )
    return num % prime * pow(den % prime, prime - 2, prime) % prime


def _reduce_vector_by_pivots(vec: dict, pivots: Mapping, prime: int) -> dict:
    """Eliminate every pivot column from ``vec`` using the (mutually reduced) RREF pivot rows."""
    w = {c: v % prime for c, v in vec.items() if v % prime}
    for col, prow in pivots.items():
        f = w.get(col, 0)
        if not f:
            continue
        for cc, vv in prow.items():  # prow has a 1 at ``col`` and only free-col entries otherwise
            nv = (w.get(cc, 0) - f * vv) % prime
            if nv:
                w[cc] = nv
            elif cc in w:
                del w[cc]
    return w


def _relation_vector_mod_p(
    target_label: Label,
    terms: Mapping[Label, object],
    lhs_terms: Mapping[Label, object] | None,
    sample: Mapping,
    prime: int,
) -> dict:
    """Assemble the relation vector ``sum L_j J_j - sum C_i J_i`` mod ``prime`` (may raise
    :class:`BadSpecialization`)."""
    if lhs_terms is None:
        lhs_terms = {target_label: 1}
    relation: dict = {}
    for label in sorted(lhs_terms):
        c = _coeff_mod_p(lhs_terms[label], sample, prime)
        relation[label] = (relation.get(label, 0) + c) % prime
    for label in sorted(terms):
        c = _coeff_mod_p(terms[label], sample, prime)
        relation[label] = (relation.get(label, 0) - c) % prime
    return {c: v for c, v in relation.items() if v}


def verify_reduction_relations_mod_p(
    family: ParametricFamily,
    rows: Sequence[Row],
    relations: Mapping[Label, Mapping[Label, object]],
    sample: Mapping,
    prime: int,
    column_order: Sequence[Label] | None = None,
    lhs_terms_map: Mapping[Label, Mapping[Label, object]] | None = None,
) -> dict[Label, CertificateResult]:
    """Perf.5: certify *several* claimed relations at ONE ``(sample, prime)`` point, sharing the
    assemble + RREF work.

    ``relations`` maps target label -> ``terms`` (as in
    :func:`verify_reduction_relation_mod_p`); ``lhs_terms_map`` optionally supplies a per-target
    ``lhs_terms``. Returns ``{target: CertificateResult}`` with each result identical to what the
    single-relation verifier would produce at the same point: the matrix assembly and RREF do not
    depend on the claimed relation, so only the cheap per-relation coefficient evaluation and
    pivot reduction are repeated. A bad specialization *in the rows* rejects the whole point (all
    targets get ``BadSpecialization``); a bad specialization in one relation's claimed
    coefficients rejects only that target.
    """
    rows = list(rows)
    targets = list(relations)
    lhs_terms_map = lhs_terms_map or {}

    def _all(status: str) -> dict[Label, CertificateResult]:
        return {
            tgt: CertificateResult(
                status=status,
                in_span=False,
                target_label=tgt,
                prime=prime,
                sample=dict(sample),
            )
            for tgt in targets
        }

    if not rows:
        return _all(STATUS_EMPTY_SYSTEM)
    try:
        matrix = assemble_rows_mod_p(family, rows, dict(sample), prime)
    except BadSpecialization:
        return _all(STATUS_BAD_SPECIALIZATION)
    if not matrix:
        return _all(STATUS_EMPTY_SYSTEM)

    res = rref_mod_p(matrix, prime, column_order=column_order)
    out: dict[Label, CertificateResult] = {}
    for tgt in targets:
        base = dict(in_span=False, target_label=tgt, prime=prime, sample=dict(sample))
        try:
            relation = _relation_vector_mod_p(
                tgt, relations[tgt], lhs_terms_map.get(tgt), sample, prime
            )
        except BadSpecialization:
            out[tgt] = CertificateResult(status=STATUS_BAD_SPECIALIZATION, **base)
            continue
        residual = _reduce_vector_by_pivots(relation, res.pivots, prime)
        in_span = residual == {}
        out[tgt] = CertificateResult(
            status=STATUS_IN_SPAN if in_span else STATUS_NOT_IN_SPAN,
            in_span=in_span,
            target_label=tgt,
            prime=prime,
            sample=dict(sample),
            relation=relation,
            residual=residual,
            nrows=len(matrix),
            rank=res.rank,
        )
    return out


def verify_reduction_relation_mod_p(
    family: ParametricFamily,
    rows: Sequence[Row],
    target_label: Label,
    terms: Mapping[Label, object],
    sample: Mapping,
    prime: int,
    column_order: Sequence[Label] | None = None,
    lhs_terms: Mapping[Label, object] | None = None,
) -> CertificateResult:
    """Certify (mod ``prime``, at ``sample``) that ``J[target] = sum terms[label]*J[label]``
    is in the span of ``rows``.

    ``terms`` maps master label -> claimed coefficient (SymPy expr, ``ParamExpr``, ``int`` or
    ``Fraction``). ``lhs_terms`` optionally generalizes the left-hand side to a linear
    combination ``sum lhs_terms[label]*J[label]`` (same coefficient types); when ``None`` the
    LHS is the classic single ``J[target]`` with unit coefficient. ``column_order`` optionally
    fixes the RREF pivot order (span membership does not depend on it; exposed for determinism
    experiments). Returns a :class:`CertificateResult`; bad specializations reject the point
    honestly (status, never a patch).
    """
    rows = list(rows)
    base = dict(
        in_span=False, target_label=target_label, prime=prime, sample=dict(sample)
    )
    if not rows:
        return CertificateResult(status=STATUS_EMPTY_SYSTEM, **base)

    try:
        matrix = assemble_rows_mod_p(family, rows, dict(sample), prime)
    except BadSpecialization:
        return CertificateResult(status=STATUS_BAD_SPECIALIZATION, **base)
    if not matrix:
        return CertificateResult(status=STATUS_EMPTY_SYSTEM, **base)

    # relation vector for sum_j L_j * J[lhs_label_j] - sum_i C_i * J[label_i] = 0
    # (default LHS: the single target with unit coefficient, i.e. J[target] - sum C_i J_i = 0)
    if lhs_terms is None:
        lhs_terms = {target_label: 1}
    try:
        relation: dict = {}
        for label in sorted(lhs_terms):
            c = _coeff_mod_p(lhs_terms[label], sample, prime)
            relation[label] = (relation.get(label, 0) + c) % prime
        for label in sorted(terms):
            c = _coeff_mod_p(terms[label], sample, prime)
            relation[label] = (relation.get(label, 0) - c) % prime
    except BadSpecialization:
        return CertificateResult(status=STATUS_BAD_SPECIALIZATION, **base)
    relation = {c: v for c, v in relation.items() if v}

    res = rref_mod_p(matrix, prime, column_order=column_order)
    residual = _reduce_vector_by_pivots(relation, res.pivots, prime)
    in_span = residual == {}
    return CertificateResult(
        status=STATUS_IN_SPAN if in_span else STATUS_NOT_IN_SPAN,
        in_span=in_span,
        target_label=target_label,
        prime=prime,
        sample=dict(sample),
        relation=relation,
        residual=residual,
        nrows=len(matrix),
        rank=res.rank,
    )
