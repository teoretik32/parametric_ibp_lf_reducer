"""Dual LF-obstruction certificate modulo a prime (Method.6, External Int2).

Companion of :mod:`lf_feasibility`. Where the span test answers *whether* the target can be
reduced through the allowed (locally finite) labels, this module produces, for an ``Obstructed``
system, an explicit **dual witness** ``w`` proving it.

Math. The projected rows span relations restricted to the forbidden columns (incl. the target).
``e_target`` lies in that row span iff the system is ``Feasible``. When it does not, there is a
vector ``w`` in the RIGHT nullspace of the projected matrix — ``<row, w> == 0`` for every
projected row — with ``w[target] != 0``. Such a ``w`` certifies obstruction: were
``e_target = sum c_i row_i`` in the row span, then ``w[target] = <e_target, w> = sum c_i <row_i, w>
= 0``, a contradiction. We normalize ``w[target] == 1``.

Construction (from the RREF pivot rows, each of the form ``e_c + sum_{f free} a_{c,f} e_f``):

* target is a FREE column (equivalently: target is not a pivot, the ``Obstructed`` case whose
  canonical residual is exactly ``{target}``) — take the nullvector indexed by ``f0 = target``:
  ``w = {target: 1}`` and ``w[c] = -a_{c,target}`` for each pivot ``c``;
* target is a PIVOT column — reducing ``e_target`` leaves ``-sum_f a_{target,f} e_f``. If every
  ``a_{target,f}`` vanishes, ``e_target`` is in the row span (``Feasible``, empty witness).
  Otherwise take the first free column ``f0`` with ``a_{target,f0} != 0``, form the nullvector
  ``v`` indexed by ``f0`` (``v[target] = -a_{target,f0} != 0``) and rescale so ``w[target] == 1``.

Codimension correction (binding wording). ``residual_support == (target,)`` says only that the
canonical residual of ``e_target`` lands on the target coordinate; it does NOT imply the quotient
dimension is one. That dimension is the nullity ``n_projected_cols - rank`` and may exceed 1 (see
``tests/test_lf_obstruction_witness.py::test_nullity_gt_one_target_only_residual``: nullity 2 with
``residual_support == (target,)``). Earlier "codimension-one obstruction" phrasing is retracted.

Determinism rules.

1. Default ``column_order`` = ``sorted({columns present in projected rows} | {target})``
   (lexicographic tuple sort), passed to :func:`rref_mod_p`.
2. Target-pivot case: free column ``f0`` = first in that column order with ``r_target[f] != 0``.
3. The ``witness`` tuple is sorted by label; coefficients are reduced into ``[0, prime)``.
4. Rows are consumed in the given order (callers supply deterministic, dedup-merged row lists).

Scope. Read-only diagnostics reusing ``assemble_rows_mod_p``/``rref_mod_p`` and the
``lf_feasibility`` primitives; no reducer state, LF gate or certificate semantics is touched. A
``Witness`` is an honest per-(sample, prime), per-label-box negative — never a global impossibility
claim.

pytest gotcha. :func:`test_rows_against_obstruction_witness` matches pytest's ``test_*`` collection
pattern. Test modules MUST import this module (``import ... as low``) and never
``from ... import test_rows_against_obstruction_witness`` nor star-import it.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

from .labels import Label
from .lf_feasibility import (
    STATUS_BAD_SPECIALIZATION,
    STATUS_FEASIBLE,
    _partition_labels,
    _sample_key,
)
from .modular_normal_form import BadSpecialization, assemble_rows_mod_p
from .row_generation import Row
from .sparse_rref import rref_mod_p

STATUS_WITNESS = "Witness"


@dataclass(frozen=True)
class LFObstructionWitness:
    """Dual certificate of LF-obstruction at one ``(sample, prime)`` point.

    ``status`` is ``Witness`` (obstructed, witness produced), ``Feasible`` (target in the
    projected row span, empty witness) or ``BadSpecialization``.
    """

    status: str
    prime: int
    sample: tuple  # sorted (name, value) pairs — deterministic, JSON-safe
    target_label: Label
    witness: tuple  # ((label, coeff), ...) sorted by label, coeff in [0, prime); empty unless Witness
    rank: int  # rank of the projected matrix
    n_projected_cols: int  # |columns present in projected rows U {target}|
    nullity: int  # n_projected_cols - rank (quotient-dimension proxy; may exceed 1)
    nrows: int  # specialized rows
    n_projected_rows: int  # nonzero rows after deleting allowed columns
    n_allowed: int  # allowed (LF-True) labels, target excluded
    n_forbidden: int  # forbidden columns actually present (target excluded)
    check_annihilation: bool  # exact: <row, w> == 0 for every projected row
    check_target_unit: bool  # exact: w[target] == 1
    detail: str = ""


@dataclass(frozen=True)
class RowPairing:
    """Pairing of one candidate row against a stored witness ``w``."""

    index: int  # index into the candidate rows list
    kind: str  # Row.kind
    pairing: int  # <row_specialized, w> mod prime, in [0, prime)
    breaks: bool  # pairing != 0 -> this row could cure THIS witness (necessary, not sufficient)


def _bad_specialization_result(
    prime: int, sample: Mapping, target_label: Label, n_allowed: int, detail: str
) -> LFObstructionWitness:
    return LFObstructionWitness(
        status=STATUS_BAD_SPECIALIZATION,
        prime=prime,
        sample=_sample_key(sample),
        target_label=target_label,
        witness=(),
        rank=0,
        n_projected_cols=0,
        nullity=0,
        nrows=0,
        n_projected_rows=0,
        n_allowed=n_allowed,
        n_forbidden=0,
        check_annihilation=False,
        check_target_unit=False,
        detail=detail,
    )


def lf_obstruction_witness_mod_p(
    rows: Sequence[Row],
    labels: Sequence[Label],
    target_label: Label,
    lf_flags: Mapping[Label, object],
    sample: Mapping,
    prime: int,
    column_order: Sequence[Label] | None = None,
) -> LFObstructionWitness:
    """Dual witness for LF-obstruction of ``target_label`` (mod ``prime`` at ``sample``).

    Returns a :class:`LFObstructionWitness`: ``Witness`` with an explicit right-nullspace vector
    ``w`` (``w[target] == 1``, ``<row, w> == 0`` for every projected row) when the target cannot be
    reduced through the allowed labels; ``Feasible`` with an empty witness when it can;
    ``BadSpecialization`` (never raised) when a coefficient denominator vanishes at the point.
    """
    allowed = _partition_labels(labels, target_label, lf_flags)
    try:
        # ``assemble_rows_mod_p`` never dereferences its family argument (kept for signature
        # compatibility with the certificate call sites), so ``None`` is safe here.
        matrix = assemble_rows_mod_p(None, rows, sample, prime)
    except BadSpecialization as exc:
        return _bad_specialization_result(prime, sample, target_label, len(allowed), str(exc))

    projected: list[dict] = []
    forbidden_present: set[Label] = set()
    for row in matrix:
        proj = {c: v for c, v in row.items() if c not in allowed}
        if proj:
            projected.append(proj)
            forbidden_present.update(proj)
    forbidden_present.discard(target_label)

    cols = sorted({c for r in projected for c in r} | {target_label})
    order = list(column_order) if column_order is not None else cols
    res = rref_mod_p(projected, prime, column_order=order)
    pivots = res.pivots
    n_projected_cols = len(cols)
    nullity = n_projected_cols - res.rank

    if target_label not in pivots:
        # target is a free column -> the nullvector indexed by the target has w[target] == 1.
        w: dict[Label, int] = {target_label: 1}
        for col, prow in pivots.items():
            coeff = prow.get(target_label)
            if coeff:
                w[col] = (-coeff) % prime
        status = STATUS_WITNESS
    else:
        r_t = pivots[target_label]
        free_support = [f for f in order if f not in pivots and r_t.get(f)]
        if not free_support:
            # e_target lies in the projected row span: feasible, no witness.
            return LFObstructionWitness(
                status=STATUS_FEASIBLE,
                prime=prime,
                sample=_sample_key(sample),
                target_label=target_label,
                witness=(),
                rank=res.rank,
                n_projected_cols=n_projected_cols,
                nullity=nullity,
                nrows=len(matrix),
                n_projected_rows=len(projected),
                n_allowed=len(allowed),
                n_forbidden=len(forbidden_present),
                check_annihilation=True,
                check_target_unit=False,
                detail="target unit vector lies in projected row span",
            )
        f0 = free_support[0]
        v: dict[Label, int] = {f0: 1}
        for col, prow in pivots.items():
            coeff = prow.get(f0)
            if coeff:
                v[col] = (-coeff) % prime
        scale = pow(v[target_label], prime - 2, prime)  # v[target] = -r_t[f0] != 0
        w = {lab: (x * scale) % prime for lab, x in v.items() if (x * scale) % prime}
        status = STATUS_WITNESS

    check_annihilation = all(
        sum(v * w.get(c, 0) for c, v in row.items()) % prime == 0 for row in projected
    )
    check_target_unit = w.get(target_label) == 1
    witness = tuple(sorted((lab, coeff % prime) for lab, coeff in w.items() if coeff % prime))
    return LFObstructionWitness(
        status=status,
        prime=prime,
        sample=_sample_key(sample),
        target_label=target_label,
        witness=witness,
        rank=res.rank,
        n_projected_cols=n_projected_cols,
        nullity=nullity,
        nrows=len(matrix),
        n_projected_rows=len(projected),
        n_allowed=len(allowed),
        n_forbidden=len(forbidden_present),
        check_annihilation=check_annihilation,
        check_target_unit=check_target_unit,
        detail="dual obstruction witness (right nullspace, w[target]=1)",
    )


def witness_to_payload(result: LFObstructionWitness) -> dict:
    """JSON-safe dict (deterministic key order; labels as lists). FULL witness preserved."""
    return {
        "status": result.status,
        "prime": result.prime,
        "sample": [list(pair) for pair in result.sample],
        "target_label": list(result.target_label),
        "witness": [[list(lab), coeff] for lab, coeff in result.witness],
        "rank": result.rank,
        "n_projected_cols": result.n_projected_cols,
        "nullity": result.nullity,
        "nrows": result.nrows,
        "n_projected_rows": result.n_projected_rows,
        "n_allowed": result.n_allowed,
        "n_forbidden": result.n_forbidden,
        "check_annihilation": result.check_annihilation,
        "check_target_unit": result.check_target_unit,
        "detail": result.detail,
    }


def witness_from_payload(payload: Mapping) -> LFObstructionWitness:
    """Inverse of :func:`witness_to_payload` (lists -> tuples), for probe-mode reload."""
    return LFObstructionWitness(
        status=payload["status"],
        prime=payload["prime"],
        sample=tuple(tuple(pair) for pair in payload["sample"]),
        target_label=tuple(payload["target_label"]),
        witness=tuple((tuple(lab), int(coeff)) for lab, coeff in payload["witness"]),
        rank=payload["rank"],
        n_projected_cols=payload["n_projected_cols"],
        nullity=payload["nullity"],
        nrows=payload["nrows"],
        n_projected_rows=payload["n_projected_rows"],
        n_allowed=payload["n_allowed"],
        n_forbidden=payload["n_forbidden"],
        check_annihilation=payload["check_annihilation"],
        check_target_unit=payload["check_target_unit"],
        detail=payload.get("detail", ""),
    )


def test_rows_against_obstruction_witness(
    rows: Sequence[Row],
    witness: LFObstructionWitness,
    sample: Mapping,
    prime: int,
) -> tuple[RowPairing, ...]:
    """Pair each candidate row against a stored witness ``w`` (NO RREF).

    ``breaks=True`` (``pairing != 0``) means the row is not orthogonal to ``w``, so adding it to
    the system would invalidate THIS obstruction witness — a NECESSARY, not sufficient, condition
    for curing the obstruction (other nullvectors may still obstruct). ``breaks=False`` means the
    row annihilates ``w`` and cannot cure THIS particular witness. ``w`` is supported only on
    forbidden columns (incl. target), so pairing against the full row equals the projected pairing.
    """
    w = dict(witness.witness)
    spec = assemble_rows_mod_p(None, rows, sample, prime)
    out: list[RowPairing] = []
    for i, full in enumerate(spec):
        pairing = sum(v * w[c] for c, v in full.items() if c in w) % prime
        out.append(RowPairing(index=i, kind=rows[i].kind, pairing=pairing, breaks=pairing != 0))
    return tuple(out)


def pairings_to_payload(pairings: Sequence[RowPairing]) -> list[dict]:
    """JSON-safe list of pairing records (deterministic key order)."""
    return [
        {"index": p.index, "kind": p.kind, "pairing": p.pairing, "breaks": p.breaks}
        for p in pairings
    ]
