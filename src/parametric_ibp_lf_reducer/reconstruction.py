"""Coefficient reconstruction from modular normal forms (spec §5.10, §7).

Turns integer normal-form values collected over many ``(prime, sample)`` points into exact
rational coefficient functions ``C_a(params)``:

1. rank-consistency record selection (Pass D4.3): only ``Reduced``/formal-success records are
   coefficient records at all, and by default only those at the *maximal observed RREF rank* are
   consumed — a rank-deficient specialization solves a smaller linear system, so its
   shrunken/shifted normal-form support must never be union-0-filled into the value table;
2. multi-prime CRT + rational reconstruction -> exact rational *value* at each sample point;
3. union support across the selected samples (a term missing at a *max-rank* point contributes
   value 0, not dropped);
4. rational-function interpolation (univariate via a degree search; multivariate via a dense
   linear-algebra ansatz), each with an INDEPENDENT holdout validation -> :class:`InterpolationFailed`
   if it does not validate.

Bad specializations are skipped (and counted), never patched. The multivariate reconstruction is
dense (not sparse/Zippel) and validated; when it cannot be pinned down and validated it raises
rather than guessing. No ``Success`` here.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from fractions import Fraction
from itertools import combinations_with_replacement
from math import gcd, isqrt

import sympy as sp
from sympy.polys.polyfuncs import rational_interpolate

from .modular_normal_form import STATUS_REDUCED, NormalFormResult


class InterpolationFailed(Exception):
    """Raised when coefficient reconstruction cannot be validated."""


def rational_reconstruction(a: int, m: int) -> Fraction | None:
    """Recover ``p/q`` with ``|p|, q <= sqrt(m/2)`` and ``p/q == a (mod m)``; ``None`` if none."""
    a %= m
    if a == 0:
        return Fraction(0)
    bound = isqrt(m // 2)
    old_r, r = m, a
    old_t, t = 0, 1
    while r > bound:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_t, t = t, old_t - q * t
    num, den = r, t
    if den == 0 or abs(den) > bound:
        return None
    if den < 0:
        num, den = -num, -den
    if gcd(abs(num), den) != 1 or (num - a * den) % m != 0:
        return None
    return Fraction(num, den)


def _crt(r1: int, m1: int, r2: int, m2: int) -> tuple[int, int]:
    inv = pow(m1 % m2, -1, m2)
    t = ((r2 - r1) * inv) % m2
    return (r1 + m1 * t) % (m1 * m2), m1 * m2


def reconstruct_rational(residues: dict[int, int]) -> Fraction | None:
    """CRT-combine ``{prime: residue}`` and rational-reconstruct, requiring stability."""
    acc_r, acc_m = 0, 1
    prev: Fraction | None = None
    for p in sorted(residues):
        acc_r, acc_m = _crt(acc_r, acc_m, residues[p] % p, p)
        cand = rational_reconstruction(acc_r, acc_m)
        if cand is not None and cand == prev:
            return cand
        prev = cand
    return None  # never stabilized -> need more primes


def _sample_key(sample: dict) -> tuple:
    return tuple(sorted((str(k), Fraction(v)) for k, v in sample.items()))


def _record_coeffs(r) -> dict:
    """Coefficient dict of a normal-form point, from either a ``NormalFormRecord`` (``coeffs``)
    or a single-sample ``NormalFormResult`` (``terms``)."""
    c = getattr(r, "coeffs", None)
    return getattr(r, "terms", {}) if c is None else c


# --- rank-consistency record selection (Pass D4.3) --------------------------------------------
RANK_POLICY_MAX_RANK = "max_rank"
RANK_POLICY_ALL = "all"

# --- Method.11b/11c: support-deviation vs unstable classification -----------------------------
SUPPORT_STABLE = "stable"
SUPPORT_DEVIATION_PENDING = "support_deviation_pending_validation"
SUPPORT_UNSTABLE = "unstable"
# Legacy Method.11b class name.  ``classify_sample_supports`` no longer emits it: a same-rank,
# prime-consistent missing-label pattern is NOT proof of a coefficient zero.  The cache sample
# ep=6, r=57/11 has the generic rank 26984, yet its three surviving coefficients all disagree
# with the generic coefficient functions (basis/pivot specialization); zero-filling it poisoned
# every dense interpolation in Stages 2..4.  Kept only for backward-compatible imports.
SUPPORT_SPECIAL_ZERO = "special_zero"
# Post-reconstruction verdicts for pending support deviations (Method.11c):
DEVIATION_GENUINE_SPECIAL_ZERO = "GenuineSpecialZero"
DEVIATION_BASIS_SPECIALIZATION = "BasisSpecialization"


def _canonical_sample_str(sample: dict) -> str:
    """Canonical ``'ep=6,r=57/11'``-style key for one exact parameter sample."""
    return ",".join(
        f"{k}={Fraction(v)}" for k, v in sorted(sample.items(), key=lambda kv: str(kv[0]))
    )


def classify_sample_supports(records: Iterable) -> dict:
    """Classify per-sample support deviations of formally reduced records (Method.11b).

    The generic rank is the maximal rank observed over valid (``Reduced`` +
    ``formal_success``) records; the reference support is the union of coefficient labels
    over generic-rank records. Per exact parameter sample:

    - ``support_deviation_pending_validation`` -- every record at the sample has the generic
      rank, the per-prime supports agree, and the same non-empty label set is missing versus
      the reference support.  This is *not* proof of a coefficient zero: the normal-form
      basis/pivot structure may specialize at the sample while keeping the generic rank
      (Method.11c, ep=6, r=57/11).  Such samples are excluded from fitting by
      :func:`collect_value_table` and must only be re-admitted as genuine special zeros by
      :func:`validate_support_deviations` after a candidate reconstruction.
    - ``unstable`` -- ranks differ at the sample, per-prime supports disagree, a record is
      not formally reduced, or labels outside the reference union appear; such samples must
      never be zero-filled and are rejected by :func:`collect_value_table`.
    - ``stable`` -- full reference support at every prime.

    Returns ``{"by_sample": {canonical_sample: entry}, "summary": {...}}``.
    """
    records = list(records)
    valid = [r for r in records if r.status == STATUS_REDUCED and r.formal_success]
    summary: dict = {
        "generic_rank": None,
        "reference_support_size": 0,
        "n_samples": 0,
        "n_stable": 0,
        "n_support_deviation_pending_validation": 0,
        "n_unstable": 0,
        "support_deviation_samples": [],
        "unstable_samples": [],
    }
    if not valid:
        return {"by_sample": {}, "summary": summary}
    generic_rank = max(r.rank for r in valid)
    ref_set = {lab for r in valid if r.rank == generic_rank for lab in _record_coeffs(r)}
    summary["generic_rank"] = generic_rank
    summary["reference_support_size"] = len(ref_set)

    grouped: dict = {}
    for r in records:
        grouped.setdefault(_canonical_sample_str(r.sample), []).append(r)

    by_sample: dict = {}
    for key in sorted(grouped):
        recs = grouped[key]
        vrecs = [r for r in recs if r.status == STATUS_REDUCED and r.formal_success]
        supports = {r.prime: frozenset(_record_coeffs(r)) for r in vrecs}
        ranks = sorted({r.rank for r in vrecs})
        seen = set().union(*supports.values()) if supports else set()
        missing = sorted(ref_set - seen)
        extra = sorted(seen - ref_set)
        if len(vrecs) != len(recs) or not vrecs:
            cls = SUPPORT_UNSTABLE  # target not pivot / bad specialization at this sample
        elif ranks != [generic_rank]:
            cls = SUPPORT_UNSTABLE  # rank drop: shrunken/shifted support, not a value zero
        elif len(set(supports.values())) > 1:
            cls = SUPPORT_UNSTABLE  # primes disagree on the missing-label pattern
        elif extra:
            cls = SUPPORT_UNSTABLE  # labels beyond the reference union (defensive)
        elif missing:
            # Same labels absent for every prime at the generic rank: a support deviation.
            # NOT auto-classified as a special zero (Method.11c) -- the basis/pivot structure
            # may specialize while keeping the rank; exclude from fitting pending validation.
            cls = SUPPORT_DEVIATION_PENDING
        else:
            cls = SUPPORT_STABLE
        by_sample[key] = {
            "classification": cls,
            "n_records": len(recs),
            "n_valid_records": len(vrecs),
            "ranks": ranks,
            "missing_labels": [list(lab) for lab in missing],
        }
        summary["n_samples"] += 1
        summary[f"n_{cls}"] += 1
        if cls == SUPPORT_DEVIATION_PENDING:
            summary["support_deviation_samples"].append(key)
        elif cls == SUPPORT_UNSTABLE:
            summary["unstable_samples"].append(key)
    return {"by_sample": by_sample, "summary": summary}


def select_records_for_reconstruction(records: Iterable, rank_policy: str = RANK_POLICY_MAX_RANK):
    """Select the normal-form records that reconstruction may consume as coefficient records.

    Only ``Reduced`` records with ``formal_success`` are coefficient records at all. Under the
    default ``"max_rank"`` policy only the valid records whose RREF ``rank`` equals the maximal
    observed rank are kept: a specialization's rank can only *drop* below the generic rank, so a
    rank-deficient point solves a smaller system and its shrunken/shifted support must not be
    union-0-filled into the value table (the spurious zeros poison interpolation). A label
    *missing* from a max-rank record still means an exact zero at that point. ``"all"`` keeps
    every valid record (the pre-D4.3 behaviour, for tests/debugging).

    Returns ``(selected, diagnostics)`` where ``diagnostics`` reports the selection honestly:
    ``rank_policy``, ``n_records``, ``n_valid_records_before_rank_filter``, ``selected_rank``
    (``None`` under ``"all"`` or with no valid records), ``n_rank_filtered_records``,
    ``rank_histogram`` and ``support_after_rank_filter``.
    """
    if rank_policy not in (RANK_POLICY_MAX_RANK, RANK_POLICY_ALL):
        raise ValueError(f"unknown rank_policy {rank_policy!r}")
    records = list(records)
    valid = [r for r in records if r.status == STATUS_REDUCED and r.formal_success]
    histogram: dict[int, int] = {}
    for r in valid:
        histogram[r.rank] = histogram.get(r.rank, 0) + 1

    selected_rank: int | None = None
    if rank_policy == RANK_POLICY_ALL or not valid:
        selected = list(valid)
    else:  # max_rank
        selected_rank = max(histogram)
        selected = [r for r in valid if r.rank == selected_rank]

    support = sorted({lab for r in selected for lab in _record_coeffs(r)})
    diagnostics = {
        "rank_policy": rank_policy,
        "n_records": len(records),
        "n_valid_records_before_rank_filter": len(valid),
        "selected_rank": selected_rank,
        "n_selected_records": len(selected),
        "n_rank_filtered_records": len(valid) - len(selected),
        "rank_histogram": dict(sorted(histogram.items())),
        "support_after_rank_filter": tuple(support),
    }
    return selected, diagnostics


def collect_value_table(
    results: Iterable[NormalFormResult], rank_policy: str = RANK_POLICY_MAX_RANK
):
    """Group reduced results by sample, reconstruct exact rational values with union support.

    Records are first passed through :func:`select_records_for_reconstruction` (default:
    ``"max_rank"``), so rank-deficient specializations never contribute zero-filled coefficients.
    Under the ``"max_rank"`` policy, samples classified ``unstable`` by
    :func:`classify_sample_supports` (rank drops, per-prime support disagreement) are rejected
    entirely, and ``support_deviation_pending_validation`` samples (same generic rank, agreed
    non-empty missing-label set) are likewise excluded from the value table -- never zero-filled
    (Method.11c); a genuine special zero can only be recognized *after* reconstruction via
    :func:`validate_support_deviations`. Returns ``(labels, table, samples, n_skipped)`` where
    ``table[label][sample_key]`` is an exact :class:`Fraction` and ``n_skipped`` counts records
    not consumed (non-reduced records plus rank/stability/deviation-filtered ones).
    """
    results = list(results)
    reduced, _ = select_records_for_reconstruction(results, rank_policy=rank_policy)
    if rank_policy == RANK_POLICY_MAX_RANK:
        excluded = {
            key
            for key, entry in classify_sample_supports(results)["by_sample"].items()
            if entry["classification"] in (SUPPORT_UNSTABLE, SUPPORT_DEVIATION_PENDING)
        }
        if excluded:
            reduced = [r for r in reduced if _canonical_sample_str(r.sample) not in excluded]
    n_skipped = len(results) - len(reduced)

    by_sample: dict = {}
    all_labels: set = set()
    for r in reduced:
        key = _sample_key(r.sample)
        entry = by_sample.setdefault(key, {"sample": dict(r.sample), "byprime": {}})
        coeffs = _record_coeffs(r)
        entry["byprime"][r.prime] = dict(coeffs)
        all_labels.update(coeffs)

    labels = sorted(all_labels)
    table: dict = {lab: {} for lab in labels}
    samples: dict = {}
    for key, entry in by_sample.items():
        samples[key] = entry["sample"]
        byprime = entry["byprime"]
        for lab in labels:
            residues = {p: terms.get(lab, 0) for p, terms in byprime.items()}
            val = reconstruct_rational(residues)
            if val is None:
                raise InterpolationFailed(f"value reconstruction failed for {lab} at {key}")
            table[lab][key] = val
    return labels, table, samples, n_skipped


def interpolate_univariate(
    values: dict[Fraction, Fraction], param: str, max_num_deg: int = 12, min_validation: int = 2
) -> sp.Expr:
    """Reconstruct a single-parameter rational function, validated on held-out points."""
    xs = sorted(values)
    if len(xs) < min_validation + 2:
        raise InterpolationFailed("insufficient sample points for univariate reconstruction")
    sym = sp.Symbol(param)
    hold = xs[-min_validation:]
    fit = xs[:-min_validation]

    def _rat(v: Fraction) -> sp.Rational:
        return sp.Rational(v.numerator, v.denominator)

    data = [(_rat(x), _rat(values[x])) for x in fit]
    for degnum in range(0, min(max_num_deg, len(data) - 1) + 1):
        try:
            expr = rational_interpolate(data, degnum, X=sym)
        except Exception:
            continue
        if expr is None:
            continue
        if all(sp.simplify(expr.subs(sym, _rat(x)) - _rat(values[x])) == 0 for x in hold):
            return sp.simplify(expr)
    raise InterpolationFailed("univariate rational interpolation did not validate")


def _rat(v: Fraction) -> sp.Rational:
    return sp.Rational(v.numerator, v.denominator)


def _monomials(syms: Sequence[sp.Symbol], max_deg: int) -> list[sp.Expr]:
    """All monomials in ``syms`` of total degree ``0..max_deg`` (constant first)."""
    mons: list[sp.Expr] = [sp.Integer(1)]
    for deg in range(1, max_deg + 1):
        for combo in combinations_with_replacement(range(len(syms)), deg):
            m = sp.Integer(1)
            for idx in combo:
                m *= syms[idx]
            mons.append(m)
    return mons


def _try_rational_fit(values, syms, fit, hold, num_deg, den_deg):
    """Fit ``N/D`` with the given numerator/denominator degrees.

    Returns ``(expr, "validated")`` iff the (unique) fit validates on every held-out point,
    else ``(None, reason)`` where ``reason`` explains the rejection (Method.11b Phase C:
    ``"underdetermined (needs N fit points, have M)"`` reasons let callers compute the exact
    point deficit for a capacity-bound failure report).
    """
    num_mons = _monomials(syms, num_deg)
    den_mons = _monomials(syms, den_deg)
    needed = len(num_mons) + len(den_mons) - 1
    if len(fit) < needed:
        return None, f"underdetermined (needs {needed} fit points, have {len(fit)})"
    matrix = []
    for pt in fit:
        subs = {s: c for s, c in zip(syms, pt)}
        v = _rat(values[pt])
        row = [m.subs(subs) for m in num_mons]
        row += [-v * m.subs(subs) for m in den_mons]
        matrix.append(row)
    nullspace = sp.Matrix(matrix).nullspace()
    if len(nullspace) != 1:  # ambiguous or overdetermined -> not this degree
        return None, f"no unique nullspace solution (dim {len(nullspace)})"
    vec = nullspace[0]
    na = len(num_mons)
    numer = sum(vec[i] * num_mons[i] for i in range(na))
    denom = sum(vec[na + j] * den_mons[j] for j in range(len(den_mons)))
    if denom == 0:
        return None, "identically zero denominator"
    expr = sp.cancel(numer / denom)
    for pt in hold:  # independent validation
        subs = {s: c for s, c in zip(syms, pt)}
        if denom.subs(subs) == 0:
            return None, "holdout point on denominator zero locus"
        if sp.simplify(expr.subs(subs) - _rat(values[pt])) != 0:
            return None, "holdout validation failed"
    return sp.simplify(expr), "validated"


def interpolate_multivariate(
    values: dict[tuple, Fraction],
    params: Sequence[str],
    max_deg: int = 6,
    min_validation: int = 2,
    info: dict | None = None,
) -> sp.Expr:
    """Reconstruct a multi-parameter rational function via a dense ansatz + degree search.

    ``values`` is keyed by ordered parameter-value tuples (matching ``params``). Numerator and
    denominator degrees are searched from low to high; the first degree pair whose (unique)
    nullspace solution validates on held-out points is returned. Raises :class:`InterpolationFailed`
    if nothing validates within ``max_deg``.

    If ``info`` is given, it is filled in place with the fit/holdout split and the per-degree
    attempt log (Method.11b Phase C), including the exact point deficits of underdetermined
    degree pairs, so callers can report the minimum extra samples needed after a failure.
    """
    syms = [sp.Symbol(p) for p in params]
    pts = sorted(values, key=lambda t: tuple(map(_rat, t)))
    if len(pts) < min_validation + 2:
        raise InterpolationFailed("insufficient sample points for multivariate reconstruction")
    subs_pts = [tuple(_rat(c) for c in pt) for pt in pts]
    val_map = {sub: values[pt] for sub, pt in zip(subs_pts, pts)}
    hold = subs_pts[-min_validation:]
    fit = subs_pts[:-min_validation]
    attempts: list = []
    if info is not None:
        info.update(n_points=len(pts), n_fit=len(fit), n_hold=len(hold), attempts=attempts)
    for total in range(0, 2 * max_deg + 1):  # simplest (lowest combined degree) first
        for num_deg in range(0, min(total, max_deg) + 1):
            den_deg = total - num_deg
            if den_deg > max_deg:
                continue
            expr, status = _try_rational_fit(val_map, syms, fit, hold, num_deg, den_deg)
            if info is not None:
                attempts.append({"num_deg": num_deg, "den_deg": den_deg, "status": status})
            if expr is not None:
                if info is not None:
                    info.update(status="validated", num_deg=num_deg, den_deg=den_deg)
                return expr
    if info is not None:
        info.update(status="failed")
    raise InterpolationFailed("multivariate rational interpolation did not validate")


def reconstruct_coefficients(
    results: Iterable[NormalFormResult],
    param_names: Iterable[str],
    rank_policy: str = RANK_POLICY_MAX_RANK,
    report: dict | None = None,
) -> dict:
    """Reconstruct ``{label: C_label(params)}`` from modular normal forms (uni- or multivariate).

    Consumes only the records selected by :func:`select_records_for_reconstruction` (default:
    max-rank records; pass ``rank_policy="all"`` to restore the pre-D4.3 behaviour).

    If ``report`` is given, it is filled in place with a per-label interpolation report
    (``{str(label): info}``; see :func:`interpolate_multivariate`), including the label being
    processed when an :class:`InterpolationFailed` escapes (Method.11b Phase C).
    """
    param_names = list(param_names)
    labels, table, samples, _ = collect_value_table(results, rank_policy=rank_policy)
    if not param_names:
        raise InterpolationFailed("no parameters given for reconstruction")
    coeffs: dict = {}
    if len(param_names) == 1:
        param = param_names[0]
        for lab in labels:
            point_values = {Fraction(samples[key][param]): val for key, val in table[lab].items()}
            coeffs[lab] = interpolate_univariate(point_values, param)
            if report is not None:
                report[str(lab)] = {"status": "validated", "mode": "univariate"}
    else:
        for lab in labels:
            point_values = {
                tuple(Fraction(samples[key][p]) for p in param_names): val
                for key, val in table[lab].items()
            }
            info: dict | None = None
            if report is not None:
                info = {}
                report[str(lab)] = info
            coeffs[lab] = interpolate_multivariate(point_values, param_names, info=info)
    return coeffs


# --- Method.11c: post-reconstruction validation of pending support deviations ------------------
def _eval_coeff_exact(expr, sample: dict) -> Fraction | None:
    """Exact rational value of a reconstructed coefficient at one sample.

    Returns ``None`` when the expression has a pole at the sample (vanishing denominator).
    """
    subs = {
        sp.Symbol(str(k)): sp.Rational(Fraction(v).numerator, Fraction(v).denominator)
        for k, v in sample.items()
    }
    val = sp.cancel(sp.sympify(expr).subs(subs))
    if not val.is_rational:
        return None
    val = sp.Rational(val)
    return Fraction(int(val.p), int(val.q))


def _eval_coeff_mod(expr, sample: dict, prime: int) -> int | None:
    """``_eval_coeff_exact`` reduced mod ``prime`` (``None`` on pole or non-invertible denom)."""
    v = _eval_coeff_exact(expr, sample)
    if v is None or v.denominator % prime == 0:
        return None
    return v.numerator * pow(v.denominator, -1, prime) % prime


def validate_support_deviations(
    coeffs: dict, records: Iterable, classification: dict | None = None
) -> dict:
    """Verdict for every ``support_deviation_pending_validation`` sample (Method.11c).

    Given a candidate reconstruction ``coeffs`` (``{label: C_label(params)}``) fitted on
    stable samples only, evaluate every reconstructed coefficient at each pending sample:

    - ``GenuineSpecialZero`` -- every label *present* in a record matches the recorded
      residue mod that record's prime, and every *missing* label evaluates exactly to zero
      at the sample: the deviation is an honest coefficient zero and nothing else changed;
    - ``BasisSpecialization`` -- anything else (a present coefficient disagrees, a missing
      label does not evaluate to zero, a pole, or support outside ``coeffs``): the
      normal-form basis/pivot structure specializes at the sample; it must stay excluded
      from fitting.

    ``unstable`` samples are never upgraded (rank drops / prime disagreement stay unusable).
    Returns ``{"by_sample": {key: entry}, "summary": {...}}``.
    """
    records = list(records)
    cls = classification if classification is not None else classify_sample_supports(records)
    pending = {
        key
        for key, entry in cls["by_sample"].items()
        if entry["classification"] == SUPPORT_DEVIATION_PENDING
    }
    grouped: dict = {}
    for r in records:
        grouped.setdefault(_canonical_sample_str(r.sample), []).append(r)

    exprs = {lab: sp.sympify(expr) for lab, expr in coeffs.items()}
    by_sample: dict = {}
    summary = {
        "n_pending": len(pending),
        "n_genuine_special_zero": 0,
        "n_basis_specialization": 0,
        "genuine_special_zero_samples": [],
        "basis_specialization_samples": [],
    }
    for key in sorted(pending):
        recs = [
            r for r in grouped.get(key, []) if r.status == STATUS_REDUCED and r.formal_success
        ]
        sample = dict(recs[0].sample) if recs else {}
        seen: set = set()
        for r in recs:
            seen.update(_record_coeffs(r))
        present_mismatches: list = []
        missing_nonzero: list = []
        foreign_labels = sorted(lab for lab in seen if lab not in exprs)
        for r in recs:
            rc = _record_coeffs(r)
            for lab, expr in exprs.items():
                if lab not in rc:
                    continue
                pred = _eval_coeff_mod(expr, r.sample, r.prime)
                if pred is None or pred != rc[lab] % r.prime:
                    present_mismatches.append(
                        {
                            "label": list(lab),
                            "prime": r.prime,
                            "predicted": pred,
                            "recorded": rc[lab] % r.prime,
                        }
                    )
        for lab, expr in exprs.items():
            if lab in seen:
                continue
            v = _eval_coeff_exact(expr, sample)
            if v is None or v != 0:
                missing_nonzero.append(
                    {"label": list(lab), "value": None if v is None else str(v)}
                )
        genuine = bool(recs) and not (present_mismatches or missing_nonzero or foreign_labels)
        verdict = (
            DEVIATION_GENUINE_SPECIAL_ZERO if genuine else DEVIATION_BASIS_SPECIALIZATION
        )
        by_sample[key] = {
            "verdict": verdict,
            "n_records": len(recs),
            "n_present_mismatches": len(present_mismatches),
            "present_mismatches": present_mismatches,
            "missing_nonzero": missing_nonzero,
            "foreign_labels": [list(lab) for lab in foreign_labels],
        }
        if genuine:
            summary["n_genuine_special_zero"] += 1
            summary["genuine_special_zero_samples"].append(key)
        else:
            summary["n_basis_specialization"] += 1
            summary["basis_specialization_samples"].append(key)
    return {"by_sample": by_sample, "summary": summary}
