#!/usr/bin/env python
"""External Int2 Method.12R: mixed-boundary contradiction audit.

Method.12N found that the certified master ``L4 = [0,0,1,-1,0,0,-1]`` is NOT locally finite:
``J_4(ep) = -1/(r*ep) + O(1)`` as ``ep -> 0^-``, from the joint boundary ``x5, x7 -> Infinity`` at
fixed ``x2``. This script audits the contradiction with the recorded ``LocallyFinite -> True`` and
with the Laurent order of the certified four-term relation. Four independent phases:

A. **Symbolic Laurent consistency** — exact leading orders of the certified ``C1..C4`` and the
   ``ep^-3`` coefficient of ``C1*L1 + ... + C4*L4``, against the independently derived target
   leading behaviour ``JTarget = -2/(3*r*ep^2) + O(1/ep)`` (Method.2 leading-pole audit).
B. **LF ray completeness** — the heuristic candidate rays vs the complete polyhedral ray set, the
   missing joint-infinity ray, per-label scores on it, and the LF verdicts it changes. Includes a
   second, unrelated family (D4) so the finding is not specific to one label.
C. **Surface-row revalidation** (heavy) — every Method.10 chamber-only row re-checked from its
   provenance against the complete ray set and against the corrected coordinate facet criterion,
   plus a re-run of the recorded Method.7 obstruction-witness screening over the survivors.
D. **Direct numeric identity check** in the convergence chamber ``ep = -3/5, r = 1``, where all
   five integrals converge: ``JTarget`` vs ``C1*L1 + C2*L2 + C3*L3 + C4*L4`` by two independent
   routes (3-D graded Gauss-Legendre in float64; exact ``x7`` preintegration + 2-D graded
   Gauss-Legendre in mpmath), with refinement ladders and residuals.

NO RREF and NO modular solve is performed anywhere: phase C only generates rows and evaluates
recorded dual witnesses by a dot product. Nothing is re-derived: coefficients, labels and the
surface policy are read from the recorded artifacts.

Usage::

    python scripts/audit_external_int2_method12r.py --stages a,b,d
    python scripts/audit_external_int2_method12r.py --stages c --allow-heavy
    python scripts/audit_external_int2_method12r.py --stages all --allow-heavy
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from fractions import Fraction
from itertools import product
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from parametric_ibp_lf_reducer import (  # noqa: E402
    SurfacePolicy,
    base_score,
    complete_polyhedral_rays,
    coordinate_primitive_surface_free,
    explain_local_finiteness,
    heuristic_candidate_rays,
    is_locally_finite,
    parse_family_text,
    vector_field_surface_free,
)
from parametric_ibp_lf_reducer.row_generation import (  # noqa: E402
    coordinate_ibp_primitive_row,
    generate_algebraic_rows,
    monomials_up_to,
    tangent_ibp_primitive_row,
)
from parametric_ibp_lf_reducer.tangent_fields import generate_tangent_fields  # noqa: E402
from parametric_ibp_lf_reducer.valuations import score_from_exponents  # noqa: E402

INPUT = REPO_ROOT / "examples" / "external_int2_dimensionless_input.wl.txt"
D4_INPUT = REPO_ROOT / "examples" / "d4_cli_example_input.wl.txt"
CERT = REPO_ROOT / "validation" / "external_int2_lf_certificate.json"
METHOD10 = REPO_ROOT / "validation" / "external_int2_method10.json"
WITNESS = REPO_ROOT / "validation" / "external_int2_method7_witness.json"
DEFAULT_OUT = REPO_ROOT / "validation" / "external_int2_method12r.json"

EP = sp.Symbol("ep")
R = sp.Symbol("r")

#: convergence chamber used by the certified system (recorded in the certificate)
EP_CHAMBER = Fraction(-3, 5)

#: Method.10 Level-0 box (n ranges x m range), 3072 labels
N_RANGES = ((-1, 1), (0, 1), (0, 1))
M_RANGE = (-3, 0)
MAX_IBP_DEGREE = 2
TANGENT_BLOCKS = ((1, 1), (2, 2))
EXTRA_TANGENT_BLOCK = (3, 3)

TARGET_LABEL = (0, 0, 0, 0, 0, 0, 0)
MASTERS = {
    "L1": (-1, 0, 0, -1, -1, 0, 0),
    "L2": (-1, 0, 0, 0, -1, 0, -1),
    "L3": (0, 0, 0, -1, 0, 0, -1),
    "L4": (0, 0, 1, -1, 0, 0, -1),
}

#: joint-infinity ray (x5, x7 -> Infinity at fixed x2) in the package's ray convention
MIXED_RAY = (0, -1, -1)

#: independently derived leading behaviour of the pure target (Method.2 leading-pole audit,
#: notes/EXTERNAL_INT2_LEADING_POLE_AUDIT.md). NOT the Laurent oracle: only this leading term.
TARGET_LEADING = -2 / (3 * R * EP**2)
TARGET_LEADING_ORDER = -2


# --------------------------------------------------------------------------------------------
# helpers


def _family() -> object:
    return parse_family_text(INPUT.read_text(encoding="utf-8"))


def certified_coefficients() -> dict[str, sp.Expr]:
    """The four certified coefficients, read verbatim from the recorded certificate."""
    claimed = json.loads(CERT.read_text(encoding="utf-8"))["claimed_coefficients"]
    out = {}
    for name, label in MASTERS.items():
        key = str(list(label))
        out[name] = sp.sympify(claimed[key], locals={"ep": EP, "r": R})
    return out


def laurent_terms(expr: sp.Expr, lo: int, hi: int) -> dict[int, sp.Expr]:
    """Exact Laurent coefficients of ``expr`` in ``ep`` for orders ``lo..hi``."""
    ser = sp.series(sp.together(expr), EP, 0, max(hi + 1, 1)).removeO()
    ser = sp.expand(sp.simplify(ser))
    return {k: sp.simplify(ser.coeff(EP, k)) for k in range(lo, hi + 1)}


def leading_order(expr: sp.Expr, lo: int = -4, hi: int = 2) -> int | None:
    terms = laurent_terms(expr, lo, hi)
    for k in range(lo, hi + 1):
        if terms[k] != 0:
            return k
    return None


def _label_exponents(family, label, at_ep=None):
    """``(e_i, f_l)`` for ``label``; substituted at ``ep = at_ep`` when given."""
    e, f = family.exponent_at_label(tuple(label))
    e = [x.to_sympy() for x in e]
    f = [x.to_sympy() for x in f]
    if at_ep is not None:
        sub = {EP: sp.Rational(at_ep.numerator, at_ep.denominator)}
        e = [x.subs(sub) for x in e]
        f = [x.subs(sub) for x in f]
    return e, f


# --------------------------------------------------------------------------------------------
# Phase A: symbolic Laurent consistency


def phase_a() -> dict:
    coeffs = certified_coefficients()
    orders = {}
    for name, c in coeffs.items():
        k = leading_order(c)
        orders[name] = {
            "expr": str(c),
            "leading_order": k,
            "leading_coefficient": str(laurent_terms(c, k, k)[k]),
        }

    # master Laurent models: L1..L3 regular at ep=0, L4 has the Method.12N simple pole.
    l1, l2, l3, l40, l41 = sp.symbols("l1 l2 l3 l40 l41")
    masters = {
        "L1": l1,
        "L2": l2,
        "L3": l3,
        "L4": -1 / (R * EP) + l40 + l41 * EP,
    }
    rhs = sum(coeffs[k] * masters[k] for k in MASTERS)
    rhs_terms = laurent_terms(rhs, -3, -1)
    c4l4_terms = laurent_terms(coeffs["L4"] * masters["L4"], -3, -3)

    pole3 = sp.simplify(rhs_terms[-3])
    consistent = bool(sp.simplify(pole3) == 0)
    return {
        "coefficients": orders,
        "master_model": {k: str(v) for k, v in masters.items()},
        "rhs_ep_minus_3": str(pole3),
        "rhs_ep_minus_3_from_C4_L4_alone": str(sp.simplify(c4l4_terms[-3])),
        "rhs_ep_minus_3_depends_on_finite_parts": bool(
            pole3.free_symbols & {l1, l2, l3, l40, l41}
        ),
        "rhs_leading_order": -3 if not consistent else None,
        "target_leading": str(TARGET_LEADING),
        "target_leading_order": TARGET_LEADING_ORDER,
        "laurent_orders_match": consistent,
        "verdict": (
            "CONSISTENT"
            if consistent
            else "INCONSISTENT: RHS carries an ep^-3 pole that the target does not have"
        ),
    }


# --------------------------------------------------------------------------------------------
# Phase B: LF ray completeness


def symbolic_score(family, label, direction) -> sp.Expr:
    """Scaling score along ``direction`` keeping ``ep`` symbolic.

    ``base_score`` evaluates the exponents at ``ep = 0`` (that is the local-finiteness question);
    here the regulator is kept so the same ray can also be read inside the convergence chamber.
    """
    e, f = _label_exponents(family, label)
    return sp.expand(score_from_exponents(e, f, family, direction))


def phase_b() -> dict:
    family = _family()
    heur = [ray.direction for ray in heuristic_candidate_rays(family)]
    poly, complete = complete_polyhedral_rays(family)
    poly_dirs = [ray.direction for ray in poly]
    missing = [d for d in poly_dirs if d not in set(heur)]

    per_label = {}
    for name, label in {"Target": TARGET_LABEL, **MASTERS}.items():
        score_mixed = sp.simplify(symbolic_score(family, label, MIXED_RAY))
        per_label[name] = {
            "label": list(label),
            "score_on_mixed_ray": str(score_mixed),
            "score_on_mixed_ray_at_ep0": str(sp.simplify(base_score(family, label, MIXED_RAY))),
            "score_on_mixed_ray_in_chamber": str(
                sp.simplify(score_mixed.subs(EP, sp.Rational(-3, 5)))
            ),
            "locally_finite": _verdict_str(is_locally_finite(family, label)),
            "failing_rays": [
                {"direction": list(rv.ray.direction), "kind": rv.ray.kind, "score": str(rv.score)}
                for rv in explain_local_finiteness(family, label).failing_rays
            ],
            # convergence of the regulated integral at the chamber point: every ray of the
            # complete set must score > 0 there (this is what makes phase D well-posed).
            "chamber_failing_rays": [
                [list(d), str(s)]
                for d, s in (
                    (d, sp.simplify(symbolic_score(family, label, d).subs(EP, sp.Rational(-3, 5))))
                    for d in dict.fromkeys(heur + poly_dirs)
                )
                if not s.free_symbols and s <= 0
            ],
        }

    # genericity: an unrelated family (D4) also has a marginal label on a missing joint ray
    d4 = parse_family_text(D4_INPUT.read_text(encoding="utf-8"))
    d4_heur = {ray.direction for ray in heuristic_candidate_rays(d4)}
    d4_label = (1, 1, 1, 0, -1, 0, -1)
    d4_ray = (-1, 0, -1, 0)
    generic = {
        "family": "examples/d4_cli_example_input.wl.txt",
        "label": list(d4_label),
        "ray": list(d4_ray),
        "ray_in_heuristic_set": d4_ray in d4_heur,
        "score_at_ep0": str(sp.simplify(base_score(d4, d4_label, d4_ray).subs(EP, 0))),
        "locally_finite": _verdict_str(is_locally_finite(d4, d4_label)),
        "n_heuristic_rays": len(d4_heur),
        "n_complete_rays": len(complete_polyhedral_rays(d4)[0]),
    }

    return {
        "n_heuristic_rays": len(heur),
        "n_complete_polyhedral_rays": len(poly_dirs),
        "enumeration_complete": complete,
        "rays_missing_from_heuristic_set": [list(d) for d in missing],
        "mixed_ray": list(MIXED_RAY),
        "mixed_ray_in_heuristic_set": MIXED_RAY in set(heur),
        "per_label": per_label,
        "genericity_check": generic,
        "verdict": (
            "L4 IS NOT LOCALLY FINITE (marginal on the joint infinity ray); the heuristic ray set "
            "never tested that ray"
            if per_label["L4"]["locally_finite"] == "False"
            else "unexpected: L4 still reported locally finite"
        ),
    }


def _verdict_str(v) -> str:
    return v if isinstance(v, str) else str(v)


# --------------------------------------------------------------------------------------------
# Phase C: surface-row revalidation (heavy)


def box_labels() -> list[tuple[int, ...]]:
    n_iters = [range(lo, hi + 1) for lo, hi in N_RANGES]
    m_iter = list(product(range(M_RANGE[0], M_RANGE[1] + 1), repeat=4))
    return [tuple(list(n) + list(m)) for n in product(*n_iters) for m in m_iter]


def corrected_coordinate_surface_free(family, label, var_index, multiplier_exps, rays, policy):
    """Facet criterion: the boundary form on ``x_i = 0`` / ``x_i = inf`` must vanish on EVERY ray.

    The surface term of ``d/dx_i (P F)`` at the ``x_i`` facets is an ``(N-1)``-dimensional integral;
    it vanishes only if the integrand-times-transverse-measure is boundary-suppressed along every
    degeneration of the *remaining* variables as well. Dropping ``dx_i`` from the scaling score
    gives ``score(d) - d_i``, which must be ``> 0`` for every candidate ray with ``d_i != 0``. On
    the pure coordinate rays ``+-e_i`` this reduces exactly to the production check
    (``exp_zero > 0`` and ``exp_inf < 0``); the mixed rays are the new content.
    """
    e_syms, f_syms = _label_exponents(family, label)
    p = tuple(int(x) for x in multiplier_exps)
    e_shift = [e_syms[k] + p[k] for k in range(family.nvars)]
    for d in rays:
        if d[var_index] == 0:
            continue
        score = sp.expand(score_from_exponents(e_shift, f_syms, family, d) - d[var_index])
        # same sign semantics as the production filter, only on more rays
        if policy.sign_at(score, family.regulators) != "pos":
            return (list(d), str(sp.simplify(score)))
    return True


class FastFluxFilter:
    """Exact integer re-implementation of :func:`vector_field_surface_free` for one policy.

    The flux score is affine in the label and in the field monomial::

        score(d; label, i, c) = base_const(d) + base_ep(d)*ep
                                + sum_k d_k n_k + sum_l m_l val_d(G_l)     (label part)
                                + sum_k d_k c_k - d_i                       (field part)

    so every verdict is a strict sign test on exact rationals. This is used only to make the
    Level-0 sweep affordable; :meth:`verify` cross-checks it against the library filter on a random
    sample of (label, field) pairs and the audit aborts on any disagreement.
    """

    def __init__(self, family, rays, policy: SurfacePolicy):
        self.family = family
        self.policy = policy
        self.nvars = family.nvars
        self.npolys = len(family.poly_names)
        zero = tuple([0] * (self.nvars + self.npolys))
        e0, f0 = _label_exponents(family, zero)
        self.rays = [tuple(int(x) for x in d) for d in rays]
        self.vals = {
            d: [family.polynomials[nm].valuation(d) for nm in family.poly_names] for d in self.rays
        }
        self.base = {}
        for d in self.rays:
            expr = sp.expand(
                sum(d[k] * (e0[k] + 1) for k in range(self.nvars))
                + sum(f0[j] * self.vals[d][j] for j in range(self.npolys))
            )
            poly = sp.Poly(expr, EP)
            coeffs = poly.all_coeffs()[::-1]
            if len(coeffs) > 2 or any(c.free_symbols for c in coeffs):
                raise RuntimeError("exponents are not affine in ep with rational coefficients")
            const = Fraction(sp.Rational(coeffs[0]).p, sp.Rational(coeffs[0]).q)
            slope = (
                Fraction(sp.Rational(coeffs[1]).p, sp.Rational(coeffs[1]).q)
                if len(coeffs) > 1
                else Fraction(0)
            )
            self.base[d] = (const, slope)
        if policy.mode == "chamber":
            self.ep_value = policy.point["ep"]
        else:
            self.ep_value = None
            self.direction = policy.direction

    def _sign(self, const: Fraction, slope: Fraction) -> str:
        if self.ep_value is not None:
            v = const + slope * self.ep_value
            return "pos" if v > 0 else ("neg" if v < 0 else "zero")
        if const:
            return "pos" if const > 0 else "neg"
        if not slope:
            return "zero"
        if self.direction == "minus":
            return "pos" if slope < 0 else "neg"
        return "pos" if slope > 0 else "neg"

    def field_terms(self, field) -> list[tuple[int, tuple[int, ...]]]:
        """``(i, c)`` pairs of the vector field, i.e. one flux monomial per component."""
        out = []
        for i in range(self.nvars):
            qi = field.components[i]
            if qi.is_zero:
                continue
            for c in qi.support():
                out.append((i, tuple(int(x) for x in c)))
        return out

    def accepts(self, label, terms) -> bool:
        n = label[: self.nvars]
        m = label[self.nvars :]
        for d in self.rays:
            const0, slope = self.base[d]
            lab_part = sum(d[k] * n[k] for k in range(self.nvars)) + sum(
                m[j] * self.vals[d][j] for j in range(self.npolys)
            )
            for i, c in terms:
                const = const0 + lab_part + sum(d[k] * c[k] for k in range(self.nvars)) - d[i]
                if self._sign(const, slope) != "pos":
                    return False
        return True

    def verify(self, labels, fields, n_samples: int = 240, seed: int = 20260728) -> dict:
        import random

        rng = random.Random(seed)
        pairs = [
            (rng.choice(labels), rng.randrange(len(fields)))
            for _ in range(min(n_samples, len(labels) * len(fields)))
        ]
        agree = 0
        for label, fi in pairs:
            tf = fields[fi]
            lib = vector_field_surface_free(
                self.family, label, list(tf.components), rays=self.rays, policy=self.policy
            )
            fast = self.accepts(label, self.field_terms(tf))
            if (lib is True) != fast:
                raise SystemExit(
                    f"fast flux filter disagrees with the library filter at label={label} "
                    f"field#{fi}: library={lib!r}, fast={fast!r}"
                )
            agree += 1
        return {"pairs_checked": agree, "disagreements": 0}


def _witness_points() -> list[dict]:
    doc = json.loads(WITNESS.read_text(encoding="utf-8"))
    return [p for p in doc["points"] if p["status"] == "Witness"]


def _pair_row_with_witness(row, witness_vec: dict, sample: dict, prime: int) -> int:
    total = 0
    for label, coeff in row.terms.items():
        w = witness_vec.get(label)
        if w is None:
            continue
        v = coeff.eval_mod_p(sample, prime)
        if v is None:
            raise RuntimeError(f"coefficient denominator vanishes mod {prime}")
        total = (total + v * w) % prime
    return total % prime


def _assemble(family, labels, policy, log, flux: FastFluxFilter | None = None) -> dict:
    """Level-0 rows under ``policy``: algebraic + coordinate + tangent (blocks + extra (3,3)).

    Returns the deduped rows keyed by ``dedup_key`` together with EVERY accepted justification per
    key, so a row can be declared invalid only when *all* of its justifications fail.
    """
    t0 = time.time()
    rows: dict[frozenset, dict] = {}

    def offer(row, kind, prov):
        if row.is_trivial():
            return
        key = row.dedup_key()
        rec = rows.get(key)
        if rec is None:
            rows[key] = {"row": row, "kind": kind, "provenances": [prov]}
        else:
            rec["provenances"].append(prov)

    for row in generate_algebraic_rows(family, labels):
        offer(row, "algebraic", {"kind": "algebraic"})
    log(f"    algebraic: {len(rows)} rows ({time.time() - t0:.1f}s)")

    t = time.time()
    n_coord = 0
    mons = list(monomials_up_to(family.nvars, MAX_IBP_DEGREE))
    for label in labels:
        for var_index in range(family.nvars):
            for p in mons:
                if coordinate_primitive_surface_free(
                    family, label, var_index, multiplier_exps=p, policy=policy
                ) is not True:
                    continue
                n_coord += 1
                offer(
                    coordinate_ibp_primitive_row(family, label, var_index, p),
                    "coordinate_ibp",
                    {"kind": "coordinate_ibp", "label": label, "var": var_index, "P": p},
                )
    log(f"    coordinate: {n_coord} accepted justifications ({time.time() - t:.1f}s)")

    t = time.time()
    fields = generate_tangent_fields(family, [*TANGENT_BLOCKS, EXTRA_TANGENT_BLOCK])
    n_tan = 0
    for fi, tf in enumerate(fields):
        if not tf.is_tangent(family):
            continue
        terms = flux.field_terms(tf) if flux is not None else None
        for label in labels:
            if flux is not None:
                if not flux.accepts(label, terms):
                    continue
            elif vector_field_surface_free(
                family, label, list(tf.components), policy=policy
            ) is not True:
                continue
            n_tan += 1
            offer(
                tangent_ibp_primitive_row(family, label, tf),
                "tangent_ibp",
                {"kind": "tangent_ibp", "label": label, "field_index": fi, "field": tf},
            )
    log(
        f"    tangent: {len(fields)} fields, {n_tan} accepted justifications "
        f"({time.time() - t:.1f}s)"
    )
    log(f"    total deduped rows: {len(rows)} ({time.time() - t0:.1f}s)")
    return rows


def phase_c(log) -> dict:
    family = _family()
    labels = box_labels()
    poly_rays, complete = complete_polyhedral_rays(family)
    ray_dirs = [ray.direction for ray in poly_rays]
    heur_dirs = [ray.direction for ray in heuristic_candidate_rays(family)]
    all_dirs = list(dict.fromkeys(heur_dirs + ray_dirs))
    chamber = SurfacePolicy.chamber({"ep": EP_CHAMBER})
    limit = SurfacePolicy.limit("minus")

    log(f"[12R-C] Level-0 box: {len(labels)} labels; chamber ep={EP_CHAMBER}")
    fields = generate_tangent_fields(family, [*TANGENT_BLOCKS, EXTRA_TANGENT_BLOCK])
    flux_limit = FastFluxFilter(family, heur_dirs, limit)
    flux_chamber = FastFluxFilter(family, heur_dirs, chamber)
    # corrected filter: same criterion, complete ray set
    flux_chamber_full = FastFluxFilter(family, all_dirs, chamber)
    flux_limit_full = FastFluxFilter(family, all_dirs, limit)
    gate = {
        "limit_heuristic_rays": flux_limit.verify(labels, fields),
        "chamber_heuristic_rays": flux_chamber.verify(labels, fields),
        "chamber_complete_rays": flux_chamber_full.verify(labels, fields),
        "limit_complete_rays": flux_limit_full.verify(labels, fields),
    }
    log(f"[12R-C] fast flux filter verified against the library filter: {gate}")

    log("[12R-C] assembling limit-policy rows")
    lim_rows = _assemble(family, labels, limit, log, flux=flux_limit)
    log("[12R-C] assembling chamber-policy rows")
    cham_rows = _assemble(family, labels, chamber, log, flux=flux_chamber)

    chamber_only = {k: v for k, v in cham_rows.items() if k not in lim_rows}
    limit_only = {k: v for k, v in lim_rows.items() if k not in cham_rows}
    by_kind: dict[str, int] = {}
    for rec in chamber_only.values():
        by_kind[rec["kind"]] = by_kind.get(rec["kind"], 0) + 1
    log(f"[12R-C] chamber_only={len(chamber_only)} {by_kind}; limit_only={len(limit_only)}")

    def revalidate(records, flux_full, policy):
        """Re-check every justification of every row under the corrected criteria."""
        still, rejected = [], []
        counters = {"coordinate_ibp": [0, 0], "tangent_ibp": [0, 0], "algebraic": [0, 0]}
        cache: dict = {}
        for key, rec in records.items():
            ok = False
            witness = None
            for prov in rec["provenances"]:
                if prov["kind"] == "algebraic":
                    ok = True  # exact identity, no surface term
                    break
                if prov["kind"] == "coordinate_ibp":
                    ck = (prov["label"], prov["var"], prov["P"])
                    if ck not in cache:
                        cache[ck] = corrected_coordinate_surface_free(
                            family, prov["label"], prov["var"], prov["P"], all_dirs, policy
                        )
                    res = cache[ck]
                else:
                    ck = (prov["label"], prov["field_index"])
                    if ck not in cache:
                        tf = prov["field"]
                        cache[ck] = (
                            True
                            if flux_full.accepts(prov["label"], flux_full.field_terms(tf))
                            else (None, "flux not suppressed on the complete ray set")
                        )
                    res = cache[ck]
                if res is True:
                    ok = True
                    break
                witness = res
            counters[rec["kind"]][0 if ok else 1] += 1
            (still if ok else rejected).append((key, rec, witness))
        return still, rejected, counters

    log("[12R-C] revalidating chamber-only rows against the complete ray set")
    still, rejected, counters = revalidate(chamber_only, flux_chamber_full, chamber)
    log(f"[12R-C] chamber-only: still valid {len(still)}, newly rejected {len(rejected)}")

    log("[12R-C] revalidating the production (limit-policy) rows")
    lim_still, lim_rejected, lim_counters = revalidate(lim_rows, flux_limit_full, limit)
    log(f"[12R-C] limit-policy: still valid {len(lim_still)}, newly rejected {len(lim_rejected)}")

    # witness re-screening over the chamber-only rows that survive the correction
    witness_report = []
    for wp in _witness_points():
        wvec = {tuple(lab): int(v) for lab, v in wp["witness"]}
        sample = {k: Fraction(v) for k, v in wp["sample"]}
        prime = int(wp["prime"])
        n_break_all = n_break_valid = 0
        for key, rec, _w in still:
            if _pair_row_with_witness(rec["row"], wvec, sample, prime):
                n_break_valid += 1
        for key, rec in ((k, r) for k, r in chamber_only.items()):
            if _pair_row_with_witness(rec["row"], wvec, sample, prime):
                n_break_all += 1
        witness_report.append(
            {
                "prime": prime,
                "sample": [list(pair) for pair in wp["sample"]],
                "n_chamber_only": len(chamber_only),
                "n_breaks_recorded_rows": n_break_all,
                "n_breaks_surviving_rows": n_break_valid,
            }
        )
        log(
            f"[12R-C] witness p={prime} sample={wp['sample']}: breaks {n_break_all} -> "
            f"{n_break_valid} after the correction"
        )

    recorded = json.loads(METHOD10.read_text(encoding="utf-8"))["rows"]
    return {
        "ray_enumeration_complete": complete,
        "fast_flux_filter_gate": gate,
        "n_rays_used": len(all_dirs),
        "rows": {
            "n_limit": len(lim_rows),
            "n_chamber": len(cham_rows),
            "n_chamber_only": len(chamber_only),
            "n_limit_only": len(limit_only),
            "chamber_only_by_kind": by_kind,
            "recorded_method10": recorded,
            "reproduces_recorded_chamber_only": len(chamber_only) == recorded["n_chamber_only"],
        },
        "chamber_only_revalidation": {
            "still_surface_valid": len(still),
            "newly_rejected": len(rejected),
            "by_kind": {
                k: {"still_valid": v[0], "newly_rejected": v[1]} for k, v in counters.items()
            },
            "examples_rejected": [
                {
                    "kind": rec["kind"],
                    "provenance": _prov_payload(rec["provenances"][0]),
                    "witness_ray": w[0] if w else None,
                    "witness_score": w[1] if w else None,
                }
                for _k, rec, w in rejected[:8]
            ],
        },
        "production_limit_policy_revalidation": {
            "still_surface_valid": len(lim_still),
            "newly_rejected": len(lim_rejected),
            "by_kind": {
                k: {"still_valid": v[0], "newly_rejected": v[1]} for k, v in lim_counters.items()
            },
        },
        "witness_rescreening": witness_report,
        "note": (
            "A row is declared invalid only when EVERY accepted justification of it fails the "
            "corrected criteria. No row is added to or removed from any production artifact here."
        ),
    }


def _prov_payload(prov: dict) -> dict:
    out = {k: v for k, v in prov.items() if k != "field"}
    if "label" in out:
        out["label"] = list(out["label"])
    if "P" in out:
        out["P"] = list(out["P"])
    return out


# --------------------------------------------------------------------------------------------
# Phase D: direct numeric identity check in the convergence chamber

EP_NUM = -0.6
R_NUM = 1.0


def _graded_nodes(order: int, depth: int):
    """Gauss-Legendre nodes on dyadically graded panels of (0,1), refined toward both ends."""
    xg, wg = np.polynomial.legendre.leggauss(order)
    cuts = sorted(
        {0.0, 1.0}
        | {0.5 ** (depth - k) for k in range(depth + 1)}
        | {1 - 0.5 ** (depth - k) for k in range(depth + 1)}
    )
    xs, ws = [], []
    for a, b in zip(cuts[:-1], cuts[1:]):
        xs.append(0.5 * (b - a) * xg + 0.5 * (a + b))
        ws.append(0.5 * (b - a) * wg)
    return np.concatenate(xs), np.concatenate(ws)


def _integrands_3d(X2, X5, X7, ep, r):
    G0, G1, G2 = 1 + X2, 1 + X5, 1 + X7
    G3 = 1 + X7 + X2 * X7 + r * X2 * X5

    def f(a2, b0, b1, b2, b3, extra=1.0):
        return extra * X2**a2 * G0**b0 * G1**b1 * G2**b2 * G3**b3

    return {
        "Target": f(1 + ep, ep, ep, -1 - ep, -1 + ep),
        "L1": f(ep, ep - 1, ep - 1, -ep - 1, ep - 1),
        "L2": f(ep, ep, ep - 1, -ep - 1, ep - 2),
        "L3": f(ep + 1, ep - 1, ep, -ep - 1, ep - 2),
        "L4": f(ep + 1, ep - 1, ep, -ep - 1, ep - 2, extra=X7),
    }


def route_3d_tensor(ladder, ep=EP_NUM, r=R_NUM) -> dict:
    """Route A: 3-D graded Gauss-Legendre tensor quadrature in float64.

    The x2 direction is looped over rather than broadcast, so memory stays at one 2-D slice
    regardless of the ladder level.
    """
    out: dict[str, list[float]] = {}
    levels = []
    for order, depth in ladder:
        u, w = _graded_nodes(order, depth)
        x = u / (1 - u)
        wj = w / (1 - u) ** 2
        X5 = x[:, None]
        X7 = x[None, :]
        W2 = wj[:, None] * wj[None, :]
        acc: dict[str, float] = {}
        for x2, w2 in zip(x, wj):
            for k, v in _integrands_3d(x2, X5, X7, ep, r).items():
                acc[k] = acc.get(k, 0.0) + w2 * float(np.sum(v * W2))
        vals = {k: float(v) for k, v in acc.items()}
        for k, v in vals.items():
            out.setdefault(k, []).append(v)
        levels.append({"order": order, "depth": depth, "n_1d_nodes": int(len(u)), "values": vals})
    return {"levels": levels, "series": out}


def _x7_closed_forms():
    """Exact ``x7`` primitives, derived once with sympy (verified numerically before use).

    ``F(A,B) = Int[(1+x7)^(-1-ep) (A+B x7)^(ep-1), {x7,0,Inf}] = (B^ep - A^ep)/(ep (B-A))``
    (Method.2 leading-pole audit, G&R). Differentiating in ``A`` gives the ``G3^(ep-2)`` cases, and
    ``x7 = ((A + B x7) - A)/B`` gives the ``x7`` numerator case — all elementary.
    """
    A, B, e = sp.symbols("A B e")
    F = (B**e - A**e) / (e * (B - A))
    F2 = sp.diff(F, A) / (e - 1)  # Int[(1+x7)^(-1-ep) (A+B x7)^(ep-2)]
    F4 = (F - A * F2) / B  # Int[x7 (1+x7)^(-1-ep) (A+B x7)^(ep-2)]
    return A, B, e, sp.simplify(F), sp.simplify(F2), sp.simplify(F4)


#: relative tolerance for the closed-form verification (adaptive quadrature on a semi-infinite
#: interval with an algebraic tail is the weaker of the two sides here)
X7_FORM_TOLERANCE = 1e-9


def _verify_x7_forms(ep=EP_NUM, r=R_NUM, dps=30) -> list[dict]:
    """Numerically verify the three closed forms against adaptive quadrature."""
    A, B, e, F, F2, F4 = _x7_closed_forms()
    closed_fns = {
        "F": sp.lambdify((A, B, e), F, "mpmath"),
        "F2": sp.lambdify((A, B, e), F2, "mpmath"),
        "F4": sp.lambdify((A, B, e), F4, "mpmath"),
    }
    checks = []
    with mp.workdps(dps):
        epm = mp.mpf(ep)
        for x2f, x5f in ((0.3, 0.7), (2.0, 0.5), (0.1, 5.0), (4.0, 1 / r + 1e-3)):
            x2, x5 = mp.mpf(x2f), mp.mpf(x5f)
            a = 1 + mp.mpf(r) * x2 * x5
            b = 1 + x2
            num = {
                "F": mp.quad(
                    lambda t: (1 + t) ** (-1 - epm) * (a + b * t) ** (epm - 1),
                    [0, 1, mp.inf],
                    maxdegree=8,
                ),
                "F2": mp.quad(
                    lambda t: (1 + t) ** (-1 - epm) * (a + b * t) ** (epm - 2),
                    [0, 1, mp.inf],
                    maxdegree=8,
                ),
                "F4": mp.quad(
                    lambda t: t * (1 + t) ** (-1 - epm) * (a + b * t) ** (epm - 2),
                    [0, 1, mp.inf],
                    maxdegree=8,
                ),
            }
            for key, fn in closed_fns.items():
                rel = abs(mp.mpf(fn(a, b, epm)) - num[key]) / abs(num[key])
                checks.append(
                    {
                        "x2": x2f,
                        "x5": x5f,
                        "form": key,
                        "relative_difference": float(rel),
                        "ok": bool(rel < X7_FORM_TOLERANCE),
                    }
                )
    return checks


def route_2d_preintegrated(ladder, ep=EP_NUM, r=R_NUM, dps=30) -> dict:
    """Route B: exact ``x7`` preintegration, then 2-D graded Gauss-Legendre in mpmath."""
    A, B, e, F, F2, F4 = _x7_closed_forms()
    fF = sp.lambdify((A, B, e), F, "mpmath")
    fF2 = sp.lambdify((A, B, e), F2, "mpmath")
    fF4 = sp.lambdify((A, B, e), F4, "mpmath")

    levels = []
    series: dict[str, list[float]] = {}
    with mp.workdps(dps):
        epm = mp.mpf(ep)
        rm = mp.mpf(r)
        for order, depth in ladder:
            u, w = _graded_nodes(order, depth)
            nodes = [(mp.mpf(ui) / (1 - mp.mpf(ui)), mp.mpf(wi) / (1 - mp.mpf(ui)) ** 2)
                     for ui, wi in zip(u, w)]
            acc = dict.fromkeys(("Target", "L1", "L2", "L3", "L4"), mp.mpf(0))
            for x2, w2 in nodes:
                g0 = 1 + x2
                for x5, w5 in nodes:
                    g1 = 1 + x5
                    a = 1 + rm * x2 * x5
                    b = g0
                    wt = w2 * w5
                    i1 = fF(a, b, epm)
                    i2 = fF2(a, b, epm)
                    i4 = fF4(a, b, epm)
                    acc["Target"] += wt * x2 ** (1 + epm) * g0**epm * g1**epm * i1
                    acc["L1"] += wt * x2**epm * g0 ** (epm - 1) * g1 ** (epm - 1) * i1
                    acc["L2"] += wt * x2**epm * g0**epm * g1 ** (epm - 1) * i2
                    acc["L3"] += wt * x2 ** (1 + epm) * g0 ** (epm - 1) * g1**epm * i2
                    acc["L4"] += wt * x2 ** (1 + epm) * g0 ** (epm - 1) * g1**epm * i4
            vals = {k: float(v) for k, v in acc.items()}
            for k, v in vals.items():
                series.setdefault(k, []).append(v)
            levels.append(
                {"order": order, "depth": depth, "n_1d_nodes": int(len(u)), "values": vals}
            )
    return {"levels": levels, "series": series, "dps": dps}


def phase_d(ladder_3d, ladder_2d, dps, log) -> dict:
    ep_exact = sp.Rational(-3, 5)
    r_exact = sp.Integer(1)
    coeffs = certified_coefficients()
    c_exact = {k: sp.nsimplify(v.subs({EP: ep_exact, R: r_exact})) for k, v in coeffs.items()}
    c_float = {k: float(v) for k, v in c_exact.items()}
    log(f"[12R-D] exact coefficients at ep=-3/5, r=1: "
        f"{ {k: str(v) for k, v in c_exact.items()} }")

    log("[12R-D] route A: 3-D graded Gauss-Legendre (float64)")
    t = time.time()
    route_a = route_3d_tensor(ladder_3d)
    log(f"[12R-D]   done ({time.time() - t:.1f}s): "
        f"{ {k: f'{v[-1]:.8f}' for k, v in route_a['series'].items()} }")

    log("[12R-D] verifying the exact x7 primitives against adaptive quadrature")
    x7_checks = _verify_x7_forms()
    if not all(c["ok"] for c in x7_checks):
        raise SystemExit("x7 closed forms failed their numeric verification; refusing to continue")

    log(f"[12R-D] route B: exact x7 preintegration + 2-D graded GL (mpmath dps={dps})")
    t = time.time()
    route_b = route_2d_preintegrated(ladder_2d, dps=dps)
    log(f"[12R-D]   done ({time.time() - t:.1f}s): "
        f"{ {k: f'{v[-1]:.8f}' for k, v in route_b['series'].items()} }")

    def summarize(series):
        out = {}
        for k, vals in series.items():
            inc = abs(vals[-1] - vals[-2]) if len(vals) > 1 else None
            out[k] = {
                "value": vals[-1],
                "ladder": vals,
                "last_increment": inc,
                "rel_last_increment": (inc / abs(vals[-1])) if inc is not None else None,
            }
        return out

    sum_a = summarize(route_a["series"])
    sum_b = summarize(route_b["series"])

    per_integral = {}
    for k in ("Target", "L1", "L2", "L3", "L4"):
        va, vb = sum_a[k]["value"], sum_b[k]["value"]
        per_integral[k] = {
            "route_3d_tensor": sum_a[k],
            "route_2d_preintegrated": sum_b[k],
            "route_agreement_abs": abs(va - vb),
            "route_agreement_rel": abs(va - vb) / max(abs(va), abs(vb)),
        }

    results = {}
    for tag, summ in (("route_3d_tensor", sum_a), ("route_2d_preintegrated", sum_b)):
        rhs = sum(c_float[k] * summ[k]["value"] for k in MASTERS)
        target = summ["Target"]["value"]
        # error propagated from the ladder increments of every term
        err = (summ["Target"]["last_increment"] or 0.0) + sum(
            abs(c_float[k]) * (summ[k]["last_increment"] or 0.0) for k in MASTERS
        )
        results[tag] = {
            "target": target,
            "rhs": rhs,
            "residual_abs": target - rhs,
            "residual_rel": (target - rhs) / target,
            "quadrature_error_estimate": err,
            "residual_over_error": abs(target - rhs) / err if err else None,
        }

    relations_hold = all(
        abs(v["residual_abs"]) <= 10 * (v["quadrature_error_estimate"] or 0.0)
        for v in results.values()
    )
    return {
        "point": {"ep": "-3/5", "r": "1"},
        "coefficients_exact": {k: str(v) for k, v in c_exact.items()},
        "coefficients_float": c_float,
        "x7_closed_form_checks": x7_checks,
        "per_integral": per_integral,
        "identity": results,
        "relation_holds_numerically": relations_hold,
        "verdict": (
            "RELATION HOLDS numerically within the quadrature error"
            if relations_hold
            else "RELATION FAILS numerically: target and RHS differ far beyond every error estimate"
        ),
    }


# --------------------------------------------------------------------------------------------
# main


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--stages", default="a,b,d", help="comma list of a,b,c,d or 'all'")
    ap.add_argument("--allow-heavy", action="store_true", help="required for stage c")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--dps", type=int, default=30, help="mpmath precision for route B")
    ap.add_argument("--fast", action="store_true", help="short quadrature ladders (smoke run)")
    args = ap.parse_args(argv)

    stages = ("a", "b", "c", "d") if args.stages == "all" else tuple(
        s.strip().lower() for s in args.stages.split(",") if s.strip()
    )
    if "c" in stages and not args.allow_heavy:
        ap.error("stage c re-assembles the Level-0 row system; pass --allow-heavy")

    def log(msg: str) -> None:
        print(msg, flush=True)

    t0 = time.time()
    report: dict = {
        "script": "scripts/audit_external_int2_method12r.py",
        "method": "External Int2 Method.12R: mixed-boundary contradiction audit",
        "scope_note": (
            "Audit only. No RREF, no modular solve, no new records, no HyperInt, no re-derivation "
            "of the reduction. Coefficients, labels, surface policy and the Method.10/Method.7 "
            "artifacts are read as recorded. Phase B changes the local-finiteness GATE only "
            "(complete polyhedral ray set); production row generation is untouched and phase C "
            "reports what correcting it would remove."
        ),
        "stages_run": list(stages),
        "chamber_ep": str(EP_CHAMBER),
    }
    if "a" in stages:
        log("[12R-A] symbolic Laurent consistency")
        report["phase_a_laurent_consistency"] = phase_a()
        log(f"[12R-A] {report['phase_a_laurent_consistency']['verdict']}")
    if "b" in stages:
        log("[12R-B] LF ray completeness")
        report["phase_b_ray_completeness"] = phase_b()
        log(f"[12R-B] {report['phase_b_ray_completeness']['verdict']}")
    if "c" in stages:
        log("[12R-C] surface-row revalidation (heavy)")
        report["phase_c_row_revalidation"] = phase_c(log)
    if "d" in stages:
        log("[12R-D] direct numeric identity check at ep=-3/5, r=1")
        ladder_3d = ((12, 6), (16, 8)) if args.fast else ((12, 6), (16, 8), (20, 10), (24, 12))
        ladder_2d = ((8, 4), (10, 6)) if args.fast else ((10, 6), (12, 8), (16, 10))
        report["phase_d_numeric_identity"] = phase_d(ladder_3d, ladder_2d, args.dps, log)
        log(f"[12R-D] {report['phase_d_numeric_identity']['verdict']}")

    verdicts = []
    if "a" in stages and not report["phase_a_laurent_consistency"]["laurent_orders_match"]:
        verdicts.append("laurent_order_contradiction")
    if "b" in stages and report["phase_b_ray_completeness"]["per_label"]["L4"][
        "locally_finite"
    ] == "False":
        verdicts.append("lf_gate_ray_gap")
    if "d" in stages and not report["phase_d_numeric_identity"]["relation_holds_numerically"]:
        verdicts.append("numeric_identity_refuted")
    report["findings"] = verdicts
    report["verdict"] = (
        "FOUR_TERM_RELATION_INVALID" if verdicts else "NO_CONTRADICTION_FOUND"
    )
    report["elapsed_seconds"] = round(time.time() - t0, 1)

    if args.out:
        merged = report
        if args.out.exists():
            try:
                old = json.loads(args.out.read_text(encoding="utf-8"))
                merged = {**old, **report}
                merged["stages_run"] = sorted(set(old.get("stages_run", [])) | set(stages))
                merged["findings"] = sorted(set(old.get("findings", [])) | set(verdicts))
                merged["verdict"] = (
                    "FOUR_TERM_RELATION_INVALID"
                    if merged["findings"]
                    else "NO_CONTRADICTION_FOUND"
                )
            except json.JSONDecodeError:
                pass
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(merged, indent=1, sort_keys=True) + "\n", encoding="utf-8")
        log(f"[12R] wrote {args.out} ({report['elapsed_seconds']}s)")
        report = merged  # report the merged verdict, not this run's stages alone
    log(f"[12R] verdict: {report['verdict']} findings={report['findings']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
