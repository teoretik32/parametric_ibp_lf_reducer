#!/usr/bin/env python
"""External Int2 Method.13: complete mixed-boundary surface validation audit.

Stages (``--stages int2,controls`` by default; ``preflight`` needs ``--allow-heavy``):

- ``int2``    — Phase D: rebuild (NOT solve) the Method.10 Level-0 candidate rows under the OLD
  heuristic surface policy (replicated in-script: component-local coordinate test = the corrected
  facet criterion restricted to ``+-e_i``; per-term all-components flux on the heuristic rays)
  and under the CORRECTED complete policy (the production filters, cross-verified against them).
  Reports totals, per-kind breakdowns, old rows invalidated, rows gained (exact flux cancellation
  admits fields the old per-term test rejected), chamber-only survival, witness-breaking rows
  remaining valid, and the comparison with the Method.12R preliminary figures.
- ``preflight`` — Phase D item 16: cheap generic LF-feasibility of the corrected chamber row set
  (1 prime x 3 samples, ``lf_reduction_feasible_mod_p``). This is the ONLY RREF in the audit and
  the gate for any future reconstruction.
- ``controls`` — Phase E: re-audit of every tracked Success artifact (D4 CLI, notebook examples
  1/2, corrected Example 4*, External Int1): master/target LF verdicts under complete rays, row
  survival under the corrected criteria, certificate applicability, independent evidence;
  writes ``validation/surface_reaudit_manifest.json``.

Criteria are the ones derived in ``docs/TORIC_SURFACE_VALIDATION.md``. Exact arithmetic only.
No new LF search, no reconstruction, no row-set changes to any committed artifact.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import random
import sys
import time
from fractions import Fraction
from itertools import product
from pathlib import Path

import sympy as sp

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

from parametric_ibp_lf_reducer import (  # noqa: E402
    SurfacePolicy,
    boundary_ray_data,
    coordinate_primitive_surface_free,
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

EP = sp.Symbol("ep")

INT2_INPUT = REPO_ROOT / "examples" / "external_int2_dimensionless_input.wl.txt"
METHOD10 = REPO_ROOT / "validation" / "external_int2_method10.json"
METHOD12R = REPO_ROOT / "validation" / "external_int2_method12r.json"
WITNESS = REPO_ROOT / "validation" / "external_int2_method7_witness.json"
AUDIT_OUT = REPO_ROOT / "validation" / "external_int2_method13_surface_audit.json"
MANIFEST_OUT = REPO_ROOT / "validation" / "surface_reaudit_manifest.json"

EP_CHAMBER = Fraction(-3, 5)
INT2_N_RANGES = ((-1, 1), (0, 1), (0, 1))
INT2_M_RANGE = (-3, 0)
MAX_IBP_DEGREE = 2
TANGENT_BLOCKS = ((1, 1), (2, 2), (3, 3))  # baseline blocks + the Method.10 extra block


# --------------------------------------------------------------------------------------------
# exact score engine (shared by old and corrected criteria)


class ScoreEngine:
    """Exact per-family score machinery: sign of ``S0(d) + integer shift`` under a policy.

    ``S0(d)`` (zero-label full-measure score) must be affine in the single regulator with
    rational coefficients — true for every family in this audit; anything else raises.
    """

    def __init__(self, family):
        self.family = family
        self.nvars = family.nvars
        self.npolys = family.npolys
        dirs, complete = boundary_ray_data(family)
        self.complete_dirs = list(dirs)
        self.complete = complete
        self.heuristic_dirs = [r.direction for r in heuristic_candidate_rays(family)]
        self.coord_dirs = [
            tuple(s if k == i else 0 for k in range(self.nvars))
            for i in range(self.nvars)
            for s in (1, -1)
        ]
        self.all_dirs = list(dict.fromkeys(self.complete_dirs + self.heuristic_dirs))
        zero = tuple([0] * (self.nvars + self.npolys))
        e0, f0 = family.exponent_at_label(zero)
        e0 = [x.to_sympy() for x in e0]
        f0 = [x.to_sympy() for x in f0]
        self.vals = {}
        self.base = {}
        for d in self.all_dirs:
            self.vals[d] = tuple(
                family.polynomials[name].valuation(d) for name in family.poly_names
            )
            expr = sp.expand(score_from_exponents(e0, f0, family, d))
            poly = sp.Poly(expr, EP)
            coeffs = poly.all_coeffs()[::-1]
            if len(coeffs) > 2 or any(c.free_symbols for c in coeffs):
                raise RuntimeError(f"S0 not affine-rational in ep on ray {d}: {expr}")
            c0 = Fraction(sp.Rational(coeffs[0]).p, sp.Rational(coeffs[0]).q)
            c1 = (
                Fraction(sp.Rational(coeffs[1]).p, sp.Rational(coeffs[1]).q)
                if len(coeffs) > 1
                else Fraction(0)
            )
            self.base[d] = (c0, c1)

    def sign(self, pol: SurfacePolicy, d, shift: int) -> str:
        c0, c1 = self.base[d]
        if pol.mode == "chamber":
            v = c0 + shift + c1 * pol.point["ep"]
            return "pos" if v > 0 else ("neg" if v < 0 else "zero")
        v0 = c0 + shift
        if v0 > 0:
            return "pos"
        if v0 < 0:
            return "neg"
        if c1 == 0:
            return "zero"
        if pol.direction == "minus":
            return "pos" if c1 < 0 else "neg"
        return "pos" if c1 > 0 else "neg"

    def label_shift(self, d, label) -> int:
        n = label[: self.nvars]
        m = label[self.nvars :]
        return sum(dk * nk for dk, nk in zip(d, n)) + sum(
            v * ml for v, ml in zip(self.vals[d], m)
        )

    # --- coordinate criteria -----------------------------------------------------------------
    def coordinate_ok(self, pol, label, var, p, dirs) -> bool:
        """Facet criterion ``score(PF, d) - d_i > 0`` over ``dirs`` with ``d_i != 0``.

        With ``dirs = coord_dirs`` this is EXACTLY the legacy component-local test
        (docs/TORIC_SURFACE_VALIDATION.md §2 sanity check); with ``complete_dirs`` it is the
        corrected Method.13 criterion.
        """
        for d in dirs:
            if d[var] == 0:
                continue
            shift = self.label_shift(d, label) + sum(dk * pk for dk, pk in zip(d, p)) - d[var]
            if self.sign(pol, d, shift) != "pos":
                return False
        return True

    # --- tangent criteria --------------------------------------------------------------------
    def field_terms(self, field):
        out = []
        for i in range(self.nvars):
            qi = field.components[i]
            if qi.is_zero:
                continue
            for c, coeff in qi.terms.items():
                out.append((i, tuple(int(x) for x in c), coeff))
        return out

    def flux_shifts_old(self, terms, d):
        """Legacy per-term shifts: every component, every monomial, no cancellation."""
        return [
            sum(dk * ck for dk, ck in zip(d, c)) - d[i] for i, c, _coeff in terms
        ]

    def flux_shifts_new(self, terms, d):
        """Corrected: exact normal component with cancellation; tangential terms dropped."""
        acc: dict[tuple, object] = {}
        for i, c, coeff in terms:
            if d[i] == 0:
                continue
            mono = tuple(ck - (1 if k == i else 0) for k, ck in enumerate(c))
            contrib = coeff.scale_int(d[i])
            acc[mono] = acc[mono] + contrib if mono in acc else contrib
        return [
            sum(dk * mk for dk, mk in zip(d, mono))
            for mono, coeff in acc.items()
            if not coeff.is_zero
        ]

    def tangent_ok(self, pol, label, shifts_by_dir) -> bool:
        for d, shifts in shifts_by_dir:
            base = self.label_shift(d, label)
            for s in shifts:
                if self.sign(pol, d, base + s) != "pos":
                    return False
        return True


# --------------------------------------------------------------------------------------------
# row assembly under a chosen criterion


def box_labels(n_ranges, m_range, npolys):
    n_iters = [range(lo, hi + 1) for lo, hi in n_ranges]
    m_iter = list(product(range(m_range[0], m_range[1] + 1), repeat=npolys))
    return [tuple(list(n) + list(m)) for n in product(*n_iters) for m in m_iter]


def assemble(family, engine: ScoreEngine, labels, pol, criterion: str, log, tangent_blocks):
    """Level rows under ``criterion`` in {"old", "new"}; returns key -> record with provenances."""
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
    n_alg = len(rows)

    coord_dirs = engine.coord_dirs if criterion == "old" else engine.complete_dirs
    n_coord = 0
    mons = list(monomials_up_to(family.nvars, MAX_IBP_DEGREE))
    for label in labels:
        for var in range(family.nvars):
            for p in mons:
                if not engine.coordinate_ok(pol, label, var, p, coord_dirs):
                    continue
                n_coord += 1
                offer(
                    coordinate_ibp_primitive_row(family, label, var, p),
                    "coordinate_ibp",
                    {"kind": "coordinate_ibp", "label": label, "var": var, "P": p},
                )

    fields = generate_tangent_fields(family, list(tangent_blocks))
    n_tan = 0
    flux_dirs = engine.heuristic_dirs if criterion == "old" else engine.complete_dirs
    for fi, tf in enumerate(fields):
        if not tf.is_tangent(family):
            continue
        terms = engine.field_terms(tf)
        shifts_fn = engine.flux_shifts_old if criterion == "old" else engine.flux_shifts_new
        shifts_by_dir = [(d, shifts_fn(terms, d)) for d in flux_dirs]
        for label in labels:
            if not engine.tangent_ok(pol, label, shifts_by_dir):
                continue
            n_tan += 1
            offer(
                tangent_ibp_primitive_row(family, label, tf),
                "tangent_ibp",
                {"kind": "tangent_ibp", "label": label, "field_index": fi},
            )
    by_kind: dict[str, int] = {}
    for rec in rows.values():
        by_kind[rec["kind"]] = by_kind.get(rec["kind"], 0) + 1
    log(
        f"    [{criterion}] rows={len(rows)} {by_kind} "
        f"(alg {n_alg}, coord just. {n_coord}, tangent just. {n_tan}; {time.time() - t0:.1f}s)"
    )
    return rows, by_kind


def verify_engine_against_library(family, engine, labels, pol, n_samples, rng):
    """The corrected in-script criteria must agree with the production filters exactly."""
    mons = list(monomials_up_to(family.nvars, MAX_IBP_DEGREE))
    fields = [tf for tf in generate_tangent_fields(family, list(TANGENT_BLOCKS))
              if tf.is_tangent(family)]
    n_c = n_t = 0
    for _ in range(n_samples):
        label = rng.choice(labels)
        var = rng.randrange(family.nvars)
        p = rng.choice(mons)
        lib = coordinate_primitive_surface_free(family, label, var, p, policy=pol)
        fast = engine.coordinate_ok(pol, label, var, p, engine.complete_dirs)
        if (lib is True) != fast:
            raise SystemExit(f"coordinate engine mismatch at {label} var={var} P={p}: "
                             f"lib={lib!r} fast={fast!r}")
        n_c += 1
        if fields:
            tf = rng.choice(fields)
            lib = vector_field_surface_free(family, label, list(tf.components), policy=pol)
            terms = engine.field_terms(tf)
            sbd = [(d, engine.flux_shifts_new(terms, d)) for d in engine.complete_dirs]
            fast = engine.tangent_ok(pol, label, sbd)
            if (lib is True) != fast:
                raise SystemExit(f"tangent engine mismatch at {label}: lib={lib!r} fast={fast!r}")
            n_t += 1
    return {"coordinate_pairs": n_c, "tangent_pairs": n_t, "disagreements": 0}


# --------------------------------------------------------------------------------------------
# Phase D: Int2 Level-0 re-audit


def _witness_points():
    doc = json.loads(WITNESS.read_text(encoding="utf-8"))
    return [p for p in doc["points"] if p["status"] == "Witness"]


def _pair(row, wvec, sample, prime) -> int:
    total = 0
    for label, coeff in row.terms.items():
        w = wvec.get(label)
        if w is None:
            continue
        v = coeff.eval_mod_p(sample, prime)
        if v is None:
            raise RuntimeError("bad specialization in witness pairing")
        total = (total + v * w) % prime
    return total


def stage_int2(log) -> dict:
    family = parse_family_text(INT2_INPUT.read_text(encoding="utf-8"))
    engine = ScoreEngine(family)
    labels = box_labels(INT2_N_RANGES, INT2_M_RANGE, family.npolys)
    chamber = SurfacePolicy.chamber({"ep": EP_CHAMBER})
    limit = SurfacePolicy.limit("minus")
    log(f"[13-D] Level-0 box: {len(labels)} labels; complete rays: "
        f"{len(engine.complete_dirs)} (complete={engine.complete})")

    rng = random.Random(20260728)
    gate = {
        "limit": verify_engine_against_library(family, engine, labels, limit, 150, rng),
        "chamber": verify_engine_against_library(family, engine, labels, chamber, 150, rng),
    }
    log(f"[13-D] engine verified against production filters: {gate}")

    out: dict = {"engine_gate": gate, "n_labels": len(labels),
                 "n_complete_rays": len(engine.complete_dirs)}
    sets: dict = {}
    for pname, pol in (("limit", limit), ("chamber", chamber)):
        for crit in ("old", "new"):
            log(f"[13-D] assembling {pname}/{crit}")
            rows, by_kind = assemble(family, engine, labels, pol, crit, log, TANGENT_BLOCKS)
            sets[(pname, crit)] = rows
            out[f"rows_{pname}_{crit}"] = {"total": len(rows), "by_kind": by_kind}

    recorded = json.loads(METHOD10.read_text(encoding="utf-8"))["rows"]
    out["method10_recorded"] = recorded
    out["old_counts_reproduce_method10"] = (
        len(sets[("limit", "old")]) == recorded["n_limit"]
        and len(sets[("chamber", "old")]) == recorded["n_chamber"]
    )

    def key_delta(a, b):
        ka, kb = set(a), set(b)
        return ka - kb, kb - ka

    inval_limit, gained_limit = key_delta(sets[("limit", "old")], sets[("limit", "new")])
    inval_cham, gained_cham = key_delta(sets[("chamber", "old")], sets[("chamber", "new")])

    def kind_count(keys, src):
        d: dict[str, int] = {}
        for k in keys:
            kind = src[k]["kind"]
            d[kind] = d.get(kind, 0) + 1
        return d

    out["limit_old_invalidated"] = {
        "n": len(inval_limit), "by_kind": kind_count(inval_limit, sets[("limit", "old")])}
    out["limit_new_gained"] = {
        "n": len(gained_limit), "by_kind": kind_count(gained_limit, sets[("limit", "new")])}
    out["chamber_old_invalidated"] = {
        "n": len(inval_cham), "by_kind": kind_count(inval_cham, sets[("chamber", "old")])}
    out["chamber_new_gained"] = {
        "n": len(gained_cham), "by_kind": kind_count(gained_cham, sets[("chamber", "new")])}

    cham_only_old = {k: v for k, v in sets[("chamber", "old")].items()
                     if k not in sets[("limit", "old")]}
    cham_only_valid = {k: v for k, v in cham_only_old.items() if k in sets[("chamber", "new")]}
    out["chamber_only"] = {
        "n_old": len(cham_only_old),
        "by_kind_old": kind_count(cham_only_old, sets[("chamber", "old")]),
        "n_remaining_valid": len(cham_only_valid),
        "by_kind_remaining": kind_count(cham_only_valid, sets[("chamber", "old")]),
    }
    log(f"[13-D] chamber-only: {len(cham_only_old)} old, {len(cham_only_valid)} remain valid")

    witness_rows = []
    for wp in _witness_points():
        wvec = {tuple(lab): int(v) for lab, v in wp["witness"]}
        sample = {k: Fraction(v) for k, v in wp["sample"]}
        prime = int(wp["prime"])
        n_old = sum(1 for k, r in cham_only_old.items() if _pair(r["row"], wvec, sample, prime))
        n_new = sum(1 for k, r in cham_only_valid.items() if _pair(r["row"], wvec, sample, prime))
        witness_rows.append({"prime": prime, "sample": [list(x) for x in wp["sample"]],
                             "breaks_old_chamber_only": n_old,
                             "breaks_remaining_valid": n_new})
        log(f"[13-D] witness p={prime}: breaks {n_old} -> {n_new}")
    out["witness_rescreening"] = witness_rows

    m12r = {}
    if METHOD12R.exists():
        c = json.loads(METHOD12R.read_text(encoding="utf-8")).get("phase_c_row_revalidation", {})
        m12r = {
            "chamber_only_rejected_12r": c.get("chamber_only_revalidation", {}).get(
                "newly_rejected"),
            "limit_rejected_12r": c.get("production_limit_policy_revalidation", {}).get(
                "newly_rejected"),
        }
    out["comparison_with_method12r"] = {
        **m12r,
        "chamber_only_rejected_13": len(cham_only_old) - len(cham_only_valid),
        "limit_rejected_13": len(inval_limit),
        "explanation": (
            "Method.12R phase C was a preliminary bound: its tangent check was per-term on the "
            "complete rays with NO exact cancellation and INCLUDING tangential (d_i = 0) "
            "components, i.e. strictly harsher than the true normal-flux criterion; its "
            "coordinate check is identical to Method.13. The Method.13 counts below therefore "
            "differ only in the tangent column (fewer tangent rows are invalidated, and some "
            "fields the OLD per-term test rejected are now legitimately gained back)."
        ),
    }

    out["int2_lf_verdicts"] = {
        str(lab): str(is_locally_finite(family, lab))
        for lab in ((0, 0, 0, 0, 0, 0, 0), (-1, 0, 0, -1, -1, 0, 0), (-1, 0, 0, 0, -1, 0, -1),
                    (0, 0, 0, -1, 0, 0, -1), (0, 0, 1, -1, 0, 0, -1))
    }
    # keep the corrected chamber row set for the preflight
    out["_sets"] = sets
    out["_family"] = family
    out["_labels"] = labels
    return out


def stage_preflight(int2_out, log) -> dict:
    """Phase D item 16: 1-prime / 3-sample generic LF-feasibility of the corrected chamber rows."""
    spec = importlib.util.spec_from_file_location(
        "t2", REPO_ROOT / "scripts" / "run_external_int2_t2_rankrepair.py")
    T2 = importlib.util.module_from_spec(spec)
    sys.modules["t2"] = T2
    spec.loader.exec_module(T2)

    family = int2_out["_family"]
    labels = int2_out["_labels"]
    rows = [rec["row"] for rec in int2_out["_sets"][("chamber", "new")].values()]
    target = tuple([0] * (family.nvars + family.npolys))
    log(f"[13-preflight] corrected chamber rows: {len(rows)}; computing LF map")
    lf_map = {lab: is_locally_finite(family, lab) for lab in labels}
    lf_map[target] = is_locally_finite(family, target)
    n_lf_true = sum(1 for v in lf_map.values() if v is True)
    log(f"[13-preflight] LF-True labels in box: {n_lf_true}/{len(labels)}")

    samples = T2.default_scattered_samples(family.parameters)
    generic = [s for i, s in enumerate(samples) if i != T2.SPECIAL_SAMPLE_INDEX][:3]
    prime = 2147483647
    points = []
    for sample in generic:
        t0 = time.time()
        res = T2.lf_reduction_feasible_mod_p(rows, labels, target, lf_map, sample, prime)
        payload = T2.feasibility_to_payload(res)
        payload["solve_seconds"] = round(time.time() - t0, 1)
        points.append(payload)
        log(f"[13-preflight] {T2._sample_repr(sample)} p={prime}: {res.status} "
            f"rank={res.rank} ({payload['solve_seconds']}s)")
        if res.status == T2.STATUS_OBSTRUCTED:
            break
    statuses = [p["status"] for p in points]
    verdict = ("Feasible(modular)" if statuses and all(s == "Feasible" for s in statuses)
               else "Obstructed" if "Obstructed" in statuses else "Mixed")
    return {
        "prime": prime,
        "n_lf_true_labels": n_lf_true,
        "points": points,
        "verdict": verdict,
        "note": (
            "1-prime/3-sample preflight on the CORRECTED chamber row set (Phase D item 16). "
            "Feasible = the target lies in the row span through LF-True columns at these "
            "points only; NOT an LF basis, NOT a reconstruction."
        ),
    }


# --------------------------------------------------------------------------------------------
# Phase E: project-wide re-audit of tracked Success artifacts


CONTROLS = {
    "d4_cli": {
        "input": "examples/d4_cli_example_input.wl.txt",
        "target": (0, 0, 0, 0, 0, 0, 0),
        "masters": [(0, 1, 1, 0, -2, -1, 0), (1, 1, 0, 0, -2, -1, 0), (0, 1, 1, 0, -3, -1, 0),
                    (1, 1, 0, 0, -3, -1, 0), (0, 1, 1, 0, -4, -1, 0)],
        "box": ((( 0, 1), (0, 1), (0, 1), (0, 0)), ((-4, 0), (-1, 0), (0, 0))),
        "tangent_blocks": ((1, 1), (2, 2)),
        "live_revalidation": (
            "tests/test_d4_vertical.py (reference relation re-certified in the CORRECTED row "
            "span at 3 (ep,r,p) points; reducer Success with SurfaceValidationStatus=Passed) "
            "and tests/test_d4_cli_e2e.py — both green in this session's suite"
        ),
        "independent_evidence": "validation/expected_d4_coefficients.json (notebook reference)",
    },
    "notebook_example2_n3_five_term": {
        "input": "examples/notebook_example2_n3_five_term_explicit_family.wl.txt",
        "target": (0, 0, 0, 0, 0, 0, 0),
        "masters": [(0, 1, 0, -1, -1, 0, 0), (0, 1, 0, -1, -2, 0, 0), (0, 1, 0, -1, -1, -1, 0),
                    (0, 1, 0, 0, -2, -1, 0), (0, 1, 0, 0, -2, -2, 0)],
        # the certifying expand-1 box of the adaptive run (shallow box m-deepened by one)
        "box": (((0, 0), (0, 1), (0, 0)), ((-1, 0), (-2, 0), (-2, 0), (-1, 0))),
        "tangent_blocks": ((1, 1),),
        "live_revalidation": (
            "tests/test_adaptive_real_family.py::test_default_schedule_certifies_real_five_term_"
            "family — adaptive Success + certificate Passed + coefficients equal the notebook "
            "fixture symbolically, re-run green under the corrected filters in this session"
        ),
        "independent_evidence": "validation/notebook_example2_n3_five_term_expected.json",
    },
    "notebook_example1_d4_alt": {
        "input": "examples/notebook_example1_d4_alt_explicit_family.wl.txt",
        "target": (0, 0, 0, 0, 0, 0, 0),
        # 6-term alternative basis from the notebook fixture (relative_factors)
        "masters": [(0, 0, 1, 1, -3, -2, 0), (0, 0, 1, 1, -2, -2, 0), (0, 0, 2, 1, -3, -2, 0),
                    (0, 1, 1, 0, -3, -2, 0), (0, 1, 1, 0, -2, -2, 0), (0, 1, 1, 0, -2, -1, 0)],
        "box": None,  # no committed Success artifact/row config: fixture-only reference
        "tangent_blocks": (),
        "live_revalidation": None,
        "independent_evidence": "validation/notebook_example1_d4_alt_expected.json",
    },
    "example4_star_corrected": {
        "input": "examples/example4_star_corrected_input.wl.txt",
        "target": (0, 0, 0, 0, 0, 0, 0),
        "masters": [(1, 1, 0, -1, 0, 0, 0), (1, 1, 0, 0, 0, -1, 0)],
        "box": (((0, 1), (0, 2), (0, 1)), ((-2, 0), (-2, 0), (-2, 0), (-2, 0))),
        "tangent_blocks": ((1, 1), (2, 2)),
        "live_revalidation": None,  # end-to-end re-run is env-gated (RUN_EXAMPLE4_STAR)
        "independent_evidence": (
            "validation/notebook_star_example4_known_value_expansion.txt (known value); "
            "validation/example4_star_corrected_diagnostics.json (legacy certificate)"
        ),
    },
    "external_int1_corrected": {
        "input": "examples/external_int1_corrected_input.wl.txt",
        "target": (0, 0, 0, 0, 0),
        "masters": [(0, 0, 0, -1, -1), (0, 0, 0, 0, -2)],
        "box": (((0, 1), (0, 1)), ((-2, 0), (-2, 0), (-2, 0))),
        "tangent_blocks": ((1, 1), (2, 2)),
        "live_revalidation": None,  # end-to-end re-run is env-gated (RUN_EXTERNAL_INT1)
        "independent_evidence": (
            "validation/external_int1_laurent_audit.json + notes/EXTERNAL_INT1_LAURENT_AUDIT.md "
            "— independent numeric Laurent comparison of the assembled result (item 19: this "
            "evidence stands regardless of row provenance)"
        ),
    },
}


def stage_controls(log) -> dict:
    manifest: dict = {}
    for name, cfg in CONTROLS.items():
        log(f"[13-E] control: {name}")
        family = parse_family_text((REPO_ROOT / cfg["input"]).read_text(encoding="utf-8"))
        entry: dict = {"input": cfg["input"]}
        lf = {}
        for lab in cfg["masters"]:
            lf[str(list(lab))] = str(is_locally_finite(family, lab))
        entry["master_lf_under_complete_rays"] = lf
        entry["target_lf"] = str(is_locally_finite(family, cfg["target"]))
        masters_ok = all(v == "True" for v in lf.values())

        rows_delta = None
        if cfg["box"] is not None:
            engine = ScoreEngine(family)
            labels = [
                tuple(list(n) + list(m))
                for n in product(*[range(lo, hi + 1) for lo, hi in cfg["box"][0]])
                for m in product(*[range(lo, hi + 1) for lo, hi in cfg["box"][1]])
            ]
            pol = SurfacePolicy.limit("minus")
            old_rows, old_kinds = assemble(
                family, engine, labels, pol, "old", log, cfg["tangent_blocks"])
            new_rows, new_kinds = assemble(
                family, engine, labels, pol, "new", log, cfg["tangent_blocks"])
            inval = set(old_rows) - set(new_rows)
            gained = set(new_rows) - set(old_rows)
            rows_delta = {
                "n_labels": len(labels),
                "old_total": len(old_rows), "old_by_kind": old_kinds,
                "new_total": len(new_rows), "new_by_kind": new_kinds,
                "old_invalidated": len(inval),
                "new_gained": len(gained),
            }
            entry["rows"] = rows_delta

        entry["live_revalidation"] = cfg["live_revalidation"]
        entry["independent_evidence"] = cfg["independent_evidence"]

        if not masters_ok:
            status = "Revoked"
        elif cfg["live_revalidation"]:
            # the relation/result was re-derived and re-certified against the CORRECTED row
            # system inside this session's test suite — the old certificate is not merely
            # inherited.
            status = "UnchangedValid"
        elif cfg["box"] is None:
            status = "LegacyInsufficientData"
        elif rows_delta and rows_delta["old_invalidated"] == 0:
            # identical row set -> the legacy certificate applies verbatim
            status = "UnchangedValid"
        else:
            status = "NeedsRecompute"
        entry["status"] = status
        log(f"[13-E] {name}: {status} (masters LF ok: {masters_ok}, "
            f"rows {rows_delta and rows_delta['old_invalidated']} invalidated)")
        manifest[name] = entry
    return manifest


# --------------------------------------------------------------------------------------------


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--stages", default="int2,controls",
                    help="comma list from {int2, preflight, controls}")
    ap.add_argument("--allow-heavy", action="store_true",
                    help="required for the preflight (the only RREF in this audit)")
    args = ap.parse_args(argv)
    stages = tuple(s.strip() for s in args.stages.split(",") if s.strip())
    if "preflight" in stages and not args.allow_heavy:
        ap.error("--stages preflight requires --allow-heavy (RREF)")
    if "preflight" in stages and "int2" not in stages:
        ap.error("preflight needs the int2 stage in the same run (it uses the corrected rows)")

    def log(msg):
        print(msg, flush=True)

    t0 = time.time()
    if "int2" in stages or "preflight" in stages:
        int2 = stage_int2(log)
        report = {k: v for k, v in int2.items() if not k.startswith("_")}
        if "preflight" in stages:
            report["lf_feasibility_preflight"] = stage_preflight(int2, log)
        report["scope_note"] = (
            "Method.13 Phase D surface re-audit: rebuild-only comparison of the Method.10 "
            "Level-0 candidate rows under the legacy vs the corrected complete-face surface "
            "criteria (docs/TORIC_SURFACE_VALIDATION.md). No production artifact is modified; "
            "the preflight (if present) is the only RREF and is span evidence only."
        )
        report["elapsed_seconds"] = round(time.time() - t0, 1)
        merged = report
        if AUDIT_OUT.exists():
            try:
                old = json.loads(AUDIT_OUT.read_text(encoding="utf-8"))
                merged = {**old, **report}
            except json.JSONDecodeError:
                pass
        AUDIT_OUT.write_text(json.dumps(merged, indent=1, sort_keys=True) + "\n",
                             encoding="utf-8")
        log(f"[13] wrote {AUDIT_OUT}")

    if "controls" in stages:
        manifest = {
            "schema": "surface-reaudit-manifest/v1",
            "method": "External Int2 Method.13 Phase E: project-wide Success re-audit",
            "criteria": "docs/TORIC_SURFACE_VALIDATION.md (complete-face surface validation)",
            "results": stage_controls(log),
        }
        MANIFEST_OUT.write_text(json.dumps(manifest, indent=1, sort_keys=True) + "\n",
                                encoding="utf-8")
        log(f"[13] wrote {MANIFEST_OUT}")
    log(f"[13] done ({round(time.time() - t0, 1)}s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
