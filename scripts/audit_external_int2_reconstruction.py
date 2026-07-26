"""External Int2 Method.11c: offline reconstruction audit of the Method.11 record cache.

Reads the existing 38-sample x 4-prime modular record cache and answers -- without
generating any new modular records and without invoking any package RREF backend or
reducer solve -- whether the four LF coefficients admit a sparse/factored/structured
reconstruction, and if not, what the minimal Stage-5 sample schedule is.

Phases:

* A: compact exact value table (per-prime residues, CRT-lifted rationals,
  special-zero zero-fills, source cache keys) + integrity verification.
* B: per-coefficient dense degree diagnostics over the Stage-4 ladder
  (monomial/unknown counts, nullspace dimensions, consistency), cross-checked
  against the Stage-4 artifact attempt log.
* C: structured experiments on cached values only: univariate fibers (A/B),
  sparse tensor supports (C), factor-aware denominator ansatz (D),
  shared-denominator ansatz (E), dense total-degree baseline (F).
* D: minimal deterministic new-sample schedule per viable ansatz.

All linear algebra below is local dense elimination on interpolation design
matrices (<= ~150 unknowns) built from cached values; the package's
``sparse_rref`` / ``reducer`` machinery is never imported or invoked, and no
modular reduction records are created.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import random
import re
import sys
from fractions import Fraction
from math import comb, gcd
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

import sympy as sp  # noqa: E402

from parametric_ibp_lf_reducer.reconstruction import (  # noqa: E402
    classify_sample_supports,
    collect_value_table,
    reconstruct_rational,
)

RECORDS_PATH = REPO_ROOT / "validation" / "external_int2_method11_records.jsonl"
ARTIFACT_PATH = REPO_ROOT / "validation" / "external_int2_method11b_stage4p4.json"
VALUE_TABLE_OUT = REPO_ROOT / "validation" / "external_int2_method11_value_table.json"
AUDIT_OUT = REPO_ROOT / "validation" / "external_int2_method11_reconstruction_audit.json"
NOTE_OUT = REPO_ROOT / "notes" / "EXTERNAL_INT2_RECONSTRUCTION_AUDIT.md"

CACHE_SCHEMA = "m11-record-cache/v2"
EXPECTED_PRIMES = (2147483579, 2147483587, 2147483629, 2147483647)
EXPECTED_RANK = 26984
EXPECTED_SAMPLES = 38
SPECIAL_ZERO_CANONICAL = "ep=6,r=57/11"
MIN_VALIDATION = 2
SCHEDULE_SEED = 20260726

# Dense ladder replicating Stage-4 (all (num_deg, den_deg) <= (6,6), by total degree),
# plus the two extended determined cells (7,0) / (0,7) reachable with 36 fit points.
LADDER = sorted(((i, j) for i in range(7) for j in range(7)), key=lambda c: (c[0] + c[1], c[0]))
LADDER_EXTENDED = [(0, 7), (7, 0)]


# --------------------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------------------
def _load_runner():
    """Import the Method.11 runner for its faithful record JSON codec."""
    name = "run_external_int2_method11"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def load_records(path: Path):
    """Read the cache: returns (header, records, raw_rows). Never writes."""
    lines = path.read_text(encoding="utf-8").splitlines()
    header = json.loads(lines[0])
    if header.get("schema") != CACHE_SCHEMA:
        raise SystemExit(f"unexpected cache schema: {header.get('schema')!r}")
    runner = _load_runner()
    raw_rows = [json.loads(x) for x in lines[1:] if x.strip()]
    records = [runner._record_from_json(row) for row in raw_rows]
    return header, records, raw_rows


def _lab_str(label) -> str:
    return str(tuple(label))


def _frac_str(v: Fraction) -> str:
    return str(Fraction(v))


def _canon(sample: dict) -> str:
    return ",".join(f"{k}={Fraction(v)}" for k, v in sorted(sample.items(), key=lambda kv: kv[0]))


# --------------------------------------------------------------------------------------
# Phase A: value table
# --------------------------------------------------------------------------------------
def build_value_table(records, raw_rows):
    """Exact value table + verification block (pure read-only reprocessing)."""
    classification = classify_sample_supports(records)
    labels, table, samples, n_skipped = collect_value_table(records)
    by_sample: dict = {}
    for rec, raw in zip(records, raw_rows):
        canon = _canon(rec.sample)
        entry = by_sample.setdefault(
            canon, {"sample": rec.sample, "keys": [], "byprime": {}, "ranks": set()}
        )
        entry["keys"].append(raw["key"])
        entry["byprime"][rec.prime] = dict(rec.coeffs)
        entry["ranks"].add(rec.rank)

    ranks_all = sorted({r.rank for r in records})
    problems: list[str] = []
    if len(samples) != EXPECTED_SAMPLES:
        problems.append(f"expected {EXPECTED_SAMPLES} samples, got {len(samples)}")
    if ranks_all != [EXPECTED_RANK]:
        problems.append(f"rank set {ranks_all} != [{EXPECTED_RANK}]")

    canon_by_key = {key: _canon(samples[key]) for key in samples}
    key_by_canon = {v: k for k, v in canon_by_key.items()}
    sz_entries = {
        canon: e
        for canon, e in classification["by_sample"].items()
        if e["classification"] == "special_zero"
    }

    out_samples = []
    for canon in sorted(
        by_sample, key=lambda c: tuple(Fraction(p.split("=")[1]) for p in c.split(","))
    ):
        entry = by_sample[canon]
        primes = sorted(entry["byprime"])
        if tuple(primes) != EXPECTED_PRIMES:
            problems.append(f"{canon}: incomplete prime group {primes}")
        supports = {p: tuple(sorted(entry["byprime"][p])) for p in primes}
        if len(set(supports.values())) != 1:
            problems.append(f"{canon}: per-prime support disagreement")
        cls_entry = classification["by_sample"].get(canon, {})
        cls = cls_entry.get("classification")
        skey = key_by_canon[canon]
        coeffs_out = {}
        for lab in labels:
            residues = {p: entry["byprime"][p].get(lab) for p in primes}
            zero_filled = all(v is None for v in residues.values())
            if not zero_filled and any(v is None for v in residues.values()):
                problems.append(f"{canon}/{_lab_str(lab)}: label present at some primes only")
            value = table[lab][skey]
            for p, res in residues.items():
                want = 0 if res is None else res
                got = (value.numerator * pow(value.denominator, -1, p)) % p
                if got != want:
                    problems.append(f"{canon}/{_lab_str(lab)}: value != residue mod {p}")
            coeffs_out[_lab_str(lab)] = {
                "residues": {str(p): (0 if r is None else r) for p, r in residues.items()},
                "zero_filled": zero_filled,
                "value": _frac_str(value),
            }
        if cls == "special_zero":
            missing = {p: tuple(sorted(set(labels) - set(entry["byprime"][p]))) for p in primes}
            if len(set(missing.values())) != 1:
                problems.append(f"{canon}: special-zero missing labels differ across primes")
        out_samples.append(
            {
                "sample": {k: _frac_str(v) for k, v in sorted(entry["sample"].items())},
                "canonical": canon,
                "classification": cls,
                "selected_rank": sorted(entry["ranks"])[0],
                "labels": [_lab_str(lab) for lab in labels],
                "coefficients": coeffs_out,
                "cache_keys": sorted(entry["keys"]),
            }
        )

    if sorted(sz_entries) != [SPECIAL_ZERO_CANONICAL]:
        problems.append(f"special-zero samples {sorted(sz_entries)} != [{SPECIAL_ZERO_CANONICAL}]")

    verification = {
        "n_samples": len(samples),
        "n_records": len(records),
        "primes": list(EXPECTED_PRIMES),
        "selected_rank": ranks_all[0] if len(ranks_all) == 1 else ranks_all,
        "n_skipped_records": n_skipped,
        "support_size": len(labels),
        "classification_summary": classification["summary"],
        "special_zero": {
            canon: {"missing_labels": [_lab_str(t) for t in e.get("missing_labels", ())]}
            for canon, e in sz_entries.items()
        },
        "problems": problems,
        "ok": not problems,
    }
    doc = {
        "schema": "m11c-value-table/v1",
        "source": {"records": str(RECORDS_PATH.name), "fingerprint": header_fp},
        "labels": [_lab_str(lab) for lab in labels],
        "samples": out_samples,
        "verification": verification,
    }
    return doc, labels, table, samples, classification


header_fp: dict = {}  # populated in main() before build_value_table is called


# --------------------------------------------------------------------------------------
# small dense linear algebra mod p (local interpolation design matrices only)
# --------------------------------------------------------------------------------------
def _row_reduce_mod(rows: list[list[int]], p: int):
    """Local elimination on a copy; returns (rank, pivot_cols, echelon matrix)."""
    m = [row[:] for row in rows]
    ncols = len(m[0]) if m else 0
    piv_cols: list[int] = []
    r = 0
    for c in range(ncols):
        pr = next((i for i in range(r, len(m)) if m[i][c] % p), None)
        if pr is None:
            continue
        m[r], m[pr] = m[pr], m[r]
        inv = pow(m[r][c], -1, p)
        m[r] = [(x * inv) % p for x in m[r]]
        for i in range(len(m)):
            if i != r and m[i][c]:
                f = m[i][c]
                m[i] = [(a - f * b) % p for a, b in zip(m[i], m[r])]
        piv_cols.append(c)
        r += 1
        if r == len(m):
            break
    return r, piv_cols, m


def _nullspace_dim_mod(rows: list[list[int]], p: int):
    """(dim, nullvector-if-dim-1). dim==0 mod p proves dim==0 over QQ."""
    ncols = len(rows[0]) if rows else 0
    rank, piv_cols, ech = _row_reduce_mod(rows, p)
    dim = ncols - rank
    if dim != 1:
        return dim, None
    free = next(c for c in range(ncols) if c not in piv_cols)
    vec = [0] * ncols
    vec[free] = 1
    for r, c in enumerate(piv_cols):
        vec[c] = (-ech[r][free]) % p
    return 1, vec


def _monoms_total(d: int):
    return sorted(
        ((a, b) for a in range(d + 1) for b in range(d + 1 - a)), key=lambda m: (sum(m), m[0])
    )


def _monoms_tensor(a: int, b: int):
    return sorted(
        ((i, j) for i in range(a + 1) for j in range(b + 1)), key=lambda m: (sum(m), m[0])
    )


def _mod_frac(v: Fraction, p: int) -> int:
    return (v.numerator * pow(v.denominator, -1, p)) % p


def _eval_monoms_mod(pt, monoms, p: int):
    ep, r = pt
    epm, rm = _mod_frac(ep, p), _mod_frac(r, p)
    return [(pow(epm, a, p) * pow(rm, b, p)) % p for a, b in monoms]


def _eval_poly_exact(coeffs, monoms, pt) -> Fraction:
    ep, r = pt
    return sum((c * ep**a * r**b for c, (a, b) in zip(coeffs, monoms)), Fraction(0))


# --------------------------------------------------------------------------------------
# Phase B: degree diagnostics
# --------------------------------------------------------------------------------------
def split_fit_hold(samples: dict):
    """Deterministic fit/holdout split mirroring the runner (last 2 points held out)."""
    pts = sorted(((s["ep"], s["r"]) for s in samples.values()), key=lambda t: (t[0], t[1]))
    return pts[:-MIN_VALIDATION], pts[-MIN_VALIDATION:]


def _values_by_point(table_lab: dict, samples: dict):
    return {(samples[k]["ep"], samples[k]["r"]): v for k, v in table_lab.items()}


def _dense_rows(vals, pts, num_monoms, den_monoms, p):
    rows = []
    for pt in pts:
        vm = _mod_frac(vals[pt], p)
        row = _eval_monoms_mod(pt, num_monoms, p)
        row += [(-vm * x) % p for x in _eval_monoms_mod(pt, den_monoms, p)]
        rows.append(row)
    return rows


def _validate_rational_exact(num_c, num_m, den_c, den_m, vals, hold, special_pt, special_zero):
    """Acceptance per Method.11c item 7: all samples exact, >=2 holdout, special zero."""
    for pt, v in vals.items():
        den = _eval_poly_exact(den_c, den_m, pt)
        if den == 0:
            return False, f"denominator vanishes at {pt}"
        if _eval_poly_exact(num_c, num_m, pt) / den != v:
            return False, f"mismatch at {pt}"
    for pt in hold:
        if pt not in vals:
            return False, "holdout point missing from value table"
    if special_pt in vals and special_zero and vals[special_pt] != 0:
        return False, "special-zero value not zero in table"
    return True, "exact on all fit+holdout points (incl. special-zero)"


def cell_diagnostics(vals, fit, cell, primes):
    i, j = cell
    num_m, den_m = _monoms_total(i), _monoms_total(j)
    unknowns = len(num_m) + len(den_m)
    entry = {
        "num_deg": i,
        "den_deg": j,
        "n_num_monomials": len(num_m),
        "n_den_monomials": len(den_m),
        "unknowns": unknowns,
        "n_fit": len(fit),
        "n_hold": MIN_VALIDATION,
    }
    if unknowns - 1 > len(fit):
        entry["status"] = f"underdetermined (needs {unknowns - 1} fit points, have {len(fit)})"
        entry["nullspace_dim"] = None
        return entry, None
    dims = {}
    nullvec = None
    for p in primes:
        dim, vec = _nullspace_dim_mod(_dense_rows(vals, fit, num_m, den_m, p), p)
        dims[str(p)] = dim
        if p == primes[0]:
            nullvec = vec
    entry["nullspace_dim"] = dims[str(primes[0])]
    entry["nullspace_dims_by_prime"] = dims
    dim0 = dims[str(primes[0])]
    if len(set(dims.values())) != 1:
        entry["status"] = f"prime-dependent dim {sorted(set(dims.values()))} (borderline)"
    elif dim0 == 0:
        entry["status"] = "inconsistent (dim 0): no rational function of this cell fits"
    elif dim0 == 1:
        entry["status"] = "candidate (dim 1)"
    else:
        entry["status"] = f"ambiguous (dim {dim0}): degenerate sample geometry for this cell"
    return entry, (nullvec if dim0 == 1 else None)


def phase_b(labels, table, samples, primes, artifact):
    fit, hold = split_fit_hold(samples)
    special_pt = (Fraction(6), Fraction(57, 11))
    report = {}
    for lab in labels:
        vals = _values_by_point(table[lab], samples)
        cells = []
        first_admissible = None
        for cell in LADDER + LADDER_EXTENDED:
            entry, vec = cell_diagnostics(vals, fit, cell, primes)
            entry["extended"] = cell in LADDER_EXTENDED
            if vec is not None:
                num_m, den_m = _monoms_total(cell[0]), _monoms_total(cell[1])
                exact = _solve_cell_exact(vals, fit, hold, num_m, den_m, primes, special_pt)
                entry["exact_validation"] = exact["status"]
                if exact["ok"]:
                    entry["validated_candidate"] = exact["expr"]
            cells.append(entry)
            if first_admissible is None and (
                entry["status"].startswith("underdetermined")
                or entry.get("validated_candidate") is not None
            ):
                first_admissible = {
                    k: entry[k] for k in ("num_deg", "den_deg", "unknowns", "status")
                }
        report[_lab_str(lab)] = {
            "ladder_cells": cells,
            "first_potentially_admissible_dense_cell": first_admissible,
        }
    star = "(-1, 0, 0, -1, -1, 0, 0)"
    report["artifact_crosscheck"] = _crosscheck_artifact(report.get(star), artifact, star)
    report[star]["failure_explanation"] = _explain_star(report[star])
    return report, fit, hold


def _crosscheck_artifact(star_entry, artifact, star):
    if not artifact or star_entry is None:
        return {"available": False}
    try:
        attempts = artifact["reconstruction"]["detail"]["coefficients"][star]["attempts"]
    except (KeyError, TypeError):
        return {"available": False}
    ours = {
        (c["num_deg"], c["den_deg"]): c for c in star_entry["ladder_cells"] if not c["extended"]
    }
    mismatches = []
    for att in attempts:
        cell = (att["num_deg"], att["den_deg"])
        mine = ours.get(cell)
        if mine is None:
            mismatches.append({"cell": cell, "why": "cell not scanned offline"})
            continue
        m_dim = re.search(r"dim (\d+)", att["status"])
        m_need = re.search(r"needs (\d+) fit points", att["status"])
        if m_dim and mine["nullspace_dim"] != int(m_dim.group(1)):
            mismatches.append(
                {"cell": cell, "artifact": att["status"], "offline_dim": mine["nullspace_dim"]}
            )
        if m_need and f"needs {m_need.group(1)} fit points" not in mine["status"]:
            mismatches.append({"cell": cell, "artifact": att["status"], "offline": mine["status"]})
    return {
        "available": True,
        "n_attempts": len(attempts),
        "mismatches": mismatches,
        "matches": not mismatches,
    }


def _explain_star(entry):
    cells = [c for c in entry["ladder_cells"] if not c["extended"]]
    det = [c for c in cells if not c["status"].startswith("underdetermined")]
    dim0 = [c for c in det if c["nullspace_dim"] == 0]
    ambig = [c for c in det if (c["nullspace_dim"] or 0) >= 2]
    under = [c for c in cells if c["status"].startswith("underdetermined")]
    return (
        f"Through (6,6): {len(dim0)} determined cells are inconsistent (nullspace dim 0 mod all "
        f"4 primes, which proves dim 0 over QQ -- the coefficient is not a rational function of "
        f"those degrees); {len(ambig)} cells "
        f"({[(c['num_deg'], c['den_deg']) for c in ambig]}) have nullspace dim >= 2 because 36 "
        f"scattered points with all-distinct ep and r are degenerate for these supports (no "
        f"unique candidate, holdout cannot arbitrate); the remaining {len(under)} cells need "
        f"37..55 fit points but only 36 exist. Hence no dense cell <= (6,6) can validate: the "
        f"true rational form lies outside the dense (6,6) ladder reachable with 38 samples."
    )


# --------------------------------------------------------------------------------------
# exact solving of a dim-1 cell (shared by phases B/C)
# --------------------------------------------------------------------------------------
def _to_frac(x) -> Fraction:
    q = sp.Rational(x)
    return Fraction(int(q.p), int(q.q))


def _solve_cell_exact(vals, fit, hold, num_m, den_m, primes, special_pt):
    """Exact QQ nullspace for one candidate cell + full acceptance validation."""
    rows = []
    for pt in fit:
        v = vals[pt]
        ep, r = pt
        row = [Fraction(ep**a * r**b) for a, b in num_m]
        row += [-v * ep**a * r**b for a, b in den_m]
        rows.append([sp.Rational(f.numerator, f.denominator) for f in row])
    ns = sp.Matrix(rows).nullspace()
    if len(ns) != 1:
        return {"ok": False, "status": f"exact nullspace dim {len(ns)} != 1"}
    vec = [_to_frac(x) for x in ns[0]]
    num_c, den_c = vec[: len(num_m)], vec[len(num_m) :]
    special_zero = special_pt in vals and vals[special_pt] == 0
    ok, why = _validate_rational_exact(
        num_c, num_m, den_c, den_m, vals, hold, special_pt, special_zero
    )
    if not ok:
        return {"ok": False, "status": f"exact candidate rejected: {why}"}
    return {
        "ok": True,
        "status": "validated: " + why,
        "expr": _expr_str(num_c, num_m, den_c, den_m),
    }


def _expr_str(num_c, num_m, den_c, den_m):
    ep, r = sp.symbols("ep r")
    num = sum(
        sp.Rational(c.numerator, c.denominator) * ep**a * r**b for c, (a, b) in zip(num_c, num_m)
    )
    den = sum(
        sp.Rational(c.numerator, c.denominator) * ep**a * r**b for c, (a, b) in zip(den_c, den_m)
    )
    return str(sp.simplify(num / den))


# --------------------------------------------------------------------------------------
# Phase C: structured experiments
# --------------------------------------------------------------------------------------
def phase_c(labels, table, samples, primes, fit, hold):
    special_pt = (Fraction(6), Fraction(57, 11))
    out = {}
    # A/B: univariate fibers
    by_r: dict = {}
    by_ep: dict = {}
    for s in samples.values():
        by_ep.setdefault(s["ep"], []).append(s["r"])
        by_r.setdefault(s["r"], []).append(s["ep"])
    usable_ep = {str(k): len(v) for k, v in by_ep.items() if len(v) >= MIN_VALIDATION + 2}
    usable_r = {str(k): len(v) for k, v in by_r.items() if len(v) >= MIN_VALIDATION + 2}
    out["A_univariate_ep_on_fixed_r"] = {
        "usable_fibers": usable_r,
        "statement": "no repeated-r fibers exist in the cached samples (all 38 r values are "
        "distinct); univariate interpolation in ep is impossible without new data"
        if not usable_r
        else "fibers available",
    }
    out["B_univariate_r_on_fixed_ep"] = {
        "usable_fibers": usable_ep,
        "statement": "no repeated-ep fibers exist in the cached samples (all 38 ep values are "
        "distinct); univariate interpolation in r is impossible without new data"
        if not usable_ep
        else "fibers available",
    }

    # C: sparse tensor-grid supports
    tensor = {}
    p0 = primes[0]
    for lab in labels:
        vals = _values_by_point(table[lab], samples)
        tried = n_dim0 = n_ambig = 0
        candidates = []
        for a1 in range(5):
            for a2 in range(5):
                for b1 in range(5):
                    for b2 in range(5):
                        num_m = _monoms_tensor(a1, a2)
                        den_m = _monoms_tensor(b1, b2)
                        unknowns = len(num_m) + len(den_m)
                        if unknowns - 1 > len(fit):
                            continue
                        tried += 1
                        dim, _ = _nullspace_dim_mod(_dense_rows(vals, fit, num_m, den_m, p0), p0)
                        if dim == 0:
                            n_dim0 += 1
                        elif dim == 1:
                            exact = _solve_cell_exact(
                                vals, fit, hold, num_m, den_m, primes, special_pt
                            )
                            candidates.append(
                                {
                                    "support": [a1, a2, b1, b2],
                                    "unknowns": unknowns,
                                    "exact": exact["status"],
                                    **({"expr": exact["expr"]} if exact["ok"] else {}),
                                }
                            )
                        else:
                            n_ambig += 1
        tensor[_lab_str(lab)] = {
            "cells_tried": tried,
            "inconsistent": n_dim0,
            "ambiguous": n_ambig,
            "dim1_candidates": candidates,
            "validated": [c for c in candidates if "expr" in c],
        }
    out["C_sparse_tensor_supports"] = tensor

    # D: factor-aware denominator ansatz
    out["D_factor_aware"] = _factor_aware(labels, table, samples, primes, fit, hold, special_pt)

    # E: shared denominator across the four coefficients
    out["E_shared_denominator"] = _shared_denominator(
        labels, table, samples, primes, fit, hold, special_pt
    )

    # G: piecewise / chamber-boundary hypothesis (low height values + dense dim-0
    # everywhere is the signature of a piecewise-rational coefficient)
    out["G_piecewise_chamber_split"] = _piecewise_split(labels, table, samples, primes)

    out["F_dense_baseline"] = {
        "statement": "dense total-degree ladder results are Phase B (all cells <= (6,6) plus "
        "extended (7,0)/(0,7) refuted or underdetermined per label)"
    }
    return out


def _uv(pt):
    ep, r = pt
    u, v = 7 * ep, 11 * r
    assert u.denominator == 1 and v.denominator == 1
    return int(u), int(v)


def _val_q(n: int, q: int) -> int:
    v = 0
    n = abs(n)
    while n and n % q == 0:
        n //= q
        v += 1
    return v


def _form_pool(points):
    """Small linear forms in (7*ep, 11*r) that vanish at NO cached sample (pole safety)."""
    forms = []
    for alpha in range(4):
        for beta in range(-3, 4):
            if alpha == 0 and beta <= 0:
                continue
            for gamma in range(-15, 16):
                if gcd(gcd(abs(alpha), abs(beta)), abs(gamma)) != 1:
                    continue
                if any(alpha * _uv(pt)[0] + beta * _uv(pt)[1] + gamma == 0 for pt in points):
                    continue
                forms.append((alpha, beta, gamma))
    return forms


def _factor_aware(labels, table, samples, primes, fit, hold, special_pt):
    """Data-driven denominator form search + polynomial fit of value * denominator."""
    forms = None
    report = {}
    for lab in labels:
        vals = _values_by_point(table[lab], samples)
        pts = [pt for pt in vals if vals[pt] != 0]
        if forms is None:
            forms = _form_pool(list(vals))
        dens = {pt: vals[pt].denominator for pt in pts}
        # residual = odd, 3-free part that a bounded coefficient denominator cannot absorb
        residual = {}
        for pt, d in dens.items():
            while d % 2 == 0:
                d //= 2
            while d % 3 == 0:
                d //= 3
            residual[pt] = d
        needed = {pt: set(sp.factorint(rd)) for pt, rd in residual.items() if rd > 1}
        chosen: list[tuple[int, int, int]] = []
        uncovered = {(pt, q) for pt, qs in needed.items() for q in qs}
        while uncovered and len(chosen) < 4:

            def score(form):
                a, b, c = form
                return sum(
                    1 for (pt, q) in uncovered if (a * _uv(pt)[0] + b * _uv(pt)[1] + c) % q == 0
                )

            best = max(forms, key=lambda f: (score(f), -abs(f[2]), f))
            if score(best) == 0:
                break
            chosen.append(best)
            uncovered = {
                (pt, q)
                for (pt, q) in uncovered
                if (best[0] * _uv(pt)[0] + best[1] * _uv(pt)[1] + best[2]) % q != 0
            }
        cover_ok = not uncovered
        entry = {
            "denominator_forms": [
                {"form": f"({a}*(7*ep) + {b}*(11*r) + {c})", "alpha_beta_gamma": [a, b, c]}
                for a, b, c in chosen
            ],
            "full_cover": cover_ok,
            "uncovered_prime_hits": len(uncovered),
        }
        if not cover_ok:
            entry["statement"] = (
                "no product of small linear forms in (7*ep, 11*r) covers all denominator "
                "primes; factor-aware ansatz refuted within the searched form family"
            )
            report[_lab_str(lab)] = entry
            continue
        # exponents: minimal demand-satisfying assignment (per-prime valuation accounting)
        base_exps = [1] * len(chosen)
        demands = [(pt, q, k) for pt in pts for q, k in sp.factorint(residual[pt]).items()]
        solvable = True
        for _ in range(80):
            viol = next(
                (
                    (pt, q, k)
                    for pt, q, k in demands
                    if sum(
                        e * _val_q(a * _uv(pt)[0] + b * _uv(pt)[1] + c, q)
                        for (a, b, c), e in zip(chosen, base_exps)
                    )
                    < k
                ),
                None,
            )
            if viol is None:
                break
            pt, q, _k = viol
            vq = [_val_q(a * _uv(pt)[0] + b * _uv(pt)[1] + c, q) for a, b, c in chosen]
            if max(vq) == 0:
                solvable = False
                break
            base_exps[max(range(len(chosen)), key=lambda i: vq[i])] += 1
        else:
            solvable = False
        entry["minimal_exponents"] = list(base_exps)
        if not solvable:
            entry["statement"] = "no bounded exponent assignment satisfies all valuations"
            report[_lab_str(lab)] = entry
            continue
        result = None
        for extra in range(3):
            exps = [e + extra for e in base_exps]
            dvals = {}
            for pt in vals:
                d = Fraction(1)
                for (a, b, c), e in zip(chosen, exps):
                    d *= Fraction(a * _uv(pt)[0] + b * _uv(pt)[1] + c) ** e
                dvals[pt] = vals[pt] * d
            fitres = _poly_fit_ladder(dvals, fit, hold, primes, special_pt)
            entry.setdefault("attempts", []).append(
                {"exponents": exps, **{k: v for k, v in fitres.items() if k != "expr"}}
            )
            if fitres["ok"]:
                den_expr = "*".join(
                    f"({a}*(7*ep) + {b}*(11*r) + {c})**{e}" for (a, b, c), e in zip(chosen, exps)
                )
                entry["validated"] = {
                    "numerator_poly": fitres["expr"],
                    "denominator": den_expr,
                    "poly_degree": fitres["degree"],
                }
                result = entry["validated"]
                break
        report[_lab_str(lab)] = entry
        if result is None and cover_ok and "validated" not in entry:
            entry["statement"] = (
                "denominator forms found, but value*denominator is not a polynomial of total "
                "degree <= 7 (the maximum determinable from 36 fit points)"
            )
    return report


def _poly_fit_ladder(wvals, fit, hold, primes, special_pt):
    """Fit w = value * D as a dense polynomial, degree ladder 0..7 (36 fit points)."""
    for d in range(8):
        monoms = _monoms_total(d)
        if len(monoms) - 1 > len(fit):
            break
        p0 = primes[0]
        rows = []
        for pt in fit:
            row = _eval_monoms_mod(pt, monoms, p0)
            row.append((-_mod_frac(wvals[pt], p0)) % p0)
            rows.append(row)
        dim, _ = _nullspace_dim_mod(rows, p0)
        if dim != 1:
            continue
        srows = [
            [
                sp.Rational(
                    (pt[0] ** a * pt[1] ** b).numerator, (pt[0] ** a * pt[1] ** b).denominator
                )
                for a, b in monoms
            ]
            + [sp.Rational(-wvals[pt].numerator, wvals[pt].denominator)]
            for pt in fit
        ]
        ns = sp.Matrix(srows).nullspace()
        if len(ns) != 1 or ns[0][-1] == 0:
            continue
        vec = [_to_frac(x / ns[0][-1]) for x in ns[0][:-1]]
        if all(_eval_poly_exact(vec, monoms, pt) == wvals[pt] for pt in wvals):
            ep, r = sp.symbols("ep r")
            expr = sum(
                sp.Rational(c.numerator, c.denominator) * ep**a * r**b
                for c, (a, b) in zip(vec, monoms)
            )
            return {
                "ok": True,
                "degree": d,
                "n_unknowns": len(monoms),
                "status": f"polynomial of total degree {d} exact on all 38 points "
                f"(36 fit + {MIN_VALIDATION} holdout)",
                "expr": str(sp.expand(expr)),
            }
    return {"ok": False, "status": "no polynomial of total degree <= 7 matches (36 fit points)"}


def _shared_denominator(labels, table, samples, primes, fit, hold, special_pt):
    """One denominator shared by all four coefficients: 4 equations per sample."""
    p0 = primes[0]
    vals = {lab: _values_by_point(table[lab], samples) for lab in labels}
    results = {"cells": [], "validated": None}
    for a in range(8):
        for b in range(13):
            num_m, den_m = _monoms_total(a), _monoms_total(b)
            unknowns = 4 * len(num_m) + len(den_m)
            n_eq = 4 * len(fit)
            if unknowns - 1 > n_eq:
                continue
            rows = []
            for li, lab in enumerate(labels):
                for pt in fit:
                    vm = _mod_frac(vals[lab][pt], p0)
                    row = [0] * (4 * len(num_m))
                    nm = _eval_monoms_mod(pt, num_m, p0)
                    row[li * len(num_m) : (li + 1) * len(num_m)] = nm
                    row += [(-vm * x) % p0 for x in _eval_monoms_mod(pt, den_m, p0)]
                    rows.append(row)
            dim, vec = _nullspace_dim_mod(rows, p0)
            cell = {
                "num_deg": a,
                "den_deg": b,
                "unknowns": unknowns,
                "n_equations": n_eq,
                "nullspace_dim_mod_p": dim,
            }
            if dim == 1:
                lift = _lift_shared(labels, vals, fit, hold, num_m, den_m, primes, special_pt)
                cell["lift"] = lift["status"]
                if lift["ok"]:
                    cell["validated"] = lift["exprs"]
                    results["validated"] = {"cell": [a, b], "exprs": lift["exprs"]}
            results["cells"].append(cell)
        if results["validated"]:
            break
    if results["validated"] is None:
        det = [c for c in results["cells"] if c["nullspace_dim_mod_p"] == 0]
        results["statement"] = (
            f"all {len(det)} determined shared-denominator cells are inconsistent (dim 0): the "
            "four coefficients do not share one denominator within the scanned degrees"
            if det and all(c["nullspace_dim_mod_p"] == 0 for c in results["cells"])
            else "no validated shared-denominator cell within scanned degrees"
        )
    return results


def _lift_shared(labels, vals, fit, hold, num_m, den_m, primes, special_pt):
    """CRT-lift the shared-denominator nullvector across all 4 primes; verify exactly."""
    vecs = {}
    for p in primes:
        rows = []
        for li, lab in enumerate(labels):
            for pt in fit:
                vm = _mod_frac(vals[lab][pt], p)
                row = [0] * (4 * len(num_m))
                row[li * len(num_m) : (li + 1) * len(num_m)] = _eval_monoms_mod(pt, num_m, p)
                row += [(-vm * x) % p for x in _eval_monoms_mod(pt, den_m, p)]
                rows.append(row)
        dim, vec = _nullspace_dim_mod(rows, p)
        if dim != 1:
            return {"ok": False, "status": f"dim {dim} != 1 mod {p}"}
        vecs[p] = vec
    pivot = next(
        (i for i in range(len(vecs[primes[0]])) if all(vecs[p][i] % p for p in primes)), None
    )
    if pivot is None:
        return {"ok": False, "status": "no common nonzero pivot coordinate"}
    lifted = []
    for i in range(len(vecs[primes[0]])):
        residues = {p: (vecs[p][i] * pow(vecs[p][pivot], -1, p)) % p for p in primes}
        f = reconstruct_rational(residues)
        if f is None:
            return {
                "ok": False,
                "status": f"CRT lift not stabilized at coordinate {i} (needs more primes)",
            }
        lifted.append(f)
    den_c = lifted[4 * len(num_m) :]
    exprs = {}
    for li, lab in enumerate(labels):
        num_c = lifted[li * len(num_m) : (li + 1) * len(num_m)]
        ok, why = _validate_rational_exact(
            num_c,
            num_m,
            den_c,
            den_m,
            vals[lab],
            hold,
            special_pt,
            special_pt in vals[lab] and vals[lab][special_pt] == 0,
        )
        if not ok:
            return {"ok": False, "status": f"{_lab_str(lab)}: {why}"}
        exprs[_lab_str(lab)] = _expr_str(num_c, num_m, den_c, den_m)
    return {"ok": True, "status": "validated on all labels/samples/holdout", "exprs": exprs}


BOUNDARY_FORMS = (
    [(1, 0, -c) for c in (45, 60, 75, 90, 105)]
    + [(0, 1, -c) for c in (32, 38, 44, 50, 56)]
    + [(1, -1, -c) for c in (-30, -15, 0, 15, 30, 45, 60)]
    + [(1, -2, -c) for c in (-60, -45, -30, -15, 0, 15)]
    + [(3, -2, -c) for c in (-12, 0, 12, 24, 36)]
    + [(1, 1, -c) for c in (60, 90, 120, 150)]
)


def _piecewise_split(labels, table, samples, primes):
    """Experiment G: split samples by a linear boundary in (7*ep, 11*r) and fit
    low-degree dense rational cells on each side separately (chamber hypothesis:
    low-height values with dense dim-0 everywhere suggest a piecewise coefficient)."""
    p0 = primes[0]
    all_pts = sorted({(s["ep"], s["r"]) for s in samples.values()})
    out = {
        "boundaries_tried": 0,
        "splits_with_enough_points": 0,
        "validated": None,
        "near_misses": [],
    }
    for a, b, c in BOUNDARY_FORMS:
        pos = [pt for pt in all_pts if a * _uv(pt)[0] + b * _uv(pt)[1] + c > 0]
        neg = [pt for pt in all_pts if a * _uv(pt)[0] + b * _uv(pt)[1] + c < 0]
        out["boundaries_tried"] += 1
        if len(pos) < 10 or len(neg) < 10:
            continue
        out["splits_with_enough_points"] += 1
        exprs = {}
        cells_used = {}
        ok_all = True
        n_sides_ok = 0
        for lab in labels:
            vals = _values_by_point(table[lab], samples)
            per_side = {}
            for side_name, side in (("pos", pos), ("neg", neg)):
                got = _fit_side(vals, side, primes, p0)
                if got is None:
                    break
                per_side[side_name] = got
                n_sides_ok += 1
            if len(per_side) != 2:
                ok_all = False
                break
            exprs[_lab_str(lab)] = {k: v["expr"] for k, v in per_side.items()}
            cells_used[_lab_str(lab)] = {k: v["cell"] for k, v in per_side.items()}
        if ok_all:
            out["validated"] = {
                "boundary": f"({a}*(7*ep) + {b}*(11*r) + {c}) > 0",
                "cells": cells_used,
                "exprs": exprs,
                "n_pos": len(pos),
                "n_neg": len(neg),
            }
            break
        if n_sides_ok >= 6:
            out["near_misses"].append(
                {"boundary": f"({a}*(7*ep) + {b}*(11*r) + {c})", "sides_fit": n_sides_ok}
            )
    if out["validated"] is None:
        out["statement"] = (
            "no tested linear chamber boundary yields exact low-degree rational fits on "
            "both sides for all four coefficients"
        )
    return out


def _fit_side(vals, side, primes, p0):
    """First dense cell <= (3,3) exact on ALL points of one side (last 2 held out of fit)."""
    side_sorted = sorted(side)
    fit_side = side_sorted[:-MIN_VALIDATION]
    svals = {pt: vals[pt] for pt in side_sorted}
    for cell in sorted(
        ((i, j) for i in range(4) for j in range(4)), key=lambda cc: (cc[0] + cc[1], cc[0])
    ):
        num_m, den_m = _monoms_total(cell[0]), _monoms_total(cell[1])
        if len(num_m) + len(den_m) - 1 > len(fit_side):
            continue
        dim, _ = _nullspace_dim_mod(_dense_rows(svals, fit_side, num_m, den_m, p0), p0)
        if dim != 1:
            continue
        exact = _solve_cell_exact(
            svals,
            fit_side,
            side_sorted[-MIN_VALIDATION:],
            num_m,
            den_m,
            primes,
            (Fraction(6), Fraction(57, 11)),
        )
        if exact["ok"]:
            return {"cell": list(cell), "expr": exact["expr"]}
    return None


# --------------------------------------------------------------------------------------
# Phase D: minimal new-sample schedule
# --------------------------------------------------------------------------------------
def phase_d(phase_b_report, phase_c_report, samples, reconstructed_all, fit):
    existing = {(s["ep"], s["r"]) for s in samples.values()}
    existing_ep = sorted({pt[0] for pt in existing})
    existing_r = sorted({pt[1] for pt in existing})

    ansatz_costs = []
    if not reconstructed_all:
        for lab, rep in phase_b_report.items():
            if not isinstance(rep, dict) or "ladder_cells" not in rep:
                continue
            fa = rep["first_potentially_admissible_dense_cell"]
            if fa and fa["status"].startswith("underdetermined"):
                add = (fa["unknowns"] - 1) - 36
                ansatz_costs.append(
                    {
                        "ansatz": f"dense ({fa['num_deg']},{fa['den_deg']}) for {lab}",
                        "equations_per_sample": 1,
                        "additional_samples": add,
                        "total_samples": 38 + add,
                    }
                )
        fa_rep = phase_c_report.get("D_factor_aware", {})
        for lab, e in fa_rep.items():
            if e.get("full_cover") and "validated" not in e:
                # next polynomial degree beyond 7: d=8 -> C(10,2)=45 unknowns
                add = (comb(10, 2) - 1) - 36
                ansatz_costs.append(
                    {
                        "ansatz": f"factor-aware poly degree 8 for {lab}",
                        "equations_per_sample": 1,
                        "additional_samples": add,
                        "total_samples": 38 + add,
                    }
                )
        sh = phase_c_report.get("E_shared_denominator", {})
        if sh.get("validated") is None and sh.get("cells"):
            frontier = next(
                cell
                for cell in sorted(
                    ((x, y) for x in range(10) for y in range(16)),
                    key=lambda cc: (cc[0] + cc[1], cc[0]),
                )
                if 4 * comb(cell[0] + 2, 2) + comb(cell[1] + 2, 2) - 1 > 4 * 36
            )
            a, b = frontier
            unknowns = 4 * comb(a + 2, 2) + comb(b + 2, 2)
            need_fit = -(-(unknowns - 1) // 4)
            add = max(1, need_fit - 36)
            ansatz_costs.append(
                {
                    "ansatz": f"shared-denominator frontier cell ({a},{b}) beyond current data",
                    "equations_per_sample": 4,
                    "additional_samples": add,
                    "total_samples": 38 + add,
                }
            )
    ansatz_costs.sort(key=lambda e: e["additional_samples"])
    minimum = ansatz_costs[0] if ansatz_costs else None

    pool = []
    fiber_eps = existing_ep[:3]
    for k, ep in enumerate(fiber_eps):
        for dm in range(1, 6):
            pool.append((ep, Fraction(60 + dm + 5 * k, 11), f"fixed-ep fiber ep={ep}"))
    fiber_rs = existing_r[:3]
    for k, r in enumerate(fiber_rs):
        for di in range(1, 6):
            pool.append((Fraction(138 + 3 * (di + 5 * k), 7), r, f"fixed-r fiber r={r}"))
    for i in range(1, 25):
        pool.append((Fraction(138 + 3 * i, 7), Fraction(60 + i + 20, 11), "generic scatter"))

    def _excluded(pt):
        ep, r = pt
        if ep == 3 or ep == 6 or r == Fraction(57, 11):
            return True
        return pt in existing

    pool = [(ep, r, why) for ep, r, why in pool if not _excluded((ep, r))]
    seen: set = set()
    pool = [(ep, r, why) for ep, r, why in pool if not ((ep, r) in seen or seen.add((ep, r)))]

    schedule = []
    target = None
    if minimum:
        # schedule the most informative structural target (factor-aware degree-8 bundle
        # if a denominator cover exists); its prefix already decides the cheaper families
        target = next((e for e in ansatz_costs if e["ansatz"].startswith("factor-aware")), minimum)
        n_new = target["additional_samples"]
        p = EXPECTED_PRIMES[-1]
        rng = random.Random(SCHEDULE_SEED)
        # greedy design-rank growth (generic-value proxy: seeded random residues stand
        # in for the unknown coefficient values at new points); existing fit rows kept
        if target["ansatz"].startswith("factor-aware"):
            monoms = _monoms_total(8)
            rows = [_eval_monoms_mod(pt, monoms, p) for pt in fit]
        else:
            m = re.search(r"dense \((\d+),(\d+)\)", target["ansatz"])
            i, j = (int(m.group(1)), int(m.group(2))) if m else (0, 7)
            num_m, den_m = _monoms_total(i), _monoms_total(j)
            rows = []
            for pt in fit:
                vm = rng.randrange(1, p)
                rows.append(
                    _eval_monoms_mod(pt, num_m, p)
                    + [(-vm * x) % p for x in _eval_monoms_mod(pt, den_m, p)]
                )
            monoms = None
        rank0, _, _ = _row_reduce_mod(rows, p)
        for ep, r, why in pool:
            if len(schedule) >= n_new:
                break
            if monoms is not None:
                row = _eval_monoms_mod((ep, r), monoms, p)
            else:
                vm = rng.randrange(1, p)
                row = _eval_monoms_mod((ep, r), num_m, p) + [
                    (-vm * x) % p for x in _eval_monoms_mod((ep, r), den_m, p)
                ]
            rank1, _, _ = _row_reduce_mod(rows + [row], p)
            if rank1 > rank0:
                rows.append(row)
                rank0 = rank1
                schedule.append(
                    {
                        "ep": _frac_str(ep),
                        "r": _frac_str(r),
                        "reason": why,
                        "design_rank_after": rank1,
                    }
                )
    return {
        "reconstructed_from_existing_data": bool(reconstructed_all),
        "viable_ansatz_costs": ansatz_costs,
        "minimum": minimum,
        "schedule_target": target,
        "recommended_new_samples": schedule,
        "notes": [
            "candidate points avoid ep=3, ep=6, r=57/11 and all existing samples",
            "fiber points additionally enable univariate cross-checks (phases C.A/C.B)",
            "greedy rank growth uses seeded random residues as generic-position proxies "
            f"(seed {SCHEDULE_SEED}); ordering is fully deterministic",
        ],
    }


# --------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------
def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--records", type=Path, default=RECORDS_PATH)
    parser.add_argument("--artifact", type=Path, default=ARTIFACT_PATH)
    parser.add_argument("--value-table-out", type=Path, default=VALUE_TABLE_OUT)
    parser.add_argument("--audit-out", type=Path, default=AUDIT_OUT)
    parser.add_argument("--note-out", type=Path, default=NOTE_OUT)
    args = parser.parse_args(argv)

    global header_fp
    header, records, raw_rows = load_records(args.records)
    header_fp = {"family_sha256": header["fingerprint"].get("family_sha256")}
    artifact = (
        json.loads(args.artifact.read_text(encoding="utf-8")) if args.artifact.exists() else None
    )

    table_doc, labels, table, samples, classification = build_value_table(records, raw_rows)
    args.value_table_out.write_text(
        json.dumps(table_doc, ensure_ascii=False, indent=1), encoding="utf-8"
    )

    primes = list(EXPECTED_PRIMES)
    b_report, fit, hold = phase_b(labels, table, samples, primes, artifact)
    c_report = phase_c(labels, table, samples, primes, fit, hold)

    reconstructed = {}
    for lab in labels:
        ls = _lab_str(lab)
        got = None
        fa = c_report["D_factor_aware"].get(ls, {})
        if "validated" in fa:
            got = {"via": "factor-aware", **fa["validated"]}
        for cand in c_report["C_sparse_tensor_supports"].get(ls, {}).get("validated", []):
            got = got or {"via": "sparse-tensor", "expr": cand["expr"]}
        sh = c_report["E_shared_denominator"].get("validated")
        if sh and ls in sh["exprs"]:
            got = got or {"via": "shared-denominator", "expr": sh["exprs"][ls]}
        for cell in b_report[ls]["ladder_cells"]:
            if cell.get("validated_candidate"):
                got = got or {"via": "dense", "expr": cell["validated_candidate"]}
        pw = c_report.get("G_piecewise_chamber_split", {}).get("validated")
        if pw and ls in pw.get("exprs", {}):
            got = got or {
                "via": "piecewise-boundary",
                "boundary": pw["boundary"],
                "expr": pw["exprs"][ls],
            }
        reconstructed[ls] = got
    all_done = all(reconstructed[_lab_str(lab)] is not None for lab in labels)

    d_report = phase_d(b_report, c_report, samples, all_done, fit)

    audit = {
        "schema": "m11c-reconstruction-audit/v1",
        "inputs": {
            "records": args.records.name,
            "artifact": args.artifact.name,
            "n_records": len(records),
            "n_samples": len(samples),
            "primes": primes,
            "selected_rank": EXPECTED_RANK,
        },
        "no_rref_or_reducer_solve": True,
        "phase_a_verification": table_doc["verification"],
        "phase_b_degree_diagnostics": b_report,
        "phase_c_structured_experiments": c_report,
        "reconstructed_coefficients": reconstructed,
        "phase_d_schedule": d_report,
    }
    args.audit_out.write_text(json.dumps(audit, ensure_ascii=False, indent=1), encoding="utf-8")
    args.note_out.write_text(_render_note(audit), encoding="utf-8")

    print(_render_summary(audit))
    return 0


def _render_note(audit) -> str:
    lines = ["# External Int2 Method.11c: offline reconstruction audit", ""]
    va = audit["phase_a_verification"]
    lines += [
        "Offline audit of the Method.11 record cache (38 samples x 4 primes); no new modular",
        "records, no RREF/reducer solve. Artifacts: `external_int2_method11_value_table.json`,",
        "`external_int2_method11_reconstruction_audit.json`.",
        "",
        "## Phase A - value table",
        f"- samples: {va['n_samples']}, records: {va['n_records']}, primes: {va['primes']}",
        f"- selected rank: {va['selected_rank']}; support size: {va['support_size']}",
        f"- special-zero: {list(va['special_zero'])} missing "
        f"{[e['missing_labels'] for e in va['special_zero'].values()]}",
        f"- verification: {'OK' if va['ok'] else va['problems']}",
        "",
        "## Phase B - dense degree diagnostics",
    ]
    for lab, rep in audit["phase_b_degree_diagnostics"].items():
        if not isinstance(rep, dict) or "ladder_cells" not in rep:
            continue
        fa = rep["first_potentially_admissible_dense_cell"]
        lines.append(f"- `{lab}`: first admissible dense cell {fa}")
        if "failure_explanation" in rep:
            lines += ["", f"  {rep['failure_explanation']}", ""]
    xc = audit["phase_b_degree_diagnostics"].get("artifact_crosscheck", {})
    lines.append(f"- Stage-4 artifact crosscheck: {xc}")
    lines += ["", "## Phase C - structured experiments"]
    c = audit["phase_c_structured_experiments"]
    lines.append(f"- A: {c['A_univariate_ep_on_fixed_r']['statement']}")
    lines.append(f"- B: {c['B_univariate_r_on_fixed_ep']['statement']}")
    for lab, t in c["C_sparse_tensor_supports"].items():
        lines.append(
            f"- C `{lab}`: {t['cells_tried']} tensor cells, {t['inconsistent']} inconsistent, "
            f"{t['ambiguous']} ambiguous, validated: {len(t['validated'])}"
        )
    for lab, e in c["D_factor_aware"].items():
        got = e.get("validated")
        lines.append(
            f"- D `{lab}`: forms {[f['form'] for f in e.get('denominator_forms', [])]} "
            f"cover={e.get('full_cover')} -> "
            + (f"VALIDATED deg {got['poly_degree']}" if got else e.get("statement", "no"))
        )
    sh = c["E_shared_denominator"]
    lines.append(f"- E: {sh.get('statement') or 'validated: ' + str(sh['validated']['cell'])}")
    g = c.get("G_piecewise_chamber_split", {})
    if g.get("validated"):
        lines.append(
            f"- G: VALIDATED piecewise split at {g['validated']['boundary']} "
            f"(pos {g['validated']['n_pos']} / neg {g['validated']['n_neg']} points; "
            f"cells {g['validated']['cells']})"
        )
    else:
        lines.append(
            f"- G: {g.get('statement', 'n/a')} ({g.get('boundaries_tried', 0)} boundaries, "
            f"{g.get('splits_with_enough_points', 0)} usable splits, "
            f"near-misses {g.get('near_misses', [])})"
        )
    lines += ["", "## Reconstructed coefficients"]
    for lab, got in audit["reconstructed_coefficients"].items():
        lines.append(f"- `{lab}`: " + (f"YES via {got['via']}" if got else "not reconstructed"))
    d = audit["phase_d_schedule"]
    lines += ["", "## Phase D - minimal new-sample schedule"]
    if d["reconstructed_from_existing_data"]:
        lines.append("- all four coefficients reconstructed from existing data; Stage-5 NOT needed")
    else:
        lines.append(f"- viable ansatz costs: {d['viable_ansatz_costs']}")
        lines.append(f"- minimum: {d['minimum']}")
        lines.append(f"- recommended points: {d['recommended_new_samples']}")
    lines += ["", "Notes: " + "; ".join(d["notes"]), ""]
    return "\n".join(lines)


def _render_summary(audit) -> str:
    va = audit["phase_a_verification"]
    d = audit["phase_d_schedule"]
    rec = audit["reconstructed_coefficients"]
    out = [
        f"samples {va['n_samples']} x primes {len(va['primes'])} = {va['n_records']} records; "
        f"rank {va['selected_rank']}; verification {'OK' if va['ok'] else 'PROBLEMS'}",
        f"reconstructed: {[lab for lab, g in rec.items() if g]}",
        f"not reconstructed: {[lab for lab, g in rec.items() if not g]}",
    ]
    if not d["reconstructed_from_existing_data"] and d["minimum"]:
        out.append(
            f"minimum new samples: {d['minimum']['additional_samples']} ({d['minimum']['ansatz']})"
        )
    return "\n".join(out)


if __name__ == "__main__":
    raise SystemExit(main())
