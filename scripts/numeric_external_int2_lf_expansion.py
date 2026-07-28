"""External Int2 Method.12N: numerical epsilon expansion of the four certified LF masters.

Computes j_{k,n} = 1/n! * Integral[ f_k(x,0) * L_k(x)^n , x in (0,inf)^3 ]
directly (no fitting in ep), assembles the certified reduction with the exact
C1..C4 and the corrected external prefactor, and freezes the resulting Laurent
coefficients BEFORE any oracle access.

Stages:
  compute  -- numeric integration + assembly; writes the frozen artifact.
  compare  -- reads the frozen artifact, verifies the freeze hash, and only
              then reads the analytic Laurent oracle for comparison.

No RREF / modular solving / record generation is imported or performed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
MPL_FILES = [
    REPO_ROOT / "validation" / "hyperint" / f"external_int2_lf_L{k}.mpl"
    for k in (1, 2, 3, 4)
]
CERTIFICATE = REPO_ROOT / "validation" / "external_int2_lf_certificate.json"
DEFAULT_OUT = REPO_ROOT / "validation" / "external_int2_lf_numeric_expansion.json"

# (s, t) kinematic points; r = s/t enters the masters, t the prefactor.
DEFAULT_POINTS = [("1", "1"), ("3/2", "1"), ("1", "3/2")]

# quadrature configs: (name, min_exp of graded mesh, Gauss-Legendre order)
CONFIGS = {
    "base": (8, 12),
    "precision_up": (8, 20),
    "subdivision_up": (12, 12),
    "final": (12, 20),
}
FAST_CONFIGS = {
    "base": (5, 6),
    "precision_up": (5, 10),
    "subdivision_up": (8, 6),
    "final": (8, 10),
}


# --------------------------------------------------------------------------
# master parsing (from the certified Method.12A .mpl audit artifacts)
# --------------------------------------------------------------------------

def parse_masters():
    """Parse the four .mpl inputs into structural exponent data.

    Every master must factor as  prod_i base_i^(a_i*ep + b_i)  with bases in
    {x2, x5, x7, 1+x2, 1+x5, 1+x7, G3}, G3 = r*x2*x5 + x2*x7 + x7 + 1.
    Returns list of dicts with label, eps_order, exponent tables.
    """
    import sympy as sp

    ep = sp.Symbol("ep")  # plain symbol: must be identical to sympify's own
    base_names = {
        sp.sympify("x2"): "x2",
        sp.sympify("x5"): "x5",
        sp.sympify("x7"): "x7",
        sp.sympify("x2 + 1"): "G0",
        sp.sympify("x5 + 1"): "G1",
        sp.sympify("x7 + 1"): "G2",
        sp.sympify("r*x2*x5 + x2*x7 + x7 + 1"): "G3",
    }
    masters = []
    for path in MPL_FILES:
        text = path.read_text(encoding="utf-8")
        mlab = re.search(r"# label \[([-0-9, ]+)\]", text)
        meps = re.search(r"epsOrder := (\d+):", text)
        mf = re.search(r"^f := (.+):$", text, re.M)
        if not (mlab and meps and mf):
            raise SystemExit(f"cannot parse {path}")
        label = [int(v) for v in mlab.group(1).split(",")]
        eps_order = int(meps.group(1))
        expr = sp.sympify(mf.group(1).replace("^", "**"))
        factors = expr.args if expr.func is sp.Mul else (expr,)
        a_tab = {}
        b_tab = {}
        for fac in factors:
            if fac.func is sp.Pow:
                base, e = fac.args
            else:
                base, e = fac, sp.Integer(1)
            name = None
            for cand, nm in base_names.items():
                if sp.simplify(base - cand) == 0:
                    name = nm
                    break
            if name is None:
                raise SystemExit(f"unexpected base {base} in {path}")
            a = e.coeff(ep)
            b = e.subs(ep, 0)
            if sp.simplify(e - (a * ep + b)) != 0:
                raise SystemExit(f"non-affine ep exponent {e} in {path}")
            a_tab[name] = a_tab.get(name, 0) + int(a)
            b_tab[name] = b_tab.get(name, 0) + int(b)
        a_tab = {k: v for k, v in a_tab.items() if v != 0}
        b_tab = {k: v for k, v in b_tab.items() if v != 0}
        masters.append(
            {
                "file": str(path.relative_to(REPO_ROOT)).replace("\\", "/"),
                "label": label,
                "eps_order": eps_order,
                "a": a_tab,
                "b": b_tab,
                "integrand": mf.group(1),
            }
        )
    # common log factor L(x): identical a-tables required
    a0 = masters[0]["a"]
    for m in masters[1:]:
        if m["a"] != a0:
            raise SystemExit("masters do not share a common log factor L(x)")
    return masters


def load_certified_coefficients(masters):
    cert = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    claimed = cert["claimed_coefficients"]
    coeffs = []
    for m in masters:
        key = str(m["label"])
        if key not in claimed:
            raise SystemExit(f"label {key} not in certificate")
        coeffs.append(claimed[key])
    return coeffs


# --------------------------------------------------------------------------
# quadrature engine (numpy, graded Gauss-Legendre tensor mesh)
# --------------------------------------------------------------------------

def graded_axis(min_exp: int, order: int):
    """Return quadrature nodes/weights on (0,1) with dyadic grading at 0 and 1."""
    lo = [2.0 ** -k for k in range(min_exp, 0, -1)]
    hi = [1.0 - 2.0 ** -k for k in range(2, min_exp + 1)]
    bps = np.array([0.0] + lo + hi + [1.0])
    xg, wg = np.polynomial.legendre.leggauss(order)
    nodes, weights = [], []
    for a, b in zip(bps[:-1], bps[1:]):
        h = 0.5 * (b - a)
        nodes.append(h * xg + 0.5 * (a + b))
        weights.append(h * wg)
    return np.concatenate(nodes), np.concatenate(weights)


def sweep(masters, r: float, min_exp: int, order: int, n_max: int = 4):
    """One full tensor sweep: returns j_raw[k][n] = Integral[f_k0 * L^n] (no 1/n!)."""
    u, w = graded_axis(min_exp, order)
    x = u / (1.0 - u)
    wt = w / (1.0 - u) ** 2  # weight * Jacobian dx/du = 1/(1-u)^2
    lnx = np.log(x)
    ln1p = np.log1p(x)
    a = masters[0]["a"]

    # per-axis pieces of L = a_x2*ln x2 + a_G0*ln(1+x2) + a_G1*ln(1+x5)
    #                        + a_x7*ln x7 + a_G2*ln(1+x7) + a_G3*ln G3
    l2 = a.get("x2", 0) * lnx + a.get("G0", 0) * ln1p
    l5 = a.get("x5", 0) * lnx + a.get("G1", 0) * ln1p
    l7 = a.get("x7", 0) * lnx + a.get("G2", 0) * ln1p
    a_g3 = a.get("G3", 0)

    # per-axis ep^0 factors (weights folded in), and G3 exponent per master
    ax2, ax5, ax7, mG3 = [], [], [], []
    for m in masters:
        b = m["b"]
        f2 = x ** b.get("x2", 0) * (1.0 + x) ** b.get("G0", 0)
        f5 = x ** b.get("x5", 0) * (1.0 + x) ** b.get("G1", 0)
        f7 = x ** b.get("x7", 0) * (1.0 + x) ** b.get("G2", 0)
        ax2.append(wt * f2)
        ax5.append(wt * f5)
        ax7.append(wt * f7)
        if b.get("G3", 0) >= 0:
            raise SystemExit("expected negative G3 exponent")
        mG3.append(-b["G3"])

    n = x.size
    acc = np.zeros((len(masters), n_max + 1))
    chunk = max(1, 6_000_000 // (n * n))
    rx5 = r * x
    for i0 in range(0, n, chunk):
        i1 = min(i0 + chunk, n)
        xb = x[i0:i1]
        g3 = (
            xb[:, None, None] * rx5[None, :, None]
            + (1.0 + xb)[:, None, None] * x[None, None, :]
            + 1.0
        )
        lng3 = np.log(g3)
        inv = 1.0 / g3
        big_l = (
            l2[i0:i1][:, None, None]
            + l5[None, :, None]
            + l7[None, None, :]
            + a_g3 * lng3
        )
        invsq = inv * inv
        for k in range(len(masters)):
            core = inv if mG3[k] == 1 else invsq
            t_k = (
                ax2[k][i0:i1][:, None, None]
                * ax5[k][None, :, None]
                * ax7[k][None, None, :]
                * core
            )
            pw = None
            for nn in range(0, n_max + 1):
                if pw is None:
                    pw = np.ones_like(t_k)
                else:
                    pw = pw * big_l
                acc[k, nn] += float(np.sum(t_k * pw))
    return acc


# --------------------------------------------------------------------------
# high-precision spot check (mpmath tanh-sinh, n=0 coefficients)
# --------------------------------------------------------------------------

def mpmath_j0(masters, k: int, r, dps: int = 25):
    """Nested tanh-sinh integration of f_k0 over the u-cube at precision dps."""
    import mpmath as mp

    m = masters[k]
    b = m["b"]
    with mp.workdps(dps):
        rr = mp.mpf(r.numerator) / mp.mpf(r.denominator)

        def fk0(x2, x5, x7):
            g3 = rr * x2 * x5 + x2 * x7 + x7 + 1
            val = (
                x2 ** b.get("x2", 0)
                * (1 + x2) ** b.get("G0", 0)
                * (1 + x5) ** b.get("G1", 0)
                * x7 ** b.get("x7", 0)
                * (1 + x7) ** b.get("G2", 0)
                * g3 ** b.get("G3", 0)
            )
            return val

        def integrand(u2, u5, u7):
            x2 = u2 / (1 - u2)
            x5 = u5 / (1 - u5)
            x7 = u7 / (1 - u7)
            jac = 1 / ((1 - u2) ** 2 * (1 - u5) ** 2 * (1 - u7) ** 2)
            return fk0(x2, x5, x7) * jac

        val = mp.quad(
            lambda u2: mp.quad(
                lambda u5: mp.quad(lambda u7: integrand(u2, u5, u7), [0, 1]),
                [0, 1],
            ),
            [0, 1],
        )
        return float(val)


# --------------------------------------------------------------------------
# exact series pieces: certified C_k(ep) and external prefactor P(ep)
# --------------------------------------------------------------------------

def coefficient_series(coeff_str: str, r: Fraction, lo: int = -2, hi: int = 4):
    """Exact Laurent coefficients of a certified C_k at rational r via sympy."""
    import sympy as sp

    ep = sp.Symbol("ep")
    expr = sp.sympify(coeff_str).subs(sp.Symbol("r"), sp.Rational(r.numerator, r.denominator))
    shifted = sp.expand(sp.cancel(expr * ep ** (-lo)))  # analytic at ep=0
    ser = sp.series(shifted, ep, 0, hi - lo + 1).removeO()
    poly = sp.Poly(ser, ep)
    out = {}
    for i in range(0, hi - lo + 1):
        c = poly.coeff_monomial(ep ** i) if i else poly.coeff_monomial(1)
        c = sp.Rational(c)
        out[i + lo] = Fraction(int(c.p), int(c.q))
    return out


def prefactor_series(t: Fraction, n_terms: int = 8):
    """Laurent coefficients of ExternalPrefactor2 via exact series arithmetic.

    P(ep) = Exp[2 ep EulerGamma] t^(-3-ep) Gamma[1-ep] Gamma[-ep]^3 Gamma[ep]
            / (Gamma[-1-3 ep] Gamma[-2 ep])
          = (6 (1+3 ep) / ep^2) * t^(-3-ep) e^{2 gamma ep}
            * Gamma[1-ep]^4 Gamma[1+ep] / (Gamma[1-3 ep] Gamma[1-2 ep]).
    Returns {m: float} for m in [-2, n_terms-3].
    """
    import mpmath as mp

    with mp.workdps(60):
        n = n_terms
        lnt = mp.log(mp.mpf(t.numerator) / mp.mpf(t.denominator))

        def lngamma1p(c):
            # series of ln Gamma(1 + c*ep)
            out = [mp.mpf(0)] * n
            out[1] = -mp.euler * c
            for m in range(2, n):
                out[m] = (-1) ** m * mp.zeta(m) / m * mp.mpf(c) ** m
            return out

        def add(s1, s2, f=1):
            return [a + f * b for a, b in zip(s1, s2)]

        # ln(1+3 ep)
        log13 = [mp.mpf(0)] * n
        for m in range(1, n):
            log13[m] = (-1) ** (m + 1) * mp.mpf(3) ** m / m
        ln_a = [mp.mpf(0)] * n
        ln_a[0] = -3 * lnt
        ln_a[1] = -lnt + 2 * mp.euler
        ln_a = add(ln_a, log13)
        ln_a = add(ln_a, lngamma1p(-1), 4)
        ln_a = add(ln_a, lngamma1p(1), 1)
        ln_a = add(ln_a, lngamma1p(-3), -1)
        ln_a = add(ln_a, lngamma1p(-2), -1)
        # exp of series
        a_ser = [mp.mpf(0)] * n
        a_ser[0] = mp.exp(ln_a[0])
        for m in range(1, n):
            s = mp.mpf(0)
            for j in range(1, m + 1):
                s += j * ln_a[j] * a_ser[m - j]
            a_ser[m] = s / m
        # P = 6/ep^2 * A(ep)
        return {m - 2: float(6 * a_ser[m]) for m in range(n)}


# --------------------------------------------------------------------------
# assembly
# --------------------------------------------------------------------------

def assemble(masters, coeff_strs, j_val, j_err, r: Fraction, t: Fraction):
    """Full(ep) = P(ep) * sum_k C_k(ep) * J_k(ep); Laurent orders -4..0.

    j_val/j_err: dict (k, n) -> value/abs-error.  Linear error propagation.
    """
    pref = prefactor_series(t)
    c_ser = [coefficient_series(cs, r) for cs in coeff_strs]
    full = {}
    full_err = {}
    for order in range(-4, 1):
        tot = 0.0
        tot_err = 0.0
        for k, m in enumerate(masters):
            for nn in range(0, m["eps_order"] + 1):
                w = 0.0
                for i, ci in c_ser[k].items():
                    pm = order - i - nn
                    if pm in pref:
                        w += float(ci) * pref[pm]
                tot += j_val[(k, nn)] * w
                tot_err += abs(w) * j_err[(k, nn)]
        full[order] = tot
        full_err[order] = tot_err
    return full, full_err


# --------------------------------------------------------------------------
# compute stage
# --------------------------------------------------------------------------

GROWTH_DEPTHS = (8, 12, 16, 20)


def classify_convergence(masters, r: Fraction, order: int = 12):
    """Term-wise convergence screen: j-integrals vs dyadic cutoff depth.

    Returns (status[k], growth) where status is 'converged' or
    'divergent_termwise'; growth holds the raw n=0..n_k tables and increments.
    """
    tables = {d: sweep(masters, float(r), d, order) for d in GROWTH_DEPTHS}
    status, growth = [], []
    d_last = GROWTH_DEPTHS[-1]
    for k, m in enumerate(masters):
        rows = {}
        divergent = False
        for nn in range(0, m["eps_order"] + 1):
            vals = [tables[d][k, nn] / math.factorial(nn) for d in GROWTH_DEPTHS]
            incs = [b - a for a, b in zip(vals[:-1], vals[1:])]
            rows[str(nn)] = {
                "depths": list(GROWTH_DEPTHS),
                "values": vals,
                "increments": incs,
            }
            # divergent <=> tail increments do NOT decay with cutoff depth
            # (convergent tails shrink geometrically; log-divergence keeps
            # them constant, log^n>1 makes them grow)
            last_inc, prev_inc = abs(incs[-1]), abs(incs[-2])
            scale = max(abs(tables[d_last][k, nn]), 1.0)
            if last_inc > 1e-9 * scale and last_inc >= 0.5 * prev_inc:
                divergent = True
        if divergent:
            # fitted log-slope of the n=0 growth, in units of ln(2) per depth
            n0 = rows["0"]
            slope = n0["increments"][-1] / ((GROWTH_DEPTHS[-1] - GROWTH_DEPTHS[-2]) * math.log(2))
            rows["log_slope_n0"] = slope
            rows["predicted_slope_1_over_r"] = float(1 / r)
        status.append("divergent_termwise" if divergent else "converged")
        growth.append(rows)
    return status, growth


def run_compute_stage(out_path: Path, fast: bool, skip_mpmath: bool):
    masters = parse_masters()
    coeff_strs = load_certified_coefficients(masters)
    cfgs = FAST_CONFIGS if fast else CONFIGS
    points = [(Fraction(s), Fraction(t)) for s, t in DEFAULT_POINTS]
    r_values = sorted({s / t for s, t in points})

    per_r: dict[Fraction, dict] = {}
    for r in r_values:
        t0 = time.time()
        status, growth = classify_convergence(masters, r)
        runs = {}
        for name, (min_exp, order) in cfgs.items():
            runs[name] = sweep(masters, float(r), min_exp, order)
        elapsed = time.time() - t0
        coeffs = {}
        for k, m in enumerate(masters):
            if status[k] != "converged":
                continue
            for nn in range(0, m["eps_order"] + 1):
                fac = math.factorial(nn)
                vals = {name: runs[name][k, nn] / fac for name in cfgs}
                value = vals["final"]
                stab_prec = abs(vals["precision_up"] - vals["base"])
                stab_sub = abs(vals["subdivision_up"] - vals["base"])
                err = max(
                    abs(vals["final"] - vals["precision_up"]),
                    abs(vals["final"] - vals["subdivision_up"]),
                    1e-15 * abs(value),
                )
                coeffs[(k, nn)] = {
                    "value": value,
                    "err_est": err,
                    "stab_precision": stab_prec,
                    "stab_subdivision": stab_sub,
                }
        per_r[r] = {
            "coeffs": coeffs,
            "status": status,
            "growth": growth,
            "elapsed_seconds": round(elapsed, 2),
        }

    if not skip_mpmath:
        for r in r_values:
            for k in range(len(masters)):
                if per_r[r]["status"][k] != "converged":
                    continue
                if r != Fraction(1) and k > 0:
                    continue  # full set at r=1; k=0 spot check at other r
                ref = mpmath_j0(masters, k, r)
                ent = per_r[r]["coeffs"][(k, 0)]
                ent["mpmath_dps25_check"] = ref
                ent["mpmath_rel_diff"] = abs(ent["value"] - ref) / abs(ref)

    all_ok = all(
        st == "converged" for r in r_values for st in per_r[r]["status"]
    )
    point_blocks = []
    for s, t in points:
        r = s / t
        coeffs = per_r[r]["coeffs"]
        pm = {}
        for k, m in enumerate(masters):
            label = str(m["label"])
            if per_r[r]["status"][k] == "converged":
                pm[label] = {
                    "status": "converged",
                    "coefficients": {
                        str(nn): coeffs[(k, nn)]
                        for nn in range(0, m["eps_order"] + 1)
                    },
                }
            else:
                pm[label] = {
                    "status": "divergent_termwise",
                    "growth_diagnostic": per_r[r]["growth"][k],
                    "analysis": (
                        "log-divergent at the joint (x5,x7)->infinity face; "
                        "J(ep) = -1/(r*ep) + O(1) as ep->0-, so a term-wise "
                        "Taylor expansion at ep=0 does not exist"
                    ),
                }
        blk = {
            "s": str(s),
            "t": str(t),
            "r": str(r),
            "per_master": pm,
            "sweep_elapsed_seconds": per_r[r]["elapsed_seconds"],
        }
        if all_ok:
            j_val = {kn: c["value"] for kn, c in coeffs.items()}
            j_err = {kn: c["err_est"] for kn, c in coeffs.items()}
            full, full_err = assemble(masters, coeff_strs, j_val, j_err, r, t)
            blk["assembled_laurent"] = {str(o): full[o] for o in sorted(full)}
            blk["assembled_err"] = {str(o): full_err[o] for o in sorted(full_err)}
        point_blocks.append(blk)

    frozen = {
        "schema": "external-int2-lf-numeric-expansion/v1",
        "script": "scripts/numeric_external_int2_lf_expansion.py",
        "arithmetic": "float64 quadrature; exact-Fraction C_k series; mpmath dps=60 prefactor series; mpmath dps=25 tanh-sinh spot checks",
        "quadrature": {
            "transform": "x = u/(1-u), Jacobian 1/(1-u)^2, per axis",
            "mesh": "dyadic graded Gauss-Legendre tensor mesh",
            "configs": {k: {"min_exp": v[0], "gl_order": v[1]} for k, v in cfgs.items()},
        },
        "masters": [
            {
                "file": m["file"],
                "label": m["label"],
                "eps_order": m["eps_order"],
                "integrand": m["integrand"],
                "coefficient": coeff_strs[k],
            }
            for k, m in enumerate(masters)
        ],
        "points": point_blocks,
        "all_masters_termwise_convergent": all_ok,
        "verdict": (
            "COMPUTED_AND_FROZEN"
            if all_ok
            else "FAILURE_TERMWISE_EXPANSION_DIVERGENT"
        ),
    }
    blob = json.dumps(frozen, sort_keys=True).encode()
    doc = {
        "frozen": frozen,
        "frozen_sha256": hashlib.sha256(blob).hexdigest(),
        "oracle_read_during_compute": False,
    }
    if not all_ok:
        doc["comparison"] = {
            "status": "not_run",
            "reason": (
                "assembly impossible: at least one master has no term-wise "
                "epsilon expansion at ep=0 (divergent coefficient integrals); "
                "comparing a partial sum against the oracle would be meaningless"
            ),
        }
    out_path.write_text(json.dumps(doc, indent=1), encoding="utf-8")
    print(f"[12N] compute stage done -> {out_path}")
    print(f"[12N] verdict: {frozen['verdict']}")
    print(f"[12N] frozen_sha256 = {doc['frozen_sha256']}")
    return doc


# --------------------------------------------------------------------------
# comparison stage (the ONLY place that reads the analytic Laurent oracle)
# --------------------------------------------------------------------------

def run_comparison_stage(artifact_path: Path):
    doc = json.loads(artifact_path.read_text(encoding="utf-8"))
    blob = json.dumps(doc["frozen"], sort_keys=True).encode()
    if hashlib.sha256(blob).hexdigest() != doc["frozen_sha256"]:
        raise SystemExit("freeze hash mismatch: refusing to compare")
    if not doc["frozen"].get("all_masters_termwise_convergent"):
        raise SystemExit(
            "refusing to compare: frozen artifact records divergent term-wise "
            "expansions; there is no assembled Laurent series to compare"
        )
    oracle_path = REPO_ROOT / "validation" / "external_int2_full_laurent_audit.json"
    oracle = json.loads(oracle_path.read_text(encoding="utf-8"))

    import sympy as sp

    s_sym, t_sym, ep = sp.symbols("s t ep", positive=True)
    coeff_exprs = load_oracle_laurent(oracle, sp, s_sym, t_sym, ep)

    comparison = {
        "oracle_file": str(oracle_path.relative_to(REPO_ROOT)).replace("\\", "/"),
        "per_point": [],
    }
    worst = 0.0
    for blk in doc["frozen"]["points"]:
        sv, tv = sp.Rational(blk["s"]), sp.Rational(blk["t"])
        rows = {}
        for order in range(-4, 1):
            oc = coeff_exprs[order].subs({s_sym: sv, t_sym: tv})
            oc = float(sp.N(oc, 40))
            cv = blk["assembled_laurent"][str(order)]
            err = blk["assembled_err"][str(order)]
            resid = cv - oc
            rel = abs(resid) / max(abs(oc), 1e-300)
            worst = max(worst, rel)
            rows[str(order)] = {
                "oracle": oc,
                "computed": cv,
                "abs_residual": resid,
                "rel_residual": rel,
                "propagated_err": err,
            }
        comparison["per_point"].append({"s": blk["s"], "t": blk["t"], "orders": rows})
    comparison["max_rel_residual"] = worst
    comparison["verdict"] = "PASS" if worst < 5e-7 else "FAIL"
    doc["comparison"] = comparison
    artifact_path.write_text(json.dumps(doc, indent=1), encoding="utf-8")
    print(f"[12N] comparison verdict: {comparison['verdict']}; max rel residual {worst:.3e}")
    return doc


def load_oracle_laurent(oracle, sp, s_sym, t_sym, ep):
    """Extract exact Laurent coefficient expressions F_-4..F_0 from the oracle JSON."""
    raise SystemExit(
        "oracle schema loader not yet implemented; inspect the oracle file "
        "AFTER the freeze and fill in this function"
    )


def run_divergence_diagnostic(r=Fraction(1), depths=(8, 12, 16, 20), order=12):
    """Cutoff-growth experiment: j_{k,0} vs dyadic grading depth.

    The u-mesh grading 2^-d places the effective x-cutoff at ~2^d, so a
    log-divergent coefficient integral grows LINEARLY in d (slope ln2 * c),
    while convergent ones plateau.  Prints and returns the table.
    """
    masters = parse_masters()
    table = {}
    for d in depths:
        acc = sweep(masters, float(r), d, order, n_max=0)
        table[d] = [acc[k, 0] for k in range(4)]
        print(f"[12N diag] depth {d:2d}: " + "  ".join(f"j{k+1},0={acc[k,0]:.6f}" for k in range(4)))
    print("[12N diag] increments per +4 depth (log-divergent => ~const 4*ln2/r):")
    ds = list(depths)
    for a, b in zip(ds[:-1], ds[1:]):
        inc = [table[b][k] - table[a][k] for k in range(4)]
        print(f"[12N diag] d {a}->{b}: " + "  ".join(f"dj{k+1}={inc[k]:+.6f}" for k in range(4)))
    return table


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("stage", choices=["compute", "compare", "diagnose"])
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--fast", action="store_true", help="small smoke configuration")
    ap.add_argument("--skip-mpmath", action="store_true")
    args = ap.parse_args(argv)
    if args.stage == "compute":
        run_compute_stage(args.out, args.fast, args.skip_mpmath)
    elif args.stage == "diagnose":
        run_divergence_diagnostic()
    else:
        run_comparison_stage(args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
