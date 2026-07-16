"""Standalone Laurent audit of External Int1 (corrected).

Pure mathematics: this script never imports parametric_ibp_lf_reducer and does
not touch the certified reduction or its coefficients.  J1 and J2 are evaluated
directly from their integral definitions.  It audits

    Full(ep) = Prefactor(ep) * (A(ep)*J1(ep) + B(ep)*J2(ep))        (x 1/(s*t^2))

against the reference target

    1/ep^4 - pi^2/(12*ep^2) - 43*zeta(3)/(6*ep) - pi^4/180          (x 1/(s*t^2)).

Method
------
1. The inner x2-integral of J1/J2 is reduced exactly to a Gauss 2F1 kernel
   (derivation in notes/EXTERNAL_INT1_LAURENT_AUDIT.md), so J1(ep), J2(ep)
   become fast 1-D quadratures valid for complex ep.  The reduction is
   validated against a direct 2-D quadrature.
2. Taylor coefficients of J1 (to ep^3) and J2 (to ep^4) at ep = 0 are
   extracted with the Cauchy formula on the circle |ep| = 1/32.
3. Every coefficient is identified with PSLQ in the weight-graded basis
   {1, pi^2, zeta(3), pi^4, zeta(5), pi^2*zeta(3)} (order-k coefficient may
   carry weight up to k+1); residuals are reported, and the analytically
   known anchors (J1(0), J1'(0), J2(0), J2'(0), J2''(0)/2) are cross-checked.
4. The exact series assembly (A*J1 + B*J2, external prefactor, full Laurent
   expansion through ep^0) and the target comparison are done in sympy.

Usage:  python scripts/audit_external_int1_laurent.py
Output: report on stdout + validation/external_int1_laurent_audit.json
"""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
import sympy as sp

REPO_ROOT = Path(__file__).resolve().parents[1]
OUT_JSON = REPO_ROOT / "validation" / "external_int1_laurent_audit.json"

mp.mp.dps = 45

EP = sp.Symbol("ep")
PI = sp.pi
Z3 = sp.zeta(3)
Z5 = sp.zeta(5)

A_EXPR = (4 * EP - 1) / (3 * (3 * EP + 1))
B_EXPR = (EP - 2) * (5 * EP - 2) / (3 * EP * (3 * EP + 1))

TARGET = {
    -4: sp.Integer(1),
    -3: sp.Integer(0),
    -2: -(PI**2) / 12,
    -1: -sp.Rational(43, 6) * Z3,
    0: -(PI**4) / 180,
}

# Analytically derived anchors (see the notes file for the derivations).
KNOWN = {
    (1, 0): sp.Integer(1),  # J1(0)
    (1, 1): sp.Integer(3) + PI**2 / 6,  # J1'(0)
    (2, 0): sp.Rational(1, 2),  # J2(0)
    (2, 1): sp.Rational(7, 4),  # J2'(0)
    (2, 2): sp.Rational(35, 8) + PI**2 / 12,  # J2''(0)/2
}

# Hand-derived prefactor Laurent coefficients (independent cross-check).
EXPECTED_PREFACTOR = {
    -3: sp.Rational(3, 2),
    -2: sp.Rational(9, 2),
    -1: -(PI**2) / 4,
    0: -3 * PI**2 / 4 - 16 * Z3,
    1: -48 * Z3 - sp.Rational(19, 80) * PI**4,
}

SYM_BASIS = [sp.Integer(1), PI**2, Z3, PI**4, Z5, PI**2 * Z3]
BASIS_WEIGHTS = (0, 2, 3, 4, 5, 5)
ZETA_SUBS = {sp.zeta(2): PI**2 / 6, sp.zeta(4): PI**4 / 90}


def to_mp(expr):
    """Evaluate a real sympy constant with mpmath at current precision."""
    return sp.lambdify((), expr, "mpmath")()


# ---------------------------------------------------------------------------
# J1 / J2: inner x-integral reduced exactly to a 2F1 kernel
# ---------------------------------------------------------------------------


def j_inner(y, ep, which):
    """int_0^oo (1+x)^a (1+x+y)^c dx  ==  int_0^1 w^(-a-c-2) (1+y*w)^c dw.

    J1: (a, c) = (ep, ep-2)  ->  2F1(2-ep, 1-2ep; 2-2ep; -y) / (1-2ep)
    J2: (a, c) = (ep, ep-3)  ->  2F1(3-ep, 2-2ep; 3-2ep; -y) / (2-2ep)
    """
    if which == 1:
        return mp.hyp2f1(2 - ep, 1 - 2 * ep, 2 - 2 * ep, -y) / (1 - 2 * ep)
    return mp.hyp2f1(3 - ep, 2 - 2 * ep, 3 - 2 * ep, -y) / (2 - 2 * ep)


def j_val(ep, which):
    """J1(ep) or J2(ep) as a 1-D quadrature over the outer variable."""
    bexp = ep - 1 if which == 1 else ep

    def f(y):
        return (1 + y) ** bexp * j_inner(y, ep, which)

    return mp.quad(f, [0, 1, 8, 64, 1024, mp.inf], error=True, maxdegree=7)


def j_val_2d(ep, which):
    """Direct 2-D quadrature of the defining integral (validation only)."""

    def f(x, y):
        if which == 1:
            return (1 + x) ** ep * (1 + y) ** (ep - 1) * (1 + x + y) ** (ep - 2)
        return (1 + x) ** ep * (1 + y) ** ep * (1 + x + y) ** (ep - 3)

    pts = [0, 1, 10, 100, mp.inf]
    with mp.workdps(20):
        return mp.quad(f, pts, pts)


# ---------------------------------------------------------------------------
# Taylor coefficients via the Cauchy formula + PSLQ identification
# ---------------------------------------------------------------------------

NPTS = 44


def taylor_coeffs(which, kmax, npts=NPTS):
    radius = mp.mpf(1) / 32
    half = npts // 2
    vals = [None] * npts
    qerr = mp.mpf(0)
    for j in range(half + 1):
        epj = radius * mp.expjpi(mp.mpf(2 * j) / npts)
        v, e = j_val(epj, which)
        vals[j] = v
        qerr = max(qerr, abs(e))
    for j in range(half + 1, npts):
        vals[j] = mp.conj(vals[npts - j])
    coeffs = []
    imag_max = mp.mpf(0)
    for k in range(kmax + 1):
        s = mp.mpc(0)
        for j in range(npts):
            s += vals[j] * mp.expjpi(mp.mpf(-2 * j * k) / npts)
        c = s / (npts * radius**k)
        imag_max = max(imag_max, abs(mp.im(c)))
        coeffs.append(mp.re(c))
    return coeffs, qerr, imag_max


def identify(cnum, order):
    """PSLQ-identify cnum in the weight-graded basis (max weight order+1).

    The order-k Taylor coefficient of J1/J2 is (1/k!) * int rational * log^k,
    and the two remaining integrations can raise the weight by one more unit
    (e.g. J1'(0) = 3 + zeta(2), weight 2 at order 1), so the cutoff is k+1.
    """
    idx = [i for i, w in enumerate(BASIS_WEIGHTS) if w <= order + 1]
    nb = [to_mp(e) for e in SYM_BASIS]
    vec = [cnum] + [nb[i] for i in idx]
    rel = mp.pslq(vec, tol=mp.mpf(10) ** -34, maxcoeff=10**8, maxsteps=10**6)
    if not rel or rel[0] == 0:
        raise RuntimeError(f"PSLQ failed at order {order}: {mp.nstr(cnum, 30)}")
    expr = sp.Integer(0)
    num = mp.mpf(0)
    for m, i in enumerate(idx):
        expr += sp.Integer(rel[m + 1]) * SYM_BASIS[i]
        num += rel[m + 1] * nb[i]
    expr = sp.simplify(-sp.Rational(1, rel[0]) * expr)
    resid = abs(cnum + num / rel[0])
    return expr, resid


# ---------------------------------------------------------------------------
# External prefactor (exact reduction and Laurent series)
# ---------------------------------------------------------------------------


def prefactor_direct(ep):
    return (
        mp.exp(2 * ep * mp.euler)
        * mp.gamma(1 - ep)
        * mp.gamma(-ep) ** 2
        * mp.gamma(ep)
        * mp.gamma(2 * ep)
        / (mp.gamma(-1 - 3 * ep) * mp.gamma(1 + ep))
    )


def prefactor_reduced_num(ep):
    return (
        mp.mpf(3)
        / 2
        * (1 + 3 * ep)
        / ep**3
        * mp.exp(2 * ep * mp.euler)
        * mp.gamma(1 - ep) ** 3
        * mp.gamma(1 + 2 * ep)
        / mp.gamma(1 - 3 * ep)
    )


def _log_gamma1p(z, kmax=4):
    """ln Gamma(1+z) as a truncated series with symbolic zeta values."""
    s = -sp.EulerGamma * z
    for k in range(2, kmax + 1):
        s += sp.Rational((-1) ** k, k) * sp.zeta(k) * z**k
    return s


def prefactor_series_coeffs():
    ln_reg = (
        2 * EP * sp.EulerGamma
        + 3 * _log_gamma1p(-EP)
        + _log_gamma1p(2 * EP)
        - _log_gamma1p(-3 * EP)
    )
    reg = sp.series(sp.exp(sp.expand(ln_reg)), EP, 0, 5).removeO()
    pre = sp.expand(sp.Rational(3, 2) * (1 + 3 * EP) * reg / EP**3)
    out = {}
    for k in range(-3, 2):
        out[k] = sp.simplify(sp.expand(pre.coeff(EP, k)).subs(ZETA_SUBS))
    return out


def full_numeric(ep):
    j1, _ = j_val(ep, 1)
    j2, _ = j_val(ep, 2)
    a = (4 * ep - 1) / (3 * (3 * ep + 1))
    b = (ep - 2) * (5 * ep - 2) / (3 * ep * (3 * ep + 1))
    return prefactor_direct(ep) * (a * j1 + b * j2)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    print("External Int1 Laurent audit (standalone; reducer code untouched)")
    print(f"mp.dps = {mp.mp.dps}")
    report = {"dps": mp.mp.dps}

    # 0. Validate the 2F1 reduction against direct 2-D quadrature.
    print("[0] validating 2F1 reduction of the inner integral ...")
    for which in (1, 2):
        v1, _ = j_val(mp.mpf("0.03"), which)
        v2 = j_val_2d(mp.mpf("0.03"), which)
        d = abs(v1 - v2)
        print(f"    J{which}(0.03): 1D-vs-2D |diff| = {mp.nstr(d, 3)}")
        if d > mp.mpf(10) ** -10:
            raise RuntimeError(f"2F1 reduction check failed for J{which}")

    # 0b. Validate the exact prefactor reduction numerically.
    for eptest in (mp.mpf("0.0173"), mp.mpc("0.011", "0.007")):
        pd = prefactor_direct(eptest)
        rel = abs(pd - prefactor_reduced_num(eptest)) / abs(pd)
        print(f"    prefactor reduced-form rel.diff at ep={eptest}: {mp.nstr(rel, 3)}")
        if rel > mp.mpf(10) ** -30:
            raise RuntimeError("prefactor reduction check failed")

    # 1. Prefactor Laurent series (exact) + cross-check vs hand derivation.
    pcoef = prefactor_series_coeffs()
    print("[1] prefactor Laurent coefficients (exact):")
    for k in range(-3, 2):
        if sp.simplify(pcoef[k] - EXPECTED_PREFACTOR[k]) != 0:
            raise RuntimeError(f"prefactor coefficient mismatch at ep^{k}")
        print(f"    ep^{k}: {pcoef[k]}")
    report["prefactor"] = {str(k): str(v) for k, v in pcoef.items()}

    # 2. Taylor coefficients of J1, J2 + PSLQ identification.
    print(f"[2] Cauchy-circle Taylor extraction (|ep| = 1/32, {NPTS} nodes) ...")
    j_exact = {1: {}, 2: {}}
    resids = {}
    for which, kmax in ((1, 3), (2, 4)):
        coeffs, qerr, imax = taylor_coeffs(which, kmax)
        print(f"    J{which}: quad err <= {mp.nstr(qerr, 3)}, max |Im c_k| = {mp.nstr(imax, 3)}")
        if imax > mp.mpf(10) ** -30:
            raise RuntimeError(f"imaginary residue too large for J{which}")
        for k, c in enumerate(coeffs):
            expr, resid = identify(c, k)
            j_exact[which][k] = expr
            resids[f"J{which}[{k}]"] = mp.nstr(resid, 3)
            print(f"    J{which}[ep^{k}] = {expr}   (PSLQ residual {mp.nstr(resid, 3)})")
            if resid > mp.mpf(10) ** -32:
                raise RuntimeError(f"PSLQ residual too large: J{which}[{k}]")
            if (which, k) in KNOWN and sp.simplify(expr - KNOWN[(which, k)]) != 0:
                raise RuntimeError(f"identified J{which}[{k}] contradicts analytic anchor")
    report["pslq_residuals"] = resids
    report["J1"] = {str(k): str(v) for k, v in j_exact[1].items()}
    report["J2"] = {str(k): str(v) for k, v in j_exact[2].items()}

    # 3. Exact assembly: C = A*J1 + B*J2, then Full = Prefactor * C.
    j1s = sum(j_exact[1][k] * EP**k for k in range(4))
    j2s = sum(j_exact[2][k] * EP**k for k in range(5))
    cfull = sp.expand(sp.series(A_EXPR * j1s + B_EXPR * j2s, EP, 0, 4).removeO())
    ccoef = {m: sp.simplify(cfull.coeff(EP, m)) for m in range(-1, 4)}
    print("[3] C(ep) = A*J1 + B*J2 series (trustworthy through ep^3):")
    for m in range(-1, 4):
        print(f"    ep^{m}: {ccoef[m]}")
    if sp.simplify(ccoef[-1] - sp.Rational(2, 3)) != 0:
        raise RuntimeError("C[1/ep] != 2/3")
    if sp.simplify(ccoef[0] + 2) != 0:
        raise RuntimeError("C[ep^0] != -2")
    report["C"] = {str(m): str(v) for m, v in ccoef.items()}

    full = {}
    for n in range(-4, 1):
        tot = sp.Integer(0)
        for k in range(-3, 2):
            m = n - k
            if -1 <= m <= 3:
                tot += pcoef[k] * ccoef[m]
        full[n] = sp.simplify(sp.expand(tot))
    report["full"] = {str(n): str(v) for n, v in full.items()}
    report["target"] = {str(n): str(v) for n, v in TARGET.items()}

    # 4. Comparison against the reference target, order by order.
    print("[4] full Laurent expansion (units of 1/(s*t^2)) vs target:")
    verdicts = {}
    for n in range(-4, 1):
        diff = sp.simplify(full[n] - TARGET[n])
        ok = diff == 0
        verdicts[str(n)] = "OK" if ok else f"MISMATCH diff={diff}"
        print(f"    ep^{n}: ours = {full[n]}")
        print(
            f"           target = {TARGET[n]}  ->  "
            f"{'OK' if ok else 'MISMATCH, diff = ' + str(diff)}"
        )
    report["verdicts"] = verdicts

    # 5. Numeric sanity of the assembled Laurent polynomial.
    ep_num = mp.mpf("0.02")
    truth = full_numeric(ep_num)
    laur = sum(to_mp(full[n]) * ep_num**n for n in range(-4, 1))
    print(f"[5] numeric sanity at ep = {ep_num}:")
    print(f"    direct   = {mp.nstr(truth, 25)}")
    print(f"    laurent  = {mp.nstr(laur, 25)}")
    print(f"    diff     = {mp.nstr(truth - laur, 5)} (expected O(ep^1) tail)")
    report["numeric_sanity"] = {
        "ep": mp.nstr(ep_num, 10),
        "direct": mp.nstr(truth, 30),
        "laurent_sum": mp.nstr(laur, 30),
        "diff": mp.nstr(truth - laur, 10),
    }

    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"wrote {OUT_JSON}")


if __name__ == "__main__":
    main()
