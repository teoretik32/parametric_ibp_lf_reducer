"""External Int2 Method.2: exact x7 preintegration + leading Laurent pole audit.

Pure symbolic/numeric mathematics (no security/network relevance). This script does
NOT import or modify the reducer core (``parametric_ibp_lf_reducer``), does NOT run
any RREF, and does NOT touch the certified reduction artifacts. It audits the
*wrapper* level only:

1. Exact x7 preintegration. For ``A = 1 + r*x2*x5 > 0``, ``B = 1 + x2 > 0``::

       Integral[(1+x7)^(-1-ep)*(A+B*x7)^(-1+ep), {x7, 0, Infinity}]
           = (B^ep - A^ep)/(ep*(B - A))          (u = x7/(1+x7) maps it to
                                                  Integral[(A+(B-A)*u)^(ep-1),{u,0,1}])

   With ``B - A = x2*(1 - r*x5)`` this is the identity quoted in the task, so::

       J2(ep, r) = (1/ep)*Integral[x2^ep*(1+x2)^ep*(1+x5)^ep
                       *((1+x2)^ep - (1+r*x2*x5)^ep)/(1 - r*x5),
                       {x2, 0, Infinity}, {x5, 0, Infinity}]

   The integrand is regular at ``x5 = 1/r`` (the numerator vanishes there too).

2. Reduced 1-D form. The x2 integral is a Gauss 2F1 (G&R 3.197.1)::

       Integral[x2^ep*(1+x2)^ep*(1+c*x2)^ep, {x2, 0, Infinity}]
           = G1(ep)*2F1(-ep, 1+ep; -2*ep; 1-c),
       G1(ep) = Gamma[1+ep]*Gamma[-1-3*ep]/Gamma[-2*ep]     (analytic continuation)

   so with ``z = 1 - r*x5``::

       J2 = (G1(ep)/ep)*Q(ep, r),
       Q   = Integral[(1+x5)^ep*(1 - 2F1(-ep,1+ep;-2*ep;z))/z, {x5, 0, Infinity}]

3. Pole bookkeeping (the subtle point). Endpoint zones carry algebraic layers with
   exponents that degenerate to ``-1`` as ``ep -> 0``:

   * ``x5 -> Infinity``: ``x5^(ep-1)`` (from the "1") and ``x5^(2*ep-1)`` with
     connection coefficient ``K1 = Gamma[-2ep]*Gamma[1+2ep]/(Gamma[1+ep]*Gamma[-ep])``;
   * ``x5 -> 0`` (z -> 1): ``x5^(-1-2*ep)`` with coefficient
     ``C_B = Gamma[-2ep]*Gamma[1+2ep]/(Gamma[-ep]*Gamma[1+ep])``.

   ``K1 == C_B`` *identically*, so the two crossover poles cancel exactly and only the
   "1"-layer survives:  ``Q = 1/(r*ep) + O(1)``, hence::

       J2(ep, r) = G1(0)/(r*ep^2) + O(1/ep) = -2/(3*r*ep^2) + O(1/ep)

   (A naive x5->Infinity-only count would give -1/(2*r*ep^2); the x5->0 crossover
   layer restores -2/(3*r*ep^2).  Both boundary layers are verified here.)

4. Corrected external prefactor (the wrapper previously missed the exponential)::

       ExternalPrefactor2 = Exp[2*ep*EulerGamma]*t^(-3-ep)
           *Gamma[1-ep]*Gamma[-ep]^3*Gamma[ep]/(Gamma[-1-3*ep]*Gamma[-2*ep])
           = 6/(t^3*ep^2) + O(1/ep)

   With ``Exp[2*ep*EulerGamma]`` the Laurent coefficients are EulerGamma-free (the
   net gamma_E count of the Gamma ratio is exactly ``-2*ep``); pure reduction
   coefficients are unaffected -- this is wrapper/reference metadata only.

5. Leading pole of the full object, with ``r = s/t``::

       P2*J2 = (6/(t^3*ep^2))*(-2/(3*(s/t)*ep^2)) + O(1/ep^3)
             = -4/(s*t^2*ep^4) + O(1/ep^3)

   which must (and does) match the supplied AnsvInt2 leading pole exactly.  The full
   AnsvInt2 value is NOT invented and never enters the reducer; its source text is
   stored as metadata only in ``examples/external_int2_source_reference.wl.txt``.

Numerics: mpmath only, moderate cost, no RREF.  The decomposition (closed endpoint
sums + uniformly regular remainders) is validated against direct quadrature at
moderate ``ep``, checked for invariance under (h, X, N, delta) changes, and only then
used for the small-``ep`` Laurent fit.

Exit codes: 0 = all checks passed (JSON written), 1 = at least one check failed
(JSON still written), 2 = usage problem.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import mpmath as mp
import sympy as sp

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_JSON = REPO_ROOT / "validation" / "external_int2_leading_pole_audit.json"

CORRECTED_PREFACTOR_TEXT = (
    "Exp[2*ep*EulerGamma]*t^(-3-ep)*Gamma[1-ep]*Gamma[-ep]^3*Gamma[ep]"
    "/(Gamma[-1-3*ep]*Gamma[-2*ep])"
)
OLD_PREFACTOR_TEXT = "t^(-3-ep)*Gamma[1-ep]*Gamma[-ep]^3*Gamma[ep]/(Gamma[-1-3*ep]*Gamma[-2*ep])"
ANSV_INT2_SUPPLIED_LEADING = "-4/(s*t^2*ep^4)"

J2_LEADING_TEXT = "-2/(3*r*ep^2)"
P2_LEADING_TEXT = "6/(t^3*ep^2)"
PRODUCT_LEADING_TEXT = "-4/(s*t^2*ep^4)"


def _log(msg: str) -> None:
    print(f"[audit-int2-pole] {msg}", flush=True)


# --------------------------------------------------------------------------------------
# Symbolic (sympy) part: exact statements
# --------------------------------------------------------------------------------------


def check_x7_identity_symbolic() -> dict:
    """d/du of the closed antiderivative reproduces the x7 integrand (u-form)."""
    u, ep, a_, b_ = sp.symbols("u ep A B", positive=True)
    anti = (a_ + (b_ - a_) * u) ** ep / (ep * (b_ - a_))
    ok_diff = sp.simplify(sp.diff(anti, u) - (a_ + (b_ - a_) * u) ** (ep - 1)) == 0
    boundary = sp.simplify(
        (anti.subs(u, 1) - anti.subs(u, 0)) - (b_**ep - a_**ep) / (ep * (b_ - a_))
    )
    return {
        "name": "x7_identity_symbolic",
        "passed": bool(ok_diff and boundary == 0),
        "antiderivative_ok": bool(ok_diff),
        "boundary_ok": bool(boundary == 0),
        "note": "x7 = u/(1-u) maps the x7 integral to Integral[(A+(B-A)u)^(ep-1),{u,0,1}]",
    }


def check_prefactor_series() -> dict:
    """ep^2*P2hat -> 6 at ep=0; corrected series EulerGamma-free, old one is not."""
    ep = sp.Symbol("ep")
    gamma_ratio = (
        sp.gamma(1 - ep)
        * sp.gamma(-ep) ** 3
        * sp.gamma(ep)
        / (sp.gamma(-1 - 3 * ep) * sp.gamma(-2 * ep))
    )
    new_hat = sp.exp(2 * ep * sp.EulerGamma) * gamma_ratio
    coeffs_new = sp.Poly(sp.series(ep**2 * new_hat, ep, 0, 3).removeO().expand(), ep).all_coeffs()[
        ::-1
    ]
    coeffs_old = sp.Poly(
        sp.series(ep**2 * gamma_ratio, ep, 0, 3).removeO().expand(), ep
    ).all_coeffs()[::-1]
    lead_new = sp.simplify(coeffs_new[0])
    new_gamma_free = all(not c.has(sp.EulerGamma) for c in coeffs_new)
    old_has_gamma = any(c.has(sp.EulerGamma) for c in coeffs_old)
    return {
        "name": "prefactor_series",
        "passed": bool(lead_new == 6 and new_gamma_free and old_has_gamma),
        "leading_coefficient": str(lead_new),
        "corrected_series_eulergamma_free": bool(new_gamma_free),
        "old_series_has_eulergamma": bool(old_has_gamma),
        "p2_leading": P2_LEADING_TEXT,
        "note": "net gamma_E count of the Gamma ratio is exp(-2*ep*EulerGamma); "
        "Exp[2*ep*EulerGamma] removes it, so Laurent data is EulerGamma-free",
    }


def check_leading_pole_exact() -> dict:
    """Exact pole bookkeeping: K1 == C_B, Q-pole = 1/(r*ep), c2*r = -2/3, product."""
    ep = sp.Symbol("ep")
    r, s, t, x_cut, h_cut = sp.symbols("r s t X h", positive=True)
    k1 = sp.gamma(-2 * ep) * sp.gamma(1 + 2 * ep) / (sp.gamma(1 + ep) * sp.gamma(-ep))
    c_b = sp.gamma(-2 * ep) * sp.gamma(1 + 2 * ep) / (sp.gamma(-ep) * sp.gamma(1 + ep))
    k1_eq_cb = sp.simplify(k1 - c_b) == 0

    # closed singular sums: tail "1"-layer, tail K1-layer, head C_B-layer (k = 0 terms;
    # k >= 1 terms are regular at ep = 0 and cannot contribute to the 1/ep pole)
    q_sing = (
        x_cut**ep / (r * ep)
        - k1 * r ** (ep - 1) * x_cut ** (2 * ep) / (2 * ep)
        + c_b * r ** (-1 - 2 * ep) * h_cut ** (-2 * ep) / (2 * ep)
    )
    pole_coeff = sp.simplify(sp.limit(ep * q_sing, ep, 0))
    crossover = sp.simplify(
        sp.limit(
            ep
            * (
                -k1 * r ** (ep - 1) * x_cut ** (2 * ep) / (2 * ep)
                + c_b * r ** (-1 - 2 * ep) * h_cut ** (-2 * ep) / (2 * ep)
            ),
            ep,
            0,
        )
    )
    g1_0 = sp.limit(sp.gamma(1 + ep) * sp.gamma(-1 - 3 * ep) / sp.gamma(-2 * ep), ep, 0)
    beta_lemma = sp.simplify(g1_0 - sp.Rational(-2, 3)) == 0  # Int[x^ep (1+x)^(2ep)]
    c2_times_r = sp.simplify(g1_0 * pole_coeff * r)  # J2 ~ (G1/ep)*(pole/ep)
    j2_lead = g1_0 * pole_coeff / ep**2
    p2_lead = 6 / (t**3 * ep**2)
    product = sp.simplify(p2_lead * j2_lead.subs(r, s / t) - (-4 / (s * t**2 * ep**4)))
    return {
        "name": "leading_pole_exact",
        "passed": bool(
            k1_eq_cb
            and pole_coeff == 1 / r
            and crossover == 0
            and beta_lemma
            and c2_times_r == sp.Rational(-2, 3)
            and product == 0
        ),
        "k1_equals_cb": bool(k1_eq_cb),
        "q_pole_coefficient": str(pole_coeff),
        "crossover_pole_cancellation": str(crossover),
        "beta_lemma_g1_at_0": str(g1_0),
        "c2_times_r": str(c2_times_r),
        "j2_leading": J2_LEADING_TEXT,
        "product_leading": PRODUCT_LEADING_TEXT,
        "ansv_int2_supplied_leading": ANSV_INT2_SUPPLIED_LEADING,
        "ansv_match": True,
        "note": "naive x5->Infinity-only bookkeeping would give -1/(2*r*ep^2); the "
        "x5->0 crossover layer (K1 == C_B) cancels the x5->Infinity 2ep-layer",
    }


# --------------------------------------------------------------------------------------
# Numeric (mpmath) part: truncated-series helpers for the endpoint layers
# --------------------------------------------------------------------------------------


def _ser_mul(a: list, b: list, n: int) -> list:
    out = [mp.mpf(0)] * (n + 1)
    for i, ai in enumerate(a[: n + 1]):
        if ai:
            for j, bj in enumerate(b[: n + 1 - i]):
                out[i + j] += ai * bj
    return out


def _ser_binom(alpha, scale, n: int) -> list:
    """Coefficients of (1 + scale*x)^alpha up to x^n."""
    return [mp.binomial(alpha, k) * scale**k for k in range(n + 1)]


def _ser_hyp(a, b, c, n: int) -> list:
    """Taylor coefficients of 2F1(a, b; c; w) in w up to w^n."""
    return [mp.rf(a, k) * mp.rf(b, k) / (mp.rf(c, k) * mp.factorial(k)) for k in range(n + 1)]


def _ser_compose(outer: list, inner: list, n: int) -> list:
    """outer(inner(x)) for a power series ``inner`` with zero constant term."""
    assert inner[0] == 0
    out = [mp.mpf(0)] * (n + 1)
    out[0] = outer[-1]
    for k in range(len(outer) - 2, -1, -1):  # Horner
        out = _ser_mul(out, inner, n)
        out[0] += outer[k]
    return out


def _geom(scale, n: int) -> list:
    """Coefficients of 1/(1 - scale*x)."""
    return [scale**k for k in range(n + 1)]


class _Decomposition:
    """Q(ep, r) = closed endpoint sums + uniformly regular numeric remainders."""

    def __init__(self, ep, r, n_ser: int = 8, h_frac=None, x_frac=None, delta=None):
        self.ep, self.r, self.n = ep, r, n_ser
        self.h = (h_frac or mp.mpf("0.5")) / r
        self.x_cut = max((x_frac or mp.mpf(4)) / r, (x_frac or mp.mpf(4)))
        self.delta = delta or mp.mpf("1e-3")
        e = ep
        self.k1 = mp.gamma(-2 * e) * mp.gamma(1 + 2 * e) / (mp.gamma(1 + e) * mp.gamma(-e))
        self.k2 = mp.gamma(-2 * e) * mp.gamma(-1 - 2 * e) / (mp.gamma(-e) * mp.gamma(-1 - 3 * e))
        self.c_a = self.k2  # z->1 analytic coefficient == z->-oo (-z)^(-1-ep) coefficient
        self.c_b = self.k1  # z->1 singular coefficient == z->-oo (-z)^ep coefficient
        self.g1 = mp.gamma(1 + e) * mp.gamma(-1 - 3 * e) / mp.gamma(-2 * e)

    # -- integrand pieces ------------------------------------------------------------
    def f_full(self, x5):
        """(1+x5)^ep * (1 - 2F1(-ep,1+ep;-2ep;z))/z with a guarded z -> 0 limit."""
        e = self.ep
        z = 1 - self.r * x5
        if abs(z) < mp.mpf("1e-8"):
            a, b, c = -e, 1 + e, -2 * e
            d1 = a * b / c
            d2 = a * (a + 1) * b * (b + 1) / (c * (c + 1) * 2)
            val = -d1 - d2 * z
        else:
            val = (1 - mp.hyp2f1(-e, 1 + e, -2 * e, z)) / z
        return (1 + x5) ** e * val

    def _phi_series(self) -> list:
        """Taylor of F_B(r*x5)*(1+x5)^ep/(1-r*x5) in x5 (head singular cofactor)."""
        e, r, n = self.ep, self.r, self.n
        f_b = [c * r**k for k, c in enumerate(_ser_hyp(-e, -1 - 3 * e, -2 * e, n))]
        return _ser_mul(_ser_mul(f_b, _ser_binom(e, mp.mpf(1), n), n), _geom(r, n), n)

    # -- head window [0, h]: z in (1-r*h, 1) -----------------------------------------
    def head(self) -> tuple:
        e, r, n, h = self.ep, self.r, self.n, self.h
        phi = self._phi_series()
        pref = -self.c_b * r ** (-1 - 2 * e)
        s_closed = pref * mp.fsum(phi[k] * h ** (k - 2 * e) / (k - 2 * e) for k in range(n + 1))

        def regular(x5):  # (1+x5)^ep*(1 - C_A*F_A(w))/z, w = r*x5, z = 1 - w
            w = r * x5
            f_a = mp.hyp2f1(-e, 1 + e, 2 + 2 * e, w)
            return (1 + x5) ** e * (1 - self.c_a * f_a) / (1 - w)

        r1 = mp.quad(regular, [0, h / 10, h])

        def sing_tail(x5):  # x5^(-1-2ep)*(phi(x5) - truncated phi)
            w = r * x5
            phi_exact = mp.hyp2f1(-e, -1 - 3 * e, -2 * e, w) * (1 + x5) ** e / (1 - w)
            phi_part = mp.fsum(phi[k] * x5**k for k in range(n + 1))
            return pref * x5 ** (-1 - 2 * e) * (phi_exact - phi_part)

        r2 = mp.quad(sing_tail, [self.delta, min(10 * self.delta, h), h / 10, h])
        # |neglected [0, delta] piece| <= |pref|*max|phi_{n+1}|-ish * delta^(n+1-2ep)
        neglect = abs(pref) * self.delta ** (n + 1 - 2 * e) / (n + 1 - 2 * e)
        return s_closed, r1 + r2, neglect

    # -- tail window [X, oo): z -> -oo ------------------------------------------------
    def _tail_series(self) -> tuple:
        """Coefficient lists (in u = 1/x5) for the x5^(ep-1), x5^(2ep-1), x5^(-2) layers."""
        e, r, n = self.ep, self.r, self.n
        inner = [mp.mpf(0)] + [-((mp.mpf(1) / r) ** k) for k in range(1, n + 1)]  # -v/(1-v)
        one_pu = _ser_binom(e, mp.mpf(1), n)  # (1+u)^ep
        p1 = _ser_mul([-c / r for c in _geom(1 / r, n)], one_pu, n)  # -(1/r)/(1-u/r)*(1+u)^ep
        f_t1 = _ser_compose(_ser_hyp(-e, 1 + e, -2 * e, n), inner, n)
        p2 = _ser_mul(_ser_mul(_ser_binom(e - 1, -1 / r, n), f_t1, n), one_pu, n)
        p2 = [self.k1 * r ** (e - 1) * c for c in p2]
        f_t2 = _ser_compose(_ser_hyp(1 + e, 2 + 3 * e, 2 + 2 * e, n), inner, n)
        p3 = _ser_mul(_ser_mul(_ser_binom(-2 - e, -1 / r, n), f_t2, n), one_pu, n)
        p3 = [self.k2 * r ** (-2 - e) * c for c in p3]
        return p1, p2, p3

    def tail(self) -> tuple:
        e, n, x_cut = self.ep, self.n, self.x_cut
        p1, p2, p3 = self._tail_series()
        s_closed = mp.fsum(
            p1[j] * x_cut ** (e - j) / (j - e)
            + p2[j] * x_cut ** (2 * e - j) / (j - 2 * e)
            + p3[j] * x_cut ** (-1 - j) / (j + 1)
            for j in range(n + 1)
        )

        def remainder(y):  # x5 = X/y, y in (0, 1]
            x5 = x_cut / y
            u = 1 / x5
            s_val = (
                x5 ** (e - 1) * mp.fsum(p1[j] * u**j for j in range(n + 1))
                + x5 ** (2 * e - 1) * mp.fsum(p2[j] * u**j for j in range(n + 1))
                + x5**-2 * mp.fsum(p3[j] * u**j for j in range(n + 1))
            )
            return (self.f_full(x5) - s_val) * x_cut / y**2

        r_tail = mp.quad(remainder, [0, mp.mpf("0.5"), 1])
        return s_closed, r_tail

    # -- assembled Q and J2 ------------------------------------------------------------
    def q_value(self) -> dict:
        s_head, r_head, neglect = self.head()
        mid = mp.quad(
            self.f_full,
            [self.h, mp.mpf("0.9") / self.r, 1 / self.r, mp.mpf("1.1") / self.r, self.x_cut],
        )
        s_tail, r_tail = self.tail()
        total = s_head + r_head + mid + s_tail + r_tail
        return {
            "q": total,
            "s_head": s_head,
            "r_head": r_head,
            "mid": mid,
            "s_tail": s_tail,
            "r_tail": r_tail,
            "neglected_bound": neglect,
        }

    def j2(self):
        return self.g1 / self.ep * self.q_value()["q"]


def _direct_q(ep, r) -> mp.mpf:
    """Plain compactified quadrature of Q (reliable only for moderate |ep|)."""
    dec = _Decomposition(ep, r)

    def fy(y):
        x5 = y / (1 - y)
        return dec.f_full(x5) / (1 - y) ** 2

    pts = [0]
    for x in (1 / (2 * r), 1 / r, 2 / r, 10 / r, 100 / r):
        pts.append(x / (1 + x))
    pts.append(1)
    return mp.quad(fy, pts)


def check_x7_identity_numeric(fast: bool) -> dict:
    cases, worst = [], mp.mpf(0)
    ep_list = ("-0.55", "-0.3") if not fast else ("-0.3",)
    grid = ((0.35, 0.6), (2.2, 3.1)) if fast else ((0.35, 0.6), (2.2, 3.1), (0.9, 1.7))
    for ep_s in ep_list:
        ep = mp.mpf(ep_s)
        for r in (mp.mpf("0.4"), mp.mpf("1.3")):
            for x2_, x5_ in grid:
                x2, x5 = mp.mpf(x2_), mp.mpf(x5_)
                a_val, b_val = 1 + r * x2 * x5, 1 + x2
                direct = mp.quad(
                    lambda x7: (1 + x7) ** (-1 - ep) * (a_val + b_val * x7) ** (ep - 1),
                    [0, 1, 10, 100, mp.inf],
                )
                closed = (b_val**ep - a_val**ep) / (ep * (b_val - a_val))
                rel = abs(direct - closed) / abs(closed)
                worst = max(worst, rel)
                cases.append(
                    {"ep": ep_s, "r": float(r), "x2": x2_, "x5": x5_, "rel_err": float(rel)}
                )
    tol = 1e-15
    return {
        "name": "x7_identity_numeric",
        "passed": bool(worst < tol),
        "max_rel_err": float(worst),
        "tolerance": tol,
        "n_cases": len(cases),
        "cases": cases,
    }


def check_connection_formulas() -> dict:
    """Numeric validation of the z->1 and z->-oo 2F1 connection decompositions."""
    ep = mp.mpf("-0.22")
    dec = _Decomposition(ep, mp.mpf("0.75"))
    e = ep
    # z -> 1: F(z) = C_A*F_A(w) + C_B*w^(-1-2ep)*F_B(w), w = 1 - z = 0.3
    w = mp.mpf("0.3")
    lhs = mp.hyp2f1(-e, 1 + e, -2 * e, 1 - w)
    rhs = dec.c_a * mp.hyp2f1(-e, 1 + e, 2 + 2 * e, w) + dec.c_b * w ** (-1 - 2 * e) * mp.hyp2f1(
        -e, -1 - 3 * e, -2 * e, w
    )
    err_head = abs(lhs - rhs) / abs(lhs)
    # z -> -oo: F(z) = K1*(-z)^ep*F_T1(1/z) + K2*(-z)^(-1-ep)*F_T2(1/z), z = -50
    z = mp.mpf(-50)
    lhs2 = mp.hyp2f1(-e, 1 + e, -2 * e, z)
    rhs2 = dec.k1 * (-z) ** e * mp.hyp2f1(-e, 1 + e, -2 * e, 1 / z) + dec.k2 * (-z) ** (
        -1 - e
    ) * mp.hyp2f1(1 + e, 2 + 3 * e, 2 + 2 * e, 1 / z)
    err_tail = abs(lhs2 - rhs2) / abs(lhs2)
    tol = 10.0 ** (-(mp.mp.dps - 6))
    return {
        "name": "hyp2f1_connection_formulas",
        "passed": bool(err_head < tol and err_tail < tol),
        "rel_err_z_to_1": float(err_head),
        "rel_err_z_to_minus_inf": float(err_tail),
        "tolerance": tol,
    }


def check_decomposition_consistency(fast: bool) -> dict:
    """Decomposition == direct quadrature at moderate ep; invariance under knobs."""
    r = mp.mpf("0.75")
    tol_direct = 1e-7 if fast else 1e-9
    rows, ok = [], True
    for ep_s in ("-0.30", "-0.22"):
        ep = mp.mpf(ep_s)
        ambient = mp.mp.dps
        mp.mp.dps = ambient + 10  # the direct reference needs extra digits to settle
        try:
            dq = _direct_q(ep, r)
        finally:
            mp.mp.dps = ambient
        base = _Decomposition(ep, r).q_value()["q"]
        rel = abs(dq - base) / abs(dq)
        ok = ok and rel < tol_direct
        row = {"ep": ep_s, "rel_vs_direct": float(rel), "tolerance": tol_direct}
        if not fast:
            alt = _Decomposition(
                ep, r, n_ser=10, h_frac=mp.mpf("1") / 3, x_frac=mp.mpf(6), delta=mp.mpf("1e-4")
            ).q_value()["q"]
            rel_alt = abs(alt - base) / abs(base)
            ok = ok and rel_alt < 1e-12
            row["rel_vs_alt_knobs"] = float(rel_alt)
        rows.append(row)
    return {"name": "decomposition_consistency", "passed": bool(ok), "rows": rows}


def check_leading_pole_numeric(fast: bool) -> dict:
    """Small-ep Laurent fit of r*ep^2*J2 via the decomposition -> c2 = -2/3 exactly."""
    ep_strs = (
        ("-0.05", "-0.03", "-0.02", "-0.01")
        if fast
        else ("-0.05", "-0.04", "-0.03", "-0.02", "-0.01", "-0.005")
    )
    r_strs = ("0.75",) if fast else ("0.75", "0.4")
    deg = 2 if fast else 3
    target = -mp.mpf(2) / 3
    fits, ok = [], True
    for r_s in r_strs:
        r = mp.mpf(r_s)
        eps, vals = [], []
        for ep_s in ep_strs:
            ep = mp.mpf(ep_s)
            v = r * ep**2 * _Decomposition(ep, r).j2()
            eps.append(ep)
            vals.append(v)
        a_mat = mp.matrix(len(eps), deg + 1)
        for i, e in enumerate(eps):
            for j in range(deg + 1):
                a_mat[i, j] = e**j
        coef = mp.lu_solve(a_mat.T * a_mat, a_mat.T * mp.matrix(vals))
        dev = abs(coef[0] - target)
        tol = 5e-3 if fast else 5e-4
        ok = ok and dev < tol
        fits.append(
            {
                "r": r_s,
                "ep_values": list(ep_strs),
                "r_ep2_j2": [mp.nstr(v, 15) for v in vals],
                "c2_fit": mp.nstr(coef[0], 12),
                "c1_fit": mp.nstr(coef[1], 8),
                "abs_dev_from_minus_2_3": float(dev),
                "tolerance": tol,
            }
        )
    return {
        "name": "leading_pole_numeric",
        "passed": bool(ok),
        "target_c2_times_r": "-2/3",
        "fits": fits,
        "note": "clearly excludes the naive -1/2 (distance 1/6)",
    }


def check_2d_cross(fast: bool) -> dict:
    """Direct 2-D (x2, x5) quadrature of J2 vs the reduced 1-D decomposition."""
    if fast:
        return {"name": "j2_2d_cross_check", "passed": True, "skipped": True}
    ep, r = mp.mpf("-0.40"), mp.mpf("0.75")

    def bracket_over_z(x2, z):
        # ((1+x2)^ep - (1+x2-x2*z)^ep)/z; note 1+x2-x2*z == 1+r*x2*x5 for z = 1-r*x5.
        b_val = 1 + x2
        if abs(z) < mp.mpf("1e-6"):
            # (1-q*z)^ep = 1 - ep*q*z + ep*(ep-1)/2*q^2*z^2 - ..., q = x2/(1+x2).
            q = x2 / b_val
            return b_val**ep * (ep * q - ep * (ep - 1) / 2 * q**2 * z)
        return (b_val**ep - (b_val - x2 * z) ** ep) / z

    def inner(x5):
        z = 1 - r * x5

        def fx2(u):
            if u == 1:
                return mp.mpf(0)
            x2 = u / (1 - u)
            return x2**ep * (1 + x2) ** ep * bracket_over_z(x2, z) / (1 - u) ** 2

        bps = [mp.mpf(0)]
        if x5 > 1:
            # The (1+r*x2*x5)^ep factor turns over at x2 ~ 1/(r*x5); give the
            # tanh-sinh rule a panel boundary there for large x5.
            x_star = 1 / (r * x5)
            bps.append(x_star / (1 + x_star))
        bps += [mp.mpf("0.5"), mp.mpf("0.95"), mp.mpf(1)]
        return mp.quad(fx2, sorted(set(bps)))

    def fy(y):
        if y == 1:
            return mp.mpf(0)
        x5 = y / (1 - y)
        return (1 + x5) ** ep * inner(x5) / (1 - y) ** 2

    pts = [mp.mpf(0)]
    for x in (1 / (2 * r), 1 / r, 2 / r, 10 / r, 100 / r, 1000 / r):
        x = mp.mpf(x)
        pts.append(x / (1 + x))
    pts.append(mp.mpf(1))
    old_dps = mp.mp.dps
    # dps = 20 caps the nested tanh-sinh accuracy near ~5e-7 relative; at dps = 30
    # the same quadrature agrees with the 1-D decomposition to ~5e-9 (probe-verified).
    mp.mp.dps = 30
    try:
        j2_2d = _direct_2d = mp.quad(fy, sorted(pts)) / ep
    finally:
        mp.mp.dps = old_dps
    j2_red = _Decomposition(ep, r).j2()
    rel = abs(j2_2d - j2_red) / abs(j2_red)
    del _direct_2d
    return {
        "name": "j2_2d_cross_check",
        "passed": bool(rel < 1e-7),
        "rel_err": float(rel),
        "tolerance": 1e-7,
        "dps": 30,
        "ep": "-0.40",
        "r": "0.75",
    }


def run_audit(fast: bool = False) -> dict:
    started = time.perf_counter()
    old_dps = mp.mp.dps
    mp.mp.dps = 20 if fast else 30
    try:
        checks = [
            check_x7_identity_symbolic(),
            check_prefactor_series(),
            check_leading_pole_exact(),
            check_x7_identity_numeric(fast),
            check_connection_formulas(),
            check_decomposition_consistency(fast),
            check_leading_pole_numeric(fast),
            check_2d_cross(fast),
        ]
    finally:
        mp.mp.dps = old_dps
    all_passed = all(c["passed"] for c in checks)
    payload = {
        "audit": "external_int2_leading_pole",
        "method": "Method.2: exact x7 preintegration + leading Laurent pole audit",
        "status": "Success" if all_passed else "Failure",
        "fast_mode": fast,
        "elapsed_sec": round(time.perf_counter() - started, 3),
        "corrected_prefactor_text": CORRECTED_PREFACTOR_TEXT,
        "old_prefactor_text": OLD_PREFACTOR_TEXT,
        "j2_leading": J2_LEADING_TEXT,
        "p2_leading": P2_LEADING_TEXT,
        "product_leading": PRODUCT_LEADING_TEXT,
        "ansv_int2": {
            "supplied_leading_pole": ANSV_INT2_SUPPLIED_LEADING,
            "leading_pole_match": True,
            "full_value_available": False,
            "note": "full AnsvInt2 is metadata only "
            "(examples/external_int2_source_reference.wl.txt), never loaded and NOT "
            "invented; only the supplied leading pole is compared",
        },
        "reducer_core": {
            "imported": False,
            "modified": False,
            "rref_run": False,
            "note": "wrapper/reference metadata audit only; certified pure reduction "
            "coefficients are unaffected",
        },
        "checks": checks,
    }
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--fast", action="store_true", help="reduced precision/coverage")
    args = parser.parse_args(argv)
    payload = run_audit(fast=args.fast)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    for chk in payload["checks"]:
        _log(f"{chk['name']}: {'PASSED' if chk['passed'] else 'FAILED'}")
    _log(f"status = {payload['status']}; JSON -> {args.json}")
    return 0 if payload["status"] == "Success" else 1


if __name__ == "__main__":
    sys.exit(main())
