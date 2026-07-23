"""External Int2 full Laurent audit through O(ep^0).

This script is intentionally independent of ``parametric_ibp_lf_reducer``.  It
implements the exact x7 preintegration, the one-dimensional hypergeometric ODE,
the epsilon-recursion for S=rQ, and the compact Laurent coefficients quoted in
``notes/EXTERNAL_INT2_FULL_LAURENT_DERIVATION.md``.

It does not prove or construct an LF basis.  It is an analytic oracle for future
LF-basis work: any certified LF decomposition of External Int2 must reproduce
these coefficients.

Exit codes:
    0  all symbolic checks passed and JSON was written;
    1  at least one check failed (JSON is still written);
    2  usage error.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import mpmath as mp
import sympy as sp

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_JSON = REPO_ROOT / "validation" / "external_int2_full_laurent_audit.json"


class HPL(sp.Function):
    """Formal harmonic/Goncharov polylogarithm with alphabet {0,-1}.

    The first argument is a symbol encoding the word, e.g. ``m1_0_0`` for
    ``(-1,0,0)``.  The class only needs derivative rules; it does not attempt
    symbolic evaluation of general HPLs.
    """

    nargs = 2

    @classmethod
    def eval(cls, key, x):  # noqa: ANN001 - SymPy protocol
        if isinstance(key, str):
            key = sp.Symbol(key)
        name = str(key)
        if name == "empty":
            return sp.Integer(1)
        word = _decode_word(name)
        if all(letter == 0 for letter in word):
            return sp.log(x) ** len(word) / sp.factorial(len(word))
        return None

    def fdiff(self, argindex=2):  # noqa: ANN001 - SymPy protocol
        if argindex != 2:
            raise sp.ArgumentIndexError(self, argindex)
        name = str(self.args[0])
        word = _decode_word(name)
        x = self.args[1]
        first, *rest = word
        if rest:
            tail = HPL(sp.Symbol(_encode_word(rest)), x)
        else:
            tail = sp.Integer(1)
        kernel = x if first == 0 else 1 + x
        return tail / kernel


def _decode_word(name: str) -> tuple[int, ...]:
    return tuple(-1 if token == "m1" else 0 for token in name.split("_"))


def _encode_word(word: Iterable[int]) -> str:
    return "_".join("m1" if letter == -1 else "0" for letter in word)


def hpl(word: tuple[int, ...], r: sp.Symbol) -> sp.Expr:
    return HPL(sp.Symbol(_encode_word(word)), r)


@dataclass(frozen=True)
class CompactCoefficients:
    f_m4: sp.Expr
    f_m3: sp.Expr
    f_m2: sp.Expr
    f_m1: sp.Expr
    f_0: sp.Expr

    def as_dict(self) -> dict[str, str]:
        return {
            "ep^-4": str(self.f_m4),
            "ep^-3": str(self.f_m3),
            "ep^-2": str(self.f_m2),
            "ep^-1": str(self.f_m1),
            "ep^0": str(self.f_0),
        }


def recurrence_solutions(r: sp.Symbol) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    """Return s0,...,s3 for S=rQ=1/ep+s0+ep*s1+... ."""
    log_r = sp.log(r)
    s0 = -sp.Rational(3, 4) * log_r
    s1 = sp.Rational(1, 8) * log_r**2 - sp.pi**2 / 4
    s2 = (
        sp.Rational(1, 8) * log_r**3
        - sp.Rational(1, 2) * hpl((-1, 0, 0), r)
        + sp.Rational(11, 24) * sp.pi**2 * log_r
        - sp.pi**2 / 4 * hpl((-1,), r)
        - sp.zeta(3) / 2
    )
    s3 = (
        -sp.Rational(11, 96) * log_r**4
        - sp.Rational(7, 16) * sp.pi**2 * log_r**2
        + log_r * sp.zeta(3)
        + sp.pi**2 * hpl((0, -1), r)
        + 2 * hpl((0, -1, 0, 0), r)
        + sp.pi**2 / 3 * hpl((-1, 0), r)
        + sp.Rational(3, 2) * hpl((-1, 0, 0, 0), r)
        - sp.pi**2 / 4 * hpl((-1, -1), r)
        - sp.Rational(1, 2) * hpl((-1, -1, 0, 0), r)
        - sp.Rational(1, 2) * hpl((-1,), r) * sp.zeta(3)
        - sp.pi**4 / 9
    )
    return s0, s1, s2, s3


def compact_laurent_coefficients(r: sp.Symbol, t: sp.Symbol) -> CompactCoefficients:
    """Compact full coefficients after factoring 1/(s*t^2), r=s/t."""
    log_r = sp.log(r)
    log_t = sp.log(t)
    hm100 = hpl((-1, 0, 0), r)
    hm1 = hpl((-1,), r)
    h0m1 = hpl((0, -1), r)
    h0m100 = hpl((0, -1, 0, 0), r)
    hm10 = hpl((-1, 0), r)
    hm1000 = hpl((-1, 0, 0, 0), r)
    hm1m1 = hpl((-1, -1), r)
    hm1m100 = hpl((-1, -1, 0, 0), r)

    f_m4 = sp.Integer(-4)
    f_m3 = 3 * log_r + 4 * log_t
    f_m2 = sp.Rational(5, 3) * sp.pi**2 - log_r**2 / 2 - 3 * log_r * log_t - 2 * log_t**2
    f_m1 = (
        2 * hm100
        - log_r**3 / 2
        + log_r**2 * log_t / 2
        + sp.Rational(3, 2) * log_r * log_t**2
        + sp.Rational(2, 3) * log_t**3
        - sp.Rational(7, 3) * sp.pi**2 * log_r
        - sp.Rational(5, 3) * sp.pi**2 * log_t
        + sp.pi**2 * hm1
        + sp.Rational(62, 3) * sp.zeta(3)
    )
    f_0 = (
        sp.Rational(11, 24) * log_r**4
        + log_r**3 * log_t / 2
        + log_r**2 * (-log_t**2 / 4 + sp.Rational(11, 6) * sp.pi**2)
        + log_r * (-log_t**3 / 2 + sp.Rational(7, 3) * sp.pi**2 * log_t - 18 * sp.zeta(3))
        - log_t**4 / 6
        + sp.Rational(5, 6) * sp.pi**2 * log_t**2
        - 2 * log_t * hm100
        - sp.Rational(62, 3) * log_t * sp.zeta(3)
        - 4 * sp.pi**2 * h0m1
        - 8 * h0m100
        - sp.Rational(4, 3) * sp.pi**2 * hm10
        - 6 * hm1000
        + sp.pi**2 * hm1m1
        + 2 * hm1m100
        + hm1 * (-sp.pi**2 * log_t + 2 * sp.zeta(3))
        + sp.Rational(23, 45) * sp.pi**4
    )
    return CompactCoefficients(f_m4, f_m3, f_m2, f_m1, f_0)


def _log_gamma_1_minus(a: int, ep: sp.Symbol, order: int) -> sp.Expr:
    return sp.EulerGamma * a * ep + sum(
        sp.zeta(n) * a**n * ep**n / sp.Integer(n) for n in range(2, order + 1)
    )


def c_epsilon_series(ep: sp.Symbol, order: int = 6) -> sp.Expr:
    """Series of C(ep) in the exact Q differential equation."""
    log_c = (
        sp.log(sp.Rational(3, 4))
        + sp.log(1 + 3 * ep)
        - sp.log(1 + 2 * ep)
        + 2 * _log_gamma_1_minus(2, ep, order)
        - _log_gamma_1_minus(1, ep, order)
        - _log_gamma_1_minus(3, ep, order)
    )
    return sp.series(sp.exp(log_c), ep, 0, order).removeO().expand()


def verify_x7_preintegration() -> dict:
    u, ep, a, b = sp.symbols("u ep A B", positive=True)
    primitive = (a + (b - a) * u) ** ep / (ep * (b - a))
    derivative_ok = sp.simplify(sp.diff(primitive, u) - (a + (b - a) * u) ** (ep - 1)) == 0
    endpoint_ok = sp.simplify(
        primitive.subs(u, 1)
        - primitive.subs(u, 0)
        - (b**ep - a**ep) / (ep * (b - a))
    ) == 0
    return {
        "name": "exact_x7_preintegration",
        "passed": bool(derivative_ok and endpoint_ok),
        "derivative_ok": bool(derivative_ok),
        "endpoint_ok": bool(endpoint_ok),
    }


def verify_recurrence() -> dict:
    ep, r = sp.symbols("ep r", positive=True)
    s0, s1, s2, s3 = recurrence_solutions(r)
    s = 1 / ep + s0 + ep * s1 + ep**2 * s2 + ep**3 * s3
    q = s / r
    c_ser = c_epsilon_series(ep, 6)
    lhs = (
        r**2 * (1 + r) ** 2 * sp.diff(q, r, 2)
        + r * (1 + r) * (3 * (1 + r) + ep * (3 - r)) * sp.diff(q, r)
        + ((1 + r) ** 2 + ep * (3 + r - r**2) + ep**2 * (2 - r)) * q
    )
    rhs = ep * r - (1 + 2 * ep) * (1 - c_ser + c_ser * sp.exp(-ep * sp.log(r)))
    residual = sp.series(sp.expand(lhs - rhs), ep, 0, 4).removeO()
    coefficients: dict[str, str] = {}
    passed = True
    for power in range(-1, 4):
        coefficient = sp.simplify(sp.together(sp.expand(residual).coeff(ep, power).doit()))
        coefficients[str(power)] = str(coefficient)
        passed = passed and coefficient == 0
    return {
        "name": "ode_epsilon_recurrence_s0_s3",
        "passed": bool(passed),
        "residual_coefficients": coefficients,
    }


def _log_rhat_series(ep: sp.Symbol, order: int = 6) -> sp.Expr:
    """Log of the Gamma product entering P2*G1 after EulerGamma cancellation."""
    return sum(
        (4 + 2 * (-1) ** n - 2 * 2**n) * sp.zeta(n) * ep**n / sp.Integer(n)
        for n in range(2, order + 1)
    )


def derive_full_coefficients() -> dict:
    """Assemble the Laurent coefficients from S and the Gamma product."""
    ep, r, t = sp.symbols("ep r t", positive=True)
    s0, s1, s2, s3 = recurrence_solutions(r)
    s_series = 1 / ep + s0 + ep * s1 + ep**2 * s2 + ep**3 * s3

    # P2*G1/ep after using 1/r=t/s and factoring 1/(s*t^2):
    #   -4/ep^3 * exp(-ep*log(t)) * Rhat(ep)
    # where Rhat contains the EulerGamma-cancelled Gamma ratio.
    rhat = sp.series(sp.exp(_log_rhat_series(ep, 6)), ep, 0, 5).removeO()
    kernel = -4 / ep**3 * sp.exp(-ep * sp.log(t)) * rhat
    assembled = sp.series(kernel * s_series, ep, 0, 1).removeO().expand()
    expected = compact_laurent_coefficients(r, t)
    expected_by_power = {
        -4: expected.f_m4,
        -3: expected.f_m3,
        -2: expected.f_m2,
        -1: expected.f_m1,
        0: expected.f_0,
    }
    residuals: dict[str, str] = {}
    passed = True
    for power, target in expected_by_power.items():
        obtained = sp.expand(assembled).coeff(ep, power)
        residual = sp.simplify(sp.together((obtained - target).doit()))
        residuals[str(power)] = str(residual)
        passed = passed and residual == 0
    return {
        "name": "full_laurent_coefficients_through_ep0",
        "passed": bool(passed),
        "residuals": residuals,
        "coefficients": expected.as_dict(),
    }


def _numeric_hpl(word: tuple[int, ...], x: mp.mpf) -> mp.mpf:
    if not word:
        return mp.mpf(1)
    if all(letter == 0 for letter in word):
        return mp.log(x) ** len(word) / mp.factorial(len(word))
    if word == (-1,):
        return mp.log(1 + x)
    first, *rest = word

    def integrand(value):
        if value == 0:
            return mp.mpf(0)
        kernel = 1 / value if first == 0 else 1 / (1 + value)
        return kernel * _numeric_hpl(tuple(rest), value)

    return mp.quad(integrand, [0, x])


def numeric_compact_coefficients(s: mp.mpf, t: mp.mpf, dps: int = 40) -> dict[str, str]:
    """Evaluate compact coefficients numerically at positive s,t.

    This is not a comparison with the notebook expression; it is a stable,
    machine-readable numerical fingerprint for regression tests and external
    cross-checks.
    """
    mp.mp.dps = dps
    r = s / t
    lr = mp.log(r)
    tt = mp.log(t)
    h = {
        "m100": _numeric_hpl((-1, 0, 0), r),
        "m1": _numeric_hpl((-1,), r),
        "0m1": _numeric_hpl((0, -1), r),
        "0m100": _numeric_hpl((0, -1, 0, 0), r),
        "m10": _numeric_hpl((-1, 0), r),
        "m1000": _numeric_hpl((-1, 0, 0, 0), r),
        "m1m1": _numeric_hpl((-1, -1), r),
        "m1m100": _numeric_hpl((-1, -1, 0, 0), r),
    }
    f_m4 = mp.mpf(-4)
    f_m3 = 3 * lr + 4 * tt
    f_m2 = 5 * mp.pi**2 / 3 - lr**2 / 2 - 3 * lr * tt - 2 * tt**2
    f_m1 = (
        2 * h["m100"]
        - lr**3 / 2
        + lr**2 * tt / 2
        + mp.mpf(3) / 2 * lr * tt**2
        + mp.mpf(2) / 3 * tt**3
        - mp.mpf(7) / 3 * mp.pi**2 * lr
        - mp.mpf(5) / 3 * mp.pi**2 * tt
        + mp.pi**2 * h["m1"]
        + mp.mpf(62) / 3 * mp.zeta(3)
    )
    f_0 = (
        mp.mpf(11) / 24 * lr**4
        + lr**3 * tt / 2
        + lr**2 * (-tt**2 / 4 + mp.mpf(11) / 6 * mp.pi**2)
        + lr * (-tt**3 / 2 + mp.mpf(7) / 3 * mp.pi**2 * tt - 18 * mp.zeta(3))
        - tt**4 / 6
        + mp.mpf(5) / 6 * mp.pi**2 * tt**2
        - 2 * tt * h["m100"]
        - mp.mpf(62) / 3 * tt * mp.zeta(3)
        - 4 * mp.pi**2 * h["0m1"]
        - 8 * h["0m100"]
        - mp.mpf(4) / 3 * mp.pi**2 * h["m10"]
        - 6 * h["m1000"]
        + mp.pi**2 * h["m1m1"]
        + 2 * h["m1m100"]
        + h["m1"] * (-mp.pi**2 * tt + 2 * mp.zeta(3))
        + mp.mpf(23) / 45 * mp.pi**4
    )
    return {
        "s": mp.nstr(s, 18),
        "t": mp.nstr(t, 18),
        "r": mp.nstr(r, 18),
        "F_-4": mp.nstr(f_m4, dps),
        "F_-3": mp.nstr(f_m3, dps),
        "F_-2": mp.nstr(f_m2, dps),
        "F_-1": mp.nstr(f_m1, dps),
        "F_0": mp.nstr(f_0, dps),
    }


def run_audit(*, numeric: bool, dps: int) -> dict:
    checks = [verify_x7_preintegration(), verify_recurrence(), derive_full_coefficients()]
    payload: dict = {
        "schema": "external-int2-full-laurent-audit-v1",
        "scope": (
            "Independent analytic oracle through ep^0; does not construct or certify an LF basis."
        ),
        "checks": checks,
        "passed": all(check["passed"] for check in checks),
        "reference": {
            "source": "ParametricInt examples 4 ChatGPT v2.nb / AnsvInt2",
            "comparison_claim": (
                "Compact coefficients map to the source GPL expression under r=s/t and GPL scaling."
            ),
        },
    }
    if numeric:
        payload["numeric_fingerprints"] = [
            numeric_compact_coefficients(mp.mpf("0.75"), mp.mpf("1"), dps),
            numeric_compact_coefficients(mp.mpf("1.7"), mp.mpf("0.8"), dps),
        ]
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--numeric", action="store_true", help="add HPL numerical fingerprints")
    parser.add_argument("--dps", type=int, default=40)
    args = parser.parse_args(argv)
    if args.dps < 20:
        parser.error("--dps must be at least 20")

    result = run_audit(numeric=args.numeric, dps=args.dps)
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
