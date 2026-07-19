"""External Int2 (dimensionless): pure-family certified reduction + external prefactor.

The full integral is ``P2 * Integral[F2, x2>0, x5>0, x7>0]`` with

    F2 = x2^(1+ep)*(1+x2)^ep*(1+x5)^ep*(1+x7)^(-1-ep)*(1+x7+x2*x7+r*x2*x5)^(-1+ep)

and ``P2 = t^(-3-ep)*Gamma[1-ep]*Gamma[-ep]^3*Gamma[ep]/(Gamma[-1-3*ep]*Gamma[-2*ep])``.
Dimensionless rewrite: ``r = s/t`` is the only kinematic parameter inside the family; the
overall ``t^(-3-ep)`` scaling lives exclusively in ``P2``. The reducer works ONLY with the
pure family ``F2`` (``TargetMultiplier = 1``, parameters ``ep, r``, assumption ``r > 0``);
the Gamma prefactor and the ``t``-scaling are NEVER multiplied into the family, the row
system, or the reduction coefficients.

Pipeline (public API only, no low-level internals):

1. parse ``examples/external_int2_dimensionless_input.wl.txt``; the initial search comes
   entirely from the document's ``Options`` (certificate gate stays default-ON,
   ``rref_backend = "auto"``, deterministic scattered (ep, r) samples by default);
2. bounded deterministic adaptive search (<= 3 levels): levels 0-1 keep the document's
   search box and only add scattered samples/primes (numeric hardening after the first
   run's InterpolationFailed); level 2 additionally deepens the m-ranges by -1;
3. success gate (never weakened): ``status == "Success"``, ``all_locally_finite is True``,
   row-span certificate ``"Passed"``;
4. numeric spot check at ``ep = -3/5, r = 1, t = 1`` (mpmath, dps = 30): the inner
   x7-integral of every master is taken exactly through a Gauss 2F1 kernel (``i7_kernel``,
   cross-validated against a direct quadrature at runtime), the remaining (x2, x5)
   integral uses the ``x = u/(1-u)`` map onto ``[0,1)^2``; LHS = P2*J[target] vs
   RHS = P2*sum(coeff_i*J[label_i]); rel tolerance ``1e-6``. A truncated Laurent series is
   never treated as an exact value at finite ``ep``.
5. write ``validation/external_int2_result.m`` (pure reduction, no prefactor),
   ``..._diagnostics.json`` (+ adaptive history + numeric_check block) and
   ``..._full_formula.m`` (prefactor strictly outside).

Reference value: ``AnsvInt2`` is NOT available in this repository and is NOT invented
here. If a source reference (e.g. a GPL ``G[...]`` expression) is added later it must be
preserved as metadata only — GPL values are never reducer coefficients.

Exit codes: 0 = certified Success + numeric check passed, 1 = honest failure (artifacts
still written), 2 = usage/input problem.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:  # allow running without installation
    sys.path.insert(0, str(REPO_ROOT / "src"))

import sympy as sp  # noqa: E402

from parametric_ibp_lf_reducer import ParserError, parse_family_text  # noqa: E402
from parametric_ibp_lf_reducer.adaptive import (  # noqa: E402
    AdaptiveSearchConfig,
    SearchLevel,
    reduce_family_adaptive,
)
from parametric_ibp_lf_reducer.api import build_reducer_config  # noqa: E402
from parametric_ibp_lf_reducer.cli import _diagnostics_payload  # noqa: E402
from parametric_ibp_lf_reducer.reducer import CERTIFICATE_PASSED  # noqa: E402
from parametric_ibp_lf_reducer.result import STATUS_SUCCESS, ReductionResult  # noqa: E402

INPUT_PATH = REPO_ROOT / "examples" / "external_int2_dimensionless_input.wl.txt"
DEFAULT_OUT = REPO_ROOT / "validation" / "external_int2_result.m"
DEFAULT_JSON = REPO_ROOT / "validation" / "external_int2_diagnostics.json"
DEFAULT_FULL = REPO_ROOT / "validation" / "external_int2_full_formula.m"

# The Gamma prefactor and the t-scaling live ONLY here and in the wrapper artifact; they
# are never multiplied into anything the reducer sees (family, rows, coefficients).
EXTERNAL_PREFACTOR_TEXT = (
    "t^(-3-ep)*Gamma[1-ep]*Gamma[-ep]^3*Gamma[ep]/(Gamma[-1-3*ep]*Gamma[-2*ep])"
)
# AnsvInt2 is not available in this repository; it is deliberately NOT invented. A GPL
# G[...] source expression, if ever added, is metadata only — never a reducer coefficient.
REFERENCE_METADATA = {
    "name": "AnsvInt2",
    "available": False,
    "compared_numerically": False,
    "note": (
        "no source reference bundled; not invented; GPL G[...] values are metadata only "
        "and are never reducer coefficients"
    ),
}

NUMERIC_EP = sp.Rational(-3, 5)  # exact substitution point for the spot check
NUMERIC_R = sp.Integer(1)  # r = s/t = 1 (t = 1)
NUMERIC_DPS = 30
NUMERIC_REL_TOL = 1e-6

_X2, _X5, _X7, _EP, _R = sp.symbols("x2 x5 x7 ep r")


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def search_levels(
    extra_samples: int = 0, extra_primes: int = 0, extra_expand: int = 0
) -> tuple[SearchLevel, ...]:
    """Bounded deterministic schedule; level 0 inherits the document's initial search as-is.

    ``default_search_levels`` is unsuitable here: its level 0 would lower ``MaxIBPDegree``
    to 1 and drop the tangent blocks, i.e. NOT run the documented initial search first.
    """

    def deepened_box(depth: int) -> tuple:
        n_range = ((0, 1), (0, 1), (0, 1))
        m_range = tuple((-2 - depth, 0) for _ in range(4))
        return (n_range, m_range)

    if extra_samples or extra_primes:
        # 3rd heavy run: run #2 failed with InterpolationFailed at every level even
        # after deepening (rank 22361, ~4.8h) — the reduction itself succeeded at both
        # depths and only coefficient reconstruction failed. With params (ep, r) and 24
        # samples (22 fit points) the degree search could not reach pairs beyond ~(3,3):
        # denser ansatz rows were skipped as underdetermined. Covering the full
        # max_deg=6 rational ansatz needs >= 57 samples (2*28-1 fit + 2 holdout), so the
        # fix is a single boosted base-box level; the document's search box, degrees and
        # certificate gate stay untouched.
        # 4th heavy run: run #3 (boost s48-p6, base box) reconstructed and certified,
        # but the normal form kept 2 non-locally-finite terms (labels (0,0,0,0,-1,0,0)
        # and (0,0,0,0,0,-1,0)); the run's own recommendation is to expand the label
        # box around them so a different reduction path can avoid those sectors. Run #2
        # showed the depth-1 box reaches the reduction target, so --extra-expand=1
        # combines the deepened box with the boosted samples/primes.
        box = deepened_box(extra_expand) if extra_expand else None
        suffix = f"-x{extra_expand}" if extra_expand else ""
        return (
            SearchLevel(
                name=f"base+boost-s{extra_samples}-p{extra_primes}{suffix}",
                label_box=box,
                extra_samples=extra_samples,
                extra_primes=extra_primes,
            ),
        )

    # Retry schedule (2nd heavy run). The 1st run failed with InterpolationFailed at
    # every level while the rank histogram showed exactly one degenerate sample point
    # (3 of 36 reduced records at lower rank; 33/36 kept), leaving too few good samples
    # to validate coefficient reconstruction. Per the run's own recommendation the fix
    # is more scattered samples/primes. The document's search box is kept as-is at
    # level 0 — extra samples/primes only harden the numeric reconstruction and cannot
    # change the certified mathematics.
    return (
        SearchLevel(name="base+samples", extra_samples=6, extra_primes=2),
        SearchLevel(name="base+samples-2", extra_samples=12, extra_primes=4),
        SearchLevel(
            name="deepen-1+samples",
            label_box=deepened_box(1),
            extra_samples=12,
            extra_primes=4,
        ),
    )


def master_integrand(label) -> sp.Expr:
    """Absolute integrand of ``J[label]``, label order ``(n2, n5, n7, m0, m1, m2, m3)``."""
    n2, n5, n7, m0, m1, m2, m3 = (int(v) for v in label)
    return (
        _X2 ** (1 + _EP + n2)
        * _X5**n5
        * _X7**n7
        * (1 + _X2) ** (_EP + m0)
        * (1 + _X5) ** (_EP + m1)
        * (1 + _X7) ** (-1 - _EP + m2)
        * (1 + _X7 + _X2 * _X7 + _R * _X2 * _X5) ** (-1 + _EP + m3)
    )


def _wl(expr: sp.Expr) -> str:
    """Wolfram-like rendering (``^`` instead of ``**``)."""
    return sp.sstr(expr).replace("**", "^")


def certified_success(result: ReductionResult) -> bool:
    """The strict gate: Success + AllLocallyFinite True + row-span certificate Passed."""
    cert = result.diagnostics.extra.get("certificate") or {}
    return (
        result.status == STATUS_SUCCESS
        and result.all_locally_finite is True
        and cert.get("certificate_status") == CERTIFICATE_PASSED
    )


def build_full_formula_text(result: ReductionResult) -> str:
    """Wrapper artifact: prefactor strictly OUTSIDE the certified pure reduction."""
    rendered = []
    ranges = "{x2, 0, Infinity}, {x5, 0, Infinity}, {x7, 0, Infinity}"
    for term in result.terms:
        integrand = _wl(master_integrand(term.label))
        rendered.append(f"({term.coefficient_text})*Int[{integrand}, {ranges}]")
    pure = " +\n    ".join(rendered) if rendered else "0"
    lines = [
        "(* External Int2 (dimensionless, r = s/t). The prefactor P2 (Gamma ratio and the",
        "   t^(-3-ep) scaling) is applied OUTSIDE the reducer: PureReduction contains ONLY",
        "   the certified pure-family reduction of Int[F2, x2>0, x5>0, x7>0] with",
        "   F2 = x2^(1+ep)*(1+x2)^ep*(1+x5)^ep*(1+x7)^(-1-ep)*(1+x7+x2*x7+r*x2*x5)^(-1+ep). *)",
        f"ExternalPrefactor2 = {EXTERNAL_PREFACTOR_TEXT};",
        "",
        f"PureReduction = {pure};",
        "",
        "FullIntegralReduction = ExternalPrefactor2*PureReduction;",
        "",
        "(* Reference value: AnsvInt2 is not available in this repository and is NOT",
        "   invented here. If a source reference (GPL G[...] expression) is added later it",
        "   must stay metadata only -- GPL values are never reducer coefficients. *)",
        "",
    ]
    return "\n".join(lines)


def _coefficient_at_point(term) -> sp.Rational:
    """Exact rational value of a term coefficient at ``ep = NUMERIC_EP, r = NUMERIC_R``."""
    if term.coefficient is not None:
        expr = sp.sympify(term.coefficient)
    else:
        expr = sp.sympify(term.coefficient_text.replace("^", "**"))
    return sp.Rational(sp.cancel(expr.subs({_EP: NUMERIC_EP, _R: NUMERIC_R})))


def i7_kernel(mp, c, b, a_exp, b_exp, k: int):
    """Exact inner x7-integral ``int_0^oo x7^k (1+x7)^A (c + b*x7)^B dx7`` for k in {0, 1}.

    Substituting ``w = 1/(1+x7)`` gives ``int_0^1 w^(s-1) (1-w)^k (b + (c-b)*w)^B dw`` with
    ``s = -A - B - 1 - k``, i.e. a Gauss 2F1 kernel (``int_0^1 w^(s-1)(1+z*w)^B dw =
    2F1(-B, s; s+1; -z)/s``). Requires ``s > 0`` — guaranteed for LF masters of this
    family (``s = 1 - m2 - m3 - k`` with ``m2, m3 <= 0``, ``k <= 1``, not both zero
    when ``k = 1``).
    """
    if k not in (0, 1):
        raise ValueError(f"i7_kernel supports x7-powers 0 and 1, got {k}")
    s = -a_exp - b_exp - 1 - k
    if not mp.re(s) > 0:
        raise ValueError(f"divergent inner x7-integral: s = {s}")
    z = (c - b) / b
    first = mp.hyp2f1(-b_exp, s, s + 1, -z) / s
    if k == 0:
        return b**b_exp * first
    second = mp.hyp2f1(-b_exp, s + 1, s + 2, -z) / (s + 1)
    return b**b_exp * (first - second)


def numeric_check(result: ReductionResult) -> dict:
    """Original-vs-reduced comparison at ``ep = -3/5, r = 1, t = 1``.

    Both sides carry the prefactor (``t = 1``). Each ``J[label]`` is computed as an exact
    2F1 kernel in x7 (validated below against a direct quadrature) integrated over
    ``(x2, x5)`` with the ``x = u/(1-u)`` map. No Laurent-series data is used anywhere.
    """
    import mpmath as mp  # hard transitive dependency of sympy; used only in script/tests

    old_dps = mp.mp.dps
    mp.mp.dps = NUMERIC_DPS
    try:
        ep = mp.mpf(NUMERIC_EP.p) / mp.mpf(NUMERIC_EP.q)
        r = mp.mpf(int(NUMERIC_R))
        prefactor = (
            mp.gamma(1 - ep)
            * mp.gamma(-ep) ** 3
            * mp.gamma(ep)
            / (mp.gamma(-1 - 3 * ep) * mp.gamma(-2 * ep))
        )  # t = 1 -> t^(-3-ep) = 1

        # Runtime self-check of the 2F1 kernel against direct 1-D quadratures (k = 0, 1).
        kernel_check = mp.mpf(0)
        for a_exp, b_exp, k in ((-1 - ep, -1 + ep, 0), (-2 - ep, -1 + ep, 1)):
            c, b = mp.mpf(3.7), mp.mpf(1.9)
            direct = mp.quad(
                lambda x: x**k * (1 + x) ** a_exp * (c + b * x) ** b_exp,
                [0, 1, 10, 100, mp.inf],
            )
            val = i7_kernel(mp, c, b, a_exp, b_exp, k)
            kernel_check = max(kernel_check, abs(val - direct) / abs(direct))

        def quad_label(label) -> object:
            n2, n5, n7, m0, m1, m2, m3 = (int(v) for v in label)
            a_exp = -1 - ep + m2
            b_exp = -1 + ep + m3

            def g(u, v):  # x = u/(1-u): maps [0,1)^2 onto the positive orthant
                x2 = u / (1 - u)
                x5 = v / (1 - v)
                c = 1 + r * x2 * x5
                b = 1 + x2
                return (
                    x2 ** (1 + ep + n2)
                    * x5**n5
                    * (1 + x2) ** (ep + m0)
                    * (1 + x5) ** (ep + m1)
                    * i7_kernel(mp, c, b, a_exp, b_exp, n7)
                    / ((1 - u) ** 2 * (1 - v) ** 2)
                )

            return mp.quad(g, [0, 1], [0, 1])

        lhs_quad = quad_label(result.target_label)
        rhs_quad = mp.mpf(0)
        for term in result.terms:
            coeff = _coefficient_at_point(term)
            weight = mp.mpf(coeff.p) / mp.mpf(coeff.q)
            rhs_quad += weight * quad_label(term.label)

        lhs = prefactor * lhs_quad
        rhs = prefactor * rhs_quad
        abs_diff = abs(lhs - rhs)
        rel_diff = abs_diff / abs(lhs)
        return {
            "ran": True,
            "ep": str(NUMERIC_EP),
            "r": str(NUMERIC_R),
            "t": 1,
            "dps": NUMERIC_DPS,
            "method": (
                "inner x7-integral via exact 2F1 kernel (runtime-validated); "
                "(x2, x5) quadrature via the x = u/(1-u) map"
            ),
            "kernel_check_rel": float(kernel_check),
            "lhs": mp.nstr(lhs, 20),
            "rhs": mp.nstr(rhs, 20),
            "abs_diff": float(abs_diff),
            "rel_diff": float(rel_diff),
            "rel_tol": NUMERIC_REL_TOL,
            "passed": bool(rel_diff < NUMERIC_REL_TOL),
        }
    finally:
        mp.mp.dps = old_dps


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="External Int2 (dimensionless): certified pure reduction + external prefactor"
    )
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--json", dest="json_path", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--full-formula", dest="full_path", type=Path, default=DEFAULT_FULL)
    parser.add_argument(
        "--extra-samples",
        type=int,
        default=0,
        help="replace the schedule with ONE boosted base-box level: extra scattered samples",
    )
    parser.add_argument(
        "--extra-primes",
        type=int,
        default=0,
        help="replace the schedule with ONE boosted base-box level: extra primes",
    )
    parser.add_argument(
        "--extra-expand",
        type=int,
        default=0,
        help="with --extra-samples/--extra-primes: deepen the boosted level's label box by this depth",
    )
    args = parser.parse_args(argv)

    try:
        text = INPUT_PATH.read_text(encoding="utf-8")
    except OSError as exc:
        print(f"cannot read {INPUT_PATH}: {exc}", file=sys.stderr)
        return 2
    try:
        family = parse_family_text(text)
    except ParserError as exc:
        print(f"cannot parse family: {exc}", file=sys.stderr)
        return 2

    target, config = build_reducer_config(family)
    if not config.require_certificate_for_success:  # default-ON; never weakened here
        print("certificate gate unexpectedly disabled; refusing to run", file=sys.stderr)
        return 2
    _log(f"family: vars={family.variables} polys={family.poly_names} target={target}")

    levels = search_levels(args.extra_samples, args.extra_primes, args.extra_expand)
    search = AdaptiveSearchConfig(levels=levels, max_levels=len(levels))
    started = time.perf_counter()
    result = reduce_family_adaptive(family, target, config, search=search)
    elapsed = time.perf_counter() - started
    _log(f"adaptive reduction finished in {elapsed:.1f}s: status={result.status}")

    payload = _diagnostics_payload(result)
    payload["example"] = "ExternalInt2Dimensionless"
    payload["elapsed_sec"] = round(elapsed, 3)
    payload["external_prefactor"] = EXTERNAL_PREFACTOR_TEXT
    payload["dimensionless_rewrite"] = "r = s/t; the t^(-3-ep) scaling lives in P2 only"
    payload["reference_value"] = dict(REFERENCE_METADATA)

    ok = certified_success(result)
    if ok:
        _log(f"numeric spot check: ep={NUMERIC_EP}, r={NUMERIC_R}, t=1, dps={NUMERIC_DPS} ...")
        check = numeric_check(result)
        payload["numeric_check"] = check
        print(
            "numeric check: lhs={lhs} rhs={rhs} rel_diff={rel_diff:.3e} passed={passed}".format(
                **check
            )
        )
        ok = ok and check["passed"]
    else:
        payload["numeric_check"] = {"ran": False, "reason": "no certified Success to check"}

    for path in (args.out, args.json_path, args.full_path):
        path.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(result.wolfram_style_text, encoding="utf-8")
    args.json_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    args.full_path.write_text(build_full_formula_text(result), encoding="utf-8")
    for path in (args.out, args.json_path, args.full_path):
        _log(f"wrote {path}")

    if ok:
        print("Success: certified reduction; prefactor applied outside; numeric check passed")
        return 0
    print(f"Failure: status={result.status} error={result.error}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
