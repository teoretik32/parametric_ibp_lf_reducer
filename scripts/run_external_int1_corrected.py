"""External Int1 (corrected): pure-family certified reduction + external Gamma prefactor.

The full integral is ``ExternalPrefactor1 * Integral[F1, x2>0, x6>0]`` with
``F1 = (1+x2)^ep * (1+x6)^ep * (1+x2+x6)^(-1+ep)``. The reducer works ONLY with the pure
family ``F1`` (``TargetMultiplier = 1``); the Gamma/EulerGamma prefactor and the kinematic
``1/(s*t^2)`` factor live exclusively in this script's wrapper artifact and are NEVER
multiplied into the family, the row system, or the reduction coefficients.

Pipeline (public API only, no low-level internals):

1. parse ``examples/external_int1_corrected_input.wl.txt``; the initial search comes
   entirely from the document's ``Options`` (certificate gate stays default-ON);
2. bounded deterministic adaptive search (<= 3 levels): level 0 is exactly the document's
   initial search; levels 1-2 deepen the m-ranges by -1/-2 (level 2 adds extra samples);
3. success gate (never weakened): ``status == "Success"``, ``all_locally_finite is True``,
   row-span certificate ``"Passed"``;
4. numeric spot check at ``ep = -3/5, s = t = 1`` (mpmath, dps = 40): prefactor * quad(F1)
   vs prefactor * sum(coeff_i * quad(master_i)); rel tolerance ``1e-6``;
5. write ``validation/external_int1_corrected_reduction.m`` (pure reduction, no prefactor),
   ``..._diagnostics.json`` (+ adaptive history + numeric_check block) and
   ``..._full_formula.m`` (prefactor strictly outside; reference Laurent series as text).

The reference Laurent series is an expansion around ``ep = 0`` and is NOT compared
numerically at finite ``ep`` — it is recorded as text only.

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

INPUT_PATH = REPO_ROOT / "examples" / "external_int1_corrected_input.wl.txt"
DEFAULT_OUT = REPO_ROOT / "validation" / "external_int1_corrected_reduction.m"
DEFAULT_JSON = REPO_ROOT / "validation" / "external_int1_corrected_diagnostics.json"
DEFAULT_FULL = REPO_ROOT / "validation" / "external_int1_corrected_full_formula.m"

# The Gamma/EulerGamma prefactor lives ONLY here and in the wrapper artifact; it is never
# multiplied into anything the reducer sees (family, rows, coefficients).
EXTERNAL_PREFACTOR_TEXT = (
    "Exp[2*ep*EulerGamma]*Gamma[1-ep]*Gamma[-ep]^2*Gamma[ep]*Gamma[2*ep]"
    "/(s*t^2*Gamma[-1-3*ep]*Gamma[1+ep])"
)
# Reference value expansion around ep=0 (+ O[ep]); text-only reference, NOT compared
# numerically at finite ep.
REFERENCE_LAURENT_TEXT = "1/(s*t^2)*(1/ep^4 - Pi^2/(12*ep^2) - 43*Zeta[3]/(6*ep) - Pi^4/180)"

NUMERIC_EP = sp.Rational(-3, 5)  # exact substitution point for the spot check (s = t = 1)
NUMERIC_DPS = 40
NUMERIC_REL_TOL = 1e-6

_X2, _X6, _EP = sp.symbols("x2 x6 ep")


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def search_levels() -> tuple[SearchLevel, ...]:
    """Bounded deterministic schedule; level 0 inherits the document's initial search as-is.

    ``default_search_levels`` is unsuitable here: its level 0 would lower ``MaxIBPDegree``
    to 1 and drop the tangent blocks, i.e. NOT run the documented initial search first.
    """

    def deepened_box(depth: int) -> tuple:
        n_range = ((0, 1), (0, 1))
        m_range = tuple((-2 - depth, 0) for _ in range(3))
        return (n_range, m_range)

    return (
        SearchLevel(name="base"),  # exactly the document's Options (box n:0..1 / m:-2..0)
        SearchLevel(name="deepen-1", label_box=deepened_box(1)),
        SearchLevel(name="deepen-2", label_box=deepened_box(2), extra_samples=2),
    )


def master_integrand(label) -> sp.Expr:
    """Absolute integrand of ``J[label]`` for THIS family, label order ``(n2,n6,m0,m1,m2)``."""
    n2, n6, m0, m1, m2 = (int(v) for v in label)
    return (
        _X2**n2
        * _X6**n6
        * (1 + _X2) ** (_EP + m0)
        * (1 + _X6) ** (_EP + m1)
        * (1 + _X2 + _X6) ** (_EP - 1 + m2)
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
    for term in result.terms:
        integrand = _wl(master_integrand(term.label))
        rendered.append(
            f"({term.coefficient_text})*Int[{integrand}, {{x2, 0, Infinity}}, {{x6, 0, Infinity}}]"
        )
    pure = " +\n    ".join(rendered) if rendered else "0"
    lines = [
        "(* External Int1 (corrected). The Gamma/EulerGamma prefactor is applied OUTSIDE the",
        "   reducer: PureReduction contains ONLY the certified pure-family reduction of",
        "   Int[(1+x2)^ep*(1+x6)^ep*(1+x2+x6)^(-1+ep), {x2, 0, Infinity}, {x6, 0, Infinity}]. *)",
        f"ExternalPrefactor1 = {EXTERNAL_PREFACTOR_TEXT};",
        "",
        f"PureReduction = {pure};",
        "",
        "FullIntegralReduction = ExternalPrefactor1*PureReduction;",
        "",
        f"ReferenceLaurentSeries = {REFERENCE_LAURENT_TEXT};",
        "(* ReferenceLaurentSeries is the expansion around ep=0 (+ O[ep]); reference text",
        "   only -- NOT compared numerically at finite ep. *)",
        "",
    ]
    return "\n".join(lines)


def _coefficient_at_ep(term) -> sp.Rational:
    """Exact rational value of a term coefficient at ``ep = NUMERIC_EP``."""
    if term.coefficient is not None:
        expr = sp.sympify(term.coefficient)
    else:
        expr = sp.sympify(term.coefficient_text.replace("^", "**"))
    return sp.Rational(sp.cancel(expr.subs(_EP, NUMERIC_EP)))


def numeric_check(result: ReductionResult) -> dict:
    """Original-vs-reduced quadrature at ``ep=-3/5, s=t=1`` (both sides carry the prefactor)."""
    import mpmath as mp  # hard transitive dependency of sympy; used only in script/tests

    old_dps = mp.mp.dps
    mp.mp.dps = NUMERIC_DPS
    try:
        ep_num = mp.mpf(NUMERIC_EP.p) / mp.mpf(NUMERIC_EP.q)
        prefactor = (
            mp.exp(2 * ep_num * mp.euler)
            * mp.gamma(1 - ep_num)
            * mp.gamma(-ep_num) ** 2
            * mp.gamma(ep_num)
            * mp.gamma(2 * ep_num)
            / (mp.gamma(-1 - 3 * ep_num) * mp.gamma(1 + ep_num))
        )  # s = t = 1

        def quad_orthant(expr: sp.Expr) -> object:
            f = sp.lambdify((_X2, _X6), expr, modules="mpmath")

            def g(u, v):  # x = u/(1-u): maps [0,1)^2 onto the positive orthant
                return f(u / (1 - u), v / (1 - v)) / ((1 - u) ** 2 * (1 - v) ** 2)

            return mp.quad(g, [0, 1], [0, 1])

        lhs_quad = quad_orthant(master_integrand(result.target_label).subs(_EP, NUMERIC_EP))
        rhs_quad = mp.mpf(0)
        for term in result.terms:
            c = _coefficient_at_ep(term)
            weight = mp.mpf(c.p) / mp.mpf(c.q)
            rhs_quad += weight * quad_orthant(master_integrand(term.label).subs(_EP, NUMERIC_EP))

        lhs = prefactor * lhs_quad
        rhs = prefactor * rhs_quad
        abs_diff = abs(lhs - rhs)
        rel_diff = abs_diff / abs(lhs)
        return {
            "ran": True,
            "ep": str(NUMERIC_EP),
            "s": 1,
            "t": 1,
            "dps": NUMERIC_DPS,
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
        description="External Int1 (corrected): certified pure reduction + external prefactor"
    )
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--json", dest="json_path", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--full-formula", dest="full_path", type=Path, default=DEFAULT_FULL)
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

    search = AdaptiveSearchConfig(levels=search_levels(), max_levels=3)
    started = time.perf_counter()
    result = reduce_family_adaptive(family, target, config, search=search)
    elapsed = time.perf_counter() - started
    _log(f"adaptive reduction finished in {elapsed:.1f}s: status={result.status}")

    payload = _diagnostics_payload(result)
    payload["example"] = "ExternalInt1Corrected"
    payload["elapsed_sec"] = round(elapsed, 3)
    payload["external_prefactor"] = EXTERNAL_PREFACTOR_TEXT
    payload["reference_laurent"] = {
        "text": REFERENCE_LAURENT_TEXT,
        "compared_numerically": False,
        "note": "expansion around ep=0; reference text only",
    }

    ok = certified_success(result)
    if ok:
        _log(f"numeric spot check: ep={NUMERIC_EP}, s=t=1, dps={NUMERIC_DPS} ...")
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
