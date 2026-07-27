"""Method.12A: HyperInt linear-reducibility audit of the four certified LF masters.

Pure symbolic audit (sympy only).  For every certified master and every ordering of
(x2, x5, x7) we run the classical Fubini polynomial reduction:

  * to integrate over ``v`` every active polynomial must be linear in ``v``;
  * writing ``p_i = a_i*v + b_i`` the reduced set is the irreducible nonconstant
    factors of ``{a_i, b_i}`` together with the pairwise resultants
    ``a_i*b_j - a_j*b_i``.

An order is *valid* when every step is linear.  The union of all polynomials met
along the chosen order (plus the integration variables) is the alphabet upper
bound for the hyperlogarithmic letters.

The script does NOT touch the certified LF basis: labels, coefficients and the
surface/LF/certificate semantics are read-only inputs.  No RREF, no modular
records, no solves.  Outputs:

  * validation/external_int2_lf_linear_reducibility.json
  * validation/hyperint/external_int2_lf_L{1..4}.mpl   (only if all four masters
    are linearly reducible; exact inputs, integration NOT started).
"""
from __future__ import annotations

import json
import sys
from itertools import combinations, permutations
from pathlib import Path

import sympy as sp

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from parametric_ibp_lf_reducer import parse_family_text  # noqa: E402

FAMILY_PATH = REPO_ROOT / "examples" / "external_int2_dimensionless_input.wl.txt"
CERT_PATH = REPO_ROOT / "validation" / "external_int2_lf_certificate.json"
JSON_OUT = REPO_ROOT / "validation" / "external_int2_lf_linear_reducibility.json"
HYPERINT_DIR = REPO_ROOT / "validation" / "hyperint"

EP, R = sp.symbols("ep r")
# required epsilon expansion order per master (L1..L4 in sorted-label order)
EPS_ORDER = {1: 2, 2: 3, 3: 4, 4: 4}


def _fail(msg: str) -> None:
    print(f"REFUSE: {msg}", file=sys.stderr)
    raise SystemExit(1)


def _parse_label(text: str) -> tuple[int, ...]:
    return tuple(int(p) for p in text.strip("[]() ").split(","))


def _poly_sympy(poly, var_syms) -> sp.Expr:
    total = sp.Integer(0)
    for expts, coeff in poly.terms.items():
        term = coeff.to_sympy()
        for v, e in zip(var_syms, expts):
            term *= v**e
        total += term
    return sp.expand(total)


def _canon(expr: sp.Expr, var_syms) -> sp.Expr | None:
    """Canonical primitive form of an irreducible factor; None if constant in vars."""
    expr = sp.expand(expr)
    if not (expr.free_symbols & set(var_syms)):
        return None
    poly = sp.Poly(expr, *var_syms, R)
    _, prim = poly.primitive()
    if prim.LC() < 0:
        prim = -prim
    return prim.as_expr()


def _factors(expr: sp.Expr, var_syms) -> set[sp.Expr]:
    out: set[sp.Expr] = set()
    if expr == 0:
        return out
    for base, _mult in sp.factor_list(sp.expand(expr), *var_syms, R)[1]:
        c = _canon(base, var_syms)
        if c is not None:
            out.add(c)
    return out


def _reduce_order(polys: set[sp.Expr], order, var_syms):
    """Fubini polynomial reduction along ``order``.  Returns (valid, steps, obstruction)."""
    active = set(polys)
    steps = []
    for v in order:
        for p in active:
            if sp.degree(p, gen=v) > 1:
                return False, steps, {"variable": str(v), "polynomial": str(p)}
        containing = [p for p in active if v in p.free_symbols]
        nxt: set[sp.Expr] = {p for p in active if v not in p.free_symbols}
        ab = [(sp.expand(p.coeff(v, 1)), sp.expand(p.coeff(v, 0))) for p in containing]
        for a, b in ab:
            nxt |= _factors(a, var_syms) | _factors(b, var_syms)
        for (a1, b1), (a2, b2) in combinations(ab, 2):
            nxt |= _factors(a1 * b2 - a2 * b1, var_syms)
        steps.append({
            "variable": str(v),
            "polynomials_before": sorted(str(p) for p in active),
            "reduction_set_after": sorted(str(p) for p in nxt),
        })
        active = nxt
    return True, steps, None


def _maple(expr: sp.Expr) -> str:
    return str(expr).replace("**", "^")


def _integrand_maple(fam, label, poly_map, var_syms) -> str:
    mono_exps, poly_exps = fam.exponent_at_label(label)
    pieces: list[str] = []
    for v, pe in zip(var_syms, mono_exps):
        e = sp.simplify(pe.to_sympy())
        if e != 0:
            pieces.append(f"{v}^({_maple(e)})")
    for name, pe in zip(fam.poly_names, poly_exps):
        e = sp.simplify(pe.to_sympy())
        if e != 0:
            pieces.append(f"({_maple(poly_map[name])})^({_maple(e)})")
    return " * ".join(pieces) if pieces else "1"


def _write_hyperint_input(path: Path, idx, label, integrand, order, eps_order) -> None:
    lines = [
        f"# External Int2 Method.12A -- exact HyperInt input for LF master L{idx}",
        f"# label {list(label)}; certified basis unchanged; audit artifact only.",
        "# REVIEW BEFORE RUNNING: prepared input; do NOT start a long integration",
        "# unattended.  Uncomment the hyperInt call after review.",
        'read "HyperInt.mpl":',
        "_hyper_verbosity := 0:",
        f"epsOrder := {eps_order}:  # expand through ep^{eps_order}",
        f"f := {integrand}:",
        "fser := convert(series(f, ep=0, epsOrder+1), polynom):",
        f"intOrder := [{', '.join(order)}]:  # valid order from linear-reducibility audit",
        "# result := hyperInt(fser, intOrder):  # <-- uncomment after review",
        "# fibrationBasis(result):",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    fam = parse_family_text(FAMILY_PATH.read_text(encoding="utf-8"))
    cert = json.loads(CERT_PATH.read_text(encoding="utf-8"))
    if cert.get("verdict") != "certified":
        _fail(f"certificate verdict is {cert.get('verdict')!r}, expected 'certified'")
    var_syms = sp.symbols(" ".join(fam.variables))
    poly_map = {name: _poly_sympy(fam.polynomial(name), var_syms) for name in fam.poly_names}
    labels = sorted(_parse_label(k) for k in cert["claimed_coefficients"])
    if len(labels) != 4:
        _fail(f"expected 4 certified masters, found {len(labels)}")

    masters = []
    for idx, label in enumerate(labels, start=1):
        _mono_exps, poly_exps = fam.exponent_at_label(label)
        initial: set[sp.Expr] = set(var_syms)
        for name, pe in zip(fam.poly_names, poly_exps):
            if sp.simplify(pe.to_sympy()) != 0:
                initial |= _factors(poly_map[name], var_syms)
        orders_report = []
        chosen = None
        chosen_steps = None
        for order in permutations(var_syms):
            valid, steps, obstruction = _reduce_order(initial, order, var_syms)
            orders_report.append({
                "order": [str(v) for v in order],
                "valid": valid,
                "steps": steps,
                "obstruction": obstruction,
            })
            if valid and chosen is None:
                chosen = [str(v) for v in order]
                chosen_steps = steps
        alphabet: set[str] = {str(v) for v in var_syms}
        if chosen_steps is not None:
            for st in chosen_steps:
                alphabet |= set(st["polynomials_before"]) | set(st["reduction_set_after"])
        integrand = _integrand_maple(fam, label, poly_map, var_syms)
        masters.append({
            "index": idx,
            "name": f"L{idx}",
            "label": list(label),
            "integrand_maple": integrand,
            "initial_polynomials": sorted(str(p) for p in initial),
            "orders": orders_report,
            "n_valid_orders": sum(1 for o in orders_report if o["valid"]),
            "chosen_order": chosen,
            "alphabet": sorted(alphabet),
            "linearly_reducible": chosen is not None,
            "required_ep_order": EPS_ORDER[idx],
        })

    all_ok = all(m["linearly_reducible"] for m in masters)
    if all_ok:
        HYPERINT_DIR.mkdir(parents=True, exist_ok=True)
        for m in masters:
            path = HYPERINT_DIR / f"external_int2_lf_L{m['index']}.mpl"
            _write_hyperint_input(
                path, m["index"], m["label"], m["integrand_maple"],
                m["chosen_order"], m["required_ep_order"],
            )
            m["hyperint_input"] = str(path.relative_to(REPO_ROOT)).replace("\\", "/")

    report = {
        "schema": "external-int2-lf-linear-reducibility/v1",
        "method": "External Int2 Method.12A",
        "family": "examples/external_int2_dimensionless_input.wl.txt",
        "certificate": "validation/external_int2_lf_certificate.json",
        "variables": [str(v) for v in var_syms],
        "parameters": ["r"],
        "masters": masters,
        "all_linearly_reducible": all_ok,
        "scope_note": (
            "Audit only: enumerates HyperInt integration orders and polynomial "
            "reduction sets for the four certified LF masters.  The certified LF "
            "basis, labels, coefficients and certificate semantics are unchanged.  "
            "No RREF, no modular records, no integration performed."
        ),
    }
    JSON_OUT.write_text(json.dumps(report, indent=1, sort_keys=False) + "\n", encoding="utf-8")

    for m in masters:
        print(f"[M12A] {m['name']} label={m['label']} reducible={m['linearly_reducible']} "
              f"valid_orders={m['n_valid_orders']}/6 order={m['chosen_order']} "
              f"ep_order={m['required_ep_order']}")
        print(f"[M12A]   alphabet: {m['alphabet']}")
    print(f"[M12A] all_linearly_reducible={all_ok}; wrote {JSON_OUT.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
