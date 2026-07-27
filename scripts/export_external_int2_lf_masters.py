"""Export the certified External Int2 LF-basis reduction as ReductionResult artifacts.

Reads ONLY existing certified artifacts (no re-derivation):
  * ``examples/external_int2_dimensionless_input.wl.txt``      -- the parametric family
  * ``validation/external_int2_method11_reconstruction_audit.json`` -- validated coefficients (152 records)
  * ``validation/external_int2_lf_certificate.json``           -- Phase C off-sample row-span certificate

Gates (any failure -> non-zero exit, nothing written):
  1. certificate verdict must be exactly "certified";
  2. audit coefficient expressions must match certificate ``claimed_coefficients`` (sympy-exact);
  3. every master label must pass the polyhedral LF check; the target must NOT be LF
     (otherwise an LF basis would be pointless / mislabeled);
  4. the assembled ReductionResult must gate to ``Success`` with AllLocallyFinite -> True.

Outputs:
  * ``validation/external_int2_lf_result.m``       -- Wolfram-style ReductionResult text
  * ``validation/external_int2_lf_full_formula.m`` -- self-contained formula with COMPLETE integrands

Scope note is propagated verbatim from the certificate: this certifies the row-span relation at
exact rational points modulo primes; it does NOT claim the analytic Laurent value.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import sympy as sp

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from parametric_ibp_lf_reducer import is_locally_finite, parse_family_text  # noqa: E402
from parametric_ibp_lf_reducer.result import (  # noqa: E402
    STATUS_SUCCESS,
    build_reduction_result_from_reconstruction,
)
from parametric_ibp_lf_reducer.wolfram_text_export import coeff_to_wolfram_text  # noqa: E402

FAMILY_PATH = REPO_ROOT / "examples" / "external_int2_dimensionless_input.wl.txt"
AUDIT_PATH = REPO_ROOT / "validation" / "external_int2_method11_reconstruction_audit.json"
CERT_PATH = REPO_ROOT / "validation" / "external_int2_lf_certificate.json"
RESULT_OUT = REPO_ROOT / "validation" / "external_int2_lf_result.m"
FORMULA_OUT = REPO_ROOT / "validation" / "external_int2_lf_formula" ".m"  # overwritten below
FORMULA_OUT = REPO_ROOT / "validation" / "external_int2_lf_full_formula.m"

EP, R = sp.symbols("ep r")


def _fail(msg: str) -> None:
    print(f"REFUSE: {msg}", file=sys.stderr)
    raise SystemExit(1)


def _parse_label(text: str) -> tuple[int, ...]:
    return tuple(int(p) for p in text.strip("[]() ").split(","))


def _sympify(expr_text: str) -> sp.Expr:
    return sp.sympify(expr_text, locals={"ep": EP, "r": R})


def _poly_to_wolfram(poly) -> str:
    """Render a SparsePoly as Wolfram text in the family variables."""
    var_syms = sp.symbols("x2 x5 x7")
    total = sp.Integer(0)
    for expts, coeff in poly.terms.items():
        term = coeff.to_sympy()
        for v, e in zip(var_syms, expts):
            term *= v**e
        total += term
    return coeff_to_wolfram_text(sp.expand(total), factor=False)


def _full_integrand_text(fam, label) -> str:
    """COMPLETE integrand at ``label``: monomial and G-polynomial parametric exponents."""
    mono_exps, poly_exps = fam.exponent_at_label(label)
    pieces: list[str] = []
    for name, pe in zip(fam.variables, mono_exps):
        e = sp.simplify(pe.to_sympy())
        if e == 0:
            continue
        pieces.append(name if e == 1 else f"{name}^({coeff_to_wolfram_text(e)})")
    for name, pe in zip(fam.poly_names, poly_exps):
        e = sp.simplify(pe.to_sympy())
        if e == 0:
            continue
        pieces.append(f"{name}^({coeff_to_wolfram_text(e)})")
    return " * ".join(pieces) if pieces else "1"


def main() -> int:
    fam = parse_family_text(FAMILY_PATH.read_text(encoding="utf-8"))
    cert = json.loads(CERT_PATH.read_text(encoding="utf-8"))
    audit = json.loads(AUDIT_PATH.read_text(encoding="utf-8"))

    # Gate 1: certificate verdict.
    if cert.get("verdict") != "certified":
        _fail(f"certificate verdict is {cert.get('verdict')!r}, not 'certified'")
    if audit.get("schema") != "external-int2-method11-reconstruction-audit/v1":
        _fail(f"unexpected audit schema {audit.get('schema')!r}")

    target_label = tuple(cert["target_label"])
    claimed = {_parse_label(k): _sympify(v) for k, v in cert["claimed_coefficients"].items()}

    # Gate 2: audit coefficients must be sympy-identical to certificate claims.
    audit_coeffs: dict[tuple[int, ...], sp.Expr] = {}
    for k, entry in audit["coefficients"].items():
        if entry.get("status") != "validated":
            _fail(f"audit coefficient for {k} has status {entry.get('status')!r}")
        audit_coeffs[_parse_label(k)] = _sympify(entry["expression"])
    if set(audit_coeffs) != set(claimed):
        _fail("audit/certificate label sets differ")
    for lab in sorted(claimed):
        if sp.simplify(claimed[lab] - audit_coeffs[lab]) != 0:
            _fail(f"coefficient mismatch (audit vs certificate) at {lab}")

    # Gate 3: LF flags from the polyhedral checker (no RREF, no cached claims).
    lf_flags = {lab: is_locally_finite(fam, lab) for lab in sorted(claimed)}
    for lab, flag in lf_flags.items():
        if flag is not True:
            _fail(f"master {lab} failed the LF check: {flag!r}")
    if is_locally_finite(fam, target_label) is not False:
        _fail("target label unexpectedly passed the LF check; LF basis claim would be vacuous")

    messages = (
        "coefficient source: " + AUDIT_PATH.name + f" (n_records={audit['n_records']}, "
        f"n_samples={audit.get('n_samples')}, n_primes={audit.get('n_primes')})",
        "independent certificate: " + CERT_PATH.name
        + f" (n_certified_points={cert['n_certified_points']}, n_in_span={cert['n_in_span']}/"
        + f"{cert['n_checks']}, generic_rank={cert['generic_rank_expected']}, "
        + f"primes={cert['primes']})",
        "scope: " + cert["scope_note"],
    )

    # Gate 4: strict Success gate of the package itself.
    result = build_reduction_result_from_reconstruction(
        fam,
        target_label,
        claimed,
        lf_flags,
        reconstruction_verified=True,
        independent_validation_passed=True,
        formal_success=True,
        n_records=int(audit["n_records"]),
        messages=messages,
    )
    if result.status != STATUS_SUCCESS:
        _fail(f"assembled result gated to {result.status!r}: {result.error!r}")
    if result.all_locally_finite is not True:
        _fail(f"AllLocallyFinite is {result.all_locally_finite!r}, expected True")

    provenance = "\n".join([
        "",
        "(* Method.12A: certificate provenance (recorded, no re-derivation). *)",
        "CertificateProvenance = <|",
        '  "CertificateStatus" -> "Passed",',
        '  "CertificateArtifact" -> "validation/external_int2_lf_certificate.json",',
        f'  "CertificatePoints" -> {int(cert["n_certified_points"])},',
        f'  "CertificateChecks" -> {int(cert["n_checks"])},',
        f'  "GenericRank" -> {int(cert["generic_rank_expected"])},',
        '  "SurfacePolicy" -> "convergence_chamber",',
        f'  "SurfaceChamber" -> "ep={cert["chamber_ep"]}"',
        "|>;",
    ])
    RESULT_OUT.write_text(result.wolfram_style_text + "\n" + provenance + "\n", encoding="utf-8")

    # --- self-contained full formula (Method.12A: integral identity, not pointwise) ---
    lines: list[str] = []
    lines.append("(* External Int2: certified four-term LF-basis reduction (Method.11/.12A). *)")
    lines.append("(* " + cert["scope_note"] + " *)")
    lines.append("(* Variables: " + ", ".join(fam.variables)
                 + "; parameters: " + ", ".join(fam.parameters) + ". *)")
    lines.append("(* The equality below holds AFTER integration over the positive orthant; *)")
    lines.append("(* no pointwise integrand identity is implied. *)")
    for name in fam.poly_names:
        lines.append(f"{name} = {_poly_to_wolfram(fam.polynomial(name))};")
    lines.append("")
    dom_txt = ", ".join("{" + v + ", 0, Infinity}" for v in fam.variables)
    lines.append(f"J[f_] := Inactive[Integrate][f, {dom_txt}];")
    lines.append("")
    terms_txt: list[str] = []
    for i, lab in enumerate(sorted(claimed), start=1):
        c_txt = coeff_to_wolfram_text(claimed[lab])
        lines.append(f"C{i} = {c_txt};  (* label {list(lab)}, LocallyFinite -> True *)")
        lines.append(f"LFIntegrand{i} = {_full_integrand_text(fam, lab)};")
        terms_txt.append(f"C{i}*J[LFIntegrand{i}]")
    lines.append("")
    lines.append(f"(* target label {list(target_label)} *)")
    lines.append(f"TargetIntegrand = {_full_integrand_text(fam, target_label)};")
    lines.append("")
    lines.append("ReductionIdentity =")
    lines.append("  J[TargetIntegrand] == " + " + ".join(terms_txt) + ";")
    lines.append("")
    lines.append("(* Method.12A physical wrapper (dimensionful External Int2 normalization). *)")
    lines.append("r = s/t;")
    lines.append("")
    lines.append("ExternalPrefactor2 =")
    lines.append("  Exp[2*ep*EulerGamma]")
    lines.append("  *t^(-3-ep)")
    lines.append("  *Gamma[1-ep]*Gamma[-ep]^3*Gamma[ep]/")
    lines.append("   (Gamma[-1-3*ep]*Gamma[-2*ep]);")
    lines.append("")
    lines.append("FullPhysicalReduction =")
    lines.append("  ExternalPrefactor2 * ReductionIdentity /. r -> s/t;")
    FORMULA_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[M11] wrote {RESULT_OUT.relative_to(REPO_ROOT)} and "
          f"{FORMULA_OUT.relative_to(REPO_ROOT)}")
    print(f"[M11] status={result.status}; AllLocallyFinite={result.all_locally_finite}; "
          f"n_terms={result.diagnostics.n_terms}; n_records={result.diagnostics.n_records}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
