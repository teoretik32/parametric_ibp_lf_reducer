"""Corrected Example 4* (known-value-only): certified LF reduction via linearity.

The corrected integrand carries an extra polynomial multiplier ``15*ep + 24*ep*x7``,
which is not parameter-only, so ``TargetMultiplier`` cannot express it. Instead of
touching the math core, this script certifies the corrected reduction by linearity
over ONE shared row system:

    I_corrected = 15*ep * J[{0,0,0,0,0,0,0}] + 24*ep * J[{0,1,0,0,0,0,0}]

(label order ``(n4, n7, n8, m0, m1, m2, m3)``; the ``x7``-monomial of the multiplier
is exactly the ``n7 = +1`` shift). Pipeline:

1. parse ``examples/example4_star_corrected_input.wl.txt`` (the LHS decomposition is
   read from the document's ``Options -> "LHSTerms"``, so the claim lives with the data);
2. enumerate ONE label box, generate ONE row system, evaluate local finiteness once;
3. run the full strict reduction gate (LF + row-span certificate, default-ON) for each
   LHS label via ``reduce_rows_once`` over the SAME rows;
4. require every sub-run to be ``Success`` with pairwise-equal ``selected_rank``;
5. combine the two reductions symbolically (exact SymPy ``cancel``; zero terms drop);
6. re-certify the COMBINED relation ``sum_j L_j*J[lhs_j] - sum_i C_i*J[master_i] = 0``
   in the row span at independent off-sample points (generic ``lhs_terms`` extension of
   the certificate; never a patched verdict);
7. write ``validation/example4_star_corrected_result.m`` and
   ``validation/example4_star_corrected_diagnostics.json``.

Exit codes: 0 = Success (all gates), 1 = honest failure (artifacts still written),
2 = usage/input problem.

Heavy (two reductions over a 972-label box; the single-target 648-label run took
~45-50 min on this machine). Known-value-only policy (docs/05): the notebook provides
only the value expansion, not a reference LF decomposition; nothing here compares
against it.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from contextlib import contextmanager
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:  # allow running without installation
    sys.path.insert(0, str(REPO_ROOT / "src"))

import sympy as sp  # noqa: E402

from parametric_ibp_lf_reducer import parse_family_text, zero_label  # noqa: E402
from parametric_ibp_lf_reducer.api import build_reducer_config  # noqa: E402
from parametric_ibp_lf_reducer.cli import _diagnostics_payload  # noqa: E402
from parametric_ibp_lf_reducer.reducer import (  # noqa: E402
    CERTIFICATE_PASSED,
    _default_certificate_points,
    _enumerate_labels,
    _generate_rows,
    _run_certificate_step,
    reduce_rows_once,
)
from parametric_ibp_lf_reducer.result import (  # noqa: E402
    FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE,
    FAILURE_VERIFICATION_FAILED,
    STATUS_SUCCESS,
)
from parametric_ibp_lf_reducer.timing import new_stage_timings  # noqa: E402
from parametric_ibp_lf_reducer.valuations import is_locally_finite  # noqa: E402

INPUT_PATH = REPO_ROOT / "examples" / "example4_star_corrected_input.wl.txt"
DEFAULT_OUT = REPO_ROOT / "validation" / "example4_star_corrected_result.m"
DEFAULT_JSON = REPO_ROOT / "validation" / "example4_star_corrected_diagnostics.json"


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


@contextmanager
def _timed(store: dict, key: str):
    """Perf.4: accumulate wall-clock seconds of the ``with`` body into ``store[key]``.

    Pure observability — no effect on math results, statuses, or gates.
    """
    t0 = time.perf_counter()
    try:
        yield
    finally:
        store[key] = store.get(key, 0.0) + (time.perf_counter() - t0)


def _target_key(label: tuple[int, ...]) -> str:
    """Stable per-target timing key: zero label -> ``target_zero``, n7=+1 shift -> ``target_x7``."""
    if all(x == 0 for x in label):
        return "target_zero"
    if label[1] == 1 and all(x == 0 for i, x in enumerate(label) if i != 1):
        return "target_x7"
    return "target_" + "_".join(str(x) for x in label)


#: Subrun stage keys reported per target (from ``diagnostics.extra["timings"]``).
_SUBRUN_KEYS = (
    "records_total",
    "ranking_once",
    "assemble_rows_mod_p",
    "rref_mod_p",
    "extract_normal_form",
    "reconstruction",
    "certificate_total",
    "row_generation_total",  # must be 0.0: rows are shared, not regenerated per target
)


def _parse_label(text) -> tuple[int, ...]:
    if isinstance(text, str):
        inner = text.strip().strip("{}")
        return tuple(int(p) for p in inner.split(",") if p.strip())
    return tuple(int(x) for x in text)


def _wl(expr) -> str:
    return sp.sstr(expr).replace("**", "^").replace(" ", "")


def lhs_terms_from_document(family) -> dict[tuple[int, ...], sp.Expr]:
    """Read the claimed LHS decomposition from ``Options -> "LHSTerms"``.

    Keys are labels, values are parameter-only coefficient strings (Wolfram-like ``^``).
    """
    raw = (family.options or {}).get("LHSTerms")
    if not isinstance(raw, dict) or not raw:
        raise ValueError('document Options must provide a non-empty "LHSTerms" association')
    width = family.nvars + family.npolys
    out: dict[tuple[int, ...], sp.Expr] = {}
    for key, value in raw.items():
        label = _parse_label(key)
        if len(label) != width:
            raise ValueError(f"LHSTerms label {key!r} has length {len(label)}, expected {width}")
        out[label] = sp.sympify(str(value).replace("^", "**"))
    return out


def result_coefficients(result) -> dict[tuple[int, ...], sp.Expr]:
    """Exact SymPy coefficients of a reduction result, keyed by master label."""
    coeffs: dict[tuple[int, ...], sp.Expr] = {}
    for term in result.terms:
        if term.coefficient is not None:
            expr = sp.sympify(sp.sstr(term.coefficient))
        else:
            expr = sp.sympify(term.coefficient_text.replace("^", "**"))
        coeffs[tuple(term.label)] = expr
    return coeffs


def combine_coefficients(parts) -> dict[tuple[int, ...], sp.Expr]:
    """Combine ``(weight, {label: coeff})`` parts linearly; exact cancel, zeros dropped."""
    combined: dict[tuple[int, ...], sp.Expr] = {}
    for weight, coeffs in parts:
        w = sp.sympify(weight)
        for label, coeff in coeffs.items():
            label = tuple(label)
            combined[label] = combined.get(label, sp.Integer(0)) + w * sp.sympify(coeff)
    out: dict[tuple[int, ...], sp.Expr] = {}
    for label in sorted(combined):
        expr = sp.cancel(sp.together(combined[label]))
        if expr.is_zero:
            continue
        out[label] = expr
    return out


def _json_safe(obj):
    if isinstance(obj, dict):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(x) for x in obj]
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj
    return str(obj)


def _selected_rank(result) -> int | None:
    extra = result.diagnostics.extra
    cert = extra.get("certificate") or {}
    if cert.get("selected_rank") is not None:
        return cert["selected_rank"]
    return (extra.get("record_selection") or {}).get("selected_rank")


def _label_text(label) -> str:
    return "{" + ",".join(str(x) for x in label) + "}"


def _write_wolfram(path: Path, payload: dict) -> None:
    lhs_items = ",\n    ".join(
        f'<| "Label" -> {_label_text(lab)}, "Coefficient" -> {txt} |>'
        for lab, txt in payload["lhs_terms"]
    )
    term_items = ",\n    ".join(
        f'<| "Label" -> {_label_text(t["label"])}, "Coefficient" -> {t["coefficient"]}, '
        f'"LocallyFinite" -> {"True" if t["locally_finite"] is True else "False"} |>'
        for t in payload["terms"]
    )
    text = (
        "<|\n"
        f'  "Example" -> "Example4StarCorrected",\n'
        f'  "Status" -> "{payload["status"]}",\n'
        f'  "Error" -> {"None" if payload["error"] is None else json.dumps(payload["error"])},\n'
        f'  "Target" -> "{payload["target_text"]}",\n'
        f'  "LHSTerms" -> {{\n    {lhs_items}\n  }},\n'
        f'  "AllLocallyFinite" -> {"True" if payload["all_locally_finite"] is True else "False"},\n'
        f'  "Terms" -> {{\n    {term_items}\n  }},\n'
        f'  "CertificateStatus" -> "{payload["certificate_status"]}",\n'
        f'  "SelectedRank" -> {payload["selected_rank"]},\n'
        f'  "Note" -> "Known-value-only example (docs/05): certified row-span reduction of the '
        "corrected integrand via linearity over one shared row system; no reference LF "
        'decomposition exists in the notebook."\n'
        "|>\n"
    )
    path.write_text(text, encoding="utf-8")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument(
        "--cert-points",
        type=int,
        default=5,
        help="off-sample points for the combined-relation certificate (default 5)",
    )
    args = parser.parse_args(argv)

    # Ensure standard output directories exist up front. ``outputs/`` is
    # untracked and can be wiped by ``git clean``; recreating it here keeps
    # log-redirect targets (e.g. ``> outputs/...log``) valid on future runs.
    (REPO_ROOT / "outputs").mkdir(parents=True, exist_ok=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.json.parent.mkdir(parents=True, exist_ok=True)

    perf: dict[str, float] = {}  # Perf.4: script-level stage timings (seconds)

    if not INPUT_PATH.is_file():
        print(f"input not found: {INPUT_PATH}", file=sys.stderr)
        return 2
    with _timed(perf, "parse_input"):
        family = parse_family_text(INPUT_PATH.read_text(encoding="utf-8"))
        _, config = build_reducer_config(family)
        try:
            lhs = lhs_terms_from_document(family)
        except ValueError as exc:
            print(f"bad LHSTerms: {exc}", file=sys.stderr)
            return 2
    _log(f"LHS decomposition: {[(k, sp.sstr(v)) for k, v in sorted(lhs.items())]}")

    with _timed(perf, "family_label_box"):
        labels = _enumerate_labels(family, config)
    _log(f"labels: {len(labels)} (box {config.label_box})")
    with _timed(perf, "lf_flags_shared"):
        lf_map = {lab: is_locally_finite(family, lab) for lab in labels}
    row_timings = new_stage_timings()
    with _timed(perf, "row_generation_shared"):
        rows, row_diag = _generate_rows(family, labels, config, row_timings)
    _log(f"rows: {len(rows)} by_kind={row_diag['by_kind']}")
    _log(
        "row generation breakdown (s): "
        + ", ".join(
            f"{k}={row_timings[k]:.1f}"
            for k in ("algebraic_rows", "coordinate_rows", "tangent_fields", "tangent_rows")
        )
    )

    results: dict[tuple[int, ...], object] = {}
    sub_timings: dict[tuple[int, ...], dict[str, float]] = {}
    for target in sorted(lhs):
        tkey = _target_key(target)
        _log(f"reducing target {_label_text(target)} ...")
        with _timed(perf, f"{tkey}_reduction"):
            res = reduce_rows_once(
                family,
                target,
                labels,
                rows,
                config.primes,
                config.samples,
                lf_flags=lf_map,
                preferred_masters=config.preferred_masters,
                min_valid_records=config.min_valid_records,
            )
        sub_t = {k: float(v) for k, v in (res.diagnostics.extra.get("timings") or {}).items()}
        sub_timings[target] = sub_t
        _log(
            f"target {_label_text(target)}: status={res.status} "
            f"terms={len(res.terms)} rank={_selected_rank(res)} "
            f"({perf[f'{tkey}_reduction']:.0f}s)"
        )
        _log(
            f"  {tkey} stages (s): "
            + ", ".join(f"{k}={sub_t.get(k, 0.0):.1f}" for k in _SUBRUN_KEYS)
        )
        if sub_t.get("row_generation_total", 0.0) > 0.0:
            _log(
                f"  WARNING: {tkey} regenerated rows internally "
                f"({sub_t['row_generation_total']:.1f}s) — rows were NOT shared!"
            )
        results[target] = res

    all_success = all(res.status == STATUS_SUCCESS for res in results.values())
    ranks = {tgt: _selected_rank(res) for tgt, res in results.items()}
    ranks_agree = len({r for r in ranks.values()}) == 1 and None not in ranks.values()

    combined: dict[tuple[int, ...], sp.Expr] = {}
    cert: dict = {"certificate_status": "NotRun"}
    all_lf: object = "Unknown"
    error: str | None = None

    if not all_success:
        failed = {(_label_text(t)): r.status for t, r in results.items() if not r.success}
        status = next(r.status for r in results.values() if not r.success)
        error = f"sub-reduction failed: {failed}"
    elif not ranks_agree:
        status = FAILURE_VERIFICATION_FAILED
        error = f"selected ranks disagree across sub-reductions: {ranks}"
    else:
        with _timed(perf, "combined_symbolic_merge"):
            combined = combine_coefficients(
                [(lhs[tgt], result_coefficients(res)) for tgt, res in results.items()]
            )
            lf_flags = {lab: lf_map.get(lab, is_locally_finite(family, lab)) for lab in combined}
            all_lf = all(v is True for v in lf_flags.values())
        _log(f"combined: {len(combined)} terms, all_locally_finite={all_lf}")

        points = _default_certificate_points(config.samples, n=max(1, args.cert_points))
        with _timed(perf, "combined_certificate"):
            cert = _run_certificate_step(
                family,
                rows,
                zero_label(family.nvars, family.npolys),  # reporting only; LHS comes from lhs_terms
                dict(combined),
                points,
                config.primes,
                selected_rank=next(iter(ranks.values())),
                min_points=max(1, config.min_certificate_points),
                lhs_terms=lhs,
            )
        _log(
            f"combined certificate: {cert['certificate_status']} "
            f"({cert['n_certificate_points_passed']}/{cert['n_certificate_points']} passed, "
            f"filtered {cert['n_certificate_rank_filtered']}, "
            f"exceeded {cert['n_certificate_rank_exceeded']}, bad {cert['n_certificate_bad_points']})"
        )

        if all_lf is not True:
            status = FAILURE_NORMAL_FORM_NOT_LOCALLY_FINITE
            error = f"combined term(s) not locally finite: {sorted(lab for lab, v in lf_flags.items() if v is not True)}"
        elif cert["certificate_status"] != CERTIFICATE_PASSED:
            status = FAILURE_VERIFICATION_FAILED
            error = f"combined-relation certificate not passed: {cert['certificate_status']}"
        else:
            status = STATUS_SUCCESS

    # --- Perf.4 stage timing summary (observability only) ---------------------------------
    rows_generated_once = all(
        st.get("row_generation_total", 0.0) == 0.0 for st in sub_timings.values()
    )
    subrun_sums = {k: sum(st.get(k, 0.0) for st in sub_timings.values()) for k in _SUBRUN_KEYS}
    sharing_note = (
        "shared across targets: labels, rows, lf_flags (computed once in this script); "
        "recomputed per target inside reduce_rows_once: ranking_once, assemble_rows_mod_p, "
        "rref_mod_p, extract_normal_form, records, reconstruction, certificate "
        "(rank_labels(..., target=...) is target-dependent, so RREF/ranking reuse is NOT a "
        "one-line change; not done here per Perf.4 scope)"
    )
    _log("=== Perf.4 stage timings (s) ===")
    for key in sorted(perf):
        _log(f"  {key}: {perf[key]:.1f}")
    _log(
        "  subrun totals across targets (s): "
        + ", ".join(f"{k}={subrun_sums[k]:.1f}" for k in _SUBRUN_KEYS)
    )
    _log(f"  rows_generated_once={rows_generated_once}")
    _log(f"  sharing: {sharing_note}")

    target_text = " + ".join(f"{_wl(lhs[tgt])}*J[{_label_text(tgt)}]" for tgt in sorted(lhs))
    payload = {
        "example": "example4_star_corrected",
        "status": status,
        "success": status == STATUS_SUCCESS,
        "error": error,
        "target_text": target_text,
        "lhs_terms": [(tuple(t), _wl(lhs[t])) for t in sorted(lhs)],
        "all_locally_finite": all_lf,
        "terms": [
            {
                "label": list(lab),
                "coefficient": _wl(coeff),
                "locally_finite": lf_map.get(lab, "Unknown"),
            }
            for lab, coeff in sorted(combined.items())
        ],
        "certificate_status": cert.get("certificate_status", "NotRun"),
        "selected_rank": next(iter(ranks.values())) if ranks_agree else None,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    _write_wolfram(args.out, payload)

    json_doc = {
        "example": "example4_star_corrected",
        "input": str(INPUT_PATH.relative_to(REPO_ROOT)),
        "status": status,
        "success": payload["success"],
        "error": error,
        "target": target_text,
        "lhs_terms": [{"label": list(t), "coefficient": _wl(lhs[t])} for t in sorted(lhs)],
        "combined": {
            "all_locally_finite": all_lf,
            "terms": payload["terms"],
            "certificate": _json_safe(cert),
        },
        "selected_ranks": {_label_text(t): r for t, r in ranks.items()},
        "subruns": {_label_text(t): _diagnostics_payload(res) for t, res in results.items()},
        "config": {
            "label_box": _json_safe(config.label_box),
            "n_labels": len(labels),
            "n_rows": len(rows),
            "rows_by_kind": row_diag["by_kind"],
            "max_ibp_degree": config.max_ibp_degree,
            "tangent_degree_blocks": _json_safe(config.tangent_degree_blocks),
            "primes": list(config.primes),
            "n_samples": len(list(config.samples)),
            "min_valid_records": config.min_valid_records,
            "n_certificate_points": max(1, args.cert_points),
        },
        "perf4_timings": {
            "script_stages_seconds": {k: float(v) for k, v in sorted(perf.items())},
            "row_generation_breakdown_seconds": {
                k: float(row_timings[k])
                for k in (
                    "row_generation_total",
                    "algebraic_rows",
                    "coordinate_rows",
                    "tangent_fields",
                    "tangent_rows",
                )
            },
            "subrun_stage_seconds": {
                _label_text(t): {k: float(v) for k, v in sorted(st.items())}
                for t, st in sub_timings.items()
            },
            "subrun_stage_sums_seconds": {k: float(v) for k, v in subrun_sums.items()},
            "rows_generated_once": rows_generated_once,
            "sharing": sharing_note,
        },
        "note": (
            "Known-value-only (docs/05): corrected integrand certified via linearity "
            "I_corrected = sum_j L_j*J[lhs_j] over one shared row system; combined relation "
            "re-certified in the row span at independent off-sample points. No reference LF "
            "decomposition exists; the notebook value expansion is not a reducer coefficient."
        ),
    }
    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(
        json.dumps(_json_safe(json_doc), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    _log(f"status={status}; wrote {args.out} and {args.json}")
    return 0 if payload["success"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
