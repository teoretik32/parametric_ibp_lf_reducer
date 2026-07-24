"""External Int2 Method.11: certified LF reconstruction under an explicit chamber surface policy.

Phase B runner. Assembles the Level-0 rank-repair system (baseline blocks + the ``(3, 3)``
extra tangent block) with ``SurfacePolicy.chamber({"ep": -3/5})`` passed through the
first-class core API (no monkey-patching), then runs one fixed certified reduction pass:
``reduce_rows_once`` with the row-span certificate required for Success.

Sample policy: the default scattered samples are kept as-is (including the ``ep = 3``
rank-deficient special sample at index 2); the record-selection rank-consistency filter is
the authoritative mechanism that discards bad specializations, and its counters are recorded
in the artifact.

Row-count crosscheck: when run at the recorded chamber point ``ep = -3/5`` the chamber-policy
row total must reproduce the Method.10 diagnostic artifact (which obtained the same rows via
a script-local sign swap). A mismatch is an integrity stop, not a warning.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from dataclasses import asdict, replace
from fractions import Fraction
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

SCRIPTS_DIR = Path(__file__).resolve().parent

from parametric_ibp_lf_reducer import SurfacePolicy, parse_family_text  # noqa: E402
from parametric_ibp_lf_reducer.api import default_scattered_samples  # noqa: E402
from parametric_ibp_lf_reducer.reducer import _generate_rows, reduce_rows_once  # noqa: E402
from parametric_ibp_lf_reducer.row_generation import generate_tangent_ibp_rows  # noqa: E402
from parametric_ibp_lf_reducer.sparse_rref import RREF_BACKEND_CHOICES  # noqa: E402
from parametric_ibp_lf_reducer.tangent_fields import generate_tangent_fields  # noqa: E402


def _load_script(stem: str):
    """Load a sibling runner as a module (``scripts/`` is not a package)."""
    if stem in sys.modules:
        return sys.modules[stem]
    spec = importlib.util.spec_from_file_location(stem, SCRIPTS_DIR / f"{stem}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[stem] = mod
    spec.loader.exec_module(mod)
    return mod


T2 = _load_script("run_external_int2_t2_rankrepair")

DEFAULT_EP = Fraction(-3, 5)
DEFAULT_PRIMES = (2_147_483_647, 2_147_483_629, 2_147_483_587)
DEFAULT_N_SAMPLES = 12
DEFAULT_MIN_VALID_RECORDS = 6
DEFAULT_MIN_CERTIFICATE_POINTS = 2
DEFAULT_OUT = REPO_ROOT / "validation" / "external_int2_method11.json"
DEFAULT_EXPORT = REPO_ROOT / "validation" / "external_int2_method11_reduction.txt"
METHOD10_ARTIFACT = REPO_ROOT / "validation" / "external_int2_method10.json"

SCOPE_NOTE = (
    "Method.11 Phase B: one fixed Level-0 reduction pass under the explicit exact chamber "
    "surface policy (policy is part of ReducerConfig; no monkey-patch). 'Success' claims only "
    "the row-span-certified reduction identity at the recorded (sample, prime) points, with "
    "per-term LF flags as reported. It does NOT claim the analytic Laurent value; the "
    "cross-check against the analytic Laurent oracle is a separate later phase. The package "
    "default surface policy (limit, eps -> 0^-) is unchanged."
)


def chamber_policy(ep_value: Fraction) -> SurfacePolicy:
    """Exact convergence-chamber policy at ``ep = ep_value`` (floats rejected on construction)."""
    return SurfacePolicy.chamber({"ep": ep_value})


def assemble_chamber_rows(family, labels, config) -> dict:
    """Level-0 rows (baseline + the ``(3, 3)`` extra tangent block) under ``config.surface_policy``.

    Mirrors ``T2.assemble_level_rows`` but forwards the explicit policy to the extra tangent
    block as well: T2 predates the first-class policy API, so its extra block is always
    limit-policy filtered. Baseline rows pick the policy up from ``config.surface_policy``
    inside ``_generate_rows``.
    """
    policy = config.surface_policy
    base_rows, row_diag = _generate_rows(family, labels, config, None)
    fields = generate_tangent_fields(family, [T2.BASELINE["extra_block"]])
    extra = generate_tangent_ibp_rows(family, labels, fields, policy=policy)
    rejected: dict[str, int] = {}
    for rej in extra.rejected:
        rejected[rej.reason] = rejected.get(rej.reason, 0) + 1
    merged, n_dup = T2._merge_extra_rows(base_rows, extra.rows)
    by_kind: dict[str, int] = {}
    for row in merged:
        by_kind[row.kind] = by_kind.get(row.kind, 0) + 1
    return {
        "merged": merged,
        "n_base_rows": len(base_rows),
        "n_extra_offered": len(extra.rows),
        "n_extra_new": len(merged) - len(base_rows),
        "n_extra_duplicate": n_dup,
        "n_rows_total": len(merged),
        "by_kind": by_kind,
        "row_diagnostics": row_diag,
        "richer_row_rejections": rejected,
    }


def crosscheck_method10(
    n_rows_total: int, ep_value: Fraction, artifact: Path = METHOD10_ARTIFACT
) -> dict:
    """Compare the chamber-policy row total against the recorded Method.10 artifact.

    Comparable only when the artifact exists, is a Level-0 run, and used the same chamber
    point. ``matches is False`` on a comparable artifact means the first-class policy path
    disagrees with the recorded diagnostic-swap path — an integrity stop for the caller.
    """
    out: dict = {
        "artifact": str(artifact),
        "available": False,
        "comparable": False,
        "recorded_n_chamber": None,
        "recorded_ep": None,
        "matches": None,
    }
    if not artifact.exists():
        return out
    data = json.loads(artifact.read_text(encoding="utf-8"))
    out["available"] = True
    recorded = data.get("rows", {}).get("n_chamber")
    out["recorded_n_chamber"] = recorded
    out["recorded_ep"] = data.get("chamber_ep")
    if (
        data.get("level") == 0
        and data.get("chamber_ep") == str(ep_value)
        and isinstance(recorded, int)
    ):
        out["comparable"] = True
        out["matches"] = recorded == n_rows_total
    return out


def result_payload(result) -> dict:
    """JSON-ready view of a :class:`ReductionResult` (coefficients as exported text)."""
    return {
        "status": result.status,
        "all_locally_finite": result.all_locally_finite,
        "formal_success": result.formal_success,
        "error": result.error,
        "n_terms": len(result.terms),
        "terms": [
            {
                "label": list(term.label),
                "coefficient_text": term.coefficient_text,
                "integrand_text": term.integrand_text,
                "locally_finite": term.locally_finite,
            }
            for term in result.terms
        ],
        "diagnostics": asdict(result.diagnostics),
    }


def _write_report(out_path: Path, report: dict) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, default=str) + "\n", encoding="utf-8")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", type=Path, default=T2.DEFAULT_INPUT)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--export-text",
        type=Path,
        default=DEFAULT_EXPORT,
        help="Wolfram-style reduction text, written only on certified Success",
    )
    parser.add_argument(
        "--ep",
        type=str,
        default=str(DEFAULT_EP),
        help="exact rational chamber point for the surface policy (default -3/5)",
    )
    parser.add_argument("--primes", type=str, default=",".join(str(p) for p in DEFAULT_PRIMES))
    parser.add_argument(
        "--n-samples",
        type=int,
        default=DEFAULT_N_SAMPLES,
        help="scattered samples (<= 37 keeps coordinates pairwise distinct)",
    )
    parser.add_argument(
        "--min-valid-records",
        type=int,
        default=DEFAULT_MIN_VALID_RECORDS,
        help="runner floor for Success (stricter than the package default 1); "
        "reconstruction + certificate remain the authoritative gates",
    )
    parser.add_argument(
        "--min-certificate-points",
        type=int,
        default=DEFAULT_MIN_CERTIFICATE_POINTS,
        help="informative off-sample certificate points required for Success",
    )
    parser.add_argument(
        "--jobs", type=int, default=1, help="worker processes for (prime, sample) record collection"
    )
    parser.add_argument(
        "--rref-backend",
        type=str,
        default=None,
        choices=list(RREF_BACKEND_CHOICES),
        help="sparse RREF backend (backend selection only; identical results)",
    )
    parser.add_argument(
        "--allow-heavy",
        action="store_true",
        help="required for the certified solve; without it only row assembly "
        "and the Method.10 crosscheck run",
    )
    parser.add_argument(
        "--ignore-rows-mismatch",
        action="store_true",
        help="proceed past a failed Method.10 row-count crosscheck (diagnosis only)",
    )
    args = parser.parse_args(argv)

    primes = tuple(int(p) for p in args.primes.split(",") if p.strip())
    if not primes or any(p < 3 for p in primes):
        parser.error("primes must be odd primes > 2")
    if not 1 <= args.n_samples <= 37:
        parser.error("--n-samples must be in [1, 37] (pairwise-distinct coordinate guarantee)")
    if args.min_valid_records < 1:
        parser.error("--min-valid-records must be >= 1")
    if args.min_certificate_points < 1:
        parser.error("--min-certificate-points must be >= 1")
    if args.jobs < 1:
        parser.error("--jobs must be >= 1")

    t0 = time.time()
    ep_value = Fraction(args.ep)
    policy = chamber_policy(ep_value)
    family = parse_family_text(args.input.read_text(encoding="utf-8"))
    _target, base_config = T2.build_reducer_config(family)
    labels, target_label, config = T2.build_level_config(family, base_config, T2.LEVELS[0])
    config = replace(config, surface_policy=policy)
    print(f"[m11] Level 0 box: {len(labels)} labels; chamber ep={ep_value}", flush=True)

    t_asm = time.time()
    asm = assemble_chamber_rows(family, labels, config)
    print(
        f"[m11] chamber-policy rows: {asm['n_rows_total']} {asm['by_kind']} "
        f"({time.time() - t_asm:.1f}s)",
        flush=True,
    )

    cross = crosscheck_method10(asm["n_rows_total"], ep_value)
    print(
        f"[m11] Method.10 crosscheck: comparable={cross['comparable']} "
        f"matches={cross['matches']} (recorded n_chamber={cross['recorded_n_chamber']})",
        flush=True,
    )

    report = {
        "script": "scripts/run_external_int2_method11.py",
        "method": "External Int2 Method.11: certified LF reconstruction under an explicit "
        "chamber surface policy",
        "scope_note": SCOPE_NOTE,
        "chamber_ep": str(ep_value),
        "level": 0,
        "n_labels": len(labels),
        "target_label": list(target_label),
        "rows": {k: v for k, v in asm.items() if k not in ("merged", "row_diagnostics")},
        "surface_policy": asm["row_diagnostics"].get("surface_policy"),
        "method10_crosscheck": cross,
    }

    if cross["comparable"] and cross["matches"] is False and not args.ignore_rows_mismatch:
        report["phase"] = "assembly-only"
        report["verdict"] = "Aborted(rows-mismatch)"
        report["elapsed_seconds"] = round(time.time() - t0, 3)
        _write_report(args.out, report)
        print(
            f"[m11] ABORT: chamber row total {asm['n_rows_total']} != recorded "
            f"{cross['recorded_n_chamber']}; wrote {args.out}",
            flush=True,
        )
        return 2

    if not args.allow_heavy:
        report["phase"] = "assembly-only"
        report["verdict"] = "NotRun(assembly-only)"
        report["elapsed_seconds"] = round(time.time() - t0, 3)
        _write_report(args.out, report)
        print(
            f"[m11] assembly-only (no --allow-heavy); verdict: {report['verdict']}; "
            f"wrote {args.out}",
            flush=True,
        )
        return 0

    samples = default_scattered_samples(family.parameters, args.n_samples)
    print(
        f"[m11] certified solve: {len(samples)} samples x {len(primes)} primes, "
        f"jobs={args.jobs}, rref_backend={args.rref_backend or 'default'}",
        flush=True,
    )
    t_solve = time.time()
    result = reduce_rows_once(
        family,
        target_label,
        labels,
        asm["merged"],
        primes,
        samples,
        min_valid_records=args.min_valid_records,
        require_certificate_for_success=True,
        min_certificate_points=args.min_certificate_points,
        jobs=args.jobs,
        rref_backend=args.rref_backend,
    )
    solve_seconds = round(time.time() - t_solve, 3)
    verdict = "Success(certified)" if result.success else f"Failure({result.status})"

    report["phase"] = "certified-solve"
    report["solve"] = {
        "primes": list(primes),
        "n_samples": len(samples),
        "min_valid_records": args.min_valid_records,
        "min_certificate_points": args.min_certificate_points,
        "jobs": args.jobs,
        "rref_backend": args.rref_backend,
        "solve_seconds": solve_seconds,
    }
    report["reduce"] = result_payload(result)
    report["verdict"] = verdict
    report["elapsed_seconds"] = round(time.time() - t0, 3)
    _write_report(args.out, report)

    if result.success:
        args.export_text.parent.mkdir(parents=True, exist_ok=True)
        args.export_text.write_text(result.wolfram_style_text + "\n", encoding="utf-8")
        print(f"[m11] exported reduction text: {args.export_text}", flush=True)
    print(
        f"[m11] verdict: {verdict}; n_terms={len(result.terms)}; "
        f"all_locally_finite={result.all_locally_finite}; wrote {args.out} "
        f"({report['elapsed_seconds']}s total)",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
