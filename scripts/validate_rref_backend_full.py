# Perf.13: full-pipeline validation of the RREF backend selection (dict / numba / auto).
#
# Runs the COMPLETE reduction pipeline from run_example4_star_corrected.py
# (parse -> labels -> LF flags -> rows -> reduce_rows_multi incl. shared-point
# certification) once per requested backend on the real corrected Example 4*
# family, then proves the certified mathematical result is bit-identical across
# backends and reports honest end-to-end wall times.
#
# Default mode is "medium": every label-box range is shaved by one unit off the
# top (same convention as profile_rref_real_matrix.py) so the run finishes in
# minutes while staying a real matrix/real ranking/real certificate. Pass
# --full for the untouched 972-label box (tens of minutes per backend on this
# machine). The mode is recorded in the JSON so numbers are never mislabeled.
#
# What is compared per target (must be EXACTLY equal across backends):
#   - status / formal_success / error / all_locally_finite
#   - every term: (label, coefficient_text, integrand_text, locally_finite)
#   - selected_rank
#   - the full certificate payload and record_selection diagnostics
#     (minus keys matching *time*/*seconds* — timings are reported, not compared)
#
# Rows/labels/LF flags are generated ONCE (they cannot depend on the RREF
# backend); only reduce_rows_multi is repeated per backend, with a fresh
# certificate cache each time so nothing leaks between runs.
#
# Read-only: no repo artifacts are written unless --out is given. No stdout
# except progress lines on stderr and the final JSON on stdout. Deterministic
# apart from wall-clock timings.
from __future__ import annotations

import argparse
import dataclasses
import json
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:  # allow running without installation
    sys.path.insert(0, str(REPO_ROOT / "src"))

from parametric_ibp_lf_reducer import parse_family_text  # noqa: E402
from parametric_ibp_lf_reducer.api import build_reducer_config  # noqa: E402
from parametric_ibp_lf_reducer.modular_normal_form import assemble_rows_mod_p  # noqa: E402
from parametric_ibp_lf_reducer.reducer import (  # noqa: E402
    _enumerate_labels,
    _generate_rows,
    reduce_rows_multi,
)
from parametric_ibp_lf_reducer.sparse_rref import (  # noqa: E402
    AUTO_RREF_BACKEND,
    DEFAULT_RREF_BACKEND,
    NUMBA_RREF_BACKEND,
    RREF_BACKEND_CHOICES,
    _numba_available,
    rref_backend_available,
    rref_mod_p,
    select_rref_backend,
)
from parametric_ibp_lf_reducer.valuations import is_locally_finite  # noqa: E402

# Sibling script module (scripts/ is on sys.path when run as a script): the LHS
# decomposition helper is deliberately script-local, not package API.
from run_example4_star_corrected import lhs_terms_from_document  # noqa: E402

INPUT_PATH = REPO_ROOT / "examples" / "example4_star_corrected_input.wl.txt"

_VOLATILE_KEY_FRAGMENTS = ("time", "seconds")


def _log(msg: str) -> None:
    print(msg, file=sys.stderr, flush=True)


def _shrink(rng):
    """Shave one unit off the top of every (lo, hi) range; handles nested range lists."""
    if isinstance(rng, (list, tuple)) and len(rng) == 2 and all(isinstance(x, int) for x in rng):
        lo, hi = rng
        return (lo, max(lo, hi - 1))
    return tuple(_shrink(r) for r in rng)


def _strip_volatile(obj):
    """Drop dict keys whose name looks like a timing; recurse into containers."""
    if isinstance(obj, dict):
        return {
            k: _strip_volatile(v)
            for k, v in obj.items()
            if not any(frag in str(k).lower() for frag in _VOLATILE_KEY_FRAGMENTS)
        }
    if isinstance(obj, (list, tuple)):
        return [_strip_volatile(v) for v in obj]
    return obj


def _label_text(label) -> str:
    return json.dumps(label, default=str)


def _snapshot(result) -> dict:
    """Backend-independent view of one ReductionResult (everything that is claimed)."""
    extra = result.diagnostics.extra
    return {
        "status": result.status,
        "formal_success": bool(result.formal_success),
        "error": result.error,
        "all_locally_finite": str(result.all_locally_finite),
        "terms": [
            {
                "label": _label_text(t.label),
                "coefficient_text": t.coefficient_text,
                "integrand_text": t.integrand_text,
                "locally_finite": str(t.locally_finite),
            }
            for t in result.terms
        ],
        "certificate": _strip_volatile(extra.get("certificate") or {}),
        "record_selection": _strip_volatile(extra.get("record_selection") or {}),
    }


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Perf.13: full-pipeline RREF backend validation (dict/numba/auto)"
    )
    ap.add_argument("--full", action="store_true", help="use the untouched label box (slow)")
    default_backends = [DEFAULT_RREF_BACKEND]
    if rref_backend_available(NUMBA_RREF_BACKEND):
        default_backends.append(NUMBA_RREF_BACKEND)
    default_backends.append(AUTO_RREF_BACKEND)
    ap.add_argument(
        "--backends",
        default=",".join(default_backends),
        help="comma-separated backend list (default: dict, numba-if-available, auto)",
    )
    ap.add_argument("--out", default=None, help="also write the JSON to this path")
    args = ap.parse_args()
    backends = [b.strip() for b in args.backends.split(",") if b.strip()]
    for b in backends:
        if b not in RREF_BACKEND_CHOICES:
            ap.error(f"unknown backend {b!r}; choices: {RREF_BACKEND_CHOICES}")
        if not rref_backend_available(b):
            ap.error(f"backend {b!r} is not available in this environment")

    if not INPUT_PATH.is_file():
        _log(f"input not found: {INPUT_PATH}")
        return 2
    family = parse_family_text(INPUT_PATH.read_text(encoding="utf-8"))
    _, config = build_reducer_config(family)
    if config.label_box is None:
        raise SystemExit("document did not configure a label box")
    lhs = lhs_terms_from_document(family)
    targets = sorted(lhs)

    mode = "full-real" if args.full else "medium-real (box shaved by 1 per range)"
    n_range, m_range = config.label_box
    if not args.full:
        n_range, m_range = _shrink(n_range), _shrink(m_range)
    cfg = dataclasses.replace(config, labels=None, label_box=(n_range, m_range))

    # Shared, backend-independent preparation (rows cannot depend on the backend).
    prep: dict[str, float] = {}
    t0 = time.perf_counter()
    labels = _enumerate_labels(family, cfg)
    prep["labels_s"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    lf_map = {lab: is_locally_finite(family, lab) for lab in labels}
    prep["lf_flags_s"] = time.perf_counter() - t0
    t0 = time.perf_counter()
    rows, row_diag = _generate_rows(family, labels, cfg)
    prep["row_generation_s"] = time.perf_counter() - t0
    _log(f"prepared: {len(labels)} labels, {len(rows)} rows, {len(targets)} targets ({mode})")

    # What would "auto" pick on the first record point? (reporting only)
    prime0 = cfg.primes[0]
    sample0 = dict(cfg.samples[0])
    matrix0 = assemble_rows_mod_p(family, rows, sample0, prime0)
    nnz0 = sum(len(r) for r in matrix0)
    cols0 = len({c for r in matrix0 for c in r})
    auto_choice, auto_reason = select_rref_backend(
        AUTO_RREF_BACKEND,
        n_rows=len(matrix0),
        n_cols=cols0,
        initial_nnz=nnz0,
        prime=prime0,
        numba_available=_numba_available(),
    )
    _log(f"auto on first point would pick: {auto_choice} ({auto_reason})")
    del matrix0

    per_backend: dict[str, dict] = {}
    snapshots: dict[str, dict] = {}
    for backend in backends:
        warmup_s = 0.0
        if backend in (NUMBA_RREF_BACKEND, AUTO_RREF_BACKEND) and _numba_available():
            # One-time JIT compile, reported separately so pipeline walls stay comparable.
            t0 = time.perf_counter()
            rref_mod_p(
                [{0: 1, 1: 2}, {1: 3}], prime0, column_order=[0, 1], backend=NUMBA_RREF_BACKEND
            )
            warmup_s = time.perf_counter() - t0
        _log(f"running full pipeline with rref_backend={backend!r} ...")
        cert_cache: dict = {}
        t0 = time.perf_counter()
        results = reduce_rows_multi(
            family,
            targets,
            labels,
            rows,
            cfg.primes,
            cfg.samples,
            lf_flags=lf_map,
            preferred_masters=cfg.preferred_masters,
            min_valid_records=cfg.min_valid_records,
            certificate_rref_cache=cert_cache,
            rref_backend=backend,
        )
        wall = time.perf_counter() - t0
        shared_stages = {
            k: round(float(v), 3)
            for k, v in (
                next(iter(results.values())).diagnostics.extra.get("timings") or {}
            ).items()
        }
        statuses = {_label_text(t): r.status for t, r in results.items()}
        _log(f"  {backend}: wall={wall:.1f}s statuses={sorted(set(statuses.values()))}")
        per_backend[backend] = {
            "wall_s": round(wall, 3),
            "numba_jit_warmup_s": round(warmup_s, 3),
            "shared_stage_timings_s": shared_stages,
            "statuses": statuses,
        }
        snapshots[backend] = {_label_text(t): _snapshot(r) for t, r in results.items()}

    ref_backend = backends[0]
    mismatches: list[str] = []
    for backend in backends[1:]:
        if snapshots[backend] != snapshots[ref_backend]:
            for tgt, snap in snapshots[backend].items():
                if snap != snapshots[ref_backend].get(tgt):
                    mismatches.append(f"{backend} vs {ref_backend}: target {tgt}")

    payload = {
        "mode": mode,
        "input": str(INPUT_PATH.relative_to(REPO_ROOT)),
        "label_box_used": [n_range, m_range],
        "n_labels": len(labels),
        "n_rows": len(rows),
        "rows_by_kind": row_diag.get("by_kind") if isinstance(row_diag, dict) else row_diag,
        "targets": [_label_text(t) for t in targets],
        "prep_timings_s": {k: round(v, 3) for k, v in prep.items()},
        "numba_available": _numba_available(),
        "auto_first_point_choice": {"backend": auto_choice, "reason": auto_reason},
        "backends": per_backend,
        "reference_backend": ref_backend,
        "identical_across_backends": not mismatches,
        "mismatches": mismatches,
        "reference_snapshot": snapshots[ref_backend],
    }
    text = json.dumps(payload, indent=2, default=str)
    print(text)
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(text + "\n", encoding="utf-8")
    if mismatches:
        _log("BACKEND MISMATCH — certified results differ; see 'mismatches' in the JSON")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
