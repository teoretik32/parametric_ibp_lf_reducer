"""External Int2 — Method Audit.1 Phase A: all-row-support LF feasibility (diagnostic).

Domain feedback (working hypothesis): an LF-basis reduction in the ORIGINAL parametric
family is expected for this MofR 2me double-box region, so the recorded `Obstructed`
verdicts are treated as a likely implementation incompleteness until audited.

Audited defect candidate: `lf_reduction_feasible_mod_p` derives its allowed set from the
seed box only, so every row-support column OUTSIDE the seed box is silently forbidden
even when its full LF verdict is True. This runner compares, at ONE generic modular
point of the medium T2 box (Level 0), on the SAME merged rows:

  seed-only mode      : allowed = LF-True labels of the seed box (production behaviour);
  all-row-support mode: allowed = every row-support column whose full LF verdict is True
                        (target excluded) — explicit named mode, no new rows.

Nothing runs by default. ``--phase a`` is HEAVY (support LF classification + 2 Level-0
RREFs, ~15-30 min) and requires ``--allow-heavy``. Diagnostic only: the production
seed-box default, reducer core, certificates and LF gates are untouched; no Level 1/2
rerun, no re-elimination. If the all-support point is Feasible the run stops the audit
line with a prominent banner (per protocol: stop immediately once Feasible).
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

T2_SCRIPT = Path(__file__).resolve().parent / "run_external_int2_t2_rankrepair.py"


def _load_t2_module():
    """Load the Method.6 runner as a module (scripts/ is not a package)."""
    name = "run_external_int2_t2_rankrepair"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, T2_SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod  # dataclass annotation resolution needs the module registered
    spec.loader.exec_module(mod)
    return mod


T2 = _load_t2_module()

from parametric_ibp_lf_reducer.lf_feasibility import (  # noqa: E402
    STATUS_FEASIBLE,
    classify_row_support_lf,
    collect_row_support_labels,
    feasibility_to_payload,
    lf_reduction_feasible_all_support_mod_p,
    lf_reduction_feasible_mod_p,
)
from parametric_ibp_lf_reducer.valuations import is_locally_finite  # noqa: E402

LEVEL = 0  # the medium T2 box; Levels 1/2 are heavy and are NEVER re-eliminated here
N_SCATTERED = 4
DEFAULT_PRIME = 2147483647
OUT_NAME = "external_int2_audit1_allsupport.json"

SCOPE_NOTE = (
    "Method Audit.1 Phase A: same merged Level-0 rows; allowed set widened from seed-box "
    "LF-True to all-row-support LF-True at one generic modular point; diagnostic only — "
    "production seed-box default unchanged, no new rows, no re-elimination. " + T2.SCOPE_NOTE
)


def _sample_short(sample: dict) -> str:
    return ",".join(f"{k}={v}" for k, v in sample.items())


def run_phase_a(
    family,
    base_config,
    *,
    prime: int,
    sample_index: int,
    spec=None,
    out_path: Path | None = None,
) -> dict:
    """Seed-only vs all-row-support feasibility at one generic Level-0 modular point."""
    t0 = time.time()
    spec = spec if spec is not None else T2.LEVELS[LEVEL]
    labels, target_label, config = T2.build_level_config(family, base_config, spec)
    lf_map_seed = {lab: is_locally_finite(family, lab) for lab in labels}
    if target_label not in lf_map_seed:
        lf_map_seed[target_label] = is_locally_finite(family, target_label)
    rows = T2.assemble_level_rows(family, labels, config)
    merged = rows["merged"]

    support = collect_row_support_labels(merged)
    box = set(labels)
    out_of_box = [lab for lab in support if lab not in box]
    print(
        f"[a1] level {LEVEL}: {len(labels)} seed labels, {rows['n_rows_total']} rows, "
        f"{len(support)} support columns ({len(out_of_box)} out-of-box); classifying LF...",
        flush=True,
    )
    t1 = time.time()
    verdicts = classify_row_support_lf(family, merged, known_lf_flags=lf_map_seed)
    classify_seconds = round(time.time() - t1, 3)

    def _count(labs, flag) -> int:
        return sum(1 for lab in labs if verdicts[lab] == flag or verdicts[lab] is flag)

    allowed_seed = {lab for lab in labels if lf_map_seed[lab] is True and lab != target_label}
    allowed_all = {
        lab for lab in support if verdicts[lab] is True and lab != target_label
    }
    newly_allowed = sorted(allowed_all - allowed_seed)
    print(
        f"[a1] support LF verdicts ({classify_seconds}s): out-of-box "
        f"True={_count(out_of_box, True)} False={_count(out_of_box, False)} "
        f"Unknown={_count(out_of_box, 'Unknown')}; allowed seed={len(allowed_seed)} "
        f"all={len(allowed_all)} newly={len(newly_allowed)}",
        flush=True,
    )

    samples = T2._level_samples(family, spec, N_SCATTERED)
    generic = [s for i, s in enumerate(samples) if i != T2.SPECIAL_SAMPLE_INDEX]
    sample = generic[sample_index]

    points: dict[str, dict] = {}
    full_map = {**lf_map_seed, **verdicts}
    for mode in ("seed_only", "all_support"):
        t2 = time.time()
        if mode == "seed_only":
            res = lf_reduction_feasible_mod_p(
                merged, labels, target_label, lf_map_seed, sample, prime
            )
        else:
            res, _ = lf_reduction_feasible_all_support_mod_p(
                None, merged, target_label, full_map, sample, prime
            )
        payload = feasibility_to_payload(res)
        payload["solve_seconds"] = round(time.time() - t2, 3)
        points[mode] = payload
        print(
            f"[a1] {mode} {_sample_short(sample)} p={prime}: {res.status} "
            f"({payload['solve_seconds']}s)",
            flush=True,
        )

    feasible_now = points["all_support"]["status"] == STATUS_FEASIBLE
    if feasible_now:
        print(
            "[a1] *** DEFECT CONFIRMED (Phase A): all-row-support mode is FEASIBLE on the "
            "same rows where the seed-box mode reports "
            f"{points['seed_only']['status']} — the seed-box allowed set is the "
            "obstruction source. STOPPING the audit line per protocol. ***",
            flush=True,
        )

    payload = {
        "script": "run_external_int2_audit1.py",
        "method": "External Int2 Method Audit.1: all-row-support LF feasibility (Phase A)",
        "scope_note": SCOPE_NOTE,
        "level": LEVEL,
        "target": list(target_label),
        "n_labels_seed": len(labels),
        "n_rows_total": rows["n_rows_total"],
        "support": {
            "n_support_columns": len(support),
            "n_out_of_box": len(out_of_box),
            "out_of_box_lf": {
                "true": _count(out_of_box, True),
                "false": _count(out_of_box, False),
                "unknown": _count(out_of_box, "Unknown"),
            },
            "n_allowed_seed": len(allowed_seed),
            "n_allowed_all": len(allowed_all),
            "n_newly_allowed": len(newly_allowed),
            "newly_allowed_head": [list(lab) for lab in newly_allowed[:32]],
            "classify_seconds": classify_seconds,
        },
        "sample": T2._sample_repr(sample),
        "prime": int(prime),
        "points": points,
        "defect_confirmed": feasible_now,
        "elapsed_seconds": round(time.time() - t0, 3),
    }
    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    return payload


def _describe(out_dir: Path) -> None:
    print(__doc__)
    print(f"phase a -> {OUT_NAME} (HEAVY: LF classification + 2 Level-0 RREFs; --allow-heavy)")
    path = out_dir / OUT_NAME
    if path.exists():
        data = json.loads(path.read_text(encoding="utf-8"))
        sup = data.get("support", {})
        print(
            f"existing artifact: defect_confirmed={data.get('defect_confirmed')} "
            f"seed={data['points']['seed_only']['status']} "
            f"all={data['points']['all_support']['status']} "
            f"newly_allowed={sup.get('n_newly_allowed')}"
        )


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=SCOPE_NOTE)
    parser.add_argument("--input", type=Path, default=T2.DEFAULT_INPUT)
    parser.add_argument("--out-dir", type=Path, default=T2.DEFAULT_OUT_DIR)
    parser.add_argument("--describe", action="store_true", help="print plan + artifact summary")
    parser.add_argument("--phase", choices=["a"], help="nothing runs by default; 'a' is HEAVY")
    parser.add_argument("--allow-heavy", action="store_true", help="required for phase a")
    parser.add_argument("--prime", type=int, default=DEFAULT_PRIME)
    parser.add_argument("--sample-index", type=int, default=0, help="index into generic samples")
    args = parser.parse_args(argv)

    if args.describe:
        _describe(args.out_dir)
        return 0
    if args.phase is None:
        parser.error("nothing runs by default; pass --phase a (with --allow-heavy) or --describe")
    if not args.allow_heavy:
        parser.error("phase a runs 2 Level-0 RREFs + LF classification (HEAVY); pass --allow-heavy")

    text = args.input.read_text(encoding="utf-8")
    family = T2.parse_family_text(text)
    _target, base_config = T2.build_reducer_config(family)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    run_phase_a(
        family,
        base_config,
        prime=args.prime,
        sample_index=args.sample_index,
        out_path=args.out_dir / OUT_NAME,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
