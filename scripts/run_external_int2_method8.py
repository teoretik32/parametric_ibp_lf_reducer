"""External Int2 Method.8: targeted Level-0 re-elimination with the ``ibp_deg3`` rows (HEAVY).

Method.7 screening found exactly ONE breaking candidate family at the medium T2 box: degree-3
coordinate-IBP rows (6972 of 15504 new-vs-baseline rows pair nontrivially against every stable
Level-0 dual obstruction witness), with the binding caveat that a breaking row is NECESSARY,
not sufficient, to cure the obstruction. This runner performs the justified follow-up probe:
the SAME merged Level-0 rows PLUS the ``ibp_deg3`` rows, re-eliminated from scratch with the
production seed-box semantics (witness-mode RREF) at the generic scattered samples x primes.

Per point, on the enlarged system:

  ``Feasible``  -> the obstruction is CURED at that point (prominent banner; the remaining
                   points still run so stability of the cure is recorded);
  ``Witness``   -> still obstructed; the NEW dual witness must annihilate every added
                   ``ibp_deg3`` row (recorded consistency check, 0 breaks expected), and the
                   rank/nullity deltas vs the recorded Method.7 baseline are recorded.

Nothing runs by default. ``--phase rerun`` is HEAVY (each Level-0 RREF on ~62k rows,
~10-15 min per point, default 6 points) and requires ``--allow-heavy``. Diagnostic only:
reducer core, certificates, LF gates and the recorded Method.6/7 artifacts are untouched;
no Level 1/2 rerun. Honest per-(sample, prime), per-label-box statements only — neither
verdict is a global (im)possibility claim.
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

from parametric_ibp_lf_reducer.lf_obstruction_witness import (  # noqa: E402
    STATUS_FEASIBLE,
    STATUS_WITNESS,
    lf_obstruction_witness_mod_p,
    test_rows_against_obstruction_witness,
    witness_to_payload,
)
from parametric_ibp_lf_reducer.row_generation import generate_coordinate_ibp_rows  # noqa: E402
from parametric_ibp_lf_reducer.valuations import is_locally_finite  # noqa: E402

LEVEL = 0  # the medium T2 box; Levels 1/2 are heavy and are NEVER re-eliminated here
N_SCATTERED = 4
MIN_GENERIC = 3
IBP_DEGREE = 3
DEFAULT_PRIMES = (2147483647, 2147483629)
OUT_NAME = "external_int2_method8_reelim.json"
WITNESS_RECORDED_NAME = "external_int2_method7_witness.json"
SCREEN_RECORDED_NAME = "external_int2_method7_screen.json"

SCOPE_NOTE = (
    "Method.8: targeted Level-0 re-elimination with the ibp_deg3 rows added (the one family "
    "that breaks the stable Method.7 witnesses); production seed-box semantics, witness-mode "
    "RREF per generic (sample, prime). Recorded Method.6/7 artifacts are never overwritten. "
    + T2.SCOPE_NOTE
)


def _sample_short(sample: dict) -> str:
    return ",".join(f"{k}={v}" for k, v in sample.items())


def _recorded_expectations(recorded_dir: Path | None) -> dict:
    """Cross-check numbers recorded by Method.7 (baseline ranks; ibp_deg3 family counts)."""
    out: dict = {"witness_artifact": None, "screen_artifact": None}
    if recorded_dir is None:
        return out
    wpath = recorded_dir / WITNESS_RECORDED_NAME
    if wpath.exists():
        data = json.loads(wpath.read_text(encoding="utf-8"))
        out["witness_artifact"] = {
            "n_labels": data.get("n_labels"),
            "n_rows_total": data.get("n_rows_total"),
            "points": [
                {
                    "sample": p.get("sample"),
                    "prime": p.get("prime"),
                    "rank": p.get("rank"),
                    "nullity": p.get("nullity"),
                    "support_size": p.get("support_size"),
                }
                for p in data.get("points", [])
            ],
        }
    spath = recorded_dir / SCREEN_RECORDED_NAME
    if spath.exists():
        data = json.loads(spath.read_text(encoding="utf-8"))
        fam = next(
            (f for f in data.get("families", []) if f.get("family") == "candidate_ibp_deg3"),
            None,
        )
        if fam is not None:
            out["screen_artifact"] = {
                "n_candidates": fam.get("n_candidates"),
                "n_breaks_total": fam.get("n_breaks_total"),
            }
    return out


def _baseline_point(expect: dict, sample_repr: list, prime: int) -> dict | None:
    wa = expect.get("witness_artifact")
    if not wa:
        return None
    for p in wa["points"]:
        if p["prime"] == prime and [list(x) for x in p["sample"]] == sample_repr:
            return p
    return None


def run_rerun_phase(
    family,
    base_config,
    *,
    primes,
    max_points: int | None = None,
    spec=None,
    out_path: Path | None = None,
    recorded_dir: Path | None = None,
) -> dict:
    """Re-eliminate the enlarged (merged + ibp_deg3) Level-0 system at generic points."""
    t0 = time.time()
    spec = spec if spec is not None else T2.LEVELS[LEVEL]
    labels, target_label, config = T2.build_level_config(family, base_config, spec)
    lf_map = {lab: is_locally_finite(family, lab) for lab in labels}
    if target_label not in lf_map:
        lf_map[target_label] = is_locally_finite(family, target_label)
    rows = T2.assemble_level_rows(family, labels, config)
    merged = rows["merged"]
    baseline_keys = {row.dedup_key() for row in merged}

    t1 = time.time()
    gen3 = generate_coordinate_ibp_rows(family, labels, IBP_DEGREE)
    rejected: dict[str, int] = {}
    for rej in gen3.rejected:
        rejected[rej.reason] = rejected.get(rej.reason, 0) + 1
    new3 = [r for r in gen3.rows if r.dedup_key() not in baseline_keys]
    generation_seconds = round(time.time() - t1, 3)
    enlarged = list(merged) + new3

    expect = _recorded_expectations(recorded_dir)
    screen_expect = expect.get("screen_artifact")
    generation_consistent = (
        None if screen_expect is None else screen_expect["n_candidates"] == len(new3)
    )
    print(
        f"[m8] level {LEVEL}: {len(labels)} seed labels, {rows['n_rows_total']} baseline rows "
        f"+ {len(new3)} ibp_deg{IBP_DEGREE} rows (generated {len(gen3.rows)}, "
        f"{generation_seconds}s) -> {len(enlarged)} total"
        + ("" if generation_consistent is None else f"; matches screen artifact: {generation_consistent}"),
        flush=True,
    )

    samples = T2._level_samples(family, spec, N_SCATTERED)
    generic = [s for i, s in enumerate(samples) if i != T2.SPECIAL_SAMPLE_INDEX]
    excluded = samples[T2.SPECIAL_SAMPLE_INDEX]
    if len(generic) < MIN_GENERIC:
        raise ValueError(f"need >= {MIN_GENERIC} generic samples, got {len(generic)}")

    pairs = [(s, p) for s in generic for p in primes]
    if max_points is not None:
        pairs = pairs[:max_points]

    points: list[dict] = []
    cured: list[str] = []
    for sample, prime in pairs:
        t2 = time.time()
        res = lf_obstruction_witness_mod_p(enlarged, labels, target_label, lf_map, sample, prime)
        payload = witness_to_payload(res)
        payload["solve_seconds"] = round(time.time() - t2, 3)
        payload["support_size"] = len(res.witness)
        base = _baseline_point(expect, payload["sample"], int(prime))
        if base is not None and res.status == STATUS_WITNESS:
            payload["baseline_rank"] = base["rank"]
            payload["baseline_nullity"] = base["nullity"]
            payload["delta_rank"] = res.rank - base["rank"]
            payload["delta_nullity"] = res.nullity - base["nullity"]
        if res.status == STATUS_WITNESS:
            pairings = test_rows_against_obstruction_witness(new3, res, sample, prime)
            n_breaks = sum(1 for p in pairings if p.breaks)
            payload["included_ibp3_breaks_new_witness"] = n_breaks
            payload["included_rows_annihilated"] = n_breaks == 0
        points.append(payload)
        tag = f"{_sample_short(sample)} p={prime}"
        print(
            f"[m8] rerun {tag}: {res.status} rank={res.rank} nullity={res.nullity} "
            f"support={len(res.witness)} ({payload['solve_seconds']}s)",
            flush=True,
        )
        if res.status == STATUS_FEASIBLE:
            cured.append(tag)
            print(
                f"[m8] *** OBSTRUCTION CURED at {tag}: the enlarged (merged + ibp_deg3) "
                "Level-0 system reduces the target through LF-True labels. Remaining points "
                "continue for stability. ***",
                flush=True,
            )

    statuses = [p["status"] for p in points]
    n_feasible = statuses.count(STATUS_FEASIBLE)
    n_witness = statuses.count(STATUS_WITNESS)
    if n_feasible == len(points):
        outcome = "cured_all_points"
    elif n_feasible > 0:
        outcome = "mixed"
    elif n_witness == len(points):
        outcome = "still_obstructed_all_points"
    else:
        outcome = "inconclusive"
    decision = {
        "outcome": outcome,
        "n_points": len(points),
        "n_feasible": n_feasible,
        "n_witness": n_witness,
        "cured_points": cured,
        "note": (
            "Feasible here means the target reduces through LF-True labels of THIS box at that "
            "(sample, prime); Witness means the enlarged system is still obstructed with a new "
            "recorded dual witness that annihilates every included ibp_deg3 row. Neither is a "
            "global claim."
        ),
    }

    payload = {
        "script": "run_external_int2_method8.py",
        "method": "External Int2 Method.8: targeted Level-0 re-elimination with ibp_deg3 rows",
        "scope_note": SCOPE_NOTE,
        "level": LEVEL,
        "target": list(target_label),
        "n_labels": len(labels),
        "n_rows_baseline": rows["n_rows_total"],
        "ibp_deg3": {
            "degree": IBP_DEGREE,
            "n_generated": len(gen3.rows),
            "n_new_vs_baseline": len(new3),
            "surface_rejected": rejected,
            "generation_seconds": generation_seconds,
            "matches_screen_artifact": generation_consistent,
        },
        "n_rows_enlarged": len(enlarged),
        "samples": [T2._sample_repr(s) for s, _p in pairs],
        "excluded_special_sample": T2._sample_repr(excluded),
        "primes": [int(p) for p in primes],
        "recorded_baseline": expect,
        "points": points,
        "decision": decision,
        "elapsed_seconds": round(time.time() - t0, 3),
    }
    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    return payload


def _describe(out_dir: Path) -> None:
    print(__doc__)
    print(
        f"phase rerun -> {OUT_NAME} (HEAVY: Level-0 witness-mode RREFs on the enlarged "
        "system, default 6 points; --allow-heavy)"
    )
    path = out_dir / OUT_NAME
    if path.exists():
        data = json.loads(path.read_text(encoding="utf-8"))
        dec = data.get("decision", {})
        print(
            f"existing artifact: outcome={dec.get('outcome')} "
            f"n_feasible={dec.get('n_feasible')}/{dec.get('n_points')} "
            f"enlarged_rows={data.get('n_rows_enlarged')}"
        )


def _parse_primes(text: str) -> tuple[int, ...]:
    primes = tuple(int(tok) for tok in text.split(",") if tok.strip())
    if not primes or any(p < 3 for p in primes):
        raise ValueError(f"bad primes: {text!r}")
    return primes


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=SCOPE_NOTE)
    parser.add_argument("--input", type=Path, default=T2.DEFAULT_INPUT)
    parser.add_argument("--out-dir", type=Path, default=T2.DEFAULT_OUT_DIR)
    parser.add_argument("--describe", action="store_true", help="print plan + artifact summary")
    parser.add_argument("--phase", choices=["rerun"], help="nothing runs by default")
    parser.add_argument("--allow-heavy", action="store_true", help="required for phase rerun")
    parser.add_argument(
        "--primes", default=",".join(str(p) for p in DEFAULT_PRIMES), help="comma-separated"
    )
    parser.add_argument(
        "--max-points", type=int, default=None, help="cap on (sample, prime) points (default all)"
    )
    args = parser.parse_args(argv)

    if args.describe:
        _describe(args.out_dir)
        return 0
    if args.phase is None:
        print(
            "error: --phase is required (nothing runs by default); see --describe",
            file=sys.stderr,
        )
        return 2
    if not args.allow_heavy:
        print(
            "error: phase rerun is HEAVY (Level-0 RREFs on ~62k rows); pass --allow-heavy",
            file=sys.stderr,
        )
        return 2
    try:
        primes = _parse_primes(args.primes)
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    if args.max_points is not None and args.max_points < 1:
        print("error: --max-points must be >= 1", file=sys.stderr)
        return 2

    text = args.input.read_text(encoding="utf-8")
    family = T2.parse_family_text(text)
    _target, base_config = T2.build_reducer_config(family)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    run_rerun_phase(
        family,
        base_config,
        primes=primes,
        max_points=args.max_points,
        out_path=args.out_dir / OUT_NAME,
        recorded_dir=args.out_dir,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
