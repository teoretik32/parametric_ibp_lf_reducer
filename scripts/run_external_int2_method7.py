#!/usr/bin/env python
"""External Int2 Method.7: dual-witness stability + pairing-only candidate-row screening.

Builds on Method.6 (``run_external_int2_t2_rankrepair.py``): instead of enlarging label boxes and
re-running heavy eliminations, characterize the Level-0 (medium box) obstruction through the
STABILITY of its dual LF-obstruction witness, then screen candidate row families by witness
pairing only.

Phase A (``--phase witness``; HEAVY — 6 Level-0 witness RREFs, ~1h; requires ``--allow-heavy``):
  witnesses at >=3 generic scattered samples x >=2 large primes on the Level-0 merged system.
  Records per-point rank/nullity/support + exact checks and a cross-point stability block
  (canonical support hashes, pairwise Jaccard, common core, rank/nullity constancy).
  Artifact: ``validation/external_int2_method7_witness.json``.

Phase B/C (``--phase screen``; cheap — NO RREF, never re-eliminates):
  pairs row families against every stored witness at its own ``(sample, prime)``:
  - existing rows (baseline algebraic / coordinate IBP deg<=2 / tangent (1,1)+(2,2)+(3,3), plus
    the Method.4/5/T2 tangent-(4,4) probe family) — expected to annihilate (witness exactness);
  - candidate families: coordinate IBP degree 3; ray-directed single-variable multipliers
    ``x_j^d`` (d=4..6, primary failing-ray variable first, see notes/EXTERNAL_INT2_AUDIT.md);
    tangent (5,5) under an explicit wall-clock budget (overruns recorded as ``aborted_budget``,
    never silent). Surface-rejected rows are counted as diagnostics only and can never enter a
    pairing family (``RejectedRow`` carries no terms).
  Artifact: ``validation/external_int2_method7_screen.json``.

Decision: any breaking candidate row -> targeted enlargement recommendation; none -> the
obstruction is witness-certified stable under these enlargements, with the standing caveat that a
breaking row is NECESSARY, not sufficient, to cure the obstruction.

Scope guard: reducer core, certificates and LF gates are untouched; recorded Method.6 artifacts
are read-only; no Level 1/2 RREF is ever run here. Nothing runs by default (``--phase`` is
required; ``--describe`` prints the plan).
"""

from __future__ import annotations

import argparse
import hashlib
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
    STATUS_WITNESS,
    lf_obstruction_witness_mod_p,
    test_rows_against_obstruction_witness,
    witness_from_payload,
    witness_to_payload,
)
from parametric_ibp_lf_reducer.row_generation import (  # noqa: E402
    Row,
    coordinate_ibp_primitive_row,
    generate_coordinate_ibp_rows,
    generate_tangent_ibp_rows,
)
from parametric_ibp_lf_reducer.surface import coordinate_primitive_surface_free  # noqa: E402
from parametric_ibp_lf_reducer.tangent_fields import generate_tangent_fields  # noqa: E402
from parametric_ibp_lf_reducer.valuations import is_locally_finite  # noqa: E402

LEVEL = 0  # the medium T2 box; Levels 1/2 are heavy and are NEVER re-eliminated here
N_SCATTERED = 4
MIN_GENERIC = 3
RAY_DEGREES = (4, 5, 6)
#: Primary failing ray ``(-1,0,0)`` (``x2 -> oo``) from notes/EXTERNAL_INT2_AUDIT.md: its
#: variable is probed first among the single-variable ("ray-directed") multiplier supports.
PRIMARY_RAY_VAR = 0
TANGENT55_BLOCK = (5, 5)
TANGENT55_CHUNK = 256
DEFAULT_TANGENT55_BUDGET = 1800.0
DEFAULT_MAX_CANDIDATES = 200_000

WITNESS_OUT_NAME = "external_int2_method7_witness.json"
SCREEN_OUT_NAME = "external_int2_method7_screen.json"

RERUN_CAVEAT = (
    "a breaking row is NECESSARY, not sufficient, to cure the obstruction "
    "(other nullvectors may still obstruct); no re-elimination is run here"
)

SCOPE_NOTE = (
    "Method.7: dual-witness stability + pairing-only candidate-row screening at the medium T2 "
    "box (Level 0); no Level 1/2 RREF rerun, no re-elimination. " + T2.SCOPE_NOTE
)


# --- canonical hashing / stability -------------------------------------------------------------


def support_hash(support) -> str:
    """sha256 of the canonically sorted label support (order-independent, deterministic)."""
    canon = json.dumps(sorted([int(x) for x in lab] for lab in support), separators=(",", ":"))
    return hashlib.sha256(canon.encode("ascii")).hexdigest()


def vector_hash(witness_entries) -> str:
    """sha256 of the canonically sorted ``(label, coeff)`` witness vector."""
    canon = json.dumps(
        sorted([[int(x) for x in lab], int(coeff)] for lab, coeff in witness_entries),
        separators=(",", ":"),
    )
    return hashlib.sha256(canon.encode("ascii")).hexdigest()


def jaccard(a: set, b: set) -> float:
    """Jaccard similarity; 1.0 for two empty sets."""
    if not a and not b:
        return 1.0
    return len(a & b) / len(a | b)


def stability_block(point_payloads: list[dict]) -> dict:
    """Cross-point stability report over witness payloads (Phase A core statement)."""
    witness_points = [p for p in point_payloads if p["status"] == STATUS_WITNESS]
    supports = [frozenset(tuple(int(x) for x in lab) for lab, _c in p["witness"]) for p in witness_points]
    patterns = set(supports)
    n = len(witness_points)
    jac = [
        [round(jaccard(set(supports[i]), set(supports[j])), 6) for j in range(n)] for i in range(n)
    ]
    common = set.intersection(*[set(s) for s in supports]) if supports else set()
    union = set.union(*[set(s) for s in supports]) if supports else set()
    ranks = sorted({p["rank"] for p in witness_points})
    nullities = sorted({p["nullity"] for p in witness_points})
    return {
        "n_points": len(point_payloads),
        "n_witness_points": n,
        "statuses": sorted({p["status"] for p in point_payloads}),
        "all_points_witness": n == len(point_payloads) and n > 0,
        "ranks": ranks,
        "nullities": nullities,
        "rank_constant": len(ranks) <= 1,
        "nullity_constant": len(nullities) <= 1,
        "support_sizes": [len(s) for s in supports],
        "n_distinct_support_patterns": len(patterns),
        "support_pattern_stable": len(patterns) <= 1 and n > 0,
        "pairwise_support_jaccard": jac,
        "common_core_size": len(common),
        "union_support_size": len(union),
        "common_core_hash": support_hash(common),
        "all_checks_pass": all(
            p["check_annihilation"] and p["check_target_unit"] for p in witness_points
        )
        and n > 0,
    }


# --- Phase A: witness stability ----------------------------------------------------------------


def _recorded_consistency(out_dir: Path, n_labels: int, n_rows_total: int, generic, excluded) -> dict:
    """RREF-free cross-check against the recorded Level-0 Method.6 artifact (if present)."""
    path = out_dir / T2.RECORDED_NAME.format(level=LEVEL)
    if not path.exists():
        return {"recorded_artifact": None}
    data = json.loads(path.read_text(encoding="utf-8"))
    rec_obstructed = sorted(
        tuple(map(tuple, p["sample"])) for p in data.get("points", []) if p.get("status") == "Obstructed"
    )
    rec_feasible = {
        tuple(map(tuple, p["sample"])) for p in data.get("points", []) if p.get("status") == "Feasible"
    }
    ours = sorted(tuple(map(tuple, T2._sample_repr(s))) for s in generic)
    return {
        "recorded_artifact": path.name,
        "n_labels_match": data.get("n_labels") == n_labels,
        "n_rows_total_match": data.get("n_rows_total") == n_rows_total,
        "generic_samples_match_recorded_obstructed": ours == rec_obstructed,
        "excluded_sample_recorded_feasible": tuple(map(tuple, T2._sample_repr(excluded)))
        in rec_feasible,
    }


def _sample_short(sample: dict) -> str:
    return ",".join(f"{k}={v}" for k, v in sorted(sample.items()))


def run_witness_phase(
    family,
    base_config,
    *,
    primes,
    spec=None,
    out_path: Path | None = None,
    recorded_dir: Path | None = None,
) -> dict:
    """Phase A: extract witnesses at all generic scattered samples x ``primes`` (Level-0 box)."""
    t0 = time.time()
    spec = spec if spec is not None else T2.LEVELS[LEVEL]
    labels, target_label, config = T2.build_level_config(family, base_config, spec)
    lf_map = {lab: is_locally_finite(family, lab) for lab in labels}
    if target_label not in lf_map:
        lf_map[target_label] = is_locally_finite(family, target_label)
    rows = T2.assemble_level_rows(family, labels, config)
    merged = rows["merged"]

    samples = T2._level_samples(family, spec, N_SCATTERED)
    generic = [s for i, s in enumerate(samples) if i != T2.SPECIAL_SAMPLE_INDEX]
    excluded = samples[T2.SPECIAL_SAMPLE_INDEX]
    if len(generic) < MIN_GENERIC:
        raise ValueError(f"need >= {MIN_GENERIC} generic samples, got {len(generic)}")

    consistency = (
        _recorded_consistency(recorded_dir, len(labels), rows["n_rows_total"], generic, excluded)
        if recorded_dir is not None
        else None
    )

    points: list[dict] = []
    for sample in generic:
        for prime in primes:
            t1 = time.time()
            res = lf_obstruction_witness_mod_p(merged, labels, target_label, lf_map, sample, prime)
            payload = witness_to_payload(res)
            payload["solve_seconds"] = round(time.time() - t1, 3)
            payload["support_size"] = len(res.witness)
            payload["support_hash"] = support_hash([lab for lab, _c in payload["witness"]])
            payload["vector_hash"] = vector_hash(payload["witness"])
            points.append(payload)
            print(
                f"[m7] witness {_sample_short(sample)} p={prime}: {res.status} "
                f"rank={res.rank} nullity={res.nullity} support={len(res.witness)} "
                f"({payload['solve_seconds']}s)",
                flush=True,
            )

    payload = {
        "script": "run_external_int2_method7.py",
        "method": "External Int2 Method.7: dual-witness stability (Phase A)",
        "scope_note": SCOPE_NOTE,
        "level": LEVEL,
        "label_box": [
            [list(r) for r in spec.n_ranges],
            [list(spec.m_range) for _ in range(family.npolys)],
        ],
        "baseline": {
            "max_ibp_degree": T2.BASELINE["max_ibp_degree"],
            "tangent_degree_blocks": [list(b) for b in T2.BASELINE["tangent_degree_blocks"]],
            "extra_block": list(T2.BASELINE["extra_block"]),
        },
        "target": list(target_label),
        "n_labels": len(labels),
        "n_lf_true": sum(1 for lab in labels if lf_map[lab] is True),
        "n_rows_total": rows["n_rows_total"],
        "samples": [T2._sample_repr(s) for s in generic],
        "excluded_special_sample": T2._sample_repr(excluded),
        "primes": [int(p) for p in primes],
        "points": points,
        "stability": stability_block(points),
        "recorded_consistency": consistency,
        "elapsed_seconds": round(time.time() - t0, 3),
    }
    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    return payload


# --- Phase C candidate generators (pairing-only; NEVER re-eliminates) --------------------------


def generate_ray_multiplier_rows(
    family, labels, *, degrees=RAY_DEGREES, max_candidates: int | None = None
):
    """Coordinate IBP rows for single-variable ("ray-directed") monomial multipliers ``x_j^d``.

    Same surface filter chain as :func:`generate_coordinate_ibp_rows` (surface check ->
    triviality -> dedup); the primary failing-ray variable is probed first. Returns
    ``(rows, rejected_counts, truncated_by_cap)``; rejected rows are counted only and can never
    be paired.
    """
    var_order = [PRIMARY_RAY_VAR] + [j for j in range(family.nvars) if j != PRIMARY_RAY_VAR]
    rows_out: list[Row] = []
    rejected: dict[str, int] = {}
    seen: set[frozenset] = set()
    truncated = False
    for label in labels:
        for var_index in range(family.nvars):
            for j in var_order:
                for d in degrees:
                    exps = tuple(d if k == j else 0 for k in range(family.nvars))
                    verdict = coordinate_primitive_surface_free(
                        family, label, var_index, multiplier_exps=exps, eps_direction="minus"
                    )
                    if verdict is not True:
                        reason = "surface_unknown" if verdict == "Unknown" else "surface_not_free"
                        rejected[reason] = rejected.get(reason, 0) + 1
                        continue
                    row = coordinate_ibp_primitive_row(family, label, var_index, exps)
                    if row.is_trivial():
                        rejected["trivial_row"] = rejected.get("trivial_row", 0) + 1
                        continue
                    key = row.dedup_key()
                    if key in seen:
                        continue
                    seen.add(key)
                    rows_out.append(row)
                    if max_candidates is not None and len(rows_out) >= max_candidates:
                        truncated = True
                        return rows_out, rejected, truncated
    return rows_out, rejected, truncated


def generate_tangent55_rows(
    family,
    labels,
    *,
    budget_seconds: float,
    max_candidates: int,
    block=TANGENT55_BLOCK,
    chunk: int = TANGENT55_CHUNK,
):
    """Tangent-``block`` rows under a strict wall-clock budget (checked between label chunks).

    Overruns are recorded as ``status="aborted_budget"`` with partial statistics — never silent.
    Returns ``(rows, meta)``.
    """
    t0 = time.time()
    meta: dict = {
        "block": list(block),
        "status": "complete",
        "budget_seconds": budget_seconds,
        "n_fields": 0,
        "n_labels_total": len(labels),
        "n_labels_processed": 0,
        "rejected": {},
        "truncated_by_cap": False,
    }
    fields = generate_tangent_fields(family, [block])
    meta["n_fields"] = len(fields)
    meta["fields_seconds"] = round(time.time() - t0, 3)

    rows_out: list[Row] = []
    seen: set[frozenset] = set()
    rejected: dict[str, int] = {}
    if time.time() - t0 > budget_seconds:
        meta["status"] = "aborted_budget"
        meta["elapsed_seconds"] = round(time.time() - t0, 3)
        return rows_out, meta

    for start in range(0, len(labels), chunk):
        gen = generate_tangent_ibp_rows(family, labels[start : start + chunk], fields)
        for rej in gen.rejected:
            rejected[rej.reason] = rejected.get(rej.reason, 0) + 1
        for row in gen.rows:
            key = row.dedup_key()
            if key in seen:
                continue
            seen.add(key)
            rows_out.append(row)
            if len(rows_out) >= max_candidates:
                meta["truncated_by_cap"] = True
                break
        meta["n_labels_processed"] = min(start + chunk, len(labels))
        if meta["truncated_by_cap"]:
            break
        if time.time() - t0 > budget_seconds and meta["n_labels_processed"] < len(labels):
            meta["status"] = "aborted_budget"
            break
    meta["rejected"] = rejected
    meta["elapsed_seconds"] = round(time.time() - t0, 3)
    return rows_out, meta


# --- screening (pairing only) ------------------------------------------------------------------


def _row_prov_repr(row: Row) -> dict:
    """JSON-safe provenance of one (breaking) row."""
    return {"kind": row.kind, "provenance": {k: str(v) for k, v in row.provenance.items()}}


def screen_family(name: str, category: str, rows, witness_payloads, *, extra: dict | None = None) -> dict:
    """Pair one row family against every stored witness at its own ``(sample, prime)``."""
    out: dict = {
        "family": name,
        "category": category,
        "n_candidates": len(rows),
        "n_breaks_total": 0,
        "per_witness": [],
    }
    for wp in witness_payloads:
        witness = witness_from_payload(wp)
        sample = T2._sample_from_repr([list(pair) for pair in witness.sample])
        pairings = test_rows_against_obstruction_witness(rows, witness, sample, witness.prime)
        breaking = [p for p in pairings if p.breaks]
        prov: dict[str, int] = {}
        for p in breaking:
            prov[p.kind] = prov.get(p.kind, 0) + 1
        out["per_witness"].append(
            {
                "witness_sample": [list(pair) for pair in witness.sample],
                "witness_prime": witness.prime,
                "n_breaks": len(breaking),
                "n_annihilate": len(pairings) - len(breaking),
                "breaking_kinds": prov,
                "breaking_examples": [_row_prov_repr(rows[p.index]) for p in breaking[:8]],
            }
        )
        out["n_breaks_total"] += len(breaking)
    if extra:
        out.update(extra)
    return out


_EXISTING_KIND_NAMES = {
    "algebraic": "existing_algebraic",
    "coordinate_ibp": "existing_coordinate_ibp_deg<=2",
    "tangent_ibp": "existing_tangent_(1,1)+(2,2)+(3,3)",
}


def run_screen_phase(
    family,
    base_config,
    witness_payloads,
    *,
    spec=None,
    out_path: Path | None = None,
    tangent55_budget: float = DEFAULT_TANGENT55_BUDGET,
    max_candidates: int = DEFAULT_MAX_CANDIDATES,
    witness_stability: dict | None = None,
) -> dict:
    """Phase B/C: pairing-only screening of existing + candidate row families (NO RREF)."""
    t0 = time.time()
    spec = spec if spec is not None else T2.LEVELS[LEVEL]
    labels, target_label, config = T2.build_level_config(family, base_config, spec)

    witness_payloads = [wp for wp in witness_payloads if wp["status"] == STATUS_WITNESS]
    if not witness_payloads:
        raise ValueError("no Witness-status points to screen against")
    for wp in witness_payloads:
        if tuple(wp["target_label"]) != target_label:
            raise ValueError(
                f"witness target {wp['target_label']} != level target {list(target_label)}"
            )

    rows = T2.assemble_level_rows(family, labels, config)
    merged = rows["merged"]
    baseline_keys = {row.dedup_key() for row in merged}
    print(f"[m7] screen: {len(merged)} baseline rows, {len(witness_payloads)} witnesses", flush=True)

    families: list[dict] = []

    # --- Phase B: rows already in the Method.4/5/T2 systems -----------------------------------
    by_kind: dict[str, list[Row]] = {}
    for row in merged:
        by_kind.setdefault(row.kind, []).append(row)
    for kind in sorted(by_kind):
        name = _EXISTING_KIND_NAMES.get(kind, f"existing_{kind}")
        families.append(screen_family(name, "existing", by_kind[kind], witness_payloads))
        print(f"[m7]   {name}: {families[-1]['n_candidates']} rows, breaks={families[-1]['n_breaks_total']}", flush=True)

    t44 = time.time()
    probe44 = T2.generate_probe_families(family, labels, [(4, 4)], baseline_keys)
    families.append(
        screen_family(
            "existing_tangent_(4,4)_probe",
            "existing",
            probe44["4,4"],
            witness_payloads,
            extra={"generation_seconds": round(time.time() - t44, 3)},
        )
    )
    print(f"[m7]   tangent_(4,4): {families[-1]['n_candidates']} rows, breaks={families[-1]['n_breaks_total']}", flush=True)

    # --- Phase C: candidate families (pairing only) --------------------------------------------
    t1 = time.time()
    gen3 = generate_coordinate_ibp_rows(family, labels, 3)
    rej3: dict[str, int] = {}
    for rej in gen3.rejected:
        rej3[rej.reason] = rej3.get(rej.reason, 0) + 1
    new3 = [r for r in gen3.rows if r.dedup_key() not in baseline_keys]
    families.append(
        screen_family(
            "candidate_ibp_deg3",
            "candidate",
            new3,
            witness_payloads,
            extra={
                "generation_seconds": round(time.time() - t1, 3),
                "n_generated": len(gen3.rows),
                "n_new_vs_baseline": len(new3),
                "surface_rejected": rej3,
            },
        )
    )
    print(f"[m7]   ibp_deg3: {families[-1]['n_candidates']} new rows, breaks={families[-1]['n_breaks_total']}", flush=True)
    keys_so_far = baseline_keys | {r.dedup_key() for r in gen3.rows}

    t2 = time.time()
    rows_ray, rej_ray, ray_truncated = generate_ray_multiplier_rows(
        family, labels, degrees=RAY_DEGREES, max_candidates=max_candidates
    )
    new_ray = [r for r in rows_ray if r.dedup_key() not in keys_so_far]
    families.append(
        screen_family(
            "candidate_ray_multipliers_deg4-6",
            "candidate",
            new_ray,
            witness_payloads,
            extra={
                "generation_seconds": round(time.time() - t2, 3),
                "ray_note": (
                    "single-variable monomial multipliers x_j^d; primary failing-ray variable "
                    f"(index {PRIMARY_RAY_VAR}, audit ray (-1,0,0)) probed first"
                ),
                "degrees": list(RAY_DEGREES),
                "n_generated": len(rows_ray),
                "n_new_vs_baseline_and_deg3": len(new_ray),
                "surface_rejected": rej_ray,
                "truncated_by_cap": ray_truncated,
            },
        )
    )
    print(f"[m7]   ray_multipliers: {families[-1]['n_candidates']} new rows, breaks={families[-1]['n_breaks_total']}", flush=True)

    rows55, meta55 = generate_tangent55_rows(
        family, labels, budget_seconds=tangent55_budget, max_candidates=max_candidates
    )
    new55 = [r for r in rows55 if r.dedup_key() not in baseline_keys]
    families.append(
        screen_family(
            "candidate_tangent_(5,5)",
            "candidate",
            new55,
            witness_payloads,
            extra=dict(meta55, n_new_vs_baseline=len(new55)),
        )
    )
    print(
        f"[m7]   tangent_(5,5): {families[-1]['n_candidates']} new rows, breaks={families[-1]['n_breaks_total']} "
        f"(status={meta55['status']})",
        flush=True,
    )

    # --- decision ------------------------------------------------------------------------------
    existing_breaks = sum(f["n_breaks_total"] for f in families if f["category"] == "existing")
    candidate_breaks = sum(f["n_breaks_total"] for f in families if f["category"] == "candidate")
    breaking_families = sorted(
        f["family"] for f in families if f["category"] == "candidate" and f["n_breaks_total"] > 0
    )
    if candidate_breaks > 0:
        recommendation = (
            "targeted enlargement justified: re-run the Level-0 elimination with the breaking "
            f"famil(ies) {breaking_families} added (see breaking_kinds/breaking_examples)"
        )
    else:
        recommendation = (
            "no candidate row pairs nontrivially with any stable witness; the Level-0 "
            "obstruction is witness-certified stable under IBP-degree-3, ray-multiplier "
            "(deg 4-6) and budgeted tangent-(5,5) enlargements; further enlargement along "
            "these directions is not justified"
        )

    payload = {
        "script": "run_external_int2_method7.py",
        "method": "External Int2 Method.7: pairing-only candidate screening (Phase B/C)",
        "scope_note": SCOPE_NOTE,
        "level": LEVEL,
        "label_box": [
            [list(r) for r in spec.n_ranges],
            [list(spec.m_range) for _ in range(family.npolys)],
        ],
        "baseline": {
            "max_ibp_degree": T2.BASELINE["max_ibp_degree"],
            "tangent_degree_blocks": [list(b) for b in T2.BASELINE["tangent_degree_blocks"]],
            "extra_block": list(T2.BASELINE["extra_block"]),
        },
        "target": list(target_label),
        "n_labels": len(labels),
        "n_rows_total": rows["n_rows_total"],
        "witness_points": [
            {
                "sample": [list(pair) for pair in wp["sample"]],
                "prime": wp["prime"],
                "support_size": len(wp["witness"]),
                "support_hash": support_hash([lab for lab, _c in wp["witness"]]),
            }
            for wp in witness_payloads
        ],
        "witness_stability_summary": witness_stability,
        "families": families,
        "surface_rejected_diagnostics": {
            "note": (
                "diagnostic only; RejectedRow carries no terms and can never enter a pairing "
                "family"
            ),
            "existing_richer_(3,3)": rows["richer_row_rejections"],
            "candidate_ibp_deg3": rej3,
            "candidate_ray_multipliers": rej_ray,
            "candidate_tangent_(5,5)": meta55["rejected"],
        },
        "decision": {
            "existing_rows_all_annihilate": existing_breaks == 0,
            "n_existing_breaks": existing_breaks,
            "n_candidate_breaks": candidate_breaks,
            "breaking_candidate_families": breaking_families,
            "rerun_justified": candidate_breaks > 0,
            "recommendation": recommendation,
            "rerun_caveat": RERUN_CAVEAT,
        },
        "elapsed_seconds": round(time.time() - t0, 3),
    }
    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    return payload


# --- describe / main ---------------------------------------------------------------------------


def _describe(out_dir: Path) -> None:
    spec = T2.LEVELS[LEVEL]
    print("External Int2 Method.7 (dual-witness stability + pairing-only screening):")
    print(f"  medium box: Level {LEVEL} n_ranges={spec.n_ranges} m_range={spec.m_range}")
    print(f"  baseline: {T2.BASELINE}")
    print(
        f"  phase witness: {MIN_GENERIC} generic samples x {len(T2.DEFAULT_WITNESS_PRIMES)} primes"
        f" -> {WITNESS_OUT_NAME} (HEAVY: ~6 Level-0 witness RREFs; needs --allow-heavy)"
    )
    print(
        f"  phase screen: pairing-only families (existing kinds + (4,4); ibp_deg3, "
        f"ray x_j^d d={list(RAY_DEGREES)}, tangent {list(TANGENT55_BLOCK)} budgeted)"
        f" -> {SCREEN_OUT_NAME} (NO RREF)"
    )
    for name in (WITNESS_OUT_NAME, SCREEN_OUT_NAME):
        path = out_dir / name
        if not path.exists():
            print(f"  artifact {name}: absent")
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        if "stability" in data:
            st = data["stability"]
            print(
                f"  artifact {name}: {st['n_witness_points']} witness points, "
                f"support_pattern_stable={st['support_pattern_stable']}, "
                f"nullities={st['nullities']}, all_checks_pass={st['all_checks_pass']}"
            )
        else:
            dec = data.get("decision", {})
            print(
                f"  artifact {name}: candidate_breaks={dec.get('n_candidate_breaks')}, "
                f"rerun_justified={dec.get('rerun_justified')}"
            )


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="run_external_int2_method7",
        description="External Int2 Method.7: dual-witness stability + pairing-only screening.",
    )
    parser.add_argument("--input", type=Path, default=T2.DEFAULT_INPUT)
    parser.add_argument("--out-dir", type=Path, default=T2.DEFAULT_OUT_DIR)
    parser.add_argument("--describe", action="store_true", help="print plan + artifact summaries")
    parser.add_argument(
        "--phase", choices=("witness", "screen", "all"), default=None,
        help="nothing runs by default; 'witness' is HEAVY and needs --allow-heavy",
    )
    parser.add_argument("--allow-heavy", action="store_true", help="required for phase witness/all")
    parser.add_argument(
        "--witness-primes", type=str, default=",".join(str(p) for p in T2.DEFAULT_WITNESS_PRIMES)
    )
    parser.add_argument("--tangent55-budget", type=float, default=DEFAULT_TANGENT55_BUDGET)
    parser.add_argument("--max-candidates", type=int, default=DEFAULT_MAX_CANDIDATES)
    parser.add_argument(
        "--witness-file", type=Path, default=None,
        help="witness artifact to screen against (default: <out-dir>/" + WITNESS_OUT_NAME + ")",
    )
    args = parser.parse_args(argv)

    if args.describe:
        _describe(args.out_dir)
        return 0

    if args.phase is None:
        print("error: --phase is required (nothing runs by default); see --describe", file=sys.stderr)
        return 2
    try:
        primes = [int(p.strip()) for p in args.witness_primes.split(",") if p.strip()]
    except ValueError as exc:
        print(f"error: bad --witness-primes: {exc}", file=sys.stderr)
        return 2
    if len(set(primes)) < 2 or any(p < 2 for p in primes):
        print("error: need >= 2 distinct primes, each >= 2", file=sys.stderr)
        return 2
    if args.tangent55_budget <= 0:
        print("error: --tangent55-budget must be > 0", file=sys.stderr)
        return 2
    if args.max_candidates < 1:
        print("error: --max-candidates must be >= 1", file=sys.stderr)
        return 2
    if args.phase in ("witness", "all") and not args.allow_heavy:
        print(
            f"error: phase {args.phase!r} runs ~{MIN_GENERIC * len(primes)} Level-0 witness "
            "RREFs (HEAVY); pass --allow-heavy",
            file=sys.stderr,
        )
        return 2

    try:
        text = args.input.read_text(encoding="utf-8")
    except OSError as exc:
        print(f"error: cannot read {args.input}: {exc}", file=sys.stderr)
        return 2
    try:
        family = T2.parse_family_text(text)
        _target, base_config = T2.build_reducer_config(family)
    except (T2.ParserError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    args.out_dir.mkdir(parents=True, exist_ok=True)
    witness_payload: dict | None = None

    if args.phase in ("witness", "all"):
        print(f"[m7] phase witness (primes={primes}) ...", flush=True)
        witness_payload = run_witness_phase(
            family,
            base_config,
            primes=primes,
            out_path=args.out_dir / WITNESS_OUT_NAME,
            recorded_dir=args.out_dir,
        )
        st = witness_payload["stability"]
        print(
            f"[m7] witness done: {st['n_witness_points']}/{st['n_points']} Witness points, "
            f"support_pattern_stable={st['support_pattern_stable']}, "
            f"all_checks_pass={st['all_checks_pass']}",
            flush=True,
        )

    if args.phase in ("screen", "all"):
        if witness_payload is None:
            witness_path = args.witness_file or (args.out_dir / WITNESS_OUT_NAME)
            if not witness_path.exists():
                print(
                    f"error: witness artifact {witness_path} absent; run --phase witness "
                    "--allow-heavy first (or pass --witness-file)",
                    file=sys.stderr,
                )
                return 2
            witness_payload = json.loads(witness_path.read_text(encoding="utf-8"))
        print("[m7] phase screen ...", flush=True)
        screen = run_screen_phase(
            family,
            base_config,
            witness_payload["points"],
            out_path=args.out_dir / SCREEN_OUT_NAME,
            tangent55_budget=args.tangent55_budget,
            max_candidates=args.max_candidates,
            witness_stability=witness_payload.get("stability"),
        )
        dec = screen["decision"]
        print(
            f"[m7] screen done: existing_all_annihilate={dec['existing_rows_all_annihilate']}, "
            f"candidate_breaks={dec['n_candidate_breaks']}, "
            f"rerun_justified={dec['rerun_justified']}",
            flush=True,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
