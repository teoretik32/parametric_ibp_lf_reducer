"""External Int2 T2 rank-repair: reproduce Levels 0-2 + dual LF-obstruction witness (Method.6).

Reproduces the three recorded "rank-repair" configurations (does enlarging the label box restore
LF feasibility?) and, on top of the read-only span test, produces an explicit DUAL WITNESS for the
``Obstructed`` points: a right-nullspace vector ``w`` with ``w[target] == 1`` and ``<row, w> == 0``
for every projected row (see :mod:`parametric_ibp_lf_reducer.lf_obstruction_witness`).

Levels (label boxes are SHIFTS against the Int2 base, offset convention):

* Level 0: ``n2 in [-1,1]``, ``n5,n7 in [0,1]``, ``m0..m3 in [-3,0]``  (3072 labels, 46737 rows)
* Level 1: ``n2 in [-1,2]``, ``n5,n7 in [0,1]``, ``m0..m3 in [-3,0]``  (4096 labels, 59605 rows) HEAVY
* Level 2: ``n2 in [-1,2]``, ``n5,n7 in [0,1]``, ``m0..m3 in [-4,0]``  (10000 labels, 155298 rows) HEAVY

Rows: baseline (coordinate IBP ``max_ibp_degree=2`` + tangent blocks ``((1,1),(2,2))``) plus a
richer tangent block ``(3,3)`` merged by ``dedup_key`` — the exact system the recorded feasibility
runs used, so the witness certifies precisely those rows. Levels 1-2 are HEAVY and only run with an
explicit ``--allow-heavy`` gate; nothing runs at all without ``--levels``.

Codimension correction (binding). ``residual_support == [target]`` says only that the canonical
residual of ``e_target`` lands on the target coordinate; it does NOT imply the quotient dimension is
one. That dimension is the nullity ``n_projected_cols - rank`` and may exceed 1. New artifacts and
prose use this corrected wording; the recorded JSON ``purpose`` strings are historical and left
byte-identical.

Scope. Read-only modular LF-feasibility + dual-witness diagnostics; the span test is constrained to
LF-True labels only (a span through divergent labels is never reported as an LF basis); reducer
core, certificates and LF gates are untouched. ``Obstructed``/``Witness`` are honest per-(sample,
prime), per-label-box negatives — never a global impossibility claim.
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from dataclasses import dataclass, replace
from fractions import Fraction
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:  # allow running without installation
    sys.path.insert(0, str(REPO_ROOT / "src"))

from parametric_ibp_lf_reducer import ParserError, parse_family_text  # noqa: E402
from parametric_ibp_lf_reducer.api import (  # noqa: E402
    build_reducer_config,
    default_scattered_samples,
)
from parametric_ibp_lf_reducer.lf_feasibility import (  # noqa: E402
    STATUS_FEASIBLE,
    STATUS_OBSTRUCTED,
    feasibility_to_payload,
    lf_reduction_feasible_mod_p,
)
from parametric_ibp_lf_reducer.lf_obstruction_witness import (  # noqa: E402
    STATUS_WITNESS,
    lf_obstruction_witness_mod_p,
    pairings_to_payload,
    test_rows_against_obstruction_witness,
    witness_from_payload,
    witness_to_payload,
)
from parametric_ibp_lf_reducer.reducer import _enumerate_labels, _generate_rows  # noqa: E402
from parametric_ibp_lf_reducer.row_generation import generate_tangent_ibp_rows  # noqa: E402
from parametric_ibp_lf_reducer.tangent_fields import generate_tangent_fields  # noqa: E402
from parametric_ibp_lf_reducer.valuations import is_locally_finite  # noqa: E402

DEFAULT_INPUT = REPO_ROOT / "examples" / "external_int2_dimensionless_input.wl.txt"
DEFAULT_OUT_DIR = REPO_ROOT / "validation"

VERDICT_FEASIBLE = "Feasible"
VERDICT_OBSTRUCTED = "Obstructed"
VERDICT_MIXED = "Mixed"
VERDICT_INCONCLUSIVE = "Inconclusive"

DEFAULT_PRIME = 2147483647
DEFAULT_WITNESS_PRIMES = (2147483647, 2147483629)

#: Corrected purpose phrasing for the reproduced artifacts. The recorded JSONs' flawed
#: "codimension-1 truncation hypothesis" text is intentionally NOT reused.
REPRO_PURPOSE = (
    "box-truncation obstruction probe; residual_support=[target] does not imply quotient "
    "dimension one (quotient dimension = nullity = projected cols - rank, may exceed 1)"
)

SCOPE_NOTE = (
    "T2 rank-repair reproduction + dual LF-obstruction witness (Method.6). Modular LF-feasibility "
    "only: the span test is constrained to LF-True labels; a span through divergent labels is "
    "never reported as an LF basis. Rank-deficient special points (ep=3) are excluded from "
    "verdicts. Per-(sample, prime), per-label-box statements only; 'Obstructed'/'Witness' are "
    "honest negatives, never a global impossibility claim. No reconstruction; reducer core, "
    "certificates and LF gates are untouched."
)

LEVEL2_EP3_NOTE = "ep=3 special sample omitted; its Feasible status established at Levels A/B/0/1"


@dataclass(frozen=True)
class LevelSpec:
    """One rank-repair label box (n-shift ranges per variable + uniform m-shift range)."""

    n_ranges: tuple[tuple[int, int], ...]
    m_range: tuple[int, int]
    heavy: bool
    omit_special: bool  # Level 2 drops the ep=3 (index-2) scattered sample


LEVELS: dict[int, LevelSpec] = {
    0: LevelSpec(n_ranges=((-1, 1), (0, 1), (0, 1)), m_range=(-3, 0), heavy=False, omit_special=False),
    1: LevelSpec(n_ranges=((-1, 2), (0, 1), (0, 1)), m_range=(-3, 0), heavy=True, omit_special=False),
    2: LevelSpec(n_ranges=((-1, 2), (0, 1), (0, 1)), m_range=(-4, 0), heavy=True, omit_special=True),
}

BASELINE = {
    "max_ibp_degree": 2,
    "tangent_degree_blocks": ((1, 1), (2, 2)),
    "extra_block": (3, 3),
}

#: The ep=3 rank-deficient special sample lives at scattered-sample index 2 (verified against the
#: recorded artifacts: default_scattered_samples(["ep","r"]) yields ep=3, r=54/11 at k=2).
SPECIAL_SAMPLE_INDEX = 2

RECORDED_NAME = "external_int2_t2_rankrepair_level{level}.json"
REPRO_NAME = "external_int2_t2_rankrepair_level{level}_repro.json"
WITNESS_NAME = "external_int2_t2_witness_level{level}.json"
ROWPROBE_NAME = "external_int2_t2_witness_rowprobe_level{level}.json"


@dataclass(frozen=True)
class WitnessConfig:
    """Phase C witness extraction settings (generic samples only)."""

    samples: int = 2
    primes: tuple[int, ...] = DEFAULT_WITNESS_PRIMES


# --- genericity classification (copied from run_external_int2_method5.py) ----------------------


def _analyze_points(points: list[dict], min_generic: int) -> dict:
    """Split points into generic vs rank-deficient special; verdict over generic only.

    (copied from run_external_int2_method5.py: scripts/ is not a package.)
    """
    max_rank_by_prime: dict[int, int] = {}
    for p in points:
        r = p["result"]["rank"]
        prime = p["prime"]
        max_rank_by_prime[prime] = max(max_rank_by_prime.get(prime, 0), r)
    generic, special = [], []
    for p in points:
        deficit = max_rank_by_prime[p["prime"]] - p["result"]["rank"]
        p["rank_deficit"] = deficit
        p["classification"] = "generic" if deficit == 0 else "rank_deficient_special"
        (generic if deficit == 0 else special).append(p)

    statuses = [p["result"]["status"] for p in generic]
    real = [s for s in statuses if s in (STATUS_FEASIBLE, STATUS_OBSTRUCTED)]
    if len(real) < min_generic:
        verdict = VERDICT_INCONCLUSIVE
    elif all(s == STATUS_FEASIBLE for s in real):
        verdict = VERDICT_FEASIBLE
    elif all(s == STATUS_OBSTRUCTED for s in real):
        verdict = VERDICT_OBSTRUCTED
    else:
        verdict = VERDICT_MIXED

    special_only_feasible = (
        any(p["result"]["status"] == STATUS_FEASIBLE for p in special)
        and not any(s == STATUS_FEASIBLE for s in real)
    )

    supports = [tuple(map(tuple, p["result"]["residual_support"])) for p in generic]
    support_sets = [set(s) for s in supports] or [set()]
    common = set.intersection(*support_sets) if supports else set()
    return {
        "n_points": len(points),
        "n_generic": len(generic),
        "n_special": len(special),
        "max_rank_by_prime": {str(k): v for k, v in sorted(max_rank_by_prime.items())},
        "generic_verdict": verdict,
        "min_generic_required": min_generic,
        "special_only_feasible_rejected": special_only_feasible,
        "residual_support_stable": bool(supports) and len(set(supports)) == 1,
        "common_residual_support_size": len(common),
        "common_residual_support": sorted(list(lab) for lab in common),
    }


# --- config / row assembly ---------------------------------------------------------------------


def build_box_labels(family, spec: LevelSpec) -> list[tuple[int, ...]]:
    """Explicit per-variable label enumeration for one level box."""
    if len(spec.n_ranges) != family.nvars:
        raise ValueError(
            f"level n_ranges has {len(spec.n_ranges)} axes but family has {family.nvars} variables"
        )
    n_axes = [range(lo, hi + 1) for lo, hi in spec.n_ranges]
    m_lo, m_hi = spec.m_range
    m_axes = [range(m_lo, m_hi + 1)] * family.npolys
    return [tuple(parts) for parts in itertools.product(*(n_axes + m_axes))]


def build_level_config(family, base_config, spec: LevelSpec):
    """(labels, target_label, config) for one level — baseline blocks fixed per the recorded runs."""
    labels = build_box_labels(family, spec)
    target_label = tuple([0] * (family.nvars + family.npolys))
    config = replace(
        base_config,
        labels=tuple(labels),
        label_box=None,
        max_ibp_degree=BASELINE["max_ibp_degree"],
        tangent_degree_blocks=BASELINE["tangent_degree_blocks"],
    )
    return labels, target_label, config


def _merge_extra_rows(base_rows, extra_rows) -> tuple[list, int]:
    """Merge ``extra_rows`` into a copy of ``base_rows`` by ``dedup_key`` (baseline wins)."""
    seen = {row.dedup_key() for row in base_rows}
    merged = list(base_rows)
    n_dup = 0
    for row in extra_rows:
        key = row.dedup_key()
        if key in seen:
            n_dup += 1
            continue
        seen.add(key)
        merged.append(row)
    return merged, n_dup


def assemble_level_rows(family, labels, config) -> dict:
    """Baseline rows + the ``(3,3)`` extra tangent block, merged. Returns rows + diagnostics."""
    base_rows, row_diag = _generate_rows(family, labels, config, None)
    fields = generate_tangent_fields(family, [BASELINE["extra_block"]])
    extra = generate_tangent_ibp_rows(family, labels, fields)
    rejected: dict[str, int] = {}
    for rej in extra.rejected:
        rejected[rej.reason] = rejected.get(rej.reason, 0) + 1
    merged, n_dup = _merge_extra_rows(base_rows, extra.rows)
    return {
        "base_rows": base_rows,
        "merged": merged,
        "n_base_rows": len(base_rows),
        "n_extra_offered": len(extra.rows),
        "n_extra_new": len(merged) - len(base_rows),
        "n_extra_duplicate": n_dup,
        "n_rows_total": len(merged),
        "row_diagnostics": row_diag,
        "richer_row_rejections": rejected,
    }


def _level_samples(family, spec: LevelSpec, n_samples: int) -> list[dict]:
    """Scattered samples; Level 2 drops the ep=3 special sample at index 2."""
    samples = default_scattered_samples(family.parameters, n_samples)
    if spec.omit_special:
        samples = [s for i, s in enumerate(samples) if i != SPECIAL_SAMPLE_INDEX]
    return samples


def _sample_repr(sample: dict) -> list:
    return sorted((k, str(v)) for k, v in sample.items())


def _sample_from_repr(pairs) -> dict:
    """Recover a sample dict (exact Fraction values) from stored ``[[name, value], ...]`` pairs."""
    return {name: Fraction(value) for name, value in pairs}


# --- witness extraction ------------------------------------------------------------------------


def _generic_samples(points: list[dict], samples: list[dict]) -> list[dict]:
    """Samples all of whose (per-prime) points are classified ``generic``."""
    special_reprs = {
        tuple(map(tuple, p["sample"]))
        for p in points
        if p.get("classification") == "rank_deficient_special"
    }
    out = []
    for s in samples:
        if tuple(map(tuple, _sample_repr(s))) not in special_reprs:
            out.append(s)
    return out


def _witness_section(
    merged, labels, target_label, lf_map, generic_samples, witness_cfg: WitnessConfig
) -> dict:
    """Dual witness at up to ``witness_cfg.samples`` generic samples x ``witness_cfg.primes``."""
    chosen = generic_samples[: witness_cfg.samples]
    entries: list[dict] = []
    for sample in chosen:
        for prime in witness_cfg.primes:
            res = lf_obstruction_witness_mod_p(
                merged, labels, target_label, lf_map, sample, prime
            )
            support = sorted(list(lab) for lab, _ in res.witness)
            entries.append(
                {
                    "sample": _sample_repr(sample),
                    "prime": prime,
                    "status": res.status,
                    "witness": witness_to_payload(res),
                    "support_labels": support,
                    "support_size": len(res.witness),
                    "rank": res.rank,
                    "nullity": res.nullity,
                    "checks_pass": bool(res.check_annihilation and res.check_target_unit),
                }
            )
    support_patterns = {
        tuple(tuple(lab) for lab in e["support_labels"]) for e in entries if e["status"] == STATUS_WITNESS
    }
    statuses = {e["status"] for e in entries}
    return {
        "n_witness_samples": len(chosen),
        "witness_primes": list(witness_cfg.primes),
        "points": entries,
        "support_pattern_stable": len(support_patterns) <= 1,
        "n_distinct_support_patterns": len(support_patterns),
        "all_checks_pass": all(e["checks_pass"] for e in entries if e["status"] == STATUS_WITNESS),
        "witness_obstruction_consistent": statuses <= {STATUS_WITNESS, STATUS_FEASIBLE}
        and STATUS_WITNESS in statuses,
    }


# --- level run ---------------------------------------------------------------------------------


def run_level(
    family,
    base_config,
    spec: LevelSpec,
    samples,
    primes,
    *,
    min_generic: int = 2,
    witness_cfg: WitnessConfig | None = None,
    out_path: Path | None = None,
) -> dict:
    """Reproduce one level's LF-feasibility (and optional dual witness). Returns the JSON payload."""
    t0 = time.time()
    labels, target_label, config = build_level_config(family, base_config, spec)
    assert list(_enumerate_labels(family, config)) == labels

    lf_map = {lab: is_locally_finite(family, lab) for lab in labels}
    if target_label not in lf_map:
        lf_map[target_label] = is_locally_finite(family, target_label)

    rows = assemble_level_rows(family, labels, config)
    merged = rows["merged"]

    points: list[dict] = []
    for sample in samples:
        for prime in primes:
            res = lf_reduction_feasible_mod_p(merged, labels, target_label, lf_map, sample, prime)
            points.append(
                {
                    "prime": prime,
                    "sample": _sample_repr(sample),
                    "classification": "generic",  # refined by _analyze_points
                    "result": feasibility_to_payload(res),
                }
            )
    analysis = _analyze_points(points, min_generic)

    payload: dict = {
        "script": "run_external_int2_t2_rankrepair.py",
        "method": "External Int2 T2 rank-repair + dual LF-obstruction witness (Method.6)",
        "scope_note": SCOPE_NOTE,
        "purpose": REPRO_PURPOSE,
        "label_box": [
            [list(r) for r in spec.n_ranges],
            [list(spec.m_range) for _ in range(family.npolys)],
        ],
        "baseline": {
            "max_ibp_degree": BASELINE["max_ibp_degree"],
            "tangent_degree_blocks": [list(b) for b in BASELINE["tangent_degree_blocks"]],
            "extra_block": list(BASELINE["extra_block"]),
        },
        "target": list(target_label),
        "target_lf_verdict": str(lf_map[target_label]),
        "n_labels": len(labels),
        "n_lf_true": sum(1 for lab in labels if lf_map[lab] is True),
        "n_base_rows": rows["n_base_rows"],
        "n_extra_new": rows["n_extra_new"],
        "n_rows_total": rows["n_rows_total"],
        "row_diagnostics": rows["row_diagnostics"],
        "richer_row_rejections": rows["richer_row_rejections"],
        "samples": [_sample_repr(s) for s in samples],
        "primes": list(primes),
        "points": points,
        "analysis": analysis,
        "verdict": analysis["generic_verdict"],
    }
    if spec.omit_special:
        payload["note"] = LEVEL2_EP3_NOTE

    if witness_cfg is not None:
        generic_samples = _generic_samples(points, samples)
        payload["witness"] = _witness_section(
            merged, labels, target_label, lf_map, generic_samples, witness_cfg
        )

    payload["elapsed_seconds"] = round(time.time() - t0, 3)
    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    return payload


# --- row probe (report-only; NEVER re-eliminates) ----------------------------------------------


def _parse_extra_blocks(text: str) -> list[tuple[int, int]]:
    blocks: list[tuple[int, int]] = []
    for part in text.split(";"):
        part = part.strip()
        if not part:
            continue
        pieces = part.split(",")
        if len(pieces) != 2:
            raise ValueError(f"bad block {part!r}: expected 'd_q,d_h'")
        d_q, d_h = int(pieces[0]), int(pieces[1])
        if d_q < 0 or d_h < 0:
            raise ValueError(f"bad block {part!r}: degrees must be >= 0")
        blocks.append((d_q, d_h))
    if not blocks:
        raise ValueError("no degree blocks given")
    return blocks


def generate_probe_families(family, labels, extra_blocks, baseline_keys) -> dict:
    """Genuinely-new candidate rows per extra tangent block (post-dedup vs the baseline system)."""
    families: dict[str, list] = {}
    for block in extra_blocks:
        fields = generate_tangent_fields(family, [block])
        gen = generate_tangent_ibp_rows(family, labels, fields)
        new_rows = [row for row in gen.rows if row.dedup_key() not in baseline_keys]
        families[f"{block[0]},{block[1]}"] = new_rows
    return families


def probe_rows(
    family,
    base_config,
    spec: LevelSpec,
    witness_payload,
    extra_blocks,
    *,
    candidate_families=None,
    out_path: Path | None = None,
) -> dict:
    """Pair candidate rows against a stored witness (NO RREF, never re-eliminates a level).

    ``witness_payload`` is a single witness payload (from :func:`witness_to_payload`); the pairing
    happens at that witness's own ``(sample, prime)``. ``candidate_families`` (test/orchestrator
    injection) overrides symbolic regeneration.
    """
    t0 = time.time()
    witness = witness_from_payload(witness_payload)
    sample = _sample_from_repr([list(pair) for pair in witness.sample])
    prime = witness.prime

    if candidate_families is None:
        labels, _target, config = build_level_config(family, base_config, spec)
        rows = assemble_level_rows(family, labels, config)
        baseline_keys = {row.dedup_key() for row in rows["merged"]}
        candidate_families = generate_probe_families(family, labels, extra_blocks, baseline_keys)

    families_out: list[dict] = []
    any_break = False
    for name, rows in candidate_families.items():
        pairings = test_rows_against_obstruction_witness(rows, witness, sample, prime)
        breaking = [p for p in pairings if p.breaks]
        provenance: dict[str, int] = {}
        for p in breaking:
            provenance[p.kind] = provenance.get(p.kind, 0) + 1
        any_break = any_break or bool(breaking)
        families_out.append(
            {
                "block": name,
                "n_candidates": len(rows),
                "n_breaks": len(breaking),
                "n_annihilate": len(pairings) - len(breaking),
                "breaking_provenance": provenance,
                "pairings_sample": pairings_to_payload(pairings[:8]),
            }
        )

    payload = {
        "script": "run_external_int2_t2_rankrepair.py (probe-rows)",
        "scope_note": SCOPE_NOTE,
        "witness_sample": [list(pair) for pair in witness.sample],
        "witness_prime": prime,
        "witness_support_size": len(witness.witness),
        "extra_blocks": [list(b) for b in extra_blocks],
        "families": families_out,
        "rerun_justified": any_break,
        "rerun_caveat": "a breaking row is NECESSARY, not sufficient, to cure the obstruction "
        "(other nullvectors may still obstruct); no re-elimination is run here",
        "elapsed_seconds": round(time.time() - t0, 3),
    }
    if out_path is not None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")
    return payload


# --- describe ----------------------------------------------------------------------------------


def _describe(out_dir: Path) -> None:
    print("External Int2 T2 rank-repair levels (Method.6 reproduction):")
    for level, spec in LEVELS.items():
        gate = " (HEAVY: needs --allow-heavy)" if spec.heavy else ""
        omit = "; ep=3 sample omitted" if spec.omit_special else ""
        print(
            f"  level {level}: n_ranges={spec.n_ranges} m_range={spec.m_range}{gate}{omit}"
        )
    print(f"  baseline: {BASELINE}")
    print("Recorded artifacts (read-only; new runs write *_repro/_witness/_rowprobe):")
    for level in LEVELS:
        path = out_dir / RECORDED_NAME.format(level=level)
        if path.exists():
            data = json.loads(path.read_text(encoding="utf-8"))
            ranks = sorted({p["rank"] for p in data.get("points", [])})
            print(
                f"  level {level}: n_labels={data.get('n_labels')} "
                f"n_rows_total={data.get('n_rows_total')} ranks={ranks} "
                f"verdict={data.get('verdict')}"
            )
        else:
            print(f"  level {level}: recorded artifact absent ({path.name})")


# --- main --------------------------------------------------------------------------------------


def _parse_levels(text: str) -> list[int]:
    levels = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        levels.append(int(part))
    return levels


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="run_external_int2_t2_rankrepair",
        description="External Int2 T2 rank-repair reproduction + dual LF-obstruction witness.",
    )
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--levels", type=str, default=None, help="comma list from {0,1,2}")
    parser.add_argument("--describe", action="store_true", help="print level table + recorded summaries")
    parser.add_argument("--allow-heavy", action="store_true", help="required for levels 1 and 2")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--samples", type=int, default=4)
    parser.add_argument("--primes", type=int, nargs="+", default=[DEFAULT_PRIME])
    parser.add_argument("--min-generic", type=int, default=2)
    parser.add_argument("--compare", action="store_true", help="report-only diff vs recorded JSON")
    parser.add_argument("--witness", action="store_true", help="Phase C: extract dual witness")
    parser.add_argument("--witness-samples", type=int, default=2)
    parser.add_argument(
        "--witness-primes", type=str, default=",".join(str(p) for p in DEFAULT_WITNESS_PRIMES)
    )
    parser.add_argument(
        "--probe-rows", type=str, nargs="?", const="extra-blocks", default=None,
        help="Phase C: pair candidate rows against a stored witness (report-only)",
    )
    parser.add_argument("--probe-extra-blocks", type=str, default="3,3;4,4")
    parser.add_argument("--witness-file", type=Path, default=None)
    args = parser.parse_args(argv)

    if args.describe:
        _describe(args.out_dir)
        return 0

    if args.samples < 1 or args.min_generic < 1 or any(p < 2 for p in args.primes):
        print("error: --samples/--min-generic >= 1 and --primes >= 2 required", file=sys.stderr)
        return 2
    if args.levels is None:
        print("error: --levels is required (nothing runs by default); e.g. --levels 0", file=sys.stderr)
        return 2
    try:
        levels = _parse_levels(args.levels)
    except ValueError as exc:
        print(f"error: bad --levels: {exc}", file=sys.stderr)
        return 2
    unknown = [n for n in levels if n not in LEVELS]
    if unknown or not levels:
        print(f"error: unknown levels {unknown!r} (choose from {sorted(LEVELS)})", file=sys.stderr)
        return 2
    heavy = [n for n in levels if LEVELS[n].heavy]
    if heavy and not args.allow_heavy:
        print(f"error: levels {heavy} are HEAVY; pass --allow-heavy to run them", file=sys.stderr)
        return 2

    try:
        witness_primes = [int(p.strip()) for p in args.witness_primes.split(",") if p.strip()]
    except ValueError as exc:
        print(f"error: bad --witness-primes: {exc}", file=sys.stderr)
        return 2

    try:
        text = args.input.read_text(encoding="utf-8")
    except OSError as exc:
        print(f"error: cannot read {args.input}: {exc}", file=sys.stderr)
        return 2
    try:
        family = parse_family_text(text)
        _target, base_config = build_reducer_config(family)
    except (ParserError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    witness_cfg = (
        WitnessConfig(samples=args.witness_samples, primes=tuple(witness_primes))
        if args.witness
        else None
    )
    args.out_dir.mkdir(parents=True, exist_ok=True)

    for level in levels:
        spec = LEVELS[level]
        samples = _level_samples(family, spec, args.samples)
        print(f"[t2] level {level} ...", flush=True)
        repro_path = args.out_dir / REPRO_NAME.format(level=level)
        witness_path = args.out_dir / WITNESS_NAME.format(level=level) if args.witness else None
        payload = run_level(
            family,
            base_config,
            spec,
            samples,
            args.primes,
            min_generic=args.min_generic,
            witness_cfg=witness_cfg,
            out_path=repro_path,
        )
        print(
            f"[t2] level {level}: labels={payload['n_labels']} rows={payload['n_rows_total']} "
            f"verdict={payload['verdict']} -> {repro_path}",
            flush=True,
        )
        if witness_path is not None:
            witness_path.write_text(
                json.dumps(payload["witness"], indent=2, default=str) + "\n", encoding="utf-8"
            )
            print(f"[t2] level {level}: witness -> {witness_path}", flush=True)
        if args.compare:
            recorded_path = args.out_dir / RECORDED_NAME.format(level=level)
            if recorded_path.exists():
                rec = json.loads(recorded_path.read_text(encoding="utf-8"))
                print(
                    f"[t2] compare level {level}: verdict fresh={payload['verdict']} "
                    f"recorded={rec.get('verdict')}; n_rows fresh={payload['n_rows_total']} "
                    f"recorded={rec.get('n_rows_total')}",
                    flush=True,
                )
            else:
                print(f"[t2] compare level {level}: recorded artifact absent", flush=True)

    if args.probe_rows is not None:
        if args.witness_file is None or not args.witness_file.exists():
            print("error: --probe-rows needs an existing --witness-file", file=sys.stderr)
            return 2
        try:
            blocks = _parse_extra_blocks(args.probe_extra_blocks)
        except ValueError as exc:
            print(f"error: bad --probe-extra-blocks: {exc}", file=sys.stderr)
            return 2
        stored = json.loads(args.witness_file.read_text(encoding="utf-8"))
        witness_points = stored.get("points", [])
        chosen = next(
            (e["witness"] for e in witness_points if e.get("status") == STATUS_WITNESS), None
        )
        if chosen is None:
            print("error: --witness-file has no Witness-status point to probe", file=sys.stderr)
            return 2
        probe_level = heavy[0] if heavy else levels[0]
        probe_path = args.out_dir / ROWPROBE_NAME.format(level=probe_level)
        probe_payload = probe_rows(
            family, base_config, LEVELS[probe_level], chosen, blocks, out_path=probe_path
        )
        print(
            f"[t2] probe-rows: rerun_justified={probe_payload['rerun_justified']} -> {probe_path}",
            flush=True,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
