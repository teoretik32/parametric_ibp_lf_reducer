"""Stage gate audit for the External Int2 Method.11a staged ladder (Phase C).

Reads the m11-record-cache/v2 JSONL plus (optionally) the stage artifact JSON
and decides whether the next, more expensive stage is justified:

  DONE     - the stage reconstruction succeeded; the ladder stops here.
  PROCEED  - honest capacity-bound failure: supports are stable across all
             max-rank records, every expected point has full prime coverage,
             and no unexpected rank drops occurred.  Buying more points is
             mathematically justified.
  STOP     - support instability, unexpected rank drop, missing prime
             coverage, or an artifact that failed for a non-capacity reason.
             Burning more compute would not be justified; diagnose first.

Analysis-only: never modifies the cache or the artifact.
Exit code: 0 for DONE/PROCEED, 2 for STOP.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

_SCRIPTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS_DIR.parent / "src"))

from parametric_ibp_lf_reducer.reconstruction import classify_sample_supports  # noqa: E402

# Substrings (lowercased) in a failed reconstruction section that mark the
# failure as capacity-bound (too few points), i.e. curable by the next stage.
CAPACITY_MARKERS = ("insufficient", "interpolationfailed", "underdetermined")


def load_m11():
    """Import the Method.11 runner module (single source of record schema)."""
    path = _SCRIPTS_DIR / "run_external_int2_method11.py"
    spec = importlib.util.spec_from_file_location("m11_gate_load", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_cache(path, m11):
    """Return (header|None, records) from a record-cache JSONL file."""
    header = None
    records = []
    with open(path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, start=1):
            line = line.strip()
            if not line:
                continue
            data = json.loads(line)
            if lineno == 1 and "fingerprint" in data:
                header = data
                continue
            records.append(m11._record_from_json(data))
    return header, records


def audit_records(records, m11, expected_points, expected_primes):
    """Pure record-side audit; returns a JSON-ready report dict."""
    ok = [r for r in records if r.formal_success]
    rank_histogram = {}
    for r in ok:
        rank_histogram[r.rank] = rank_histogram.get(r.rank, 0) + 1
    max_rank = max(rank_histogram) if rank_histogram else None
    top = [r for r in ok if r.rank == max_rank]

    reference = top[0].support if top else ()
    ref_set = set(reference)
    deviations = []
    for r in top:
        got = set(r.support)
        if got != ref_set:
            deviations.append(
                {
                    "key": m11.sample_prime_key(r.sample, r.prime),
                    "missing": sorted(map(list, ref_set - got)),
                    "extra": sorted(map(list, got - ref_set)),
                }
            )

    coverage = {}
    for r in top:
        # sample_prime_key -> 'p=<prime>;<canonical sample>'; strip the prime part
        point = m11.sample_prime_key(r.sample, 0).split(";", 1)[1]
        coverage.setdefault(point, set()).add(r.prime)
    incomplete = sorted(k for k, v in coverage.items() if len(v) < expected_primes)

    # Method.11b: shared library semantics for special-zero vs unstable supports.
    classification = classify_sample_supports(records)

    return {
        "n_records": len(records),
        "n_formal_success": len(ok),
        "rank_histogram": {str(k): v for k, v in sorted(rank_histogram.items())},
        "max_rank": max_rank,
        "n_max_rank_records": len(top),
        "n_rank_drops": len(ok) - len(top),
        "support_size": len(ref_set),
        "support_deviations": deviations,
        "support_classification": classification["by_sample"],
        "support_summary": classification["summary"],
        "n_points_seen": len(coverage),
        "expected_points": expected_points,
        "points_incomplete_primes": incomplete,
    }


def artifact_reconstruction(artifact):
    """Extract (success, capacity_bound, note) from a stage artifact dict."""
    rec = artifact.get("reconstruction")
    if not isinstance(rec, dict):
        return None, False, "artifact lacks a reconstruction section"
    if rec.get("formal_success") is True:
        return True, False, "reconstruction succeeded"
    blob = json.dumps(rec).lower()
    capacity = any(m in blob for m in CAPACITY_MARKERS)
    return False, capacity, "capacity-bound failure" if capacity else "non-capacity failure"


def gate_verdict(report, artifact):
    """Combine record audit and artifact into (verdict, reasons)."""
    reasons = []
    if artifact is not None:
        success, capacity, note = artifact_reconstruction(artifact)
        if success:
            return "DONE", [note]
        reasons.append(note)
        if success is False and not capacity:
            reasons.append("reconstruction failed for a non-capacity reason")
    summary = report.get("support_summary", {})
    if summary.get("n_special_zero"):
        reasons.append(
            f"special-zero support at {summary['n_special_zero']} sample(s) (retained, zero-filled)"
        )
    if summary.get("n_unstable"):
        reasons.append(f"unstable support at {summary['n_unstable']} sample(s)")
    if report["n_rank_drops"]:
        reasons.append(f"{report['n_rank_drops']} formal-success record(s) below max rank")
    if report["n_points_seen"] < report["expected_points"]:
        reasons.append(
            f"only {report['n_points_seen']}/{report['expected_points']} points at max rank"
        )
    if report["points_incomplete_primes"]:
        reasons.append(
            f"{len(report['points_incomplete_primes'])} point(s) with incomplete prime coverage"
        )
    bad = [
        r
        for r in reasons
        if r not in ("capacity-bound failure",) and not r.startswith("special-zero support")
    ]
    return ("PROCEED", reasons or ["stable supports, full coverage"]) if not bad else ("STOP", bad)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cache", default=None, help="record cache JSONL (default: runner's)")
    parser.add_argument("--artifact", default=None, help="stage artifact JSON (optional)")
    parser.add_argument("--expected-points", type=int, required=True)
    parser.add_argument("--expected-primes", type=int, default=3)
    parser.add_argument("--json-out", default=None)
    args = parser.parse_args(argv)

    m11 = load_m11()
    cache = args.cache or str(m11.DEFAULT_RECORD_CACHE)
    header, records = load_cache(cache, m11)
    artifact = None
    if args.artifact:
        artifact = json.loads(Path(args.artifact).read_text(encoding="utf-8"))

    report = audit_records(records, m11, args.expected_points, args.expected_primes)
    verdict, reasons = gate_verdict(report, artifact)
    report["cache_has_v2_header"] = header is not None
    report["verdict"] = verdict
    report["reasons"] = reasons

    if args.json_out:
        Path(args.json_out).write_text(
            json.dumps(report, indent=1, default=str) + "\n", encoding="utf-8"
        )
    print(json.dumps(report, indent=1, default=str))
    return 0 if verdict in ("DONE", "PROCEED") else 2


if __name__ == "__main__":
    sys.exit(main())
