"""External Int2 finite-numerator LF basis search (Method.2 follow-up, task #37).

Honest per-sector search for numerator-decorated masters ``N(x) * F_S`` whose
COMPLETE integrand is locally finite (single-integrand semantics: no
cancellation after integration is ever assumed).  Sectors: the six certified
Method.2 normal-form sectors plus ``1/(G1*G3)``; numerator degrees 0..2 by
default.  The offset convention is respected throughout: total exponent =
family base exponent + label shift (all labels below are SHIFTS relative to
the Int2 base ``x2^(1+ep) * G0^ep * G1^ep * G2^(-1-ep) * G3^(-1+ep)``).

If (and only if) genuinely polynomial new LF masters are found, their
defining monomial-expansion labels (individually LF by Lemma 1) are exported
so the Method.1 modular span test (``lf_reduction_feasible_mod_p``) can mark
them allowed.  Otherwise the artifact records the honest overall verdict
``NoFiniteNumeratorBasisWithinAnsatz`` together with the failing-ray evidence.

Read-only with respect to reducer state, certificates and LF gates.

Exit codes: 0 = search completed and artifact written (``NoFiniteNumerator...``
is a finding, not a failure); 2 = usage/input problem or cap violation.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:  # allow running without installation
    sys.path.insert(0, str(REPO_ROOT / "src"))

from parametric_ibp_lf_reducer import ParserError, parse_family_text  # noqa: E402
from parametric_ibp_lf_reducer.finite_numerator import (  # noqa: E402
    STATUS_ALREADY_LF,
    finite_numerator_search,
)

DEFAULT_INPUT = REPO_ROOT / "examples" / "external_int2_dimensionless_input.wl.txt"
DEFAULT_JSON = REPO_ROOT / "validation" / "external_int2_finite_numerator.json"

MAX_RANDOM_TRIALS = 256

# Label SHIFTS (offset convention) relative to the Int2 base exponents.
# Six certified Method.2 normal-form sectors + the near-normal-form 1/(G1*G3).
SECTORS: dict[str, tuple[int, ...]] = {
    "1/(x2*G0*G1)": (-1, 0, 0, -1, -1, 0, 0),
    "1/(x2*G1*G3)": (-1, 0, 0, 0, -1, 0, -1),
    "1/(G0*G3)": (0, 0, 0, -1, 0, 0, -1),
    "1/G1": (0, 0, 0, 0, -1, 0, 0),
    "1/G2": (0, 0, 0, 0, 0, -1, 0),
    "x7/(G0*G3)": (0, 0, 1, -1, 0, 0, -1),
    "1/(G1*G3)": (0, 0, 0, 0, -1, 0, -1),
}

VERDICT_FOUND = "FiniteNumeratorCandidatesFound"
VERDICT_NONE = "NoFiniteNumeratorBasisWithinAnsatz"

RECOMMENDATION = (
    "Numerator decoration cannot cure the stable remnants 1/G1, 1/G2 and the "
    "probe sector 1/(G1*G3): every failing ray is componentwise <= 0 (x -> oo "
    "type) and polynomial numerators only increase those layer scores "
    "(Lemma 2, numerator_cure_impossible_any_degree). Next directions: label "
    "shifts with POSITIVE polynomial components (raised denominator powers / "
    "dimension-shift-type masters), the Method.3 composite basis change, or "
    "an analytic treatment of the remnants outside the reducer."
)


def numerator_extended_labels(payload: dict) -> list[list[int]]:
    """Defining-expansion labels of accepted candidates in non-bare-LF sectors.

    Each such label is individually LF (Lemma 1), so the Method.1 span test may
    add it to the allowed set; the defining rows are exactly the monomial
    expansions ``M - sum_alpha c_alpha * J(sector + (alpha|0)) = 0``.
    """
    labels: set[tuple[int, ...]] = set()
    for rep in payload["reports"]:
        if rep["status"] == STATUS_ALREADY_LF:
            continue
        for cand in rep["candidates"]:
            for term in cand["defining_expansion"]:
                labels.add(tuple(term["label"]))
    return [list(lab) for lab in sorted(labels)]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--json", type=Path, default=DEFAULT_JSON)
    parser.add_argument(
        "--degrees",
        type=str,
        default="0,1,2",
        help="comma-separated numerator total degrees (default: 0,1,2)",
    )
    parser.add_argument("--random-trials", type=int, default=64)
    parser.add_argument("--seed", type=int, default=20260706)
    args = parser.parse_args(argv)

    if not (1 <= args.random_trials <= MAX_RANDOM_TRIALS):
        print(f"--random-trials must be in [1, {MAX_RANDOM_TRIALS}]", file=sys.stderr)
        return 2
    try:
        degrees = tuple(int(d) for d in args.degrees.split(","))
    except ValueError:
        print(f"cannot parse --degrees {args.degrees!r}", file=sys.stderr)
        return 2
    if not degrees or any(d < 0 for d in degrees):
        print("--degrees must be non-negative integers", file=sys.stderr)
        return 2

    try:
        family = parse_family_text(args.input.read_text(encoding="utf-8"))
    except (OSError, ParserError) as exc:
        print(f"cannot load family: {exc}", file=sys.stderr)
        return 2

    t0 = time.time()
    payload = finite_numerator_search(
        family, SECTORS, degrees, args.random_trials, args.seed
    )
    elapsed = time.time() - t0

    found = payload["new_lf_masters_found"]
    bridge_labels = numerator_extended_labels(payload)
    payload.update(
        {
            "input": str(args.input),
            "sector_labels_are_shifts": True,
            "base_exponents_note": (
                "offset convention: total exponent = base exponent + label "
                "shift; Int2 base is x2^(1+ep)*G0^ep*G1^ep*G2^(-1-ep)*G3^(-1+ep)"
            ),
            "random_trials": args.random_trials,
            "seed": args.seed,
            "elapsed_seconds": round(elapsed, 3),
            "overall_status": VERDICT_FOUND if found else VERDICT_NONE,
            "numerator_extended_labels": bridge_labels,
            "feasibility": (
                {
                    "status": "BridgeLabelsExported",
                    "note": (
                        "each accepted candidate expands into individually-LF "
                        "monomial labels (Lemma 1); feed these labels as "
                        "allowed into lf_reduction_feasible_mod_p with the "
                        "defining rows to test span feasibility"
                    ),
                }
                if found
                else {
                    "status": "SkippedNoCandidates",
                    "note": (
                        "no genuinely polynomial LF candidate in any searched "
                        "sector/degree; the Method.1 modular span test has "
                        "nothing to extend"
                    ),
                }
            ),
        }
    )
    if not found:
        payload["recommendation"] = RECOMMENDATION

    args.json.parent.mkdir(parents=True, exist_ok=True)
    args.json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    by_name: dict[str, list[str]] = {}
    for rep in payload["reports"]:
        by_name.setdefault(rep["sector_name"], []).append(rep["status"])
    print(f"sectors={len(SECTORS)} degrees={list(degrees)} elapsed={elapsed:.1f}s")
    for name, statuses in by_name.items():
        print(f"  {name:14s} {sorted(set(statuses))}")
    print(f"overall: {payload['overall_status']}")
    print(f"lemma_consistent_everywhere: {payload['lemma_consistent_everywhere']}")
    print(f"artifact: {args.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
