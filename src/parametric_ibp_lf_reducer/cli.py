"""Thin CLI wrapper over :func:`reduce_wolfram_style_input` (CLI.1).

Pure plumbing: read a Wolfram-like explicit-family document from a file, run one
reduction pass through the ordinary public API (certificate gate fully intact),
write the Wolfram-like result text and, optionally, machine-readable JSON
diagnostics. No math lives here and no Mathematica/Wolfram runtime is invoked.

Exit codes:
    0  reduction ``Status -> "Success"``
    1  reduction ran but did not reach ``Success`` (concrete reason in output/JSON)
    2  usage / I/O / malformed-document errors (nothing was reduced)
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .input_parser import ParserError
from .result import ReductionResult
from .sparse_rref import RREF_BACKEND_CHOICES
from . import api

EXIT_SUCCESS = 0
EXIT_FAILURE = 1
EXIT_USAGE = 2


def _diagnostics_payload(result: ReductionResult) -> dict:
    """JSON-safe snapshot of the result; labels become lists, tribools stay as-is."""
    d = result.diagnostics
    cert = d.extra.get("certificate") or {}
    return {
        "status": result.status,
        "exported_status": result.exported_status,
        "success": result.success,
        "error": result.error,
        # Row-span certificate verdict (gate is default-ON; "NotRun" only when the
        # reduction never reached the certificate stage).
        "certificate_status": cert.get("certificate_status"),
        "certificate": {k: v for k, v in cert.items() if isinstance(v, (str, int, type(None)))},
        "target_label": list(result.target_label),
        "all_locally_finite": result.all_locally_finite,
        "terms": [
            {
                "label": list(t.label),
                "coefficient": t.coefficient_text,
                "integrand": t.integrand_text,
                "locally_finite": t.locally_finite,
            }
            for t in result.terms
        ],
        "diagnostics": {
            "formal_success": bool(d.formal_success),
            "reconstruction_verified": bool(d.reconstruction_verified),
            "independent_validation_passed": bool(d.independent_validation_passed),
            "n_terms": d.n_terms,
            "non_lf_terms": [list(x) for x in d.non_lf_terms],
            "unknown_lf_terms": [list(x) for x in d.unknown_lf_terms],
            "n_records": d.n_records,
            "n_skipped_records": d.n_skipped_records,
            "zero_reduction": bool(d.zero_reduction),
            "messages": list(d.messages),
            # Perf.0: stage -> seconds wall-clock snapshot (observability only).
            "timings": {k: float(v) for k, v in (d.extra.get("timings") or {}).items()},
        },
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m parametric_ibp_lf_reducer",
        description="Locally-finite parametric IBP reducer (pure Python, no Wolfram runtime).",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    reduce_p = sub.add_parser(
        "reduce",
        help="run one reduction pass on a Wolfram-like explicit-family document",
    )
    reduce_p.add_argument("input", help="path to the Wolfram-like input document (e.g. input.m)")
    reduce_p.add_argument(
        "--out",
        metavar="PATH",
        help="write the Wolfram-like result text here (default: stdout)",
    )
    reduce_p.add_argument(
        "--diagnostics-json",
        metavar="PATH",
        help="also write machine-readable JSON diagnostics to this path",
    )
    reduce_p.add_argument(
        "--max-ibp-degree",
        type=int,
        metavar="N",
        help="override MaxIBPDegree (row-generation degree cap)",
    )
    reduce_p.add_argument(
        "--min-valid-records",
        type=int,
        metavar="N",
        help="override MinValidRecords (reconstruction evidence floor)",
    )
    reduce_p.add_argument(
        "--jobs",
        type=int,
        metavar="N",
        help="worker processes for (prime, sample) record collection (default: 1 = serial)",
    )
    reduce_p.add_argument(
        "--rref-backend",
        choices=RREF_BACKEND_CHOICES,
        metavar="NAME",
        help=(
            "RREF implementation for records + certificate points "
            f"(one of {', '.join(RREF_BACKEND_CHOICES)}; default: dict). "
            "'auto' (experimental, Perf.12) picks numba per matrix for large systems "
            "when numba is available and prime < 2**31, else dict. "
            "Backend selection only — all backends return identical results."
        ),
    )
    return parser


def _cmd_reduce(args: argparse.Namespace) -> int:
    try:
        input_text = Path(args.input).read_text(encoding="utf-8")
    except OSError as exc:
        print(f"error: cannot read input file: {exc}", file=sys.stderr)
        return EXIT_USAGE

    overrides = {}
    if args.max_ibp_degree is not None:
        overrides["max_ibp_degree"] = args.max_ibp_degree
    if args.min_valid_records is not None:
        overrides["min_valid_records"] = args.min_valid_records
    if args.jobs is not None:
        overrides["jobs"] = args.jobs
    if args.rref_backend is not None:
        overrides["rref_backend"] = args.rref_backend

    try:
        result = api.reduce_wolfram_style_input(input_text, overrides)
    except ParserError as exc:
        # Malformed document: nothing was reduced, so this is a usage error, not a
        # reduction failure. (ParserNeedsExplicitFamily is NOT raised here — the API
        # returns it as an honest typed failure and it goes through the normal path.)
        print(f"error: malformed input document: {exc}", file=sys.stderr)
        return EXIT_USAGE

    text = result.wolfram_style_text
    if args.out:
        Path(args.out).write_text(text + "\n", encoding="utf-8")
    else:
        print(text)

    if args.diagnostics_json:
        payload = _diagnostics_payload(result)
        Path(args.diagnostics_json).write_text(
            json.dumps(payload, indent=2) + "\n", encoding="utf-8"
        )

    if result.success:
        return EXIT_SUCCESS
    print(f"reduction failed: {result.failure_reason}", file=sys.stderr)
    if result.error:
        print(f"detail: {result.error}", file=sys.stderr)
    return EXIT_FAILURE


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "reduce":
        return _cmd_reduce(args)
    raise AssertionError(f"unhandled command {args.command!r}")  # pragma: no cover
