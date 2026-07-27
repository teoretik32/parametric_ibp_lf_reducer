"""Run a prepared External Int2 LF master through HyperInt (Maple).

Method.12A executor.  This script only *drives* the already-prepared and
reviewed inputs in ``validation/hyperint/``.  It does not derive, re-derive or
re-certify anything: the four LF labels, their coefficients, the surface policy
and the certificate are inputs, never outputs, of this script.

Deliberately stdlib-only.  It must never import RREF, modular-record or
normal-form machinery -- see ``tests/test_external_int2_hyperint_setup.py``.

Usage::

    python scripts/run_external_int2_hyperint.py --master L1 --dry-run
    python scripts/run_external_int2_hyperint.py --master L1
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
INPUT_DIR = REPO_ROOT / "validation" / "hyperint"
RESULT_DIR = INPUT_DIR / "results"
LOG_DIR = REPO_ROOT / "outputs" / "hyperint"

MASTERS = ("L1", "L2", "L3", "L4")

#: Windows subdirectory of a Maple install that holds the command-line kernel.
MAPLE_BIN_SUBDIRS = ("bin.X86_64_WINDOWS", "bin.X86_64_LINUX", "bin.APPLE_UNIVERSAL_OSX")
CMAPLE_NAMES = ("cmaple.exe", "cmaple")

#: Files that must be present for a directory to count as a HyperInt install.
HYPERINT_REQUIRED = ("HyperInt.mpl",)
HYPERINT_OPTIONAL = ("periodLookups.m", "HyperTests.mpl")


class SetupError(RuntimeError):
    """Environment or input could not be resolved."""


# --------------------------------------------------------------------------
# environment resolution
# --------------------------------------------------------------------------


def resolve_cmaple(env: dict[str, str] | None = None) -> Path:
    """Locate the Maple command-line kernel: MAPLE_CLI, then MAPLE_HOME, then PATH."""
    env = os.environ if env is None else env

    cli = env.get("MAPLE_CLI", "").strip()
    if cli:
        path = Path(cli)
        if not path.is_file():
            raise SetupError(f"MAPLE_CLI={cli!r} does not point at an existing file")
        return path

    home = env.get("MAPLE_HOME", "").strip()
    if home:
        root = Path(home)
        for sub in MAPLE_BIN_SUBDIRS:
            for name in CMAPLE_NAMES:
                cand = root / sub / name
                if cand.is_file():
                    return cand
        raise SetupError(f"MAPLE_HOME={home!r} contains no cmaple under {MAPLE_BIN_SUBDIRS}")

    for name in CMAPLE_NAMES:
        found = shutil.which(name, path=env.get("PATH"))
        if found:
            return Path(found)

    raise SetupError(
        "cannot locate cmaple; set MAPLE_CLI to the executable, MAPLE_HOME to the "
        "Maple install root, or put cmaple on PATH"
    )


def resolve_hyperint(env: dict[str, str] | None = None) -> Path:
    """Locate the HyperInt source tree via HYPERINT_HOME."""
    env = os.environ if env is None else env

    home = env.get("HYPERINT_HOME", "").strip()
    if not home:
        raise SetupError("HYPERINT_HOME is not set; point it at the HyperInt source tree")
    root = Path(home)
    if not root.is_dir():
        raise SetupError(f"HYPERINT_HOME={home!r} is not a directory")
    missing = [n for n in HYPERINT_REQUIRED if not (root / n).is_file()]
    if missing:
        raise SetupError(f"HYPERINT_HOME={home!r} is missing {', '.join(missing)}")
    return root


# --------------------------------------------------------------------------
# deterministic paths
# --------------------------------------------------------------------------


def input_path(master: str) -> Path:
    return INPUT_DIR / f"external_int2_lf_{master}.mpl"


def result_path(master: str) -> Path:
    return RESULT_DIR / f"external_int2_lf_{master}_result.mpl"


def meta_path(master: str) -> Path:
    return RESULT_DIR / f"external_int2_lf_{master}_meta.json"


def driver_path(master: str) -> Path:
    return LOG_DIR / master / f"driver_{master}.mpl"


def log_path(master: str) -> Path:
    return LOG_DIR / master / f"{master}_run.log"


def is_complete(master: str) -> bool:
    """True when a finished result for ``master`` is already on disk."""
    meta = meta_path(master)
    if not (result_path(master).is_file() and meta.is_file()):
        return False
    try:
        return json.loads(meta.read_text(encoding="utf-8")).get("status") == "complete"
    except (OSError, json.JSONDecodeError):
        return False


# --------------------------------------------------------------------------
# driver generation
# --------------------------------------------------------------------------


def _maple_str(path: Path | str) -> str:
    """Maple string literal for a path (forward slashes, no escaping needed)."""
    return '"' + str(path).replace("\\", "/") + '"'


def build_driver(master: str, hyperint_home: Path) -> str:
    """Maple source that runs the reviewed master input and saves its result.

    ``currentdir`` is set to the HyperInt tree so that the prepared input's own
    ``read "HyperInt.mpl"`` and the ``periodLookups.m`` autoload resolve without
    editing the certified input file.  The integration call is appended here --
    the prepared file leaves it commented out behind the human-review gate --
    and uses exactly the ``fser`` and ``intOrder`` that file defines.
    """
    if master not in MASTERS:
        raise SetupError(f"unknown master {master!r}; expected one of {', '.join(MASTERS)}")
    return "\n".join(
        [
            f"# External Int2 Method.12A HyperInt driver for master {master}",
            "# Generated by scripts/run_external_int2_hyperint.py -- do not edit by hand.",
            f"_hyper_autoload_periods := [{_maple_str(hyperint_home / 'periodLookups.m')}]:",
            "_hyper_verbosity := 1:",
            f"currentdir({_maple_str(hyperint_home)}):",
            f"read {_maple_str(input_path(master))}:",
            "st := time[real]():",
            "result := hyperInt(fser, intOrder):",
            "fb := fibrationBasis(result):",
            "elapsed := time[real]() - st:",
            f'printf("HYPERINT_MASTER={master}\\n");',
            'printf("HYPERINT_EPSORDER=%a\\n", epsOrder);',
            'printf("HYPERINT_INTORDER=%a\\n", intOrder);',
            'printf("HYPERINT_ELAPSED_SECONDS=%.3f\\n", elapsed);',
            f"writeto({_maple_str(result_path(master))}):",
            'printf("# External Int2 Method.12A HyperInt result, master '
            f'{master}\\n");',
            'printf("%a;\\n", fb);',
            "writeto(terminal):",
            'printf("HYPERINT_STATUS=complete\\n");',
            "quit;",
            "",
        ]
    )


def build_command(cmaple: Path, master: str) -> list[str]:
    """The exact argv used to run ``master``."""
    return [str(cmaple), "-q", str(driver_path(master))]


# --------------------------------------------------------------------------
# run
# --------------------------------------------------------------------------


def _plan(master: str, cmaple: Path, hyperint_home: Path) -> dict:
    return {
        "master": master,
        "cmaple": str(cmaple),
        "hyperint_home": str(hyperint_home),
        "input": str(input_path(master)),
        "driver": str(driver_path(master)),
        "result": str(result_path(master)),
        "meta": str(meta_path(master)),
        "log": str(log_path(master)),
        "command": build_command(cmaple, master),
        "already_complete": is_complete(master),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="run_external_int2_hyperint",
        description="Run a prepared External Int2 LF master through HyperInt.",
    )
    parser.add_argument("--master", required=True, choices=MASTERS)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="resolve the environment and print the plan; touch nothing on disk",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="overwrite an already-completed result for this master",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=None,
        help="wall-clock limit in seconds for the Maple process (default: none)",
    )
    args = parser.parse_args(argv)

    master = args.master
    try:
        cmaple = resolve_cmaple()
        hyperint_home = resolve_hyperint()
    except SetupError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    src = input_path(master)
    if not src.is_file():
        print(f"ERROR: prepared input {src} not found", file=sys.stderr)
        return 2

    plan = _plan(master, cmaple, hyperint_home)

    if args.dry_run:
        print(json.dumps(plan, indent=2))
        print(f"\nDRY RUN: would execute:\n  {' '.join(plan['command'])}")
        return 0

    if plan["already_complete"] and not args.force:
        print(
            f"ERROR: {master} already has a completed result at {result_path(master)};\n"
            "       pass --force to overwrite it.",
            file=sys.stderr,
        )
        return 3

    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    driver_path(master).parent.mkdir(parents=True, exist_ok=True)
    driver_path(master).write_text(build_driver(master, hyperint_home), encoding="utf-8")

    started = _dt.datetime.now(_dt.timezone.utc)
    print(f"running {master}: {' '.join(plan['command'])}")
    try:
        proc = subprocess.run(  # noqa: S603 - argv is fully constructed above
            plan["command"],
            capture_output=True,
            text=True,
            timeout=args.timeout,
            check=False,
        )
        stdout, stderr, returncode, timed_out = proc.stdout, proc.stderr, proc.returncode, False
    except subprocess.TimeoutExpired as exc:
        stdout = exc.stdout.decode() if isinstance(exc.stdout, bytes) else (exc.stdout or "")
        stderr = exc.stderr.decode() if isinstance(exc.stderr, bytes) else (exc.stderr or "")
        returncode, timed_out = -1, True
    finished = _dt.datetime.now(_dt.timezone.utc)

    log_path(master).write_text(
        f"# command: {' '.join(plan['command'])}\n"
        f"# started: {started.isoformat()}\n"
        f"# finished: {finished.isoformat()}\n"
        f"# returncode: {returncode}\n"
        f"# timed_out: {timed_out}\n"
        "--- stdout ---\n" + stdout + "\n--- stderr ---\n" + stderr + "\n",
        encoding="utf-8",
    )

    ok = (
        not timed_out
        and returncode == 0
        and "HYPERINT_STATUS=complete" in stdout
        and result_path(master).is_file()
    )
    meta = {
        "master": master,
        "status": "complete" if ok else "failed",
        "returncode": returncode,
        "timed_out": timed_out,
        "cmaple": str(cmaple),
        "hyperint_home": str(hyperint_home),
        "input": str(input_path(master)),
        "result": str(result_path(master)),
        "log": str(log_path(master)),
        "started_utc": started.isoformat(),
        "finished_utc": finished.isoformat(),
        "wall_seconds": round((finished - started).total_seconds(), 3),
    }
    meta_path(master).write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")

    if not ok:
        print(f"FAILED: {master}; see {log_path(master)}", file=sys.stderr)
        return 1
    print(f"OK: {master} -> {result_path(master)}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
