"""Fast gates for the External Int2 HyperInt executor (Method.12A setup).

No Maple and no HyperInt install is required: every test here works on
resolution logic, generated command strings and repository hygiene.
"""

from __future__ import annotations

import ast
import importlib.util
import json
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
RUNNER = REPO_ROOT / "scripts" / "run_external_int2_hyperint.py"
SMOKE_JSON = REPO_ROOT / "validation" / "hyperint" / "hyperint_smoke.json"
SETUP_NOTES = REPO_ROOT / "notes" / "EXTERNAL_INT2_HYPERINT_SETUP.md"


def _load_runner():
    spec = importlib.util.spec_from_file_location("_int2_hyperint_runner", RUNNER)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


runner = _load_runner()


# ---------------------------------------------------------------- environment


def test_cmaple_from_maple_cli(tmp_path: Path) -> None:
    exe = tmp_path / "cmaple.exe"
    exe.write_text("", encoding="utf-8")
    assert runner.resolve_cmaple({"MAPLE_CLI": str(exe)}) == exe


def test_cmaple_from_maple_home(tmp_path: Path) -> None:
    bindir = tmp_path / "bin.X86_64_WINDOWS"
    bindir.mkdir()
    exe = bindir / "cmaple.exe"
    exe.write_text("", encoding="utf-8")
    assert runner.resolve_cmaple({"MAPLE_HOME": str(tmp_path)}) == exe


def test_cmaple_from_path(tmp_path: Path) -> None:
    exe = tmp_path / "cmaple.exe"
    exe.write_text("", encoding="utf-8")
    found = runner.resolve_cmaple({"PATH": str(tmp_path)})
    assert found.name.startswith("cmaple")


def test_maple_cli_wins_over_maple_home(tmp_path: Path) -> None:
    """An explicit MAPLE_CLI must not be shadowed by a MAPLE_HOME guess."""
    bindir = tmp_path / "home" / "bin.X86_64_WINDOWS"
    bindir.mkdir(parents=True)
    (bindir / "cmaple.exe").write_text("", encoding="utf-8")
    explicit = tmp_path / "cmaple.exe"
    explicit.write_text("", encoding="utf-8")
    resolved = runner.resolve_cmaple(
        {"MAPLE_CLI": str(explicit), "MAPLE_HOME": str(tmp_path / "home")}
    )
    assert resolved == explicit


def test_cmaple_missing_raises() -> None:
    with pytest.raises(runner.SetupError):
        runner.resolve_cmaple({"PATH": ""})


def test_bad_maple_cli_raises(tmp_path: Path) -> None:
    with pytest.raises(runner.SetupError):
        runner.resolve_cmaple({"MAPLE_CLI": str(tmp_path / "nope.exe")})


def test_hyperint_home_resolution(tmp_path: Path) -> None:
    (tmp_path / "HyperInt.mpl").write_text("", encoding="utf-8")
    assert runner.resolve_hyperint({"HYPERINT_HOME": str(tmp_path)}) == tmp_path


def test_hyperint_home_unset_raises() -> None:
    with pytest.raises(runner.SetupError):
        runner.resolve_hyperint({})


def test_hyperint_home_without_source_raises(tmp_path: Path) -> None:
    with pytest.raises(runner.SetupError):
        runner.resolve_hyperint({"HYPERINT_HOME": str(tmp_path)})


# ------------------------------------------------------- deterministic paths


@pytest.mark.parametrize("master", runner.MASTERS)
def test_paths_are_deterministic(master: str) -> None:
    for fn in (runner.input_path, runner.result_path, runner.meta_path,
               runner.driver_path, runner.log_path):
        assert fn(master) == fn(master)
        assert master in fn(master).name


@pytest.mark.parametrize("master", runner.MASTERS)
def test_prepared_input_exists(master: str) -> None:
    assert runner.input_path(master).is_file()


def test_paths_are_distinct_per_master() -> None:
    results = {runner.result_path(m) for m in runner.MASTERS}
    assert len(results) == len(runner.MASTERS)


def test_logs_go_to_outputs_and_are_gitignored() -> None:
    log = runner.log_path("L1")
    assert (REPO_ROOT / "outputs") in log.parents
    assert "outputs/" in (REPO_ROOT / ".gitignore").read_text(encoding="utf-8").split()


# ------------------------------------------------------------ dry-run / cmd


def test_build_command_shape(tmp_path: Path) -> None:
    cmaple = tmp_path / "cmaple.exe"
    cmd = runner.build_command(cmaple, "L1")
    assert cmd[0] == str(cmaple)
    assert cmd[1] == "-q"
    assert cmd[2] == str(runner.driver_path("L1"))


def test_driver_references_absolute_paths(tmp_path: Path) -> None:
    src = runner.build_driver("L1", Path("K:/_TOOLS/HyperInt"))
    assert "K:/_TOOLS/HyperInt/periodLookups.m" in src
    assert 'currentdir("K:/_TOOLS/HyperInt")' in src
    assert str(runner.input_path("L1")).replace("\\", "/") in src
    assert "hyperInt(fser, intOrder)" in src
    assert "fibrationBasis(result)" in src


def test_build_driver_rejects_unknown_master() -> None:
    with pytest.raises(runner.SetupError):
        runner.build_driver("L9", Path("K:/_TOOLS/HyperInt"))


def test_dry_run_touches_nothing(tmp_path: Path, monkeypatch, capsys) -> None:
    cmaple = tmp_path / "cmaple.exe"
    cmaple.write_text("", encoding="utf-8")
    hyper = tmp_path / "HyperInt"
    hyper.mkdir()
    (hyper / "HyperInt.mpl").write_text("", encoding="utf-8")
    monkeypatch.setenv("MAPLE_CLI", str(cmaple))
    monkeypatch.setenv("HYPERINT_HOME", str(hyper))

    before = runner.driver_path("L1").exists()
    assert runner.main(["--master", "L1", "--dry-run"]) == 0
    assert runner.driver_path("L1").exists() == before

    plan = json.loads(capsys.readouterr().out.split("\nDRY RUN")[0])
    assert plan["master"] == "L1"
    assert plan["command"][0] == str(cmaple)


def test_dry_run_fails_cleanly_without_environment(monkeypatch, capsys) -> None:
    monkeypatch.delenv("MAPLE_CLI", raising=False)
    monkeypatch.delenv("MAPLE_HOME", raising=False)
    monkeypatch.setenv("PATH", "")
    assert runner.main(["--master", "L1", "--dry-run"]) == 2


def test_completed_result_is_not_overwritten(tmp_path: Path, monkeypatch) -> None:
    """is_complete drives the refuse-to-overwrite guard."""
    monkeypatch.setattr(runner, "RESULT_DIR", tmp_path)
    monkeypatch.setattr(runner, "result_path", lambda m: tmp_path / f"{m}_result.mpl")
    monkeypatch.setattr(runner, "meta_path", lambda m: tmp_path / f"{m}_meta.json")
    assert runner.is_complete("L1") is False
    (tmp_path / "L1_result.mpl").write_text("1;\n", encoding="utf-8")
    (tmp_path / "L1_meta.json").write_text('{"status": "complete"}', encoding="utf-8")
    assert runner.is_complete("L1") is True
    (tmp_path / "L1_meta.json").write_text('{"status": "failed"}', encoding="utf-8")
    assert runner.is_complete("L1") is False


def test_force_flag_exists() -> None:
    assert "--force" in RUNNER.read_text(encoding="utf-8")


# ------------------------------------------------------------------ hygiene


FORBIDDEN_IMPORT_SUBSTRINGS = ("rref", "modular", "records", "normal_form", "finite_field")


def test_runner_imports_no_rref_or_modular_code() -> None:
    tree = ast.parse(RUNNER.read_text(encoding="utf-8"))
    names: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names += [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            names.append(node.module or "")
            names += [a.name for a in node.names]
    lowered = [n.lower() for n in names]
    for bad in FORBIDDEN_IMPORT_SUBSTRINGS:
        assert not any(bad in n for n in lowered), f"runner imports {bad!r}: {names}"
    assert not any(n.startswith("parametric_ibp_lf_reducer") for n in names)


def test_runner_is_stdlib_only() -> None:
    tree = ast.parse(RUNNER.read_text(encoding="utf-8"))
    roots = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            roots |= {a.name.split(".")[0] for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            roots.add(node.module.split(".")[0])
    assert not roots & {"sympy", "numpy", "numba"}


def test_no_hyperint_source_committed() -> None:
    """HyperInt lives outside the repo; its sources must not be tracked here."""
    for name in ("HyperInt.mpl", "periodLookups.m", "HyperTests.mpl", "Manual.mw"):
        assert not list(REPO_ROOT.rglob(name)), f"{name} must not live in the repository"


def test_no_license_files_committed() -> None:
    for pattern in ("license.dat", "*.lic", "maple.lic", "licenseserver*"):
        hits = [p for p in REPO_ROOT.rglob(pattern) if ".git" not in p.parts]
        assert not hits, f"license artifact committed: {hits}"


def test_smoke_json_is_well_formed() -> None:
    data = json.loads(SMOKE_JSON.read_text(encoding="utf-8"))
    assert "2025" in data["maple"]["version"]
    assert data["maple"]["smoke_result"] == "MAPLE_SMOKE=2"
    assert data["hyperint"]["smoke_expression"] == "hyperInt(1/(1+x)^2, x=0..infinity)"
    assert data["hyperint"]["smoke_result"] == data["hyperint"]["smoke_expected"] == "1"
    assert len(data["hyperint"]["revision"]) == 40
    assert data["verdict"] == "PASS"
    assert data["timestamp_utc"].startswith("2026-")
    assert data["driver_probe"]["integration_executed"] is False


def test_setup_notes_exist_and_record_provenance() -> None:
    text = SETUP_NOTES.read_text(encoding="utf-8")
    rev = json.loads(SMOKE_JSON.read_text(encoding="utf-8"))["hyperint"]["revision"]
    assert rev in text
    assert "bitbucket.org/PanzerErik/hyperint" in text


def test_certified_inputs_untouched_by_this_session() -> None:
    """The prepared master inputs still carry their review gate, uncommented."""
    for master in runner.MASTERS:
        text = runner.input_path(master).read_text(encoding="utf-8")
        assert "# result := hyperInt(" in text
        assert "intOrder := [x2, x5, x7]" in text
