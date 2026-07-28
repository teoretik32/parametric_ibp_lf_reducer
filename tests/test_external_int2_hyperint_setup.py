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


# ------------------------------------------------------- basis-status gate
#
# Method.12R revoked the Method.11c four-master basis.  No real master run may
# start unless a basis-status artifact confirms all four required fields.

REVOKED_ARTIFACT = """\
(* REVOKED by External Int2 Method.12R (mixed-boundary contradiction audit). *)
<|
  "Status" -> "Revoked(Method.12R)",
  "AllLocallyFinite" -> False,
  "Diagnostics" -> <| "NonLFTerms" -> {{0,0,1,-1,0,0,-1}} |>
|>
Method12RInvalidation = <| "Method12RStatus" -> "Invalidated" |>;
"""

LEGACY_ARTIFACT = """\
<|
  "Status" -> "Success",
  "AllLocallyFinite" -> True,
  "Diagnostics" -> <| "NonLFTerms" -> {} |>
|>
"""

VALID_ARTIFACT = """\
<|
  "Status" -> "Success",
  "AllLocallyFinite" -> True,
  "IntegralIdentityStatus" -> "Verified",
  "SurfaceValidationStatus" -> "Passed"
|>
"""


def _write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def test_gate_refuses_revoked_basis(tmp_path: Path) -> None:
    allowed, report = runner.check_basis_status(_write(tmp_path, "revoked.m", REVOKED_ARTIFACT))
    assert allowed is False
    assert report["fields"]["Status"] == "Revoked(Method.12R)"
    assert report["fields"]["AllLocallyFinite"] == "False"
    assert any("revocation marker" in f for f in report["failures"])


def test_gate_refuses_legacy_unchecked_basis(tmp_path: Path) -> None:
    """Status=Success + AllLocallyFinite=True is NOT enough: the pre-12R case."""
    allowed, report = runner.check_basis_status(_write(tmp_path, "legacy.m", LEGACY_ARTIFACT))
    assert allowed is False
    assert report["fields"]["IntegralIdentityStatus"] == runner.MISSING
    assert report["fields"]["SurfaceValidationStatus"] == runner.MISSING


def test_gate_refuses_missing_artifact(tmp_path: Path) -> None:
    allowed, report = runner.check_basis_status(tmp_path / "absent.m")
    assert allowed is False
    assert any("not found" in f for f in report["failures"])


@pytest.mark.parametrize(
    "field,bad",
    [
        ("Status", '"Status" -> "Partial"'),
        ("AllLocallyFinite", '"AllLocallyFinite" -> False'),
        ("IntegralIdentityStatus", '"IntegralIdentityStatus" -> "Unknown"'),
        ("SurfaceValidationStatus", '"SurfaceValidationStatus" -> "Failed"'),
    ],
)
def test_gate_requires_every_field(tmp_path: Path, field: str, bad: str) -> None:
    """Breaking any single required field alone closes the gate."""
    original = [ln for ln in VALID_ARTIFACT.splitlines() if f'"{field}"' in ln][0]
    text = VALID_ARTIFACT.replace(original.strip().rstrip(","), bad)
    allowed, _ = runner.check_basis_status(_write(tmp_path, f"{field}.m", text))
    assert allowed is False


def test_gate_accepts_valid_synthetic_provenance(tmp_path: Path) -> None:
    allowed, report = runner.check_basis_status(_write(tmp_path, "ok.m", VALID_ARTIFACT))
    assert allowed is True, report["failures"]
    assert report["failures"] == []


def test_gate_accepts_valid_json_provenance(tmp_path: Path) -> None:
    artifact = _write(
        tmp_path,
        "ok.json",
        json.dumps(
            {
                "Status": "Success",
                "AllLocallyFinite": True,
                "IntegralIdentityStatus": "Valid",
                "SurfaceValidationStatus": "Passed",
            }
        ),
    )
    allowed, report = runner.check_basis_status(artifact)
    assert allowed is True, report["failures"]


def test_committed_basis_artifact_is_currently_blocked() -> None:
    """The repository's own artifact must not open the gate today."""
    allowed, report = runner.check_basis_status()
    assert allowed is False
    assert report["failures"]


def _env(monkeypatch, tmp_path: Path) -> None:
    cmaple = tmp_path / "cmaple.exe"
    cmaple.write_text("", encoding="utf-8")
    hyper = tmp_path / "HyperInt"
    hyper.mkdir(exist_ok=True)
    (hyper / "HyperInt.mpl").write_text("", encoding="utf-8")
    monkeypatch.setenv("MAPLE_CLI", str(cmaple))
    monkeypatch.setenv("HYPERINT_HOME", str(hyper))


def test_real_run_refused_before_subprocess_launch(tmp_path, monkeypatch, capsys) -> None:
    """A revoked basis must stop the run *before* cmaple is ever invoked."""
    _env(monkeypatch, tmp_path)
    artifact = _write(tmp_path, "revoked.m", REVOKED_ARTIFACT)

    def _boom(*args, **kwargs):  # pragma: no cover - must never run
        raise AssertionError("subprocess.run was called despite a revoked basis")

    monkeypatch.setattr(runner.subprocess, "run", _boom)
    rc = runner.main(["--master", "L1", "--basis-status", str(artifact)])
    assert rc != 0
    err = capsys.readouterr().err
    assert "REFUSED" in err
    assert "No cmaple process was started." in err
    assert not runner.driver_path("L1").exists()


def test_real_run_refused_when_provenance_missing(tmp_path, monkeypatch, capsys) -> None:
    _env(monkeypatch, tmp_path)

    def _boom(*args, **kwargs):  # pragma: no cover - must never run
        raise AssertionError("subprocess.run was called with no provenance artifact")

    monkeypatch.setattr(runner.subprocess, "run", _boom)
    rc = runner.main(["--master", "L1", "--basis-status", str(tmp_path / "absent.m")])
    assert rc != 0
    assert "REFUSED" in capsys.readouterr().err


def test_force_does_not_bypass_the_gate(tmp_path, monkeypatch) -> None:
    _env(monkeypatch, tmp_path)
    artifact = _write(tmp_path, "revoked.m", REVOKED_ARTIFACT)

    def _boom(*args, **kwargs):  # pragma: no cover - must never run
        raise AssertionError("--force bypassed the basis-status gate")

    monkeypatch.setattr(runner.subprocess, "run", _boom)
    assert runner.main(["--master", "L1", "--force", "--basis-status", str(artifact)]) != 0


def test_no_override_flag_exists() -> None:
    """No CLI switch may skip the gate."""
    text = RUNNER.read_text(encoding="utf-8").lower()
    for flag in ("--skip-basis", "--no-gate", "--ignore-revoked", "--allow-revoked",
                 "--bypass", "--i-know-what-i-am-doing"):
        assert flag not in text


def test_valid_provenance_lets_the_command_be_prepared(tmp_path, monkeypatch) -> None:
    """With a valid artifact the gate opens and the driver/command are built."""
    _env(monkeypatch, tmp_path)
    artifact = _write(tmp_path, "ok.m", VALID_ARTIFACT)
    calls: list[list[str]] = []

    class _Proc:
        stdout = "HYPERINT_STATUS=complete\n"
        stderr = ""
        returncode = 0

    def _fake_run(cmd, **kwargs):
        calls.append(cmd)
        return _Proc()

    monkeypatch.setattr(runner.subprocess, "run", _fake_run)
    monkeypatch.setattr(runner, "RESULT_DIR", tmp_path / "res")
    monkeypatch.setattr(runner, "result_path", lambda m: tmp_path / "res" / f"{m}.mpl")
    monkeypatch.setattr(runner, "meta_path", lambda m: tmp_path / "res" / f"{m}.json")
    monkeypatch.setattr(runner, "driver_path", lambda m: tmp_path / "drv" / f"{m}.mpl")
    monkeypatch.setattr(runner, "log_path", lambda m: tmp_path / "drv" / f"{m}.log")
    (tmp_path / "res").mkdir()
    (tmp_path / "res" / "L1.mpl").write_text("1;\n", encoding="utf-8")

    assert runner.main(["--master", "L1", "--basis-status", str(artifact)]) == 0
    assert len(calls) == 1
    assert calls[0][1] == "-q"
    assert (tmp_path / "drv" / "L1.mpl").is_file()


def test_refusal_message_names_every_failed_field(tmp_path: Path) -> None:
    _, report = runner.check_basis_status(_write(tmp_path, "legacy.m", LEGACY_ARTIFACT))
    msg = runner.format_refusal("L1", report)
    for key in ("Status", "AllLocallyFinite", "IntegralIdentityStatus",
                "SurfaceValidationStatus"):
        assert key in msg
    assert "REFUSED" in msg
    assert "Method.12R" in msg
    assert "no override flag" in msg
    assert "No cmaple process was started." in msg


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


def test_hyperint_smoke_metadata_unchanged() -> None:
    """Method.12R invalidates the basis, not the runtime: smoke facts must persist."""
    data = json.loads(SMOKE_JSON.read_text(encoding="utf-8"))
    assert data["maple"]["version"].startswith("Maple 2025.1")
    assert data["maple"]["smoke_result"] == "MAPLE_SMOKE=2"
    assert data["maple"]["verdict"] == "PASS"
    assert data["hyperint"]["revision"] == "ce15b287022e698d3e3b884d5c827620d20499bd"
    assert data["hyperint"]["smoke_expression"] == "hyperInt(1/(1+x)^2, x=0..infinity)"
    assert data["hyperint"]["smoke_result"] == "1"
    assert data["hyperint"]["verdict"] == "PASS"


def test_master_inputs_untouched_by_this_session() -> None:
    """The prepared master inputs still carry their review gate, uncommented."""
    for master in runner.MASTERS:
        text = runner.input_path(master).read_text(encoding="utf-8")
        assert "# result := hyperInt(" in text
        assert "intOrder := [x2, x5, x7]" in text


# ----------------------------------------------- stale-claim audit (Method.12R)


def test_notes_record_the_method12r_revocation() -> None:
    text = SETUP_NOTES.read_text(encoding="utf-8")
    assert "Method.12R" in text
    assert "REVOKED" in text
    assert "historical" in text.lower() and "fixture" in text.lower()
    # runtime setup is explicitly still usable
    assert "reusable" in text.lower()


def test_notes_contain_no_stale_readiness_claim() -> None:
    """The old 'L1 is ready, just needs approval' framing must be gone."""
    text = SETUP_NOTES.read_text(encoding="utf-8")
    assert "Blocked **only** on the Method.12N" not in text
    assert "Readiness for the real L1 run" not in text
    assert "Not ready" in text


def test_notes_do_not_call_the_current_basis_certified() -> None:
    """No surviving sentence may present the revoked basis as certified."""
    text = SETUP_NOTES.read_text(encoding="utf-8")
    for stale in (
        "the certified LF basis (4 labels, 4 coefficients), the surface\npolicy and the "
        "certificate are unchanged",
        "certified L1 input",
        "without editing the certified\ninput file",
        "Input (certified, unchanged)",
    ):
        assert stale not in text, f"stale claim still present: {stale!r}"
