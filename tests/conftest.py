"""Shared pytest fixtures for building the Fortran solver and running it
exactly the way api/main.py invokes it (same CWD convention, same PATH/
RUN_FOLDER wiring), so integration tests exercise the real code path.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src"
DATAFILES_DIR = REPO_ROOT / "datafiles"
API_DIR = REPO_ROOT / "api"

# Allow overriding the Fortran compiler. Defaults to ifx, matching
# production (see Dockerfile) so tests exercise the real compiler behavior.
FC = os.environ.get("TEST_FORTRAN_COMPILER", "ifx")
FORTRAN_FLAGS = os.environ.get("TEST_FORTRAN_FLAGS", "")
LINKER_FLAGS = os.environ.get("TEST_FORTRAN_LINKER_FLAGS", "-qmkl")


@pytest.fixture(scope="session")
def fortran_compiler() -> str:
    if not shutil.which(FC):
        pytest.skip(f"No usable Fortran compiler found (tried '{FC}')")
    return FC


@pytest.fixture(scope="session")
def build_dir(tmp_path_factory, fortran_compiler) -> Path:
    """Configure and build the CO2Transport executable once per test session."""
    build = tmp_path_factory.mktemp("build")

    cmake_args = [
        "cmake",
        f"-DCMAKE_Fortran_COMPILER={fortran_compiler}",
    ]
    if FORTRAN_FLAGS:
        cmake_args.append(f"-DCMAKE_Fortran_FLAGS={FORTRAN_FLAGS}")
    if LINKER_FLAGS:
        cmake_args.append(f"-DCMAKE_EXE_LINKER_FLAGS={LINKER_FLAGS}")
    cmake_args.append(str(REPO_ROOT))

    configure = subprocess.run(
        cmake_args,
        cwd=build,
        capture_output=True,
        text=True,
    )
    assert configure.returncode == 0, (
        f"cmake configure failed:\n{configure.stdout}\n{configure.stderr}"
    )

    build_proc = subprocess.run(
        ["cmake", "--build", ".", "--", "-j", str(os.cpu_count() or 2)],
        cwd=build,
        capture_output=True,
        text=True,
    )
    assert build_proc.returncode == 0, (
        f"build failed:\n{build_proc.stdout}\n{build_proc.stderr}"
    )

    exe = build / "CO2Transport"
    assert exe.is_file(), f"Expected executable at {exe}, build output:\n{build_proc.stdout}"
    return build


@pytest.fixture(scope="session")
def run_env(tmp_path_factory, build_dir) -> dict:
    """Stage a run directory laid out like the Dockerfile does:
    <run>/src/CO2Transport next to a sibling <run>/datafiles/, matching
    MainPath="../datafiles/" in CO2Transport_exe.f90 and RUN_FOLDER in
    api/main.py.
    """
    run_root = tmp_path_factory.mktemp("run")
    run_src = run_root / "src"
    run_datafiles = run_root / "datafiles"
    run_src.mkdir()
    run_datafiles.mkdir()

    shutil.copy(build_dir / "CO2Transport", run_src / "CO2Transport")
    (run_src / "CO2Transport").chmod(0o755)

    for f in DATAFILES_DIR.iterdir():
        if f.is_file():
            shutil.copy(f, run_datafiles / f.name)

    return {
        "build_dir": build_dir,
        "run_root": run_root,
        "run_src": run_src,
        "run_datafiles": run_datafiles,
        "exe": run_src / "CO2Transport",
    }


@pytest.fixture()
def api_module(run_env, monkeypatch):
    """Import api.main with PATH/RUN_FOLDER set up so calculate() finds and
    runs the freshly-built CO2Transport binary, same as production."""
    monkeypatch.setenv("PATH", f"{run_env['run_src']}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("RUN_FOLDER", str(run_env["run_src"]))

    sys.path.insert(0, str(API_DIR))
    try:
        import main as api_main  # type: ignore

        import importlib

        importlib.reload(api_main)
        # RUNFOLDER is read at import time, so re-read it post-reload.
        api_main.RUNFOLDER = str(run_env["run_src"])
        yield api_main
    finally:
        sys.path.remove(str(API_DIR))
        sys.modules.pop("main", None)


@pytest.fixture()
def run_calc(api_module):
    """Helper: build a SimulationRequest from kwargs and run calculate()."""
    import asyncio

    def _run(**kwargs):
        # Trace-gas defaults chosen so an aqueous phase actually forms
        # (pure CO2/H2O alone does not condense an acid phase) — matches
        # the scenario captured in tests/fixtures/baseline_default_case.json.
        default_concentrations = dict(
            h2o=5, o2=5, h2s=5, no=5, no2=5, so2=5, h2so4=5, hno3=5, s8=1
        )
        defaults = dict(
            temperature=25.0,
            pressure=50.0,
            concentrations=api_module.Concentrations(**default_concentrations),
            reactive=True,
            solubility=True,
        )
        defaults.update(kwargs)
        request = api_module.SimulationRequest(**defaults)
        return asyncio.run(api_module.calculate(request))

    return _run
