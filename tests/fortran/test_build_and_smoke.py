"""Low-level checks on the compiled Fortran binary itself, independent of
the Python API wrapper. These are smoke tests: they don't assert specific
physical results (that's the job of the integration tests), only that the
solver builds and runs to completion without crashing.

Note: the repo's checked-in datafiles/CalcInput.dat is a stale example that
predates the current ReadInput() format in CO2Transport_exe.f90 — in
production it is always overwritten by api/main.py (TEMPLATE.format(...))
before the binary runs. So these smoke tests generate a valid CalcInput.dat
via that same TEMPLATE, rather than asserting on the stale checked-in file.
"""

import subprocess
import sys


def test_binary_builds(run_env):
    assert run_env["exe"].is_file()


def _write_valid_calc_input(run_env, api_module, **overrides):
    request_kwargs = dict(
        temperature=25.0,
        pressure=50.0,
        concentrations=api_module.Concentrations(),
        reactive=True,
        solubility=True,
    )
    request_kwargs.update(overrides)
    request = api_module.SimulationRequest(**request_kwargs)
    data = request.model_dump()
    concentrations = data["concentrations"]
    concentrations["co2"] = 1_000_000 - sum(concentrations.values())

    calc_input = run_env["run_datafiles"] / "CalcInput.dat"
    calc_input.write_text(
        api_module.TEMPLATE.format(
            **concentrations,
            pressure=data["pressure"],
            temperature=data["temperature"],
            reactive=int(data["reactive"]),
            solubility=int(data["solubility"]),
            cond_factor=data["cond_factor"],
        )
    )


def test_binary_runs_to_completion_with_valid_input(run_env, monkeypatch):
    """Generates a CalcInput.dat the same way api/main.py does, then runs
    the solver directly (bypassing the API layer). Exit code 0 and the
    absence of a Fortran runtime error is the bar for this smoke test.
    """
    sys.path.insert(0, str(run_env["run_src"].parent.parent / "api"))
    try:
        import main as api_main  # type: ignore
    finally:
        sys.path.pop(0)

    _write_valid_calc_input(run_env, api_main)

    proc = subprocess.run(
        [str(run_env["exe"])],
        cwd=run_env["run_src"],
        capture_output=True,
        text=True,
        timeout=120,
    )

    assert proc.returncode == 0, (
        f"solver exited with {proc.returncode}\nstdout:\n{proc.stdout}\n"
        f"stderr:\n{proc.stderr}"
    )
    assert "runtime error" not in proc.stderr.lower()
    assert "Fortran runtime error" not in proc.stdout
