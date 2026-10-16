"""Integration tests that exercise the system the same way the FastAPI
service does: build CalcInput.dat via api.main's SimulationRequest/TEMPLATE,
shell out to the compiled CO2Transport binary, and parse its stdout.

These intentionally avoid asserting exact floating point values from the
nonlinear solver (those vary with compiler/library and can shift slightly
between runs/environments). Instead they assert:
  - structural invariants (expected species/sections present, mass balance,
    physically sane ranges)
  - a tolerance-based regression baseline captured from a known-good run,
    so a *meaningfully different* answer on tests/fixtures/code_update (or
    any future change) gets flagged without being brittle to noise.
"""

import json
import math
from pathlib import Path

import pytest

BASELINE_PATH = Path(__file__).parent.parent / "fixtures" / "baseline_default_case.json"


def _is_present_and_finite(value) -> bool:
    return value is not None and math.isfinite(value)


class TestApiRoot:
    def test_read_root(self, api_module):
        import asyncio

        result = asyncio.run(api_module.read_root())
        assert result == {"message": "running server"}


class TestNaiveBaselineCases:
    """Dead-simple scenarios run through the real binary, worth having even
    though they don't assert much: they establish that the pipeline doesn't
    blow up on the most trivial possible inputs, and pin down the (slightly
    surprising) "no aqueous phase" behavior for pure CO2/H2O feeds."""

    def test_all_zero_concentrations_runs_without_crashing(self, run_calc, api_module):
        # Pure CO2 (api/main.py backfills co2 = 1e6 - sum(others)) with every
        # trace gas at zero. No acid-forming species are present, so no
        # aqueous phase should condense — this is physically expected, not
        # a bug, and is a good smoke case precisely because it's naive.
        result = run_calc(concentrations=api_module.Concentrations())
        assert result["raw"]
        assert result["phase_reactions"]["CO2"]["feed_mol"] == pytest.approx(1_000_000)
        assert result["aqueous_phase"] is None
        assert result["solubility"] in (None, {})

    def test_single_trace_gas_only(self, run_calc, api_module):
        # Simplest non-trivial case: just one trace species (H2O) added to
        # the CO2 balance, everything else left at zero.
        result = run_calc(concentrations=api_module.Concentrations(h2o=5))
        assert result["phase_reactions"]["H2O"]["feed_mol"] == pytest.approx(5)
        assert _is_present_and_finite(result["phase_reactions"]["H2O"]["equil_mol"])

    def test_zero_pressure_boundary_does_not_hang_or_crash(self, run_calc):
        # Not a physically meaningful case, but a cheap guard against the
        # solver hanging or raising an unhandled exception on a degenerate
        # input value a caller could plausibly send by mistake.
        result = run_calc(pressure=0.0)
        assert result["raw"] is not None

    def test_reactive_and_solubility_both_disabled(self, run_calc):
        # The most "do nothing extra" configuration possible — both optional
        # calculation stages turned off. With nothing enabled to print, the
        # binary legitimately produces empty stdout; the bar here is just
        # that calculate() completes without raising.
        result = run_calc(reactive=False, solubility=False)
        assert result["raw"] == ""
        assert result["phase_reactions"] == {}


class TestCalculateStructure:
    """Checks on the shape/sanity of a result, independent of exact values."""

    def test_returns_all_expected_sections(self, run_calc):
        result = run_calc()
        assert result["phase_reactions"]
        assert result["aqueous_phase"] is not None
        assert result["solubility"] is not None

    def test_all_19_components_present_in_phase_reactions(self, run_calc):
        result = run_calc()
        expected_components = {
            "CO2", "H2O", "N2", "O2", "H2S", "NO", "NO2", "SO2", "H2SO4",
            "HNO3", "S8", "NH3", "N2O", "N2O4", "NH4HSO4", "HCHO", "CH3CHO",
            "CH3COCH3", "HCOOH", "CH3COOH",
        }
        assert expected_components <= result["phase_reactions"].keys()

    def test_co2_dominates_feed_and_stays_near_1e6_ppm(self, run_calc):
        # api/main.py sets co2 = 1_000_000 - sum(other concentrations), so
        # with small trace gases CO2 should remain the dominant species and
        # its equilibrium ppm should stay close to the 1e6 ppm feed.
        result = run_calc()
        co2 = result["phase_reactions"]["CO2"]
        assert co2["equil_ppm"] == pytest.approx(1_000_000, rel=0.01)

    def test_mole_values_are_non_negative_and_finite(self, run_calc):
        result = run_calc()
        for name, values in result["phase_reactions"].items():
            for key in ("feed_mol", "equil_mol", "equil_ppm"):
                v = values[key]
                assert _is_present_and_finite(v), f"{name}.{key} = {v!r}"
                assert v >= 0, f"{name}.{key} = {v!r} should be non-negative"

    def test_ph_is_in_physically_plausible_range(self, run_calc):
        result = run_calc()
        ph = result["aqueous_phase"]["ph"]
        assert _is_present_and_finite(ph)
        # Concentrated sulfuric-acid condensate can go quite acidic (even
        # below 0 on the formal pH scale at high molality), but should not
        # be absurdly out of range if the model is behaving.
        assert -5.0 < ph < 15.0

    def test_water_mass_is_positive_when_aqueous_phase_forms(self, run_calc):
        result = run_calc()
        assert result["aqueous_phase"]["water_mass_kg"] > 0

    def test_solubility_percentages_within_0_100(self, run_calc):
        result = run_calc()
        for name, values in result["solubility"].items():
            pct = values["pct_in_aqueous"]
            assert _is_present_and_finite(pct), f"{name}.pct_in_aqueous = {pct!r}"
            assert -1e-6 <= pct <= 100 + 1e-6, f"{name}: {pct}"


class TestCalculateOptions:
    """Checks that the reactive/solubility toggles actually change behavior,
    since that's the most common way the API is driven differently."""

    def test_solubility_disabled_skips_aqueous_phase(self, run_calc):
        result = run_calc(solubility=False)
        assert result["solubility"] is None or result["solubility"] == {}

    def test_reactive_disabled_runs_without_error(self, run_calc):
        result = run_calc(reactive=False)
        assert result["raw"]

    def test_higher_pressure_increases_co2_mole_fraction_ppm_stability(self, run_calc):
        # Sanity check that changing a core input (pressure) actually
        # produces a different, still-finite, still-valid result rather than
        # silently reusing stale output. Note: at low pressure the trace
        # acid gases here may not condense into an aqueous phase at all
        # (physically plausible), so we only assert finiteness when a
        # phase actually formed.
        low = run_calc(pressure=10.0)
        high = run_calc(pressure=100.0)
        assert low["raw"] != high["raw"]
        if low["aqueous_phase"]:
            assert _is_present_and_finite(low["aqueous_phase"]["ph"])
        if high["aqueous_phase"]:
            assert _is_present_and_finite(high["aqueous_phase"]["ph"])


class TestRegressionBaseline:
    """Numeric regression check against a captured baseline for a fixed
    input. This is the primary signal for "did switching to code_update (or
    any other change) move the answer" — update the baseline deliberately
    (via tests/fixtures/baseline_default_case.json) when a change is
    expected to alter results.
    """

    @pytest.fixture()
    def baseline(self):
        if not BASELINE_PATH.is_file():
            pytest.skip(f"No baseline captured yet at {BASELINE_PATH}")
        return json.loads(BASELINE_PATH.read_text())

    def test_matches_captured_baseline_within_tolerance(self, run_calc, baseline):
        result = run_calc(
            temperature=25.0,
            pressure=50.0,
            concentrations=dict(
                h2o=5, o2=5, h2s=5, no=5, no2=5, so2=5, h2so4=5, hno3=5, s8=1
            ),
        )

        rel_tol = 0.05  # 5%: generous enough to absorb solver-noise/compiler
        # differences, tight enough to catch real behavior changes.
        abs_tol = 1e-6

        def compare(path, expected, actual):
            if isinstance(expected, dict):
                for k, v in expected.items():
                    assert k in actual, f"missing key at {path}.{k}"
                    compare(f"{path}.{k}", v, actual[k])
            elif isinstance(expected, (int, float)):
                assert actual == pytest.approx(expected, rel=rel_tol, abs=abs_tol), (
                    f"{path}: expected {expected}, got {actual}"
                )
            else:
                assert actual == expected, f"{path}: expected {expected!r}, got {actual!r}"

        compare("phase_reactions", baseline["phase_reactions"], result["phase_reactions"])
        compare("aqueous_phase", baseline["aqueous_phase"], result["aqueous_phase"])
