"""Fast, naive unit tests for the pure-Python helpers in api/main.py.

These do not touch the Fortran binary at all — they run in milliseconds
and exercise the text-parsing/model logic directly with hand-crafted
strings, which is the cheapest way to pin down behavior (including some
rough edges) of parse_float/parse_output.
"""

import pytest


@pytest.fixture()
def api_main():
    import sys
    from pathlib import Path

    api_dir = Path(__file__).resolve().parents[1] / "api"
    sys.path.insert(0, str(api_dir))
    try:
        import main as api_main  # type: ignore

        yield api_main
    finally:
        sys.path.remove(str(api_dir))
        sys.modules.pop("main", None)


class TestParseFloat:
    def test_plain_float(self, api_main):
        assert api_main.parse_float("1.2345") == 1.2345

    def test_integer_string(self, api_main):
        assert api_main.parse_float("5") == 5.0

    def test_negative_number(self, api_main):
        assert api_main.parse_float("-3.14") == -3.14

    def test_zero(self, api_main):
        assert api_main.parse_float("0.0") == 0.0

    def test_standard_scientific_notation(self, api_main):
        assert api_main.parse_float("1.5E-05") == pytest.approx(1.5e-05)

    def test_fortran_style_missing_exponent_marker(self, api_main):
        # Some Fortran compilers emit exponents without the 'E'/'D' marker
        # when the exponent is 3 digits, e.g. "1.234-100" meaning 1.234E-100.
        # parse_float's fallback branch is expected to handle this.
        assert api_main.parse_float("1.234-100") == pytest.approx(1.234e-100)

    def test_positive_exponent_without_marker_is_not_recovered(self, api_main):
        # The fallback only looks for a '-' to re-insert an 'e', so a
        # missing-marker *positive* exponent (no '-' present) still raises.
        # This documents existing (somewhat naive) behavior rather than
        # asserting it's ideal.
        with pytest.raises(ValueError):
            api_main.parse_float("1.234+100")

    def test_garbage_input_raises(self, api_main):
        with pytest.raises(ValueError):
            api_main.parse_float("not-a-number")

    def test_leading_minus_is_not_mistaken_for_missing_exponent(self, api_main):
        # value[1:] excludes the leading sign so a simple negative float
        # must still parse via the fast path, not the fallback.
        assert api_main.parse_float("-0.0042") == pytest.approx(-0.0042)


class TestParseOutputEmptyAndNoMatch:
    def test_empty_string(self, api_main):
        result = api_main.parse_output("")
        assert result["raw"] == ""
        assert result["phase_reactions"] == {}
        assert result["aqueous_phase"] is None
        assert result["solubility"] is None

    def test_irrelevant_text_leaves_defaults(self, api_main):
        result = api_main.parse_output("hello\nworld\nno markers here\n")
        assert result["phase_reactions"] == {}
        assert result["aqueous_phase"] is None
        assert result["solubility"] is None

    def test_raw_is_preserved_verbatim(self, api_main):
        text = "some\nmultiline\ntext\n"
        result = api_main.parse_output(text)
        assert result["raw"] == text


class TestParseOutputNoAqueousPhase:
    def test_no_aqueous_phase_present_sets_empty_solubility_dict(self, api_main):
        text = "blah\nNo aquous phase present\nblah\n"
        result = api_main.parse_output(text)
        assert result["solubility"] == {}
        assert result["aqueous_phase"] is None


class TestParseOutputPhaseReactions:
    SAMPLE = (
        "=================================================================================\n"
        " CO2 PHASE REACTIONS \n"
        "=================================================================================\n"
        "  Component    Feed (mol)     Equil (mol)       Equil (ppm)\n"
        "---------------------------------------------------------------------------------\n"
        "  CO2          9.999590E+05   9.999590E+05    999956.3987\n"
        "  H2O          5.000000E+00   1.261578E+01        12.6157\n"
        "===\n"
    )

    def test_parses_each_component_row(self, api_main):
        result = api_main.parse_output(self.SAMPLE)
        assert set(result["phase_reactions"]) == {"CO2", "H2O"}

    def test_feed_equil_and_ppm_values(self, api_main):
        result = api_main.parse_output(self.SAMPLE)
        co2 = result["phase_reactions"]["CO2"]
        assert co2["feed_mol"] == pytest.approx(999959.0)
        assert co2["equil_mol"] == pytest.approx(999959.0)
        assert co2["equil_ppm"] == pytest.approx(999956.3987)

    def test_stops_at_triple_equals_terminator(self, api_main):
        text = self.SAMPLE + "  EXTRA        1.0            1.0              1.0\n"
        result = api_main.parse_output(text)
        # The extra row after "===" must not be picked up.
        assert "EXTRA" not in result["phase_reactions"]

    def test_rows_with_too_few_columns_are_skipped(self, api_main):
        text = (
            " CO2 PHASE REACTIONS \n"
            "---\n"
            "  BAD   1.0\n"
            "  CO2   1.0   2.0   3.0\n"
            "===\n"
        )
        result = api_main.parse_output(text)
        assert "BAD" not in result["phase_reactions"]
        assert "CO2" in result["phase_reactions"]


class TestParseOutputSolubility:
    SAMPLE = (
        " MOLE RATES AFTER ...\n"
        "---\n"
        "---\n"
        "  CO2   1000000.0   1000000.0   0.0311   0.0\n"
        "  H2O   12.62        11.2         1.416    11.2254\n"
        "---\n"
    )

    def test_parses_solubility_rows(self, api_main):
        result = api_main.parse_output(self.SAMPLE)
        assert set(result["solubility"]) == {"CO2", "H2O"}

    def test_solubility_fields(self, api_main):
        result = api_main.parse_output(self.SAMPLE)
        h2o = result["solubility"]["H2O"]
        assert h2o["total_mol"] == pytest.approx(12.62)
        assert h2o["gas_mol"] == pytest.approx(11.2)
        assert h2o["water_mol"] == pytest.approx(1.416)
        assert h2o["pct_in_aqueous"] == pytest.approx(11.2254)

    def test_missing_pct_column_defaults_to_zero(self, api_main):
        text = (
            " MOLE RATES AFTER ...\n"
            "---\n"
            "---\n"
            "  CO2   1000000.0   1000000.0   0.0311\n"
            "---\n"
        )
        result = api_main.parse_output(text)
        assert result["solubility"]["CO2"]["pct_in_aqueous"] == 0.0


class TestParseOutputAqueousPhase:
    SAMPLE = (
        " AQUEOUS PHASE RESULTS\n"
        " Water pH               : -1.09\n"
        " Water phase            : 0.4946 kg\n"
        " Water activity         : 0.008107\n"
        " H2SO4 in water         : 94.5592 wt%\n"
        " HNO3  in water         : 0.0044 wt%\n"
        " NH3   in water         : 0.0017 wt%\n"
        " CO2   in water         : 0.2767 wt%\n"
        "------\n"
    )

    def test_all_fields_parsed(self, api_main):
        result = api_main.parse_output(self.SAMPLE)
        aq = result["aqueous_phase"]
        assert aq["ph"] == pytest.approx(-1.09)
        assert aq["water_mass_kg"] == pytest.approx(0.4946)
        assert aq["water_activity"] == pytest.approx(0.008107)
        assert aq["h2so4_wt_pct"] == pytest.approx(94.5592)
        assert aq["hno3_wt_pct"] == pytest.approx(0.0044)
        assert aq["nh3_wt_pct"] == pytest.approx(0.0017)
        assert aq["co2_wt_pct"] == pytest.approx(0.2767)

    def test_acid_phase_results_header_variant_also_matches(self, api_main):
        text = self.SAMPLE.replace("AQUEOUS PHASE RESULTS", "ACID PHASE RESULTS")
        result = api_main.parse_output(text)
        assert result["aqueous_phase"]["ph"] == pytest.approx(-1.09)

    def test_stops_at_dash_separator_once_populated(self, api_main):
        # A dash line before anything has been recorded must not
        # prematurely terminate the block (guarded by `result["aqueous_phase"]`
        # truthiness check), but one after real fields are parsed should end it.
        text = (
            " AQUEOUS PHASE RESULTS\n"
            "------\n"  # no fields yet -> aqueous_phase is still {} (falsy)
            " Water pH               : 2.5\n"
            "------\n"  # now it should stop
            " Water phase            : 999.0 kg\n"  # must NOT be picked up
        )
        result = api_main.parse_output(text)
        assert result["aqueous_phase"]["ph"] == pytest.approx(2.5)
        assert "water_mass_kg" not in result["aqueous_phase"]


class TestConcentrationsModel:
    def test_all_fields_default_to_zero(self, api_main):
        c = api_main.Concentrations()
        assert c.h2o == 0
        assert c.hno3 == 0
        assert c.ch3cooh == 0

    def test_accepts_explicit_values(self, api_main):
        c = api_main.Concentrations(h2o=5, hno3=2.5)
        assert c.h2o == 5
        assert c.hno3 == 2.5

    def test_model_dump_round_trips_all_keys(self, api_main):
        c = api_main.Concentrations()
        dumped = c.model_dump()
        expected_keys = {
            "h2o", "o2", "so2", "no2", "h2s", "no", "h2so4", "hno3", "s8",
            "n2", "nh3", "n2o", "n2o4", "nh4hso4", "hcho", "ch3cho",
            "ch3coch3", "hcooh", "ch3cooh",
        }
        assert set(dumped) == expected_keys


class TestSimulationRequestModel:
    def test_requires_temperature_and_pressure(self, api_main):
        with pytest.raises(Exception):
            api_main.SimulationRequest(concentrations=api_main.Concentrations())

    def test_reactive_and_solubility_default_true(self, api_main):
        req = api_main.SimulationRequest(
            temperature=25.0, pressure=50.0, concentrations=api_main.Concentrations()
        )
        assert req.reactive is True
        assert req.solubility is True

    def test_booleans_can_be_overridden(self, api_main):
        req = api_main.SimulationRequest(
            temperature=25.0,
            pressure=50.0,
            concentrations=api_main.Concentrations(),
            reactive=False,
            solubility=False,
        )
        assert req.reactive is False
        assert req.solubility is False

    def test_temperature_and_pressure_accept_floats(self, api_main):
        req = api_main.SimulationRequest(
            temperature=12.5, pressure=99.9, concentrations=api_main.Concentrations()
        )
        assert req.temperature == 12.5
        assert req.pressure == 99.9


class TestTemplateRendering:
    """Naive checks that TEMPLATE.format(...) produces the expected number
    of lines/fields for the Fortran ReadInput routine, without running the
    binary. This would catch e.g. an accidentally added/removed template
    line throwing off the fixed-format parser on the Fortran side.
    """

    def _render(self, api_main, **overrides):
        concentrations = api_main.Concentrations().model_dump()
        concentrations["co2"] = 1_000_000 - sum(concentrations.values())
        concentrations.update(overrides.pop("concentrations", {}))
        return api_main.TEMPLATE.format(
            **concentrations,
            pressure=overrides.pop("pressure", 50.0),
            temperature=overrides.pop("temperature", 25.0),
            reactive=overrides.pop("reactive", 1),
            solubility=overrides.pop("solubility", 1),
            cond_factor=overrides.pop("cond_factor", 1),
        )

    def test_renders_without_error(self, api_main):
        rendered = self._render(api_main)
        assert rendered

    def test_contains_pressure_and_temperature_line(self, api_main):
        rendered = self._render(api_main, pressure=50.0, temperature=25.0)
        assert " 50.0           25.0\n" in rendered

    def test_co2_is_balance_of_1e6(self, api_main):
        concentrations = api_main.Concentrations(h2o=10, o2=20).model_dump()
        concentrations["co2"] = 1_000_000 - sum(concentrations.values())
        assert concentrations["co2"] == 1_000_000 - 30

    def test_reactive_and_solubility_render_as_bare_ints(self, api_main):
        rendered = self._render(api_main, reactive=0, solubility=1)
        lines = [line.strip() for line in rendered.splitlines()]
        assert "0" in lines
        assert "1" in lines
