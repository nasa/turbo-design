"""Tests for turbodesign.shaft_match: non-dimensional-MFP-based compressor/turbine matching."""

import numpy as np
import pytest

from turbodesign import Inlet, Outlet, Passage, PassageType, TurbineSpool
from turbodesign.compressor_spool import CompressorSpool
from turbodesign.loss.fixedpressureloss import FixedPressureLoss
from turbodesign.row_factory import make_rotor_row, make_stator_row
from turbodesign.shaft_match import ComponentSizingSpec, ShaftMatch, turbine_massflow


def test_turbine_massflow_continuity():
    assert turbine_massflow(100.0, fuel_air_ratio=0.02, bleed_fraction=0.05) == pytest.approx(97.0)
    assert turbine_massflow(100.0) == pytest.approx(100.0)


def test_shaft_match_rejects_invalid_fuel_bleed_and_efficiency():
    known = _solved_compressor()
    spec = ComponentSizingSpec(inlet_T0=1200.0, exit_static_pressure=1e5)
    with pytest.raises(ValueError):
        ShaftMatch(known, spec, reference=None, fuel_air_ratio=-0.1)
    with pytest.raises(ValueError):
        ShaftMatch(known, spec, reference=None, bleed_fraction=1.0)
    with pytest.raises(ValueError):
        ShaftMatch(known, spec, reference=None, mechanical_efficiency=0.0)


def test_shaft_match_requires_exactly_one_exit_condition():
    known = _solved_compressor()
    with pytest.raises(ValueError):
        ShaftMatch(known, ComponentSizingSpec(inlet_T0=1200.0), reference=None)
    with pytest.raises(ValueError):
        ShaftMatch(
            known,
            ComponentSizingSpec(inlet_T0=1200.0, exit_static_pressure=1e5, expansion_pressure_ratio=0.4),
            reference=None,
        )


def test_shaft_match_rejects_a_known_turbine_direction():
    """Only the compressor-known direction (compressor -> combustor -> turbine) is
    implemented; sizing a compressor from a known turbine is explicit future work."""
    reference = _reference_turbine(rpm=9549.0)
    spec = ComponentSizingSpec(inlet_T0=1200.0, exit_static_pressure=1e5)
    with pytest.raises(NotImplementedError):
        ShaftMatch(reference, spec, reference=None)


def _solved_compressor() -> CompressorSpool:
    T01 = 518.7 / 1.8
    P01 = 14.7 * 6894.76
    massflow_kg_s = 50 * 0.453592
    M1 = 0.7
    rmean = 12 * 0.0254

    area1 = 207.2 / 39.3701**2
    area2 = 179.1 / 39.3701**2
    area3 = 165.3 / 39.3701**2
    h1 = area1 / (np.pi * 4 * rmean)
    h2 = area2 / (np.pi * 4 * rmean)
    h3 = area3 / (np.pi * 4 * rmean)
    cax = 1 * 0.0254

    xhub_arr = [0, cax, 2 * cax]
    xshroud_arr = [0, cax, 2 * cax]
    rhub_arr = [rmean - h1, rmean - h2, rmean - h3]
    rshroud_arr = [rmean + h1, rmean + h2, rmean + h3]

    passage = Passage(xhub_arr, rhub_arr, xshroud_arr, rshroud_arr, passageType=PassageType.Axial, zero_phi=True)
    inlet = Inlet(hub_location=0)
    inlet.alpha2 = [40]
    inlet.init_total(P01, T01, M=M1)

    rotor = make_rotor_row(hub_location=cax / max(xhub_arr), metal_exit_angle_deg=[-23.87], loss_function=FixedPressureLoss(0))
    stator = make_stator_row(hub_location=2 * cax / max(xhub_arr), metal_exit_angle_deg=[40], loss_function=FixedPressureLoss(0), P0_ratio=1.3)

    outlet = Outlet()
    outlet.init_total(1.3 * P01, 0.5)

    rpm = 1000 * 30 / np.pi
    spool = CompressorSpool(passage, massflow_kg_s, inlet, outlet, [rotor, stator], num_streamlines=1, rpm=rpm)
    spool.solve()
    return spool


def _reference_turbine(rpm: float) -> TurbineSpool:
    T0, P0, P_exit, rmean = 1200.0, 126680.0, 90000.0, 0.15
    area1, area2, area3 = 0.05, 0.045, 0.05
    h1 = area1 / (np.pi * 4 * rmean)
    h2 = area2 / (np.pi * 4 * rmean)
    h3 = area3 / (np.pi * 4 * rmean)
    cax = 0.02

    xhub_arr = [0, cax, 2 * cax]
    xshroud_arr = [0, cax, 2 * cax]
    rhub_arr = [rmean - h1, rmean - h2, rmean - h3]
    rshroud_arr = [rmean + h1, rmean + h2, rmean + h3]
    passage = Passage(xhub_arr, rhub_arr, xshroud_arr, rshroud_arr, passageType=PassageType.Axial, zero_phi=True)

    inlet = Inlet(hub_location=0, alpha=[0])
    inlet.init_total(P0=P0, T0=T0, M=0.15)
    outlet = Outlet(num_streamlines=1)
    outlet.init_static(P=P_exit, percent_radii=[0.5])

    stator = make_stator_row(hub_location=cax / max(xhub_arr), metal_exit_angle_deg=[65.0], loss_function=FixedPressureLoss(0.03))
    rotor = make_rotor_row(hub_location=2 * cax / max(xhub_arr), metal_exit_angle_deg=[-55.0], loss_function=FixedPressureLoss(0.05))

    turbine = TurbineSpool(passage, 5.0, inlet, outlet, [stator, rotor], num_streamlines=1, rpm=rpm)
    turbine.adjust_streamlines = False
    return turbine


@pytest.fixture(scope="module")
def solved_compressor() -> CompressorSpool:
    return _solved_compressor()


class TestSizingIsPureAlgebra:
    """size() does no solving of its own - every check here should be exact algebra."""

    def test_power_balance_closes_exactly(self, solved_compressor):
        spec = ComponentSizingSpec(inlet_T0=1200.0, eta_total_guess=0.88, exit_static_pressure=101325.0)
        sizing = ShaftMatch(solved_compressor, spec, reference=None).size()
        assert sizing.required_power == pytest.approx(solved_compressor.total_power(), rel=1e-9)

    def test_continuity_with_fuel_and_bleed(self, solved_compressor):
        spec = ComponentSizingSpec(inlet_T0=1200.0, eta_total_guess=0.88, exit_static_pressure=101325.0)
        sm = ShaftMatch(solved_compressor, spec, reference=None, fuel_air_ratio=0.02, bleed_fraction=0.05)
        sizing = sm.size()
        mdot_c = float(solved_compressor._all_rows()[1].total_massflow_no_coolant)
        assert sizing.massflow == pytest.approx(turbine_massflow(mdot_c, 0.02, 0.05), rel=1e-9)

    def test_station_areas_round_trip_through_massflow(self, solved_compressor):
        """The real check: feeding each sized station's area back through the MFP
        relation at its target Mach reproduces the target massflow - this is what
        ties the sizing pass to the isentropic.py utilities under test."""
        from turbodesign.isentropic import Massflow

        spec = ComponentSizingSpec(inlet_T0=1200.0, eta_total_guess=0.88, exit_static_pressure=101325.0)
        sm = ShaftMatch(solved_compressor, spec, reference=None)
        sizing = sm.size()
        assert sizing.station_M[0] == pytest.approx(1.0)  # ngv_choked defaults to True
        station_P0 = np.array([sizing.inlet_P0, np.sqrt(sizing.inlet_P0 * sizing.exit_P0), sizing.exit_P0])
        station_T0 = np.array([sizing.inlet_T0, 0.5 * (sizing.inlet_T0 + sizing.exit_T0), sizing.exit_T0])
        recovered = Massflow(station_P0, station_T0, sizing.station_area, sizing.station_M, sizing.gamma, sizing.R)
        assert recovered == pytest.approx(sizing.massflow, rel=1e-6)

    def test_higher_inlet_temperature_raises_required_area_at_fixed_power(self, solved_compressor):
        spec_cool = ComponentSizingSpec(inlet_T0=1000.0, eta_total_guess=0.88, exit_static_pressure=101325.0)
        spec_hot = ComponentSizingSpec(inlet_T0=2000.0, eta_total_guess=0.88, exit_static_pressure=101325.0)
        sizing_cool = ShaftMatch(solved_compressor, spec_cool, reference=None).size()
        sizing_hot = ShaftMatch(solved_compressor, spec_hot, reference=None).size()
        # Hotter gas carries the same power with a smaller required temperature
        # drop fraction and higher volumetric flow -> larger sized areas.
        assert np.mean(sizing_hot.station_area) > np.mean(sizing_cool.station_area)

    def test_ngv_choked_false_falls_back_to_inlet_mach(self, solved_compressor):
        spec = ComponentSizingSpec(
            inlet_T0=1200.0, eta_total_guess=0.88, exit_static_pressure=101325.0,
            ngv_choked=False, inlet_mach=0.3,
        )
        sizing = ShaftMatch(solved_compressor, spec, reference=None).size()
        assert sizing.station_M[0] == pytest.approx(0.3)

    def test_ngv_choked_area_equals_min_area_for_massflow(self, solved_compressor):
        from turbodesign.isentropic import min_area_for_massflow

        spec = ComponentSizingSpec(inlet_T0=1200.0, eta_total_guess=0.88, exit_static_pressure=101325.0)
        sizing = ShaftMatch(solved_compressor, spec, reference=None).size()
        expected = min_area_for_massflow(sizing.massflow, sizing.inlet_P0, sizing.inlet_T0, sizing.gamma, sizing.R)
        assert sizing.station_area[0] == pytest.approx(expected, rel=1e-9)

    def test_negative_or_zero_compressor_power_raises(self):
        class _StubCompressor:
            rpm = 9549.0
            fluid = None

            def _all_rows(self):
                class Row:
                    total_massflow_no_coolant = 10.0
                    P0 = np.array([2e5])
                    T0 = np.array([400.0])

                return [None, Row(), Row()]

            def total_power(self):
                return -1.0

        spec = ComponentSizingSpec(inlet_T0=1200.0, exit_static_pressure=1e5)
        sm = ShaftMatch.__new__(ShaftMatch)
        sm.known = _StubCompressor()
        sm.spec = spec
        sm.reference = None
        sm.fuel_air_ratio = 0.0
        sm.bleed_fraction = 0.0
        sm.mechanical_efficiency = 1.0
        sm.rpm = 9549.0
        with pytest.raises(ValueError):
            sm.size()


@pytest.mark.slow
def test_match_converges_to_a_self_consistent_operating_point(solved_compressor):
    """Integration test: size + build + match a small reference turbine against
    example-9.1-style compressor. Uses a looser tolerance than the unit-level
    algebra tests, since this is a genuinely nonlinear two-knob solve against an
    arbitrarily-chosen reference turbine geometry (ScaledReferenceBuilder reuses
    its angles verbatim) - a few percent residual is the expected first-cut
    quality for this approach, not a bug. `MeanlineSynthesisBuilder`, which would
    derive angles instead of reusing a reference's, is explicit follow-up work.
    """
    rpm = float(solved_compressor.rpm)
    reference = _reference_turbine(rpm)

    spec = ComponentSizingSpec(inlet_T0=1200.0, eta_total_guess=0.88, exit_static_pressure=101325.0, num_stages=1)
    sm = ShaftMatch(solved_compressor, spec, reference=reference)

    result = sm.match(tol_rel=5e-3, max_iter=20)

    assert result.turbine is not None
    assert result.turbine is not result.compressor
    assert abs(result.massflow_residual) < 0.02
    assert abs(result.power_residual) < 0.05

    # Re-running match() must give the same answer - proves build() constructs
    # fresh state each call rather than mutating and re-solving stale rows.
    result2 = sm.match(tol_rel=5e-3, max_iter=20)
    assert result.area_scale == pytest.approx(result2.area_scale, rel=1e-3)
    assert result.power_turbine == pytest.approx(result2.power_turbine, rel=1e-3)

    for row in result.turbine._all_rows():
        if row.row_type.name in ("Rotor", "Stator", "IGV"):
            assert row.choke_margin_min > 0
