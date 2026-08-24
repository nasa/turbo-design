"""Tests for turbodesign.operating_map.sweep_operating_points."""

import numpy as np
import pytest

from turbodesign import Inlet, Outlet, Passage, PassageType
from turbodesign.compressor_spool import CompressorSpool
from turbodesign.loss.fixedpressureloss import FixedPressureLoss
from turbodesign.operating_map import sweep_operating_points
from turbodesign.row_factory import make_rotor_row, make_stator_row

from tests.test_shaft_match import _solved_compressor


@pytest.mark.slow
def test_sweep_is_ordered_rpm_major_massflow_minor():
    spool = _solved_compressor()
    base_mdot = spool._all_rows()[1].total_massflow_no_coolant
    massflow_points = [base_mdot * 0.9, base_mdot, base_mdot * 1.1]
    rpm_points = [spool.rpm, spool.rpm * 1.05]

    results = sweep_operating_points(spool, massflow_points, rpm_points)

    assert len(results) == len(massflow_points) * len(rpm_points)
    # First len(massflow_points) results are all at rpm_points[0], in massflow order.
    first_block = results[: len(massflow_points)]
    assert [r.rpm for r in first_block] == [rpm_points[0]] * len(massflow_points)
    assert [r.massflow for r in first_block] == pytest.approx(massflow_points)
    second_block = results[len(massflow_points):]
    assert [r.rpm for r in second_block] == [rpm_points[1]] * len(massflow_points)


@pytest.mark.slow
def test_choke_margin_falls_monotonically_as_massflow_rises_at_fixed_speed():
    """A pure continuity/area relationship - true regardless of how this
    particular fixed-angle geometry's pressure ratio happens to respond
    (verified separately: PR is monotonic here too, but not assumed to be
    falling - that's not universally true for a fixed-angle solve, only the
    choke-margin trend is safe to assert in general)."""
    spool = _solved_compressor()
    base_mdot = spool._all_rows()[1].total_massflow_no_coolant
    massflow_points = np.linspace(base_mdot * 0.85, base_mdot * 1.5, 6)

    results = sweep_operating_points(spool, massflow_points, [spool.rpm])

    assert all(r.converged for r in results)
    margins = [r.choke_margin_min for r in results]
    assert all(a > b for a, b in zip(margins, margins[1:])), margins
    # And the pressure ratio / efficiency / power came back as finite, sane numbers.
    for r in results:
        assert np.isfinite(r.pressure_ratio) and r.pressure_ratio > 0
        assert np.isfinite(r.power)
        assert np.isfinite(r.corrected_massflow) and r.corrected_massflow > 0
        assert np.isfinite(r.corrected_speed) and r.corrected_speed > 0


def _tiny_annulus_compressor() -> CompressorSpool:
    """A deliberately undersized single-stage compressor - infeasible at its own
    design massflow, so a sweep through it is guaranteed to hit the guard."""
    T01_K = 518.7 / 1.8
    P01_Pa = 14.7 * 6894.76
    massflow_kg_s = 50 * 0.453592
    M1 = 0.7
    rmean = 12 * 0.0254

    h = 1e-4  # deliberately tiny annulus height - see tests/unit/test_flow_capacity_guard.py
    cax = 1 * 0.0254
    xhub_arr = [0, cax, 2 * cax]
    xshroud_arr = [0, cax, 2 * cax]
    rhub_arr = [rmean - h, rmean - h, rmean - h]
    rshroud_arr = [rmean + h, rmean + h, rmean + h]

    passage = Passage(xhub_arr, rhub_arr, xshroud_arr, rshroud_arr, passageType=PassageType.Axial, zero_phi=True)
    inlet = Inlet(hub_location=0)
    inlet.alpha2 = [40]
    inlet.init_total(P01_Pa, T01_K, M=M1)

    rotor = make_rotor_row(hub_location=cax / max(xhub_arr), metal_exit_angle_deg=[-23.87], loss_function=FixedPressureLoss(0))
    stator = make_stator_row(hub_location=2 * cax / max(xhub_arr), metal_exit_angle_deg=[40], loss_function=FixedPressureLoss(0), P0_ratio=1.3)

    outlet = Outlet()
    outlet.init_total(1.3 * P01_Pa, 0.5)

    return CompressorSpool(passage, massflow_kg_s, inlet, outlet, [rotor, stator], num_streamlines=1)


def test_an_infeasible_point_is_recorded_not_raised_and_the_sweep_continues():
    spool = _tiny_annulus_compressor()

    results = sweep_operating_points(spool, [50 * 0.453592], [spool.rpm])

    assert len(results) == 1
    assert results[0].converged is False
    assert results[0].error is not None
    assert "choked" in results[0].error or "capacity" in results[0].error
    assert np.isnan(results[0].pressure_ratio)


def test_one_failed_point_does_not_abort_the_rest_of_the_sweep():
    spool = _tiny_annulus_compressor()
    # Same tiny geometry at three different massflows - all infeasible, but the
    # sweep must produce all three OperatingPoints rather than raising on the first.
    results = sweep_operating_points(spool, [5.0, 10.0, 15.0], [spool.rpm])

    assert len(results) == 3
    assert all(r.converged is False for r in results)
    assert [r.massflow for r in results] == [5.0, 10.0, 15.0]
