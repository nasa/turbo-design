"""Tests for the upfront feasibility guard (flow_math.assert_flow_capacity)."""

import numpy as np
import pytest

from turbodesign.bladerow import BladeRow
from turbodesign.enums import RowType
from turbodesign.flow_math import assert_flow_capacity, row_flow_capacity
from turbodesign.isentropic import mass_flow_function_max, min_area_for_massflow


def _row(row_type, total_area, P0, T0, massflow, gamma=1.4, R=287.0, blockage=0.0, id_=0):
    row = BladeRow(row_type=row_type)
    row.id = id_
    row.total_area = total_area
    row.P0 = np.array([P0])
    row.T0 = np.array([T0])
    row.gamma = gamma
    row.R = R
    row.blockage = blockage
    row.total_massflow = massflow
    return row


def test_feasible_row_does_not_raise():
    row = _row(RowType.Stator, total_area=0.1, P0=2e5, T0=400.0, massflow=5.0)
    assert_flow_capacity([row], "test")  # must not raise


def test_infeasible_row_raises_value_error_with_actionable_message():
    # A tiny area cannot pass this massflow at any Mach - required m~ > m~_max.
    row = _row(RowType.Stator, total_area=1e-5, P0=1e5, T0=300.0, massflow=50.0, id_=3)
    with pytest.raises(ValueError) as excinfo:
        assert_flow_capacity([row], "test-context")
    msg = str(excinfo.value)
    assert "test-context" in msg
    assert "blade row 3" in msg
    assert "50.0000" in msg or "50.00" in msg
    assert "minimum area" in msg


def test_multiple_infeasible_rows_all_listed_in_one_exception():
    row1 = _row(RowType.Stator, total_area=1e-5, P0=1e5, T0=300.0, massflow=50.0, id_=1)
    row2 = _row(RowType.Rotor, total_area=1e-5, P0=1e5, T0=300.0, massflow=50.0, id_=2)
    row2.P0R = np.array([1e5])
    row2.T0R = np.array([300.0])
    with pytest.raises(ValueError) as excinfo:
        assert_flow_capacity([row1, row2], "test")
    msg = str(excinfo.value)
    assert "blade row 1" in msg
    assert "blade row 2" in msg


def test_non_blade_rows_are_skipped():
    inlet = _row(RowType.Inlet, total_area=1e-8, P0=1e5, T0=300.0, massflow=1e6)
    assert_flow_capacity([inlet], "test")  # skipped entirely, must not raise


def test_rotor_uses_relative_frame_when_available():
    row = _row(RowType.Rotor, total_area=0.1, P0=1e5, T0=300.0, massflow=5.0)
    row.P0R = np.array([2e5])  # much more capacity in the relative frame
    row.T0R = np.array([250.0])
    cap = row_flow_capacity(row)
    assert cap["frame"] == "relative"
    assert cap["P0"] == pytest.approx(2e5)


def test_rotor_falls_back_to_absolute_frame_when_relative_unset():
    row = _row(RowType.Rotor, total_area=0.1, P0=1e5, T0=300.0, massflow=5.0)
    # P0R/T0R left at BladeRow's default of [0] - relative frame "unavailable"
    cap = row_flow_capacity(row)
    assert "absolute" in cap["frame"]
    assert cap["P0"] == pytest.approx(1e5)


def test_row_flow_capacity_feasible_flag_matches_min_area_predicate():
    P0, T0, gamma, R, massflow = 1.5e5, 350.0, 1.4, 287.0, 8.0
    min_area = min_area_for_massflow(massflow, P0, T0, gamma, R)

    just_ok = _row(RowType.Stator, total_area=min_area * 1.01, P0=P0, T0=T0, massflow=massflow, gamma=gamma, R=R)
    just_bad = _row(RowType.Stator, total_area=min_area * 0.99, P0=P0, T0=T0, massflow=massflow, gamma=gamma, R=R)

    assert row_flow_capacity(just_ok)["feasible"] is True
    assert row_flow_capacity(just_bad)["feasible"] is False


def test_explicit_targets_override_row_total_massflow():
    row = _row(RowType.Stator, total_area=0.1, P0=1e5, T0=300.0, massflow=1.0)
    # row.total_massflow says 1.0 kg/s (feasible); explicit target says 1e6 (infeasible).
    with pytest.raises(ValueError):
        assert_flow_capacity([row], "test", targets=[1e6])
    assert_flow_capacity([row], "test", targets=[None])  # falls back to row.total_massflow


def test_shrinking_example9_1_annulus_raises_value_error_before_any_solve_output(capsys):
    """Regression: shrink example9.1.py's annulus until infeasible and confirm the
    NEW ValueError guard fires - not the old, uninformative RuntimeError from
    _solve_bounded pinning at the Mach bracket edge - and fires before any
    'Loop N massflow convergence error' output (i.e. before the outer loop even starts)."""
    import numpy as np

    from turbodesign import Inlet, Outlet, Passage, PassageType
    from turbodesign.compressor_spool import CompressorSpool
    from turbodesign.loss.fixedpressureloss import FixedPressureLoss
    from turbodesign.row_factory import make_rotor_row, make_stator_row

    T01_K, P01_Pa = 518.7 / 1.8, 14.7 * 6894.76
    massflow_kg_s = 50 * 0.453592
    M1 = 0.7
    rmean = 12 * 0.0254

    # Deliberately tiny areas: far too small to pass 50 lbm/s at this rmean.
    h = 1e-4
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

    spool = CompressorSpool(passage, massflow_kg_s, inlet, outlet, [rotor, stator], num_streamlines=1)

    with pytest.raises(ValueError) as excinfo:
        spool.solve()
    assert "RuntimeError" not in type(excinfo.value).__name__
    assert "choked" in str(excinfo.value) or "capacity" in str(excinfo.value)

    captured = capsys.readouterr()
    assert "Loop 1 massflow convergence error" not in captured.out
