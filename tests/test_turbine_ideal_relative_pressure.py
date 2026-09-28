# tests/test_turbine_ideal_relative_pressure.py
"""Turbine counterpart of tests/test_ideal_relative_pressure.py (issue #40).

turbine_math.rotor_calc computes T0R2 from rothalpy conservation, so T0R moves with U^2
whenever the rotor exit radius differs from the inlet radius. It used to take P0R2 straight
from the inlet (P0R2 = P0R1 - Yp*(P0R1 - P2)), which on a radial-inflow turbine leaves P0R2
too high for T0R2: a loss-free rotor destroyed entropy.

With Yp = 0 the rotor is adiabatic and loss-free, so:

    P0R2 / P0R1 == (T0R2 / T0R1) ** (gamma / (gamma - 1))   and   ds == 0
"""

import numpy as np
import pytest

from turbodesign.bladerow import BladeRow
from turbodesign.enums import RowType
from turbodesign.turbine_math import rotor_calc

GAMMA = 1.33
R_GAS = 287.15
CP = GAMMA * R_GAS / (GAMMA - 1)
OMEGA = 5000.0  # rad/s


def _make_upstream() -> BladeRow:
    """Stator exit / rotor inlet of a radial-inflow turbine, U1 = 500 m/s."""
    row = BladeRow(row_type=RowType.Stator, r=np.array([0.10]), R=R_GAS, gamma=GAMMA, Cp=CP)
    row.rpm = OMEGA * 30 / np.pi
    row.Vm = np.array([100.0])
    row.Vx = np.array([0.0])
    row.Vr = np.array([100.0])
    row.Vt = np.array([450.0])
    row.T = np.array([1000.0])
    row.P = np.array([300000.0])
    row.P0 = row.P * (1 + (row.Vm**2 + row.Vt**2) / (2 * CP * row.T)) ** (GAMMA / (GAMMA - 1))
    row.P0_stator_inlet = row.P0
    return row


def _make_rotor(yp: float) -> BladeRow:
    """Rotor exit at half the inlet radius, so U2 = U1/2 and T0R drops ~80 K."""
    row = BladeRow(row_type=RowType.Rotor, r=np.array([0.05]), R=R_GAS, gamma=GAMMA, Cp=CP, omega=OMEGA)
    row.P = np.array([150000.0])
    row.beta2 = np.radians(np.array([-50.0]))
    row.phi = np.array([0.0])
    row.Yp = np.array([yp])
    return row


def test_loss_free_radial_rotor_ideal_exit_relative_pressure_is_isentropic():
    upstream = _make_upstream()
    row = _make_rotor(yp=0.0)

    rotor_calc(row, upstream, calculate_vm=True)

    assert row.T0R == pytest.approx(upstream.T0R - 0.75 * (OMEGA * 0.10) ** 2 / (2 * CP), rel=1e-12)
    expected_ratio = (row.T0R / upstream.T0R) ** (GAMMA / (GAMMA - 1))
    assert row.P0R / upstream.P0R == pytest.approx(expected_ratio, rel=1e-12)


def test_loss_free_radial_rotor_is_isentropic():
    upstream = _make_upstream()
    row = _make_rotor(yp=0.0)

    rotor_calc(row, upstream, calculate_vm=True)

    assert row.entropy_rise == pytest.approx([0.0], abs=1e-9)


def test_lossy_radial_rotor_generates_entropy():
    upstream = _make_upstream()
    row = _make_rotor(yp=0.1)

    rotor_calc(row, upstream, calculate_vm=True)

    assert np.all(row.entropy_rise > 0)
