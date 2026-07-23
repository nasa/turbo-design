# tests/test_ideal_relative_pressure.py
"""The lossless exit relative total pressure must be isentropic, not a copy of the inlet.

compressor_math.py sets P0R_local = upstream.P0R - Yp * (upstream.P0R - upstream.P) in
rotor_calc. At Yp = 0 that reads P0R2_ideal = P0R1, which holds only when U2 == U1. The
same function computes T0R2 from rothalpy conservation using the local U = omega*r --
and its own comment says so. Whenever the rotor's radius differs from its upstream radius,
U2 != U1, so T0R2 != T0R1, and a P0R2_ideal that does not move with T0R2 cannot be
isentropic.

The identity checked here is the same one derived independently, without importing this
library, in tests/test_relative_frame_oracle.py:

    P0R2 / P0R1 == (T0R2 / T0R1) ** (gamma / (gamma - 1))

Row setup follows tests/test_radial_velocity.py; the only difference is that the rotor's
radius is built away from the upstream radius so U2 != U1 actually occurs.
"""

import numpy as np
import pytest

from turbodesign.bladerow import BladeRow
from turbodesign.enums import RowType
from turbodesign.loss.fixedpressureloss import FixedPressureLoss
from turbodesign.compressor_math import rotor_calc


def _make_upstream() -> BladeRow:
    row = BladeRow(
        row_type=RowType.Stator,
        r=np.array([0.28]),
        R=287.15,
        gamma=1.4,
        Cp=1004.5,
    )
    row.rpm = 9549.297  # omega = rpm*pi/30 = 1000 rad/s
    row.Vx = np.array([200.0])
    row.Vt = np.array([150.0])
    row.Vr = np.array([0.0])
    row.Vm = np.array([200.0])
    row.T = np.array([288.0])
    row.P = np.array([90000.0])
    row.P0 = np.array([101325.0])
    row.total_massflow = 8.0
    return row


def _make_rotor() -> BladeRow:
    """A rotor whose radius genuinely differs from its upstream row's, so U2 != U1."""
    row = BladeRow(
        row_type=RowType.Rotor,
        r=np.array([0.32]),
        R=287.15,
        gamma=1.4,
        Cp=1004.5,
        omega=1000.0,
        total_area=0.05,
        P0_ratio=1.15,
        loss_function=FixedPressureLoss(0.0),
    )
    row.metal_exit_angle = [-30.0]
    return row


def test_ideal_exit_relative_pressure_is_isentropic_not_a_copy_of_inlet():
    upstream = _make_upstream()
    row = _make_rotor()

    rotor_calc(row, upstream, calculate_vm=True)

    assert row.r[0] != upstream.r[0]  # sanity: this row actually exercises U2 != U1
    assert row.T0R != pytest.approx(
        upstream.T0R, rel=1e-9
    )  # sanity: rothalpy moved T0R

    expected_ratio = (row.T0R / upstream.T0R) ** (row.gamma / (row.gamma - 1.0))
    assert row.P0R / upstream.P0R == pytest.approx(expected_ratio, rel=1e-9)
