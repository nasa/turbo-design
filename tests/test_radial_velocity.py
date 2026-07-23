# tests/test_radial_velocity.py
"""Drives rotor_calc directly with a non-zero meridional flow angle (phi).

Every shipped example (both Mattingly cases) constructs its Passage with
zero_phi=True, so phi is identically 0 there and sin(phi) == 0 makes both the
correct and the buggy form of Vr agree. That means the golden-file suite can
never catch a Vr regression at compressor_math.py:261. This test builds a
rotor BladeRow with phi != 0 directly and calls rotor_calc on it so the buggy
line actually executes.

Closure identity (see tests/test_velocity_triangle_oracle.py and
radial-equilibrium-derivation.md, Section 11 / Section 13's code-mapping
table): Vr = Vm*sin(phi), Vx = Vm*cos(phi), so Vm**2 == Vx**2 + Vr**2.
"""

import numpy as np
import pytest

from turbodesign.bladerow import BladeRow
from turbodesign.enums import RowType
from turbodesign.loss.fixedpressureloss import FixedPressureLoss
from turbodesign.compressor_math import rotor_calc


def _make_upstream() -> BladeRow:
    """A converged upstream row feeding the rotor (meanline, one streamline)."""
    row = BladeRow(
        row_type=RowType.Stator,
        r=np.array([0.3048]),
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


def _make_rotor(phi_deg: float) -> BladeRow:
    """A rotor row with a genuinely non-zero meridional flow angle."""
    row = BladeRow(
        row_type=RowType.Rotor,
        r=np.array([0.3048]),
        phi=np.array([np.radians(phi_deg)]),
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


def test_rotor_radial_velocity_uses_meridional_not_relative_velocity():
    """compressor_math.py:261 must build Vr from Vm, not W.

    Three other sites in the same file already do this correctly (:85, :177,
    :365), and radial-equilibrium-derivation.md documents it explicitly:
    Vr = Vm*sin(phi) (Section 11; Section 13's code-mapping table maps it to
    `Vr = Vm*np.sin(phi)` literally). Line 261 is the one site that instead
    writes `Vr = W*np.sin(phi)`, a frame error: W is the relative velocity,
    a different vector with a different magnitude than Vm.
    """
    upstream = _make_upstream()
    row = _make_rotor(phi_deg=20.0)

    rotor_calc(row, upstream, calculate_vm=True)

    assert row.phi[0] != 0.0  # sanity: the bug is invisible at phi == 0
    assert row.Vr == pytest.approx(row.Vm * np.sin(row.phi), rel=1e-9)


def test_rotor_velocity_triangle_closes_off_axis():
    """Vm**2 == Vx**2 + Vr**2 -- an identity, not a correlation.

    See tests/test_velocity_triangle_oracle.py for the independent oracle
    that derives this same closure without importing turbodesign.
    """
    upstream = _make_upstream()
    row = _make_rotor(phi_deg=20.0)

    rotor_calc(row, upstream, calculate_vm=True)

    assert row.Vm**2 == pytest.approx(row.Vx**2 + row.Vr**2, rel=1e-9)
