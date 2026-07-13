# tests/test_bracket_convergence.py
"""compressor_math.py solves for the rotor exit relative Mach number with

    minimize_scalar(calculate_vm_func, bounds=[0.01, 1], method="bounded")

`method="bounded"` cannot leave its bracket. A transonic axial fan or
front-stage compressor rotor can legitimately demand more massflow at a
given area/speed/pressure than is achievable at any relative Mach <= 1 --
the mass-flux-vs-Mach relation peaks at M=1 and falls off beyond it, so the
bounded search has nowhere left to go but the M=1 edge. It parks there,
`res.success` is `True`, and the unmatched massflow residual is silently
swallowed: the caller applies M_rel = 1 as though it were the converged
answer.

This file proves both directions of the fix:
  * a rotor whose true mass-flow-matching Mach sits well inside [0.01, 1]
    still solves normally, unchanged;
  * a rotor driven with a massflow demand the passage cannot support at any
    Mach <= 1 must raise instead of silently returning the M=1 edge.
"""

import numpy as np
import pytest

from turbodesign.bladerow import BladeRow
from turbodesign.enums import RowType
from turbodesign.loss.fixedpressureloss import FixedPressureLoss
from turbodesign.compressor_math import rotor_calc


def _make_upstream(total_massflow: float) -> BladeRow:
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
    row.total_massflow = total_massflow
    return row


def _make_rotor(total_area: float) -> BladeRow:
    """A rotor row with a small fixed exit area, so massflow demand can be
    pushed past what is achievable at relative Mach <= 1."""
    row = BladeRow(
        row_type=RowType.Rotor,
        r=np.array([0.3048]),
        phi=np.array([0.0]),
        R=287.15,
        gamma=1.4,
        Cp=1004.5,
        omega=1000.0,
        total_area=total_area,
        P0_ratio=1.15,
        loss_function=FixedPressureLoss(0.0),
    )
    row.metal_exit_angle = [-30.0]
    return row


def test_rotor_relative_mach_inside_bracket_solves_unchanged():
    """A modest massflow demand has its matching Mach well inside [0.01, 1]."""
    upstream = _make_upstream(total_massflow=8.0)
    row = _make_rotor(total_area=0.05)

    rotor_calc(row, upstream, calculate_vm=True)

    assert 0.01 < row.M_rel[0] < 0.99
    assert row.M_rel[0] == pytest.approx(0.4006, abs=1e-3)


def test_rotor_relative_mach_outside_bracket_raises_instead_of_pinning():
    """Push the massflow demand past what this area can pass at M_rel <= 1.

    Pre-fix, this silently returns M_rel pinned at the 1.0 edge of the
    bounded search, `res.success == True`, with a several-kg/s unmatched
    massflow residual that nothing surfaces. That is exactly the transonic
    condition described in the module docstring above: the true
    mass-flow-matching state is not reachable inside the bracket, so the
    bounded optimizer has nowhere to go but the edge. Post-fix this must
    raise instead of returning that edge value as though it were converged.
    """
    upstream = _make_upstream(total_massflow=15.0)
    row = _make_rotor(total_area=0.05)

    with pytest.raises(RuntimeError, match="relative Mach"):
        rotor_calc(row, upstream, calculate_vm=True)
