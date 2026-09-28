# tests/test_enthalpy_loss.py
"""compressor_math.py handles LossType.Pressure / Polytropic / Entropy but has no
implementation for LossType.Enthalpy. It used to fall through to the bare
`else: row.Yp[:] = 0`, so an enthalpy model was never called and the row came out loss-free.

It now raises NotImplementedError instead. A model that carries parasitic work (disc friction,
recirculation, leakage) gets wording that says so, because a parasitic loss has no Yp
representation at all: it belongs in the denominator of efficiency, not in a pressure-loss
coefficient.

Rows are built directly, following the pattern in test_radial_velocity.py and
test_bracket_convergence.py. The enthalpy models are minimal stubs.
"""

import numpy as np
import pytest

from turbodesign.bladerow import BladeRow
from turbodesign.enums import RowType
from turbodesign.compressor_math import rotor_calc
from turbodesign.enums import LossType
from turbodesign.loss.losstype import LossBaseClass


def _make_upstream() -> BladeRow:
    """A converged upstream row feeding the rotor (meanline, one streamline).

    Same radius as the rotor exit (see _make_rotor): U1 == U2, the axial-safe
    case this library's core assumes, so rotor_calc's own Yp = 0 baseline is
    genuinely isentropic here (row.entropy_rise == 0) and eta_tt derived from
    row.T0_is/row.P0_is is meaningful. row.P0 is derived from T, P, T0 via the
    isentropic stagnation relation rather than picked independently -- an
    upstream state where P0 is inconsistent with T/P/V corrupts every
    downstream entropy/efficiency check without ever raising.
    """
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
    V = np.sqrt(row.Vx**2 + row.Vt**2 + row.Vr**2)
    row.T0 = row.T + V**2 / (2 * row.Cp)
    row.P0 = row.P * (row.T0 / row.T) ** (row.gamma / (row.gamma - 1))
    row.total_massflow = 8.0
    row.alpha1 = np.array([np.radians(30.0)])
    row.beta1 = np.array([np.radians(-20.0)])
    row.mu = 1.8e-5
    row.rho = np.array([1.1])
    return row


def _make_rotor(loss_function, phi_deg: float = 0.0) -> BladeRow:
    """A rotor row carrying the loss model under test."""
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
        loss_function=loss_function,
    )
    row.metal_exit_angle = [-30.0]
    row.mu = 1.8e-5
    row.tip_clearance = 0.02
    # A finite chord (chord defaults to axial_chord/cos(stagger); axial_chord
    # defaults to -1 and stagger to an array): give both an explicit, physically
    # ordinary value.
    row.axial_chord = 0.05
    row.stagger = 0.0
    return row


# --- LossType.Enthalpy models must raise -------------------------------------
#
# The raise fires in compressor_math's up-front loss-type dispatch, before the loss model is
# ever called, so the stubs do not need to compute anything.


class _InternalEnthalpyLoss(LossBaseClass):
    def __init__(self):
        super().__init__(LossType.Enthalpy)

    def __call__(self, row, upstream):
        return np.zeros_like(row.r)


class _ParasiticEnthalpyLoss(LossBaseClass):
    def __init__(self):
        super().__init__(LossType.Enthalpy, is_parasitic=True)

    def __call__(self, row, upstream):
        return np.zeros_like(row.r)


@pytest.mark.parametrize(
    "model_cls,expected_phrase",
    [
        (_InternalEnthalpyLoss, "not implemented here"),
        (_ParasiticEnthalpyLoss, "parasitic"),
    ],
)
def test_enthalpy_model_raises_not_implemented(model_cls, expected_phrase):
    upstream = _make_upstream()
    row = _make_rotor(model_cls())

    with pytest.raises(NotImplementedError, match=expected_phrase) as excinfo:
        rotor_calc(row, upstream, calculate_vm=True)
    assert model_cls.__name__ in str(excinfo.value)
