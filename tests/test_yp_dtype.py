# tests/test_yp_dtype.py
"""BladeRow.Yp defaults to np.array([0]) -- an int64 array, since the literal
has no decimal point. compressor_math.py writes every trial Yp *in place*
(`row.Yp[:] = y`, in both the Polytropic and Entropy root-finds in
stator_calc/rotor_calc), and an in-place write preserves the array's existing
dtype: a fractional trial value assigned into an int64 array is silently
floor/truncated to 0.

    >>> row = BladeRow(...)
    >>> row.Yp.dtype
    dtype('int64')
    >>> row.Yp[:] = 0.35
    >>> row.Yp
    array([0])

The consequence: every `minimize_scalar` search over Yp in the Polytropic and
Entropy branches is minimizing a CONSTANT objective (Yp is always 0 inside the
solve, no matter what trial value the optimizer proposes), so the solved
state is silently lossless regardless of the requested target. This is
independent of the LossType.Enthalpy work in tests/test_enthalpy_loss.py --
which happens to depend on the fix (its own root-find writes `row.Yp[:] = y`
too) -- but the defect itself sits one level below that, in the field
default, and would have broken Polytropic/Entropy on its own.

turbine_spool.py's own Enthalpy branch is unaffected: it *rebinds* the
attribute (`row.Yp = Yp`) rather than writing into the existing array, so it
never inherits the stale int64 dtype.
"""

import numpy as np
import pytest

from turbodesign.bladerow import BladeRow
from turbodesign.enums import RowType
from turbodesign.loss.fixedpolytropic import FixedPolytropicEfficiency
from turbodesign.compressor_math import rotor_calc


def test_yp_defaults_to_a_floating_dtype():
    row = BladeRow(row_type=RowType.Rotor, r=np.array([0.3048]))
    assert np.issubdtype(row.Yp.dtype, np.floating)


def test_inplace_yp_write_does_not_truncate_a_fractional_value():
    """This is exactly the write compressor_math.py's root-finds perform
    (`row.Yp[:] = y`) on every trial value the optimizer proposes."""
    row = BladeRow(row_type=RowType.Rotor, r=np.array([0.3048]))
    row.Yp[:] = 0.35
    assert row.Yp == pytest.approx(0.35)


def _make_upstream() -> BladeRow:
    """A thermodynamically self-consistent upstream row (meanline, one
    streamline): P0 is derived from T, P, T0 via the isentropic stagnation
    relation rather than picked independently, so eta_poly computed from the
    solved pi/tau ratios is meaningful."""
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
    return row


def test_polytropic_rotor_solves_a_genuinely_nonzero_yp():
    """End-to-end proof that the Polytropic root-find is not optimizing a
    constant objective.

    row.eta_poly is set directly (rather than relying on
    FixedPolytropicEfficiency.__call__'s return value) only to route around
    an unrelated pre-existing issue in compressor_math.py's eager top-level
    `float(loss_fn(row, upstream))` call, which cannot convert a
    multi-element array to a scalar; row.eta_poly being already truthy makes
    that branch read `float(np.mean(row.eta_poly))` instead. This does not
    touch the Yp root-find itself, which is what is under test here.

    A polytropic efficiency of 0.9 is well below what this rotor achieves
    loss-free, so it can only be reached with Yp > 0. Pre-fix (Yp truncating
    to 0 on every trial), the search below raises RuntimeError: the row's
    eta_poly is constant regardless of the proposed Yp, so the optimizer
    cannot approach 0.9 and parks on the bracket's upper edge instead.
    """
    upstream = _make_upstream()
    row = BladeRow(
        row_type=RowType.Rotor,
        r=np.array([0.3048]),
        phi=np.array([0.0]),
        R=287.15,
        gamma=1.4,
        Cp=1004.5,
        omega=1000.0,
        total_area=0.05,
        P0_ratio=1.15,
        loss_function=FixedPolytropicEfficiency(eta_poly=0.9),
    )
    row.metal_exit_angle = [-30.0]
    row.eta_poly = np.array([0.9])

    rotor_calc(row, upstream, calculate_vm=True)  # must not raise

    assert np.all(row.Yp > 0)
    assert float(row.eta_poly) == pytest.approx(0.9, abs=1e-3)
