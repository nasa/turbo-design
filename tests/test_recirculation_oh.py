# tests/test_recirculation_oh.py
"""turbodesign/loss/compressor/otac.py's ImpellerRecirculationOh.__call__
wraps row.alpha2 in np.radians(...) before using it in the correlation.
row.alpha2 is already in radians -- see turbodesign/bladerow.py's "Flow
Angles (radians)" docstring (around line 99), and the metal_exit_angle
setter (bladerow.py ~line 567) that converts the input degrees to radians
once via np.radians and assigns the result straight to alpha2/beta2. The
only place alpha2 is converted back to degrees is the display boundary
(bladerow.py ~line 705, np.degrees(self.alpha2)).

So the model applies a second, spurious degrees->radians conversion on top
of an already-radians value. Since the argument is cubed inside
sinh(3.5 * x**3), the suppression compounds to roughly (pi/180)**3, and the
correlation is numerically flat regardless of flow angle for every
physically ordinary value of alpha2.

ImpellerRecirculationOh is parasitic (is_parasitic=True); compressor_math
raises NotImplementedError if it is wired up as a row's loss_function (see
test_enthalpy_loss.py). So this file drives the model directly rather than
through rotor_calc/stator_calc, and supplies the row state (U, T0) that the
solver would ordinarily have already set up before calling a loss_function.

Rows are built with the same fixtures used in test_enthalpy_loss.py.
"""

import numpy as np
import pytest

from turbodesign.loss.compressor import otac

from .test_enthalpy_loss import _make_rotor, _make_upstream


def _build(alpha2_deg: float):
    """upstream/row pair, driven directly (see module docstring)."""
    model = otac.ImpellerRecirculationOh()
    upstream = _make_upstream()
    row = _make_rotor(model)
    row.U = row.omega * row.r
    row.alpha2 = np.array([np.radians(alpha2_deg)])
    # T0 rise large enough that the model's internal cap never binds across
    # the angle range exercised below (see test_loss_grows_with_alpha2).
    row.T0 = upstream.T0 + 200.0
    return row, upstream, model


def _expected_dh(row, upstream, model, alpha2_rad: float) -> float:
    """Independent re-derivation of the Oh recirculation correlation -- not
    a call into the model under test -- with alpha2 supplied explicitly in
    radians by the caller.
    """
    r_tip_in = (
        model.radius_tip_inlet
        if model.radius_tip_inlet is not None
        else float(np.max(upstream.r))
    )
    r_exit = (
        model.radius_exit if model.radius_exit is not None else float(np.max(row.r))
    )
    W_row = float(np.linalg.norm(np.asarray(getattr(row, "W", row.V))))
    Kbl = 0.75 if model.splitter_le == 0 else 0.6
    Df = (
        1
        - W_row / max(model.surge_vrel, 1e-6)
        + Kbl
        * model.loading_coefficient
        / max(
            (model.surge_vrel / max(W_row, 1e-6))
            * (
                (model.number_of_blades / np.pi) * (1 - r_tip_in / max(r_exit, 1e-6))
                + 2 * r_tip_in / max(r_exit, 1e-6)
            ),
            1e-6,
        )
    )

    U = float(np.mean(row.U))
    dh = 8e-5 * np.sinh(3.5 * alpha2_rad**3) * Df**2 * U**2

    cp = float(np.mean([row.Cp, upstream.Cp]))
    cap = 0.5 * cp * (float(np.mean(row.T0)) - float(np.mean(upstream.T0)))
    return float(np.clip(dh * model.loss_modifier, 0.0, cap))


def test_matches_independent_correlation_with_alpha2_in_radians():
    """The model must reproduce the correlation evaluated with alpha2 taken
    as-is in radians -- no second conversion.
    """
    alpha2_deg = 60.0
    row, upstream, model = _build(alpha2_deg)

    actual = float(np.mean(model(row, upstream)))
    expected = _expected_dh(row, upstream, model, np.radians(alpha2_deg))

    assert actual == pytest.approx(expected, rel=1e-9)


def test_generated_suppression_factor_matches_pi_over_180_cubed():
    """Characterizes the magnitude of the defect via the correlation's own
    math, independent of the model under test: at a small angle,
    sinh(3.5*x**3) ~ 3.5*x**3 on both branches, so wrapping an
    already-radians x in a second np.radians(...) suppresses the result by
    ~(pi/180)**3. The ratio is generated here, not typed in from a brief.
    """
    alpha2_rad = np.radians(5.0)

    row, upstream, model = _build(5.0)
    correct = _expected_dh(row, upstream, model, alpha2_rad)
    doubly_converted = _expected_dh(row, upstream, model, np.radians(alpha2_rad))

    suppression_ratio = doubly_converted / correct
    assert suppression_ratio == pytest.approx((np.pi / 180) ** 3, rel=1e-3)


def test_loss_grows_with_alpha2():
    """The correlation must respond to the exit flow angle: sinh(3.5*x**3)
    is steeply increasing in x for the range exercised here. Pre-fix, the
    doubly-suppressed argument keeps sinh(...) pinned in its near-zero
    linear regime across this whole range, so the growth from a shallow to
    a steep exit angle is nowhere near this large.
    """
    angles_deg = [20.0, 40.0, 60.0, 80.0]
    dh_values = []
    for a in angles_deg:
        row, upstream, model = _build(a)
        dh_values.append(float(np.mean(model(row, upstream))))

    assert all(b > a for a, b in zip(dh_values, dh_values[1:]))
    assert dh_values[-1] / dh_values[0] > 1000
