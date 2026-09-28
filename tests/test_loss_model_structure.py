# tests/test_loss_model_structure.py
"""Structural checks on the turbine loss models and shared loss code.

The digitised chart lookups are replaced by constants, so these tests pin how the pieces are
assembled and called, not the chart values.

Covered:
  * On numpy >= 2.5, ``float()`` of a 1-D array raises, and ``if row.M < 0.5`` is ambiguous
    with more than one streamline. Ainley-Mathieson, Craig-Cox and Kacker-Okapuu each crashed
    on a multi-streamline row, as did the Polytropic / Entropy target in compressor_math.
  * Kacker-Okapuu took the tip-clearance gap as the blade height.
  * TD2_Reynolds_Correction divided by the cumulative (0 at the hub) massflow.
  * FixedPressureLoss returned a scalar and LossInterp modified the caller's array.
  * TurbineSpool silently ignored Polytropic / Entropy loss models.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from turbodesign.bladerow import BladeRow
from turbodesign.compressor_math import rotor_calc
from turbodesign.enums import LossType, RowType
from turbodesign.loss.fixedpolytropic import FixedPolytropicEfficiency
from turbodesign.loss.fixedpressureloss import FixedPressureLoss
from turbodesign.loss.losstype import mean_value
from turbodesign.loss.turbine.ainleymathieson import AinleyMathieson
from turbodesign.loss.turbine.craigcox import CraigCox
from turbodesign.loss.turbine.kackerokapuu import KackerOkapuu
from turbodesign.loss.turbine.TD2 import TD2_Reynolds_Correction

N = 5


def _model(cls, data):
    """A loss model with stubbed chart lookups (no TD3_HOME, no downloads)."""
    m = object.__new__(cls)
    m.data = data
    m._loss_type = LossType.Pressure
    return m


def _const(value):
    return lambda *args: value


def test_mean_value_accepts_scalars_and_arrays():
    assert mean_value(2.0) == 2.0
    assert mean_value(np.array([1.0, 2.0, 3.0])) == 2.0
    assert isinstance(mean_value(np.array([1.0])), float)


# ---------------------------------------------------------------- Kacker-Okapuu

KO_DATA = {
    "Fig01_beta0": _const(0.03),
    "Fig02": _const(0.05),
    "Fig04": _const(0.2),
    "Fig14_Impulse": _const(0.02),
    "Fig14_Axial_Entry": _const(0.02),
}


def _full(v):
    return np.full(N, float(v))


def _ko_rows(row_type, tip_clearance=0.02, r_tip=0.35):
    r = np.linspace(0.30, r_tip, N)
    common = dict(
        r=r,
        chord=_full(0.05),
        axial_chord=0.04,
        pitch=_full(0.04),
        pitch_to_chord=_full(0.8),
        te_pitch=_full(0.02),
        throat=_full(0.02),
        rho=_full(5.0),
        T=_full(1400.0),
        gamma=1.33,
        tip_clearance=tip_clearance,
    )
    if row_type == RowType.Stator:
        row = SimpleNamespace(
            row_type=RowType.Stator,
            beta1_metal=_full(10.0),
            alpha1=_full(np.radians(10.0)),
            alpha2=_full(np.radians(70.0)),
            M=_full(0.8),
            V=_full(300.0),
            **common,
        )
    else:
        row = SimpleNamespace(
            row_type=RowType.Rotor,
            beta1=_full(0.2),
            beta1_metal=_full(10.0),
            beta2_metal=_full(-60.0),
            beta2=_full(-1.0),
            M_rel=_full(0.6),
            W=_full(350.0),
            **common,
        )
    upstream = SimpleNamespace(M=_full(0.5), M_rel=_full(0.5), gamma=1.33)
    return row, upstream


@pytest.mark.parametrize("row_type", [RowType.Stator, RowType.Rotor])
def test_kacker_okapuu_runs_on_a_multi_streamline_row(row_type):
    row, upstream = _ko_rows(row_type)
    out = _model(KackerOkapuu, KO_DATA)(row, upstream)
    assert out.shape == (N,)
    assert np.all(np.isfinite(out))


def test_kacker_okapuu_stator_has_a_secondary_loss():
    """The stator used h = 0 for the blade height, which switched its secondary loss off."""
    model = _model(KackerOkapuu, KO_DATA)
    short = model(*_ko_rows(RowType.Stator, r_tip=0.35))
    tall = model(*_ko_rows(RowType.Stator, r_tip=0.40))
    assert not np.allclose(short, tall)


def test_kacker_okapuu_tip_clearance_only_moves_the_clearance_loss():
    """Ytc scales as tip_clearance**0.78 (kprime**0.78); the secondary loss must not depend
    on the clearance. It did, because the clearance gap was used as the blade height."""
    model = _model(KackerOkapuu, KO_DATA)
    y0 = model(*_ko_rows(RowType.Rotor, tip_clearance=0.0))
    y1 = model(*_ko_rows(RowType.Rotor, tip_clearance=0.02))
    y2 = model(*_ko_rows(RowType.Rotor, tip_clearance=0.04))
    assert np.all(y1 > y0)
    assert (y2 - y0) / (y1 - y0) == pytest.approx(np.full(N, 2**0.78), rel=1e-9)


# ------------------------------------------------------------ Ainley-Mathieson

AM_DATA = {
    "Fig05": _const(0.9),
    "Fig08": _const(0.02),
    "Fig04a": _const(0.04),
    "Fig04b": _const(0.06),
}


def _am_rows(row_type, mach):
    r = np.linspace(0.30, 0.35, N)
    common = dict(
        r=r,
        chord=_full(0.05),
        pitch=_full(0.04),
        throat=_full(0.02),
        tip_clearance=0.02,
        beta1_metal=_full(10.0),
    )
    if row_type == RowType.Stator:
        row = SimpleNamespace(
            row_type=RowType.Stator,
            alpha1=_full(0.1),
            alpha2=_full(1.2),
            M=_full(mach),
            **common,
        )
    else:
        row = SimpleNamespace(
            row_type=RowType.Rotor,
            beta1=_full(0.2),
            beta2=_full(-1.1),
            M=_full(mach),
            M_rel=_full(mach),
            **common,
        )
    upstream = SimpleNamespace(r=np.linspace(0.30, 0.34, N))
    return row, upstream


@pytest.mark.parametrize("row_type", [RowType.Stator, RowType.Rotor])
@pytest.mark.parametrize("mach", [0.3, 0.7])
def test_ainley_mathieson_runs_on_a_multi_streamline_row(row_type, mach):
    row, upstream = _am_rows(row_type, mach)
    out = _model(AinleyMathieson, AM_DATA)(row, upstream)
    assert out.shape == (N,)
    assert np.all(np.isfinite(out))


# ------------------------------------------------------------------- Craig-Cox

CC_DATA = {
    "Fig03": _const(1.0),
    "Fig04": _const(20.0),
    "Fig05": _const(0.1),
    "Fig06_delta_Xpt": _const(0.01),
    "Fig06_Npt": _const(1.0),
    "Fig07": _const(1.2),
    "Fig08": _const(0.01),
    "Fig09": _const(0.01),
    "Fig15": _const(1.0),
    "Fig17": _const(1.0),
    "Fig18": _const(0.05),
}


def _cc_rows():
    r = np.linspace(0.30, 0.35, N)

    def common():
        return dict(
            r=r,
            chord=_full(0.05),
            camber=_full(0.06),
            pitch=_full(0.04),
            throat=_full(0.02),
            te_pitch=_full(0.02),
            aspect_ratio=_full(1.0),
            rho=_full(5.0),
            T=_full(1400.0),
            T0=_full(1500.0),
            beta1_metal=_full(20.0),
            beta1_fixed=False,
        )

    upstream = SimpleNamespace(
        row_type=RowType.Stator,
        alpha1=_full(0.1),
        alpha2=_full(1.2),
        beta1=_full(0.1),
        M=_full(0.8),
        V=_full(300.0),
        **common(),
    )
    row = SimpleNamespace(
        row_type=RowType.Rotor,
        beta1=_full(0.3),
        beta2=_full(-1.1),
        M_rel=_full(0.6),
        W=_full(350.0),
        Cp=1150.0,
        **common(),
    )
    row.T0 = _full(1450.0)
    return row, upstream


def test_craig_cox_runs_on_a_multi_streamline_rotor():
    m = _model(CraigCox, CC_DATA)
    m.C = 1 / (200 * 32.2 * 778.16)
    out = m(*_cc_rows())
    assert out.shape == (N,)
    assert np.all(np.isfinite(out))


# ------------------------------------------------------------------------- TD2


def test_td2_reynolds_correction_has_no_hub_spike():
    """row.massflow is cumulative from the hub, so it is 0 at the hub streamline."""
    model = object.__new__(TD2_Reynolds_Correction)
    model.TD2 = lambda row, upstream: np.full(N, 0.05)
    row = SimpleNamespace(
        r=np.linspace(0.30, 0.35, N),
        massflow=np.linspace(0.0, 2.0, N),
        total_massflow=2.0,
        mu=_full(5e-5),
    )
    out = model(row, None)
    assert np.all(np.isfinite(out))
    assert out == pytest.approx(np.full(N, out[0]), rel=1e-12)
    assert 0.5 * 0.05 < out[0] < 2.0 * 0.05


# ---------------------------------------------------------- shared array contract


def test_fixed_pressure_loss_returns_an_array_for_one_streamline():
    row = SimpleNamespace(r=np.array([0.3]), percent_hub_shroud=np.array([0.5]))
    model = FixedPressureLoss(np.array([0.05, 0.07, 0.09]))
    out = model(row, None)
    assert isinstance(out, np.ndarray) and out.shape == (1,)
    assert out == pytest.approx([0.07])


def _stub_interp(is_xy):
    """A LossInterp with the spline fits replaced by simple functions."""
    from turbodesign.lossinterp import LossInterp

    interp = object.__new__(LossInterp)
    interp._name = "stub"
    interp.logX10 = False
    interp.is_xy = is_xy
    interp.c_min, interp.c_max = 1.0, 4.0
    interp.fxc_min = lambda c: 0.0
    interp.fxc_max = lambda c: 1.0
    # The scalar path indexes the bisplev-style 2-D result as [0][0].
    interp.func = (lambda x: x) if is_xy else (
        lambda x, c: np.array([[x * c]]) if np.ndim(x) == 0 else x * c
    )
    return interp


def test_loss_interp_does_not_modify_the_callers_array():
    x = np.array([-1.0, 0.5, 2.0])  # below and above the chart's x range
    before = x.copy()
    with pytest.warns(UserWarning, match="outside chart range"):
        y = _stub_interp(is_xy=False)(x, 2.0)
    assert np.array_equal(x, before)
    assert y == pytest.approx([0.0, 1.0, 2.0])  # clamped to [0, 1], times c = 2


def test_loss_interp_accepts_an_integer_x():
    assert _stub_interp(is_xy=True)(1) == pytest.approx(1.0)
    assert _stub_interp(is_xy=False)(1, 2.0) == pytest.approx(2.0)


# ------------------------------------------------------- compressor_math targets


def test_polytropic_target_from_an_array_valued_loss_model_does_not_crash():
    """compressor_math did float(loss_fn(row, upstream)), which raises for an array on
    numpy >= 2.5. test_yp_dtype.py presets row.eta_poly to avoid exactly this."""
    from tests.test_yp_dtype import _make_upstream

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

    rotor_calc(row, upstream, calculate_vm=True)  # eta_poly deliberately not preset

    assert float(row.eta_poly) == pytest.approx(0.9, abs=1e-3)


# ---------------------------------------------------------------- turbine solver


def test_turbine_spool_rejects_polytropic_loss_models():
    from tests.conftest import EXAMPLES, _run

    ns = _run(EXAMPLES / "radial-turbine" / "radial_turbine-1D.py")
    spool, rotor = ns["spool"], ns["rotor"]
    rotor.loss_function = FixedPolytropicEfficiency(eta_poly=0.9)
    with pytest.raises(NotImplementedError, match="FixedPolytropicEfficiency"):
        spool.initialize()
