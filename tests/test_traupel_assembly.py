# tests/test_traupel_assembly.py
"""Structural checks on Traupel's efficiency assembly (issue: sum of two efficiencies).

The digitised Traupel figures are replaced by constants so the expected efficiency can be
worked out by hand. This pins the assembly, not the figures.
"""

from types import SimpleNamespace

import numpy as np
import pytest

from turbodesign.enums import RowType
from turbodesign.loss.turbine.traupel import Traupel

F = 0.02
ZETA_P = 0.03
ZETA_DELTA = 0.005
ZETA_CL = 0.01


def _traupel() -> Traupel:
    t = object.__new__(Traupel)
    # No 'Fig07': a KeyError would show it being read again.
    t.data = {
        "Fig01": lambda a, b: 1.0,
        "Fig02": lambda a, b: ZETA_P,
        "Fig03_0": lambda m: 1.0,
        "Fig04": lambda s, a: ZETA_DELTA,
        "Fig05": lambda s, a: 1.0,
        "Fig06": lambda ratio, turning: F,
        "Fig08": lambda c: ZETA_CL,
    }
    return t


def _stator() -> SimpleNamespace:
    return SimpleNamespace(
        row_type=RowType.Stator,
        r=np.array([0.10, 0.15]),
        pitch=0.02,
        te_pitch=0.02,
        throat=0.01,
        alpha1=np.radians([0.0, 0.0]),
        alpha2=np.radians([70.0, 70.0]),
        beta2=np.radians([-30.0, -30.0]),
        M=np.array([0.7, 0.7]),
        V=np.array([300.0, 300.0]),
        W=np.array([250.0, 250.0]),
    )


def _rotor() -> SimpleNamespace:
    """Rotor with a different pitch and span from the stator, so a mixed-up row shows."""
    return SimpleNamespace(
        row_type=RowType.Rotor,
        r=np.array([0.10, 0.14]),
        pitch=0.015,
        te_pitch=0.02,
        throat=0.01,
        beta1=np.radians([-30.0, -30.0]),
        beta2=np.radians([-60.0, -60.0]),
        M_rel=np.array([0.6, 0.6]),
        W=np.array([400.0, 400.0]),
        tip_clearance=0.0005,
    )


def test_stator_returns_zeros():
    row = _stator()
    out = _traupel()(row, _stator())
    assert out == pytest.approx(np.zeros(2), abs=0)


def test_rotor_returns_one_efficiency_from_combined_losses():
    out = _traupel()(_rotor(), _stator())

    zeta_f_stator = 0.5 * (0.05 / 0.25) ** 2
    zeta_f_rotor = 0.5 * (0.04 / 0.24) ** 2
    zeta_s = F * 0.02 / 0.05  # stator pitch / stator span
    zeta_r = F * 0.015 / 0.04  # rotor pitch / rotor span
    zeta_stator = ZETA_P + ZETA_DELTA + zeta_f_stator + zeta_s
    zeta_rotor = ZETA_P + ZETA_DELTA + zeta_f_rotor + zeta_r + ZETA_CL

    assert out == pytest.approx(np.full(2, 1.0 - zeta_stator - zeta_rotor), rel=1e-12)
    assert np.all((out > 0.0) & (out < 1.0))
