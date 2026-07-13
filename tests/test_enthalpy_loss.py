# tests/test_enthalpy_loss.py
"""compressor_math.py handles LossType.Pressure / Polytropic / Entropy but had
no branch at all for LossType.Enthalpy -- every model in
turbodesign/loss/compressor/otac.py that declares LossType.Enthalpy fell
through to the bare `else: row.Yp[:] = 0`, and the model itself was never even
called. Sixteen models were dead code.

A later change attempted to fix this by root-finding a Yp that reproduces the
efficiency an internal loss implies. That attempt was itself wrong: it
compared the row's IDEAL exit state (row.T0_is, which derives from row.P0_is)
against its ACTUAL exit state as if that were a total-to-total efficiency.
The two are equal only at Yp = 0 (which is why a handful of smoke tests
passed), and diverge sharply away from it -- delivering several times the
loss the correlation actually asked for, with no way for either guard in that
loop to detect it (both read the same wrong quantity). This file's previous
centrepiece asserted `_eta_tt(...) == approx(0.95)`, which is tautological:
it checks that the root-find converged using the very formula the root-find
minimizes, so it cannot fail even when the formula itself is wrong.

The fix here is not a better conversion -- it is to stop converting. Every
LossType.Enthalpy model now raises NotImplementedError instead of silently
producing Yp = 0 (the original bug) or a silently-wrong Yp (the regression
this file now guards against). The five models that carry parasitic work
(disc friction, recirculation, leakage, or an aggregate of these) raise with
wording that says so explicitly, because a parasitic loss has no Yp
representation at all -- it belongs in the denominator of efficiency, not in
a pressure-loss coefficient. The other eleven raise because the internal ->
Yp conversion is simply not implemented.

Rows are built directly, following the pattern in test_radial_velocity.py and
test_bracket_convergence.py.
"""

import numpy as np
import pytest

from turbodesign.bladerow import BladeRow
from turbodesign.enums import RowType
from turbodesign.compressor_math import rotor_calc
from turbodesign.loss.compressor import otac


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
    # defaults to -1 and stagger to an array, which several otac correlations
    # cannot turn into a plain float): give both an explicit, physically
    # ordinary value.
    row.axial_chord = 0.05
    row.stagger = 0.0
    return row


# --- all sixteen LossType.Enthalpy models: every one must raise ------------
#
# Named explicitly, by class, with the parasitic/internal split hand-labelled
# here rather than read off the model's own is_parasitic flag or discovered
# by scanning the module -- otherwise this test would still pass even if
# otac.py's flags were wrong, which is exactly the failure mode being guarded
# against (see test_impeller_various_is_parasitic below for the concrete case
# that was wrong until this change).
#
# The raise fires in compressor_math's up-front loss-type dispatch, before
# the loss model is ever called, so none of these need the specific geometry
# some of their correlations would otherwise require (e.g. a genuine
# hub-to-shroud radius change for ImpellerClearanceJansen) -- solving never
# gets that far.

_ALL_ENTHALPY_MODELS = [
    ("diffuser_vaneless_stanitz", otac.DiffuserVanelessStanitz, False),
    ("impeller_blade_loading_aungier", otac.ImpellerBladeLoadingAungier, False),
    ("impeller_blade_loading_coppage", otac.ImpellerBladeLoadingCoppage, False),
    ("impeller_clearance_jansen", otac.ImpellerClearanceJansen, False),
    ("impeller_disc_friction_daily", otac.ImpellerDiscFrictionDaily, True),
    ("impeller_incidence_aungier", otac.ImpellerIncidenceAungier, False),
    ("impeller_incidence_conrad", otac.ImpellerIncidenceConrad, False),
    ("impeller_leakage_aungier", otac.ImpellerLeakageAungier, True),
    ("impeller_mixing_aungier", otac.ImpellerMixingAungier, False),
    ("impeller_mixing_johnston", otac.ImpellerMixingJohnston, False),
    ("impeller_prescribed", otac.ImpellerPrescribed, False),
    ("impeller_recirculation_aungier", otac.ImpellerRecirculationAungier, True),
    ("impeller_recirculation_oh", otac.ImpellerRecirculationOh, True),
    ("impeller_skin_friction_coppage", otac.ImpellerSkinFrictionCoppage, False),
    ("impeller_skin_friction_jansen", otac.ImpellerSkinFrictionJansen, False),
    ("impeller_various", otac.ImpellerVarious, True),
]


def test_enumerates_exactly_sixteen_models_five_parasitic():
    assert len(_ALL_ENTHALPY_MODELS) == 16
    assert sum(1 for _, _, is_parasitic in _ALL_ENTHALPY_MODELS if is_parasitic) == 5
    assert (
        sum(1 for _, _, is_parasitic in _ALL_ENTHALPY_MODELS if not is_parasitic) == 11
    )


@pytest.mark.parametrize(
    "name,model_cls,is_parasitic",
    _ALL_ENTHALPY_MODELS,
    ids=[m[0] for m in _ALL_ENTHALPY_MODELS],
)
def test_enthalpy_model_raises_not_implemented(name, model_cls, is_parasitic):
    upstream = _make_upstream()
    row = _make_rotor(model_cls())

    expected_phrase = "parasitic" if is_parasitic else "not implemented here"
    with pytest.raises(NotImplementedError, match=expected_phrase):
        rotor_calc(row, upstream, calculate_vm=True)


# --- regression: ImpellerVarious must be flagged parasitic ------------------
#
# ImpellerVarious.__call__ sums disc_friction + leakage + recirculation (all
# three parasitic) together with several internal terms into a single
# Enthalpy number. Before this change it was is_parasitic=False, so the (now
# removed) Yp root-find would have silently converted parasitic work into a
# pressure-loss coefficient. Pinning the flag directly, independent of the
# raise-wording test above, so a future edit cannot flip it back unnoticed.


def test_impeller_various_is_parasitic():
    assert otac.ImpellerVarious().is_parasitic is True
