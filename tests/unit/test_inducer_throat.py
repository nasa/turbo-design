"""Slice S3 (docs/centrifugal/17-tdd-plan.md) -- the inducer throat goes live.

THE DEFECT. ``ImpellerEntranceDiffusionAungier.delta_h``
(``turbodesign/centrifugal/losses.py``) used to set ``Wth := W1xi`` -- "consistent
with :class:`ImpellerChokeAungier`'s SAME approximation (throat == LE annulus)". That
justification expired the moment ``ImpellerChokeAungier`` was given a real throat
(``Impeller.throat_area = 0.020525 m2``, NASA Appendix C) and started using it: this
class was left behind, with ``state.A_th`` populated and sitting unused right next to
it. ``W1xi - Wth`` was identically 0, so this loss was identically 0 at EVERY
operating point, forever -- a term that cannot fire cannot be caught being wrong.

THE CLOSURE. Meroni, Zuhlsdorf, Elmegaard & Haglind (2018), *Applied Energy*
232:139-156, Eqs. (3)-(5) -- continuity + rothalpy + isentropic LE -> throat, needing
only ``A_th`` -- the same system Kosuge, Ito & Nakanishi (1982) pose as "three
nonlinear algebraic equations" (Kosuge's own primary is PAYWALLED and NOT HELD; cited
here VIA Meroni 2018 / Kovář et al. 2021 Eqs. (23)-(24)).
:func:`turbodesign.centrifugal.losses.throat_relative_velocity` implements it as a 1-D
subsonic-branch root solve, same bisection discipline as
``solve_vm_from_continuity``/``_solve_for_meridional_velocity``.

****************************************************************************
PRE-REGISTRATION OUTCOME (docs/centrifugal/15-experiments-prereg.md E2): REFUTED.
****************************************************************************

The plan pre-registered Wth = 217.4 +- 2 m/s and Delta h_edf = 0 +- 50 J/kg at the
HECC design point ("the throat is matched to 0.3%"). The HONEST closure, using
``state.U1`` -- documented as "blade speed at the RMS streamline" and the SAME
quantity :class:`~turbodesign.centrifugal.losses.ImpellerIncidenceConrad`,
:class:`~turbodesign.centrifugal.losses.ImpellerMixingAungier` and
:class:`~turbodesign.centrifugal.losses.ImpellerChokeAungier` all already use for
"W1xi" -- gives Wth = 202.3 m/s and Delta h_edf = 248.3 J/kg at design. NOT a null.

Tracked down, not adjusted to reach a target (docs/PHYSICS-RULES.md: "the prediction
can be wrong; the geometry cannot"). Reconstructing the pre-registered number with the
ARITHMETIC MEAN of the hub/shroud blade speeds instead of the RMS radius --
``U1_mean = (U1h+U1s)/2 = 165.4 m/s`` in place of ``state.U1 = 183.3 m/s`` -- gives
W1xi = 216.1 m/s (matches the pre-registered 216.4 to 0.1%) and, carrying that same
mean-radius U1 through the rothalpy closure, Wth = 211.4 m/s and a Kosuge ratio of
1.330 (matches the pre-registered 1.33 essentially exactly). All three of the
plan's headline numbers reproduce under ONE consistent explanation: the "independent
solve" that produced them used the mean-radius blade speed, not the RMS
("area-representative for an annulus") radius this codebase has used for the identical
symbol since slice S1.

Using a THIRD U1 convention for just this one closure -- while the guard comparison
(``W1xi > Wth``) and the ``- Delta h_inc`` term inside the SAME ``delta_h`` keep using
the RMS ``state.U1`` -- would be a worse defect than either convention alone: three
different "W1" values feeding one loss term is exactly the kind of cross-station
proxy conflation this project keeps finding (docs/PHYSICS-RULES.md). So this
implementation keeps the codebase's existing, established RMS convention, and the
disagreement with the pre-registered null is reported here instead of hidden.

CONSEQUENCE. The design point genuinely moves: psi 0.8268 -> 0.8259 (-0.11 pt),
stage PR 5.0040 -> 4.9914 (-0.02, i.e. +6.82% -> +6.55% vs NASA), stage eta_poly,realgas
0.8670 -> 0.8665. tests/unit/test_ledger_truth.py is re-pinned to these values in this
slice's commit, with this file as the explanation. Note the direction: this real,
sourced internal loss makes agreement marginally BETTER, not worse -- the opposite of
slices S1/S2. That is not tuning (nothing here was chosen to produce this direction);
it is what a previously-dead term does once it can finally fire.

THE SIGN GUARD, `# UNVERIFIED`. Kovář et al. Eq. (22) as printed is unguarded on sign and
fires on ACCELERATION too (Wth > W1xi, toward choke) -- a double count against
:class:`ImpellerChokeAungier`'s already-owned near-choke physics. This class applies
the loss ONLY when ``W1xi > Wth``. Whether Aungier's own primary (Aungier 1995, *J.
Turbomach.* 117(3):360-366 -- PAYWALLED, NOT HELD) restricts it so is UNVERIFIED.
"""

from __future__ import annotations

import math
import sys
import warnings
from pathlib import Path

import pytest
from scipy.optimize import brentq

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests" / "fixtures"))

import turbodesign.centrifugal.losses as L  # noqa: E402
import turbodesign.centrifugal.solver as S  # noqa: E402
from turbodesign.centrifugal import InletState  # noqa: E402
from hecc_stage import BACKSWEEP, MDOT_DESIGN, P01, RPM, T01, build  # noqa: E402


def _captured_state(mdot: float) -> L.ImpellerLossState:
    """Capture the REAL ``ImpellerLossState`` ``Stage.solve`` builds internally at
    ``mdot`` -- same monkeypatch technique documented in
    tests/unit/test_diffusion_factor_provenance.py's recipe comment. Imports the
    driver (``tests/fixtures/hecc_stage.py``) rather than duplicating its geometry --
    docs/centrifugal/17-tdd-plan.md's own lesson (S1 postmortem): a fixture that
    copies configuration cannot fail for the reason it exists.

    Returns the LAST constructed state (the parasitic-block site, ``state2`` in
    ``solver.py``) -- bit-identical to the internal-loss site's own state for every
    field this file reads (Vm1, U1, T1, rho1, A_th are all LE-only quantities, fixed
    before the impeller-exit Cm2 trial loop even starts, so which construction site
    supplies them is immaterial).
    """
    states: list = []
    old_init = L.ImpellerLossState.__init__

    def capturing_init(self, *args, **kwargs):
        old_init(self, *args, **kwargs)
        states.append(self)

    L.ImpellerLossState.__init__ = capturing_init
    S.ImpellerLossState = L.ImpellerLossState
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            build(BACKSWEEP).solve(mdot=mdot, rpm=RPM, inlet=InletState(P0=P01, T0=T01))
    finally:
        L.ImpellerLossState.__init__ = old_init
        S.ImpellerLossState = L.ImpellerLossState
    return states[-1]


def _dh_edf_at(mdot: float) -> float:
    """The entrance-diffusion loss Stage.solve actually reports at ``mdot`` --
    imports the driver rather than reconstructing an ``ImpellerLossState`` by hand, so
    this exercises the SAME code path (``EvaluatedLosses.internal``) every other
    slice's tests use.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        op = build(BACKSWEEP).solve(
            mdot=mdot, rpm=RPM, inlet=InletState(P0=P01, T0=T01)
        )
    return op.losses.internal["ImpellerEntranceDiffusionAungier"]


def _independent_Wth(state: L.ImpellerLossState) -> float:
    """A SECOND, independently-coded solve of the identical Meroni 2018 Eqs. (3)-(5)
    system (continuity + rothalpy + isentropic, constant-cp ideal gas), using scipy's
    ``brentq`` directly rather than :func:`turbodesign.centrifugal.losses.throat_relative_velocity`'s
    own bisection -- so the RED test below is not simply re-asserting that function's
    arithmetic back at itself.
    """
    fluid = state.fluid
    gamma, cp, R = fluid.gamma, fluid.cp, fluid.R
    W1xi = math.hypot(state.Vm1, state.U1)
    T0R1 = state.T1 + W1xi**2 / (2.0 * cp)

    def residual(W: float) -> float:
        T_th = T0R1 - W * W / (2.0 * cp)
        rho_th = state.rho1 * (T_th / state.T1) ** (1.0 / (gamma - 1.0))
        return rho_th * W * state.A_th - state.mdot

    W_star = math.sqrt(2.0 * gamma * R * T0R1 / (gamma + 1.0))
    return brentq(residual, 1.0, W_star - 1e-6, xtol=1e-10)


def test_throat_relative_velocity_at_hecc_design_point():
    """Wth = 202.3 m/s at the HECC design point -- NOT the pre-registered 217.4 +- 2.

    See this file's module docstring for the full accounting: the pre-registered
    figure reconstructs almost exactly (216.1 vs 216.4 quoted) if the rothalpy's U1
    term uses the arithmetic MEAN of the hub/shroud blade speeds instead of the RMS
    ("area-representative") radius ``state.U1`` -- which is what this codebase has
    used for the identical "W1xi" symbol since slice S1
    (``ImpellerIncidenceConrad``/``ImpellerMixingAungier``/``ImpellerChokeAungier``
    all read ``state.U1``). Verified TWICE: once against this module's own
    bisection, once against an independently-coded ``brentq`` solve of the same
    system.
    """
    state = _captured_state(MDOT_DESIGN)
    Wth = L.throat_relative_velocity(state)
    Wth_independent = _independent_Wth(state)

    assert Wth == pytest.approx(Wth_independent, rel=1e-6), (
        "this module's bisection must agree with an independently-coded brentq solve "
        "of the identical Meroni system"
    )
    assert Wth == pytest.approx(202.3, abs=1.0)

    # The mean-radius reconstruction of the pre-registered number -- recorded so the
    # disagreement above is traceable to an actual, checkable cause rather than
    # asserted from nowhere.
    U1_mean = 0.5 * (state.U1h + state.U1s)
    W1xi_mean_radius = math.hypot(state.Vm1, U1_mean)
    W1xi_rms = math.hypot(state.Vm1, state.U1)
    assert W1xi_mean_radius == pytest.approx(216.1, abs=0.5)
    assert W1xi_rms != pytest.approx(W1xi_mean_radius, rel=1e-2), (
        "sanity check: the RMS and mean-radius W1xi must genuinely differ -- if they "
        "did not, the disagreement analysis above would be vacuous"
    )


def test_entrance_diffusion_is_zero_at_design():
    """PRE-REGISTERED NULL, REFUTED. The plan predicted Delta h_edf = 0 +- 50 J/kg at
    design ("the throat is matched to 0.3%"); the honest, RMS-consistent closure
    gives 248.3 J/kg -- comfortably over the plan's own 200 J/kg refutation
    threshold.

    This is reported, not hidden, and not adjusted to reach the predicted null (see
    module docstring). The mechanism is genuine: Wth (202.3 m/s) is real deceleration
    below W1xi (230.1 m/s), because A_th (0.020525 m2) is ~7.9% LARGER than the LE
    relative-flow area A1*cos(beta1) (0.01903 m2) under the RMS convention -- so the
    passage genuinely widens, not narrows, from the LE to the throat at this
    operating point.

    E-R11/C3 (docs/centrifugal/42-r11-result.md): 248.3 -> 319.8 J/kg. This is a PURE
    f_c SCALING and nothing else. The RAW loss is UNCHANGED to all printed digits
    (OhLossSet(apply_head_loss_correction=False) gives 248.256 J/kg, and
    248.256 * f_c(HECC) = 248.256 * 1.2883 = 319.83) -- as it must be: every input this
    closure reads (Vm1, U1, T1, rho1, A_th) is an LE/throat quantity, fixed before the
    impeller-exit trial loop C1/C2 act on even starts. The refutation stands either way:
    the loss is far above the plan's 200 J/kg bar at f_c = 1 as well as at f_c = 1.2883.
    """
    dh_edf = _dh_edf_at(MDOT_DESIGN)
    assert dh_edf == pytest.approx(319.8, abs=5.0)  # E-R11/C3: f_c * 248.3; was 248.3
    assert dh_edf > 200.0, (
        "the pre-registered null (<= 200 J/kg) is REFUTED at the HECC design point -- "
        "see this file's module docstring for why, and that this is a closure/"
        "convention finding, not a coding defect (verified against an independent "
        "brentq solve in test_throat_relative_velocity_at_hecc_design_point)"
    )


def test_entrance_diffusion_fires_off_design():
    """IDENTICALLY ZERO before this slice, at every operating point. Now fires, and
    grows monotonically as flow falls below design -- the qualitative prediction the
    plan made (docs/centrifugal/15-experiments-prereg.md E2) DOES hold, even though
    the design-point null does not: 849.0 J/kg at 90% flow, 1461.1 at 80%, 2146.6 at
    70% -- all comfortably above the plan's ">300 J/kg at 80% flow" bar.

    E-R11/C3 (docs/centrifugal/42-r11-result.md): every pin below is the OLD pin times
    f_c, and nothing else. VERIFIED by re-solving with
    ``OhLossSet(apply_head_loss_correction=False)``, which returns the old numbers to
    four significant figures at all three points:

        flow   raw (f_c OFF)   scaled (shipped)   implied f_c
        90%       665.735          849.023           1.2753
        80%      1155.242         1461.108           1.2648
        70%      1708.907         2146.647           1.2562

    f_c is a STATE-DEPENDENT correction (Aungier Eqs. 31-33), so it is NOT one constant:
    it falls from 1.2883 at design to 1.2562 at 70% flow. The raw closure is untouched by
    C1/C2 -- it reads only LE/throat quantities.
    """
    dh_90 = _dh_edf_at(MDOT_DESIGN * 0.90)
    dh_80 = _dh_edf_at(MDOT_DESIGN * 0.80)
    dh_70 = _dh_edf_at(MDOT_DESIGN * 0.70)

    assert dh_80 > 300.0, "the plan's own refutation bar: >300 J/kg at 80% flow"
    assert dh_90 == pytest.approx(849.0, rel=0.05)  # E-R11/C3: f_c * 665.7
    assert dh_80 == pytest.approx(1461.1, rel=0.05)  # E-R11/C3: f_c * 1155.2
    assert dh_70 == pytest.approx(2146.6, rel=0.05)  # E-R11/C3: f_c * 1708.9
    assert 0.0 < dh_90 < dh_80 < dh_70, (
        "the loss must grow monotonically as flow falls further below design -- more "
        "positive incidence, more LE->throat deceleration"
    )


def test_entrance_diffusion_does_not_fire_on_acceleration():
    """THE GUARD. Kovář et al. Eq. (22) as printed is unguarded on sign and fires on
    ACCELERATION too (Wth > W1xi, i.e. toward choke) -- a double count against
    ImpellerChokeAungier's already-owned near-choke physics. `# UNVERIFIED` against
    Aungier's own primary (paywalled) whether HIS equation restricts it so; this
    class applies the loss ONLY when W1xi > Wth.

    At 113% of design flow the throat genuinely ACCELERATES the relative flow
    (Wth > W1xi -- checked directly below, not just inferred from the zero loss):
    the passage's fixed A_th becomes the tighter constraint as mdot rises, reversing
    the design-point deceleration. Unguarded, ``0.4*(W1xi-Wth)**2`` would return a
    real, non-trivial number here; guarded, it is exactly 0.0.

    ⚠️ THE PROBE MOVED 1.15 -> 1.13, AND THE SANITY BAR 50 -> 30 J/kg. NEITHER IS A
    CONVENIENCE, AND BOTH ARE REPORTED (docs/centrifugal/45-jansen-skinfriction-fix.md).
    The skin-friction W-bar fix raised the internal loss sum ~3.5x on this term, which lowers
    rho2, raises the Cm2 continuity must find, and so LOWERS THE FLOW AT WHICH THE MODEL
    CHOKES -- the same mechanism E-R11/C3 documented, but larger:

        pre-W-bar-fix (shipped):  choke onset at 5.7899 kg/s  (1.1752 x design)
        post-W-bar-fix:           choke onset at 5.6101 kg/s  (1.1387 x design)
                                  -------------------------------------------------
                                  shift: -0.1798 kg/s, -3.11%

    1.15 x design (5.6659 kg/s) now sits PAST the onset and raises "no subsonic solution".
    The largest flow the model can reach AT ALL is 1.1387 x design, where the unguarded form
    returns only ~38 J/kg -- so the old 50 J/kg bar is now UNREACHABLE ANYWHERE in the
    solvable range, and holding it would only assert that the model chokes. The bar is
    therefore re-set to 30 J/kg at the 1.13 probe (34.1 J/kg there), and pinned additionally
    to a PHYSICAL scale that cannot drift: it must be a non-trivial fraction of the
    design-point entrance-diffusion loss. The GUARD ITSELF IS UNTOUCHED; what shrank is the
    accelerating margin the model has before it chokes, which is a result, not a test artefact.

    NOTE THE DIRECTION, as E-R11 did: the model chokes ~10% HIGH against NASA's measured
    5.24 kg/s, so this moves choke flow TOWARD the measurement (5.79 -> 5.61 kg/s, still
    +7.1% high). Choke flow is set by throat area, P02/sqrt(T02) and gamma -- a loss
    correction moving it 3% does NOT close a 10% structural gap and must not be read as one.
    """
    mdot = MDOT_DESIGN * 1.13
    state = _captured_state(mdot)
    W1xi = math.hypot(state.Vm1, state.U1)
    Wth = L.throat_relative_velocity(state)
    assert Wth > W1xi, (
        "this operating point must be on the ACCELERATING side to test the guard"
    )

    dh_edf = _dh_edf_at(mdot)
    assert dh_edf == 0.0, (
        "the guard must suppress the loss exactly, not merely shrink it"
    )

    unguarded = 0.4 * (W1xi - Wth) ** 2
    assert unguarded > 30.0, (
        "sanity check: the UNGUARDED Kovář et al. Eq. (22) form must be materially non-zero "
        "here, or this test would not be exercising the guard at all"
    )
    # ... and "materially" is anchored to PHYSICS, not to a bare J/kg that can quietly rot:
    # the suppressed loss must be a non-trivial fraction of what this same closure returns at
    # the design point (~320 J/kg), or the guard is being tested on noise.
    dh_edf_design = _dh_edf_at(MDOT_DESIGN)
    assert unguarded > 0.10 * dh_edf_design, (
        f"the UNGUARDED loss here ({unguarded:.1f} J/kg) is under 10% of the design-point "
        f"entrance-diffusion loss ({dh_edf_design:.1f} J/kg) -- too small for this test to be "
        "exercising the guard against anything that would matter"
    )


def test_kosuge_stall_ratio_at_design():
    """Kosuge, Ito & Nakanishi (1982) inducer-stall criterion W1s/Wth -- cited VIA
    Meroni 2018 / Kovář et al. 2021 Eqs. (23)-(24) (Kosuge's own primary is
    PAYWALLED and NOT HELD; the 0.5/1.75 constants are CITED-ONLY, not verified).
    Diagnostic only -- NOT wired into any loss.

    1.390 at design, using the RMS-consistent Wth -- not the pre-registered 1.33
    +- 0.05. Reproduces almost exactly (1.330) if Wth is instead computed with the
    mean-radius U1 (211.4 m/s) -- the SAME convention difference documented in
    test_throat_relative_velocity_at_hecc_design_point, now propagated through this
    ratio. W1s itself (the shroud LE relative velocity) is an EXACT geometric
    quantity, not a convention choice -- only Wth's denominator moves.

    Comfortably below the 1.75 stall threshold either way: this is a design-point
    finding, not a stall.
    """
    state = _captured_state(MDOT_DESIGN)
    ratio = L.kosuge_stall_ratio(state)

    assert ratio == pytest.approx(1.390, abs=0.02)
    assert ratio == pytest.approx(
        state.W1s / L.throat_relative_velocity(state), rel=1e-12
    )
    assert ratio < 1.75, (
        "HECC's design point must sit below Kosuge's own stall threshold"
    )


def test_inducer_loss_is_continuous_through_the_guard():
    """Sweep mass flow across the deceleration/acceleration switch (design -> 115%):
    Delta h_edf must fall monotonically to exactly 0.0 with NO discontinuous jump and
    NO reappearance past the crossover -- the guard must not double-count against
    ImpellerChokeAungier on one side of the switch while leaving a step on the other.

    ⭐ THE SWEEP WAS SHORTENED AGAIN, 1.15 -> 1.13, BY THE SKIN-FRICTION (W-bar) FIX --
    SAME MECHANISM, LARGER: choke onset 5.7899 -> 5.6101 kg/s (1.1752 -> 1.1387 x design,
    -3.11%), so the old last point (1.15 x design = 5.6659 kg/s) is now past the onset.
    See test_entrance_diffusion_does_not_fire_on_acceleration's docstring for the full
    numbers and docs/centrifugal/45-jansen-skinfriction-fix.md. The guard's crossover
    (~1.105) and the zero plateau past it are BOTH still inside the trimmed sweep, so nothing
    this test exists to check has been dropped. The history below is the E-R11/C3 precedent.

    THE SWEEP WAS SHORTENED, 1.18 -> 1.15, AND THAT IS A REPORTED PHYSICS RESULT, NOT A
    CONVENIENCE (E-R11/C3, docs/centrifugal/42-r11-result.md). Aungier's head-loss
    correction f_c raises the internal loss sum ~29%, which lowers rho2, which raises the
    Cm2 the continuity solve must find -- so THE MODEL NOW CHOKES AT A LOWER MASS FLOW.
    Measured by bisecting the largest mdot that still returns a solution:

        f_c OFF (pre-E-R11):  choke onset at 5.8644 kg/s  (1.1903 x design)
        f_c ON  (shipped):    choke onset at 5.7899 kg/s  (1.1752 x design)
                              ------------------------------------------------
                              shift: -0.0746 kg/s, -1.27%

    The old sweep's last point (1.18 x design = 5.8137 kg/s) now sits PAST the onset and
    raises "no subsonic solution ... the row is CHOKED". The sweep is therefore trimmed to
    the range where the model has a solution at all; the guard's crossover (~1.105) and the
    zero plateau past it are both still inside it, so nothing this test exists to check has
    been dropped. NOTE the DIRECTION: the model already choked ~9-10% HIGH against NASA's
    measured 5.24 kg/s, so this shift moves choke flow marginally TOWARD the measurement --
    but choke flow is set by throat area, P02/sqrt(T02) and gamma, so a ~1% move from a
    loss correction does not close a ~10% structural gap and must not be read as one.
    """
    fractions = [1.00, 1.02, 1.04, 1.06, 1.08, 1.09, 1.10, 1.11, 1.12, 1.13]
    values = [_dh_edf_at(MDOT_DESIGN * f) for f in fractions]

    assert all(v >= 0.0 for v in values), "the loss must never go negative"
    # Monotonically non-increasing as flow rises (design -> acceleration side).
    for prev, nxt in zip(values, values[1:]):
        assert nxt <= prev + 1e-9, (
            f"the inducer loss must fall monotonically as flow rises through the "
            f"guard switch: {values}"
        )
    # No jump larger than one sweep step's worth of the design-point scale -- the
    # loss is quadratic in (W1xi - Wth), so steps near the design point are
    # naturally the largest in the sweep; the bound below is set from that shape,
    # not tuned to this particular data.
    #
    # E-R11/C3: the bound is DELIBERATELY LEFT AT 90.0, not widened. f_c scales the whole
    # curve up ~28%, so the largest step (design -> 1.02) grew 70.1 -> 89.8 J/kg and now
    # sits just inside it. If a future change pushes it over, the honest move is to add a
    # sweep point at 1.01 (the curvature is highest here), NOT to raise this number.
    for prev, nxt in zip(values, values[1:]):
        assert prev - nxt < 90.0, f"discontinuous jump detected in sweep: {values}"
    # Once it reaches zero, it must STAY zero (the guard, not a coincidental root).
    first_zero = next(i for i, v in enumerate(values) if v == 0.0)
    assert all(v == 0.0 for v in values[first_zero:]), (
        f"the loss must stay exactly 0.0 for every point past the guard switch: {values}"
    )
