"""Item 5e -- the incidence loss, resurrected.

WHAT WAS WRONG
--------------
``ImpellerIncidenceConrad`` returned EXACTLY 0.0 J/kg at every operating point, forever.
Not approximately zero: identically zero, by construction. With no inducer blade-metal-
angle field on ``Impeller``, the model had to define its optimum-incidence angle as the
FLOW angle (``beta_opt := beta_1xi``), which makes ``W* == 0`` for any flow whatsoever.
``f_inc`` was a phantom parameter. The slice-5 tests "verified" the model was present by
checking its name was a key in a dict -- so they passed on a corpse.

That matters beyond bookkeeping: INCIDENCE IS WHAT MAKES EFFICIENCY ROLL OVER on a
speedline. Without it, slices 9-11 cannot produce a correct map -- eta just falls
monotonically, which is the symptom the whole centrifugal effort started from.

WHAT FIXED IT
-------------
NASA publishes no inducer blade metal angle (Table 2 gives LEAN angles, a different
quantity). But Appendix C tabulates the blade SURFACE at 11 spans, so the angle is
DERIVABLE: extract_hecc_blade_angles.py -> data/hecc/blade_angles.csv.
At the RMS inlet streamline (60.8% span), beta1b = 45.46 deg.

WHAT IS STILL WRONG -- READ THIS BEFORE TRUSTING THE MAGNITUDE
--------------------------------------------------------------
The loss is now alive and correctly V-shaped, and eta now peaks and rolls over. But the V
is centred at mdot ~ 6.0 kg/s, which is ~22% ABOVE the 4.927 design point. A real
impeller's shockless (minimum-incidence) point sits AT its design point -- that is what
"design" means. So the model believes HECC runs at ~8 deg of positive incidence at
design, which is not credible.

This was ONCE thought to be the same discrepancy as the trailing-edge angle miss
(derived 28.6-37.7 deg vs NASA's stated 32-42). It is NOT: docs/centrifugal/11-blade-
angle-discrepancy.md digitized NASA's own Figure 18 and RESOLVED the TE miss -- it was
a fillet-corruption artefact in the Appendix C extraction, the TE metal angle is really
-36.0 deg (span-average), and Figure 18 also CONFIRMS the LE angles used here to ~1 deg.
So this incidence-minimum-off-design finding is a SEPARATE, still-open item: per
docs/centrifugal/11-blade-angle-discrepancy.md section 5, "not in the blade coordinates"
-- the candidate causes are the LE velocity triangle (we solve with ZERO inlet blockage,
which understates Vm1 and so overstates the flow angle) or Conrad's beta_opt, which
collapsed to the flow angle's identity because LE blade thickness t1 was 0. t1 is now
sourced (le_blade_thickness=0.00418982 m, pinned to the inducer throat station -- see
extract_hecc_le_thickness.py), which lifts beta_opt from 45.46 to 49.19 deg.
Until the shockless point is verified to land at design, the incidence loss MAGNITUDE at
design is not trustworthy -- only its SHAPE is.

Reported, not tuned. Do not "fix" this by adjusting beta1b until the V lands on the
design point: that would be fitting the geometry to the answer, which is exactly the
error that made slice 5's psi gate meaningless.
"""

from __future__ import annotations

import math

import pytest

from turbodesign.centrifugal import (
    Air,
    Impeller,
    InletState,
    MeridionalPath,
    OhLossSet,
    Stage,
    WiesnerSlip,
)
from turbodesign.centrifugal.losses import ImpellerIncidenceConrad, ImpellerLossState

R2 = 0.21581
RPM = 21789.0
MDOT_DESIGN = 4.9269
BETA1B_RMS = 45.46  # deg, data/hecc/blade_angles.csv, RMS streamline (60.8% span)
L_MAIN_M = 0.209875
L_SPLITTER_M = 0.145580

# Exit blade metal angle: NASA Figure 18 (CR-2014-218114/REV1, PDF p.36), digitized,
# span-averaged. Supersedes -32.95 (fillet-corrupted Appendix C extraction) and -37.5
# (unsourced midpoint of NASA's prose band). See docs/centrifugal/11-blade-angle-
# discrepancy.md.
BACKSWEEP_DEG = -36.0

# LE blade thickness, pinned to the inducer throat station (same cut that gives
# throat_area). Conrad's beta_opt collapses to the flow angle's identity when t1 = 0;
# this value lifts it from 45.46 to 49.19 deg. extract_hecc_le_thickness.py.
LE_BLADE_THICKNESS = 0.00418982

# MAIN blades only. The splitters start at r = 0.0675 m, downstream of the inducer eye
# (r_hub = 0.0380 m), so THEY CANNOT BLOCK THE INLET. Conrad's inlet blockage must use
# this, not Z_eff = 25.4 (ImpellerLossState.Z_le).
N_BLADES = 15


def _stage(path, beta1b=BETA1B_RMS, t1=LE_BLADE_THICKNESS):
    return Stage(
        path=path,
        impeller=Impeller(
            n_blades=N_BLADES,
            n_splitters=15,
            backsweep_deg=BACKSWEEP_DEG,
            r_te=R2,
            splitter_le_r=0.0675,
            tip_clearance=0.000305,
            inducer_blade_angle_deg=beta1b,
            le_blade_thickness=t1,
            l_main_m=L_MAIN_M,
            l_splitter_m=L_SPLITTER_M,
        ),
        diffuser=None,
        slip=WiesnerSlip(),
        losses=OhLossSet(),
        fluid=Air(),
    )


@pytest.fixture(scope="module")
def path(data_dir):
    return MeridionalPath.from_csv(
        data_dir / "hecc" / "flowpath_hub.csv",
        data_dir / "hecc" / "flowpath_shroud.csv",
    )


def _dh_inc(path, mdot):
    op = _stage(path).solve(
        mdot=mdot, rpm=RPM, inlet=InletState(P0=101325.0, T0=288.15)
    )
    return op.losses.internal["ImpellerIncidenceConrad"]


def test_incidence_loss_is_not_structurally_zero(path):
    """The regression that matters: it used to be 0.0 at EVERY point, by construction."""
    assert _dh_inc(path, MDOT_DESIGN) > 0.0


def test_incidence_loss_rises_on_BOTH_sides_of_its_minimum(path):
    """A V, not a slope.

    This is the property that makes eta roll over. A loss that only ever rises with
    falling flow would deepen the surge side without ever producing a peak. The minimum
    must be interior, with the loss climbing away from it in both directions.
    """
    flows = [4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
    dh = [_dh_inc(path, m) for m in flows]
    i_min = min(range(len(dh)), key=lambda i: dh[i])

    assert 0 < i_min < len(dh) - 1, (
        f"the incidence minimum is at the edge of the sweep (index {i_min}), so no V is "
        f"demonstrated: {list(zip(flows, [round(x, 1) for x in dh]))}"
    )
    assert dh[i_min - 1] > dh[i_min] < dh[i_min + 1], "not a strict minimum"
    # and it must climb substantially, not just wobble
    assert dh[0] > 10.0 * dh[i_min] or dh[0] > 100.0
    assert dh[-1] > 10.0 * dh[i_min] or dh[-1] > 100.0


def test_incidence_is_zero_when_the_blade_matches_the_flow_AT_ZERO_BLOCKAGE(path):
    """The model's own boundary case -- but ONLY at zero blade blockage.

    Set the blade angle TO the relative flow angle, with t1 = 0, and the loss must vanish.
    This pins the sign/datum convention: if beta1b were measured from tangential rather
    than meridional, or with the opposite sign, this would fail.

    THE ``t1 = 0`` IS LOad-BEARING, not incidental. Conrad's optimum-incidence angle is
    BLOCKAGE-CORRECTED: beta_opt = arctan[tan(beta1b) / (1 - Z*t1/(pi*D1))]. The blades
    block part of the inlet circumference, so the flow speeds up meridionally on entering
    the passage, and a stream that arrives ALIGNED WITH THE METAL therefore ends up at
    negative incidence inside it. The zero-loss direction is beta_opt, NOT beta1b -- they
    coincide only when t1 = 0.

    This test asserted the t1 = 0 identity while the fixture silently supplied a real
    thickness, and it FAILED (190 J/kg) the moment a measured t1 was wired in. The test
    was wrong, not the model: "blade matches flow => no incidence loss" is a zero-blockage
    statement being smuggled in as a general one.
    """
    op = _stage(path, t1=0.0).solve(
        mdot=MDOT_DESIGN, rpm=RPM, inlet=InletState(P0=101325.0, T0=288.15)
    )
    beta_flow_deg = math.degrees(
        math.atan2(op.stations.impeller_le.U, op.stations.impeller_le.Vm)
    )

    op2 = _stage(path, beta1b=beta_flow_deg, t1=0.0).solve(
        mdot=MDOT_DESIGN, rpm=RPM, inlet=InletState(P0=101325.0, T0=288.15)
    )
    assert op2.losses.internal["ImpellerIncidenceConrad"] == pytest.approx(
        0.0, abs=1e-6
    )


def test_with_blockage_the_zero_loss_direction_is_beta_opt_not_the_metal_angle(path):
    """The general statement the test above is the t1 = 0 corner of.

    With a real LE thickness the loss vanishes at ``beta_opt``, which sits ABOVE the metal
    angle. Aligning the blade with the flow instead leaves a real, non-zero loss -- and
    that is physically correct, not a defect.
    """
    op = _stage(path).solve(
        mdot=MDOT_DESIGN, rpm=RPM, inlet=InletState(P0=101325.0, T0=288.15)
    )
    beta_flow_deg = math.degrees(
        math.atan2(op.stations.impeller_le.U, op.stations.impeller_le.Vm)
    )

    # blade aligned with the flow, but WITH blockage -> NOT optimum, so NOT zero
    aligned = _stage(path, beta1b=beta_flow_deg).solve(
        mdot=MDOT_DESIGN, rpm=RPM, inlet=InletState(P0=101325.0, T0=288.15)
    )
    # ~62 J/kg at the throat-derived t1 (12.3% inlet blockage on the 15 MAIN blades).
    # NB: this read 190 J/kg while the loss blocked the inlet with Z_eff = 25.4 -- counting
    # splitters that do not reach the inducer eye. See ImpellerLossState.Z_le.
    assert aligned.losses.internal["ImpellerIncidenceConrad"] > 20.0, (
        "with inlet blade blockage, a blade aligned with the FLOW is not at optimum "
        "incidence -- Conrad's beta_opt is blockage-corrected above the metal angle"
    )

    # the blade angle whose beta_opt EQUALS the flow angle -> loss must vanish
    eye = path.inducer_eye()
    D1 = 2.0 * math.hypot(eye.r_hub, eye.r_shroud) / math.sqrt(2.0)
    blockage = (
        N_BLADES * LE_BLADE_THICKNESS / (math.pi * D1)
    )  # MAIN blades only at the LE
    beta1b_star = math.degrees(
        math.atan(math.tan(math.radians(beta_flow_deg)) * (1.0 - blockage))
    )

    at_opt = _stage(path, beta1b=beta1b_star).solve(
        mdot=MDOT_DESIGN, rpm=RPM, inlet=InletState(P0=101325.0, T0=288.15)
    )
    assert at_opt.losses.internal["ImpellerIncidenceConrad"] == pytest.approx(
        0.0, abs=1e-6
    )


def test_a_missing_blade_angle_RAISES_rather_than_returning_zero(path):
    """The anti-defect-D guard.

    A loss model with no blade angle cannot compute incidence. It must SAY SO, not
    quietly return 0.0 -- a silent zero is indistinguishable from "no loss here", and
    that is exactly how this model stayed dead through five slices with green tests.
    """
    state = ImpellerLossState(
        fluid=Air(),
        mdot=4.9,
        omega=2281.6,
        Z=25.4,
        backsweep_deg=BACKSWEEP_DEG,
        r1h=0.0405,
        r1s=0.108,
        U1=186.0,
        Vm1=139.0,
        T1=278.5,
        rho1=1.2,
        r2=R2,
        b2=0.0155,
        tip_clearance=0.000305,
        U2=492.0,
        Cm2=110.0,
        Vt2=430.0,
        W2=200.0,
        T2=430.0,
        beta1b_deg=None,  # <-- the whole point
    )
    with pytest.raises(ValueError, match="inducer blade metal angle"):
        ImpellerIncidenceConrad().delta_h(state)


@pytest.mark.xfail(
    strict=True,
    reason=(
        "KNOWN, UNRESOLVED (item 5e, open). The incidence minimum sits at mdot ~ 6.0 "
        "kg/s, ~22% ABOVE the 4.927 design point -- i.e. the model believes HECC runs at "
        "~8 deg positive incidence at its own design point, which is not credible. A real "
        "impeller is shockless AT design. Same signature as the unresolved ~5 deg gap "
        "between the derived TE metal angles (28.6-37.7) and NASA's stated band (32-42). "
        "Candidate causes: (a) the derived beta1b is systematically low; (b) we solve with "
        "ZERO inlet blockage, understating Vm1 and so overstating the flow angle. "
        "DO NOT close this by tuning beta1b until the V lands on the design point -- that "
        "is fitting geometry to the answer, the exact error that made the psi gate "
        "meaningless. Resolve the discrepancy, then delete this marker."
    ),
)
def test_the_impeller_is_shockless_at_its_design_point(path):
    """A well-designed inducer has its minimum incidence AT the design point."""
    flows = [4.5, 4.75, MDOT_DESIGN, 5.25, 5.5, 6.0, 6.5]
    dh = [_dh_inc(path, m) for m in flows]
    i_min = min(range(len(dh)), key=lambda i: dh[i])
    assert flows[i_min] == pytest.approx(MDOT_DESIGN, rel=0.05), (
        f"incidence minimum at mdot = {flows[i_min]} vs design {MDOT_DESIGN} "
        f"({100 * (flows[i_min] / MDOT_DESIGN - 1):+.0f}%)"
    )
