"""The tests that would have caught it.

Slice 5 shipped with THREE of its ten loss models contributing IDENTICALLY ZERO, and the
tests that "verified" them passed anyway -- because they only asserted that the model's
name was a key in the results dict::

    assert "ImpellerIncidenceConrad" in op.losses.internal      # passes on 0.0 J/kg

That is the upstream ``LossType.Enthalpy -> Yp = 0`` silent-zero failure mode
(``docs/PHYSICS-RULES.md`` defect D, ``docs/centrifugal/00-root-cause-analysis.md``) reproduced
INSIDE the module written to prevent it -- with a test certifying the dead loss as
present. A model that contributes nothing is worse than one that errors: it is a knob the
tests believe in.

Two guards, both of which FAIL on the code as it stands. That is the point: they are
written to expose the current state, not to describe it.

1. NON-VACUOUS: every loss named as present must actually contribute.
   Fails today for ImpellerIncidenceConrad (beta_opt := beta_1xi makes W* identically 0
   -- no inducer blade-metal-angle field existed), ImpellerChokeAungier and
   ImpellerEntranceDiffusionAungier (both keyed off a throat station that does not exist,
   A_th := A1, so their driving difference is 0 by construction).

2. NO PARASITIC DOMINATOR: the existing ``test_no_single_loss_dominates`` guards ONLY the
   internal set. Its own docstring says a >60% share "is probably mis-scaled -- a units
   error in a correlation looks exactly like this." The PARASITIC side had no such guard,
   which is how ImpellerRecirculationOh reached 90.1% of all parasitic work unremarked --
   and the entire psi = 0.81 "pass" rests on it (delete it: psi -> 0.776, FAIL).
   This is debt item D7.
"""

from __future__ import annotations

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

R2 = 0.21581
RPM = 21789.0
MDOT = 4.9269
# Inducer blade metal angle at the RMS streamline, from NASA's Appendix C blade
# coordinates (data/hecc/blade_angles.csv). Item 5e -- this is what makes
# ImpellerIncidenceConrad non-zero, and therefore what this file is checking.
BETA1B_RMS = 45.46

# LE blade thickness, pinned to the inducer throat station (same cut that gives
# throat_area). Conrad's beta_opt collapses to the flow angle's identity when t1 = 0;
# this value lifts it from 45.46 to 49.19 deg. extract_hecc_le_thickness.py.
LE_BLADE_THICKNESS = 0.00418982


@pytest.fixture(scope="module")
def op(data_dir):
    """HECC design point, full Oh-1997 set, published coefficients."""
    stage = Stage(
        path=MeridionalPath.from_csv(
            data_dir / "hecc" / "flowpath_hub.csv",
            data_dir / "hecc" / "flowpath_shroud.csv",
        ),
        impeller=Impeller(
            n_blades=15,
            n_splitters=15,
            # NASA Figure 18 (CR-2014-218114/REV1, PDF p.36), digitized, span-averaged.
            # Supersedes -37.5 (unsourced midpoint of NASA's prose band). See
            # docs/centrifugal/11-blade-angle-discrepancy.md.
            backsweep_deg=-36.0,
            r_te=R2,
            splitter_le_r=0.0675,
            tip_clearance=0.000305,
            inducer_blade_angle_deg=BETA1B_RMS,
            le_blade_thickness=LE_BLADE_THICKNESS,
        ),
        diffuser=None,
        slip=WiesnerSlip(),
        losses=OhLossSet(),
        fluid=Air(),
    )
    return stage.solve(mdot=MDOT, rpm=RPM, inlet=InletState(P0=101325.0, T0=288.15))


@pytest.mark.xfail(
    strict=True,
    reason=(
        "KNOWN DEFECT, not a flake. 3 of 10 losses are identically zero: "
        "ImpellerIncidenceConrad (no inducer blade-metal-angle field -> beta_opt := "
        "beta_1xi -> W* == 0; angle now derived in data/hecc/blade_angles.csv, wiring "
        "is item 5e), ImpellerChokeAungier + ImpellerEntranceDiffusionAungier (no throat "
        "station -> A_th := A1; debt D5, item 10). strict=True: when these are wired, "
        "this test MUST start passing and the marker must be removed."
    ),
)
def test_every_configured_loss_actually_contributes(op):
    """A loss model that is "present" but contributes 0.0 J/kg is a DEAD MODEL.

    Presence in the results dict is not evidence of physics -- it is evidence of a dict
    key. This asserts the thing the other tests only imply.
    """
    all_losses = {**op.losses.internal, **op.losses.parasitic}
    dead = sorted(name for name, dh in all_losses.items() if dh <= 0.0)
    assert not dead, (
        f"{len(dead)} configured loss model(s) contribute EXACTLY ZERO: {dead}. "
        "This is the silent-zero failure mode (docs/PHYSICS-RULES.md defect D). A model that "
        "cannot affect the answer must not be advertised as part of the loss set."
    )


@pytest.mark.xfail(
    strict=True,
    reason=(
        "KNOWN DEFECT (debt D7), not a flake. ImpellerRecirculationOh is 90.1% of all "
        "parasitic work, and the psi = 0.81 'pass' rests entirely on it (remove it: psi "
        "-> 0.776, FAIL). NOT expected to resolve via the Japikse two-zone model (item "
        "5c): docs/centrifugal/10-two-zone.md investigated it and its verdict is 'do not "
        "implement' -- two-zone mixing is INTERNAL (destroys pressure, no work term), "
        "recirculation is PARASITIC (adds work, no pressure); per PHYSICS-RULES rule 6 an "
        "internal model cannot substitute for a parasitic one at any coefficient value. "
        "D7 has NO known resolution in hand -- it is an open, reported defect. "
        "strict=True: when a real fix lands, remove the marker."
    ),
)
def test_no_parasitic_loss_dominates(op):
    """D7: the parasitic side needs the same dominance guard as the internal side.

    Same threshold and same rationale as test_no_single_loss_dominates: a term taking
    >60% of its category is probably mis-scaled, and a units error looks exactly like
    this. ImpellerRecirculationOh is at ~90.1% -- and it is an exponential amplifier,
    sinh(3.5*alpha2**3), evaluated at alpha2 = 76.4 deg, sinh(8.30) = 2015, where a
    1-degree error in a COMPUTED angle swings it by ~33%.
    """
    parasitic = op.losses.parasitic
    total = op.losses.parasitic_total
    assert total > 0.0, "no parasitic loss at all -- the set is not configured"
    worst, share = max(parasitic.items(), key=lambda kv: kv[1])
    frac = share / total
    assert frac <= 0.60, (
        f"{worst} is {100 * frac:.1f}% of all parasitic work ({share:.0f} of {total:.0f} "
        f"J/kg). A single term dominating its category is a mis-scaling signature, not "
        f"a physical result. Breakdown: "
        f"{ {k: round(v) for k, v in sorted(parasitic.items(), key=lambda kv: -kv[1])} }"
    )
