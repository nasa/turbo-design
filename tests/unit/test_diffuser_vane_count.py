"""Slice S5 (docs/centrifugal/17-tdd-plan.md) -- the vaned diffuser's station-dependent
vane count.

THE BUG. ``VanedDiffuser._cascade`` builds ONE :class:`_Cascade` with
``n_vanes = self.n_vanes + self.n_splitters`` (20 + 20 = 40 for HECC) and uses it for
EVERYTHING: deviation, incidence, friction, solidity, diffusion factor.

But HECC's diffuser SPLITTER VANE LEADING EDGE sits at r = 0.2488 m, while the vaned
diffuser spans r3 = 0.23133 m -> r4 = 0.28459 m (data/hecc/vane_angles.csv, derived by
extract_hecc_vane_angles.py from NASA/CR-2014-218114/REV1 Appendix C Tables
C.25-C.28). Over the FIRST THIRD of the passage only 20 vanes exist -- exactly where the
circulation per vane is largest.

Solidity sigma = chord/pitch sits in Lieblein's diffusion-factor DENOMINATOR:

    Df = 1 - V_out/V_in + (Vt_in - Vt_out) / (2 * sigma * V_in)

Overstating sigma (40 vanes instead of ~20 at the LE) LOWERS Df, which LOWERS the loading
loss (dh_load = 0.05*Df**2*V_in**2). The shipped code UNDERSTATES the diffuser loss.

THE FIX -- station-dependent, same Z_le/Z_eff discipline as the impeller's Z_eff, one
component downstream:
  - DEVIATION (``_Cascade.exit_flow_angle_deg``, Carter's rule): keeps ALL 40 vanes. This
    is a TRAILING-EDGE quantity and the splitters DO reach the TE. LEGAL. Must not change.
  - FRICTION (the hydraulic diameter used by ``friction_loss``): keeps 40 -- the wetted
    passages exist over most of the chord.
  - DIFFUSION FACTOR / LOADING LOSS: spans LE -> TE, where the splitter is absent over the
    first third. Uses a LENGTH-WEIGHTED effective vane count (``VanedDiffuser.
    Z_eff_loading``), the same Aungier formula already used for the impeller (Yang, Liu &
    Zhao 2023, *Machines* 11(1):118, Eq. (1) -- NOT "Li et al.", see data/coefficients.md):

        Z_eff,vd = n_vanes + n_splitters * (L_splitter / L_main)

    with the vane chords from NASA Appendix C (extract_hecc_vane_angles.py):
    L_main = 0.0532632 m (2.097 in, diffuser main vane), L_splitter = 0.0358396 m
    (1.411 in, diffuser splitter vane) -> Z_eff,vd = 33.4575.

THE HARD CAP -- NOT A TOLERANCE. NASA MEASURED the diffusion system's total-pressure
loss: 7.17 % (CR-2014-218114/REV1, Table A.5, design-flow row: "Total pressure loss from
the diffuser to exit" = 7.168, at 10.860 lbm/s corrected -- the HECC design point). Here
"the diffuser to exit" means from the vaned-diffuser INLET (post-vaneless-space, r3) to
the STAGE EXIT (post-EGV) -- verified against the model's OWN pre-fix number: computed
as 100*(1 - P0_stage_exit/P0_diffuser_inlet), the pre-fix model gives 6.320 %, matching
the "~6.31 %" this slice's brief pre-registered almost exactly, which is what identifies
this as the right metric definition (docs/centrifugal/15-experiments-prereg.md E3,
docs/centrifugal/16-code-update-plan.md S5).

If the converged model's diffuser P0 loss EXCEEDS ~7.2 %, that is NOT success -- it is a
compensating-error signature (a fix manufacturing an error to absorb the impeller's
excess), the exact pattern that hid the -37.5 deg backsweep. This test is that cap.

ALSO REFUTED IF the vane-LE incidence moves: it is set entirely by the (unchanged) 40-vane
cascade and the (unchanged) vaneless-space exit state, so it must be BIT-IDENTICAL.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests" / "fixtures"))

from hecc_stage import (  # noqa: E402
    BACKSWEEP,
    MDOT_DESIGN,
    P01,
    RPM,
    T01,
    build,
    diffusion_system,
)

from turbodesign.centrifugal import InletState  # noqa: E402
from turbodesign.centrifugal.diffusion import VanedDiffuser  # noqa: E402

# NASA/CR-2014-218114/REV1 Appendix C, via extract_hecc_vane_angles.py --
# reproduce with:
#   uv run python extract_hecc_vane_angles.py
CHORD_MAIN_M = 0.0532632412  # Table C.25/C.26, diffuser main vane
CHORD_SPLITTER_M = 0.0358396286  # Table C.27/C.28, diffuser splitter vane
SPLITTER_LE_R_M = 0.248797826  # diffuser splitter vane LE radius

NASA_DIFFUSER_P0_LOSS_PCT = 7.168  # Table A.5, design-flow row (10.860 lbm/s corrected)


def _hecc_diffuser(**overrides) -> VanedDiffuser:
    """The real HECC vaned diffuser (``hecc_stage.diffusion_system()``'s second
    component), with overrides -- so the test imports the shipped geometry rather than
    duplicating it (docs/centrifugal/17-tdd-plan.md rule (b))."""
    base = diffusion_system()[1]
    import dataclasses

    return dataclasses.replace(base, **overrides)


def _solve_design_point():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return build(BACKSWEEP).solve(
            mdot=MDOT_DESIGN, rpm=RPM, inlet=InletState(P0=P01, T0=T01)
        )


def _diffuser_p0_loss_pct(op) -> float:
    """100 * (1 - P0_stage_exit / P0_diffuser_inlet). ``stage_states[0]`` is the
    vaneless-space exit (== vaned-diffuser inlet); ``stage_states[-1]`` is the EGV exit
    (== stage exit). This is the metric that reproduces NASA's Table A.5 "total pressure
    loss from the diffuser to exit" column -- verified above against the pre-fix 6.32 %."""
    inlet_P0 = op.stage_states[0].P0
    exit_P0 = op.stage_states[-1].P0
    return 100.0 * (1.0 - exit_P0 / inlet_P0)


def _vane_le_incidence_deg(op) -> float:
    diffuser = diffusion_system()[1]
    inlet = op.stage_states[0]
    return inlet.alpha_deg - diffuser.beta_le_deg


# --------------------------------------------------------------------- deviation: 40, legal


def test_deviation_uses_all_40_vanes_at_the_TE():
    """Carter's rule (a TE quantity) must keep using n_vanes + n_splitters = 40, even
    when the loading-loss count (Z_eff_loading) differs. Splitters DO reach the TE --
    this is legal and must NOT change."""
    d_with_chords = _hecc_diffuser(
        chord_splitter=CHORD_SPLITTER_M, splitter_le_r=SPLITTER_LE_R_M
    )
    d_without_chords = _hecc_diffuser()  # old call site: no splitter-chord data at all

    assert d_with_chords._cascade.n_vanes == 40
    assert d_without_chords._cascade.n_vanes == 40

    r_mean = 0.5 * (d_with_chords.r3 + d_with_chords.r4)
    angle_with = d_with_chords._cascade.exit_flow_angle_deg(r_mean)
    angle_without = d_without_chords._cascade.exit_flow_angle_deg(r_mean)

    assert angle_with == pytest.approx(angle_without), (
        "the exit flow angle (Carter deviation off the 40-vane TE cascade) must be "
        "INDEPENDENT of whether splitter-chord data is supplied for the loading-loss count"
    )
    # And Z_eff_loading (33.4575, NOT 40) must differ from what deviation uses.
    assert d_with_chords.Z_eff_loading != pytest.approx(40.0)
    assert d_with_chords.Z_eff_loading == pytest.approx(33.4575, abs=0.01)


# --------------------------------------------------------------------- loading: length-weighted


def test_loading_loss_uses_the_length_weighted_count():
    """Z_eff,vd = n_vanes + n_splitters*(L_splitter/L_main) ~= 33.4575, NOT 40.

    FAILS TODAY: the shipped ``_cascade`` (and the diffusion factor / loading loss built
    from it) has no length-weighted view at all -- everything reads n_vanes+n_splitters.
    """
    d = _hecc_diffuser(chord_splitter=CHORD_SPLITTER_M, splitter_le_r=SPLITTER_LE_R_M)

    assert d.Z_eff_loading == pytest.approx(33.4575, abs=0.01)
    assert d._cascade_loading.n_vanes == pytest.approx(d.Z_eff_loading)
    assert d._cascade_loading.n_vanes != d._cascade.n_vanes, (
        "the loading cascade's vane count must differ from the deviation/friction "
        "cascade's -- that is the entire point of this slice"
    )

    # Solidity in the DENOMINATOR: fewer effective vanes -> HIGHER Df -> MORE loading loss.
    r_mean = 0.5 * (d.r3 + d.r4)
    sigma_loading = d._cascade_loading.solidity(r_mean)
    sigma_deviation = d._cascade.solidity(r_mean)
    assert sigma_loading < sigma_deviation, (
        "the length-weighted (33.4575-vane) solidity must be LOWER than the 40-vane "
        "solidity -- fewer effective vanes at the LE, same chord and radius"
    )


# --------------------------------------------------------------------- THE CAP


def test_diffuser_p0_loss_does_not_exceed_the_measurement():
    """NASA MEASURED the diffusion system's total-pressure loss: 7.17 % (Table A.5,
    design-flow row). That is NOT a tolerance -- it is the measurement.

    REFUTED IF the converged model's diffuser P0 loss exceeds ~7.2 %: a fix that pushes
    past NASA's own number is not fixing the diffuser, it is manufacturing a compensating
    error to absorb the impeller's excess -- the exact pattern that hid the -37.5 deg
    backsweep for the life of this project. The correct outcome is a PARTIAL improvement
    that leaves most of the gap open.
    """
    op = _solve_design_point()
    loss_pct = _diffuser_p0_loss_pct(op)

    assert loss_pct <= 7.2, (
        f"diffuser P0 loss = {loss_pct:.3f}% EXCEEDS NASA's measured {NASA_DIFFUSER_P0_LOSS_PCT}% "
        "(Table A.5) -- this is an OVERSHOOT FINDING, not a success. Report it; do not tune it down."
    )
    # And it must actually have MOVED toward the measurement (not be a no-op slice).
    PRE_FIX_LOSS_PCT = 6.320
    assert loss_pct > PRE_FIX_LOSS_PCT, (
        "the length-weighted vane count must RAISE the diffuser loss relative to the "
        "pre-fix (40-vane-everywhere) baseline of 6.32% -- Df rises when sigma falls"
    )


# --------------------------------------------------------------------- splitterless: bit-identical


def test_splitterless_diffuser_is_bit_identical():
    """n_splitters=0 (Eckardt-scale, unsplittered) -> Z_eff_loading == n_vanes exactly,
    regardless of whether chord_splitter/splitter_le_r are supplied -- there is no
    splitter to length-weight. The two cascades collapse to the same vane count, and the
    solved output must be bit-identical whether or not splitter-chord data is present."""
    d_no_data = VanedDiffuser(
        r3=0.231331,
        r4=0.284455,
        b=0.014199,
        n_vanes=20,
        n_splitters=0,
        beta_le_deg=74.0,
        beta_te_deg=50.0,
    )
    d_with_stray_data = VanedDiffuser(
        r3=0.231331,
        r4=0.284455,
        b=0.014199,
        n_vanes=20,
        n_splitters=0,
        beta_le_deg=74.0,
        beta_te_deg=50.0,
        chord_splitter=0.01,  # must be IGNORED -- no splitters exist
        splitter_le_r=0.25,
    )

    assert d_no_data.Z_eff_loading == pytest.approx(20.0)
    assert d_with_stray_data.Z_eff_loading == pytest.approx(20.0)
    assert d_no_data._cascade_loading.n_vanes == d_no_data._cascade.n_vanes

    from turbodesign.centrifugal.diffusion import state_from_totals
    from turbodesign.centrifugal.state import Air

    fluid = Air()
    inlet = state_from_totals(
        P0=520000.0,
        T0=475.0,
        Vm=110.0,
        Vt=440.0,
        r=0.21581,
        b=0.015467,
        fluid=fluid,
        s=0.0,
    )
    out_no_data = d_no_data.solve(inlet, 4.9269, fluid)
    out_with_stray_data = d_with_stray_data.solve(inlet, 4.9269, fluid)

    assert out_no_data.P0 == pytest.approx(out_with_stray_data.P0, rel=1e-12)
    assert out_no_data.T0 == pytest.approx(out_with_stray_data.T0, rel=1e-12)
    assert out_no_data.Vt == pytest.approx(out_with_stray_data.Vt, rel=1e-12)


# --------------------------------------------------------------------- collateral damage checks


def test_vane_le_incidence_stays_near_nasas_measurement():
    """The vane-LE incidence is set by the (unchanged) 40-vane cascade's beta_le_deg and
    the vaneless-space exit flow angle -- THIS slice (S5, diffuser vane count) touches
    only the diffusion-factor/loading-loss cascade, so incidence was BIT-IDENTICAL
    across S5 (-3.854 deg).

    R9 (docs/centrifugal/20-review-fixes.md) moved it AGAIN, to -3.844 deg: R9 fixes
    ``ImpellerClearanceJansen``/``ImpellerMixingAungier`` at the IMPELLER exit (feeding
    them the true blade count Z_exit=30 instead of the fractional Z_eff=25.4048), which
    is UPSTREAM of the vaneless space and diffuser this test is about -- less internal
    loss there genuinely shifts the impeller-exit velocity triangle the diffuser inlet
    inherits. That is legitimate re-pinning, not this slice's own drift.

    NASA: vanes staggered 10 deg from tangential, measured inflow 12-14.5 deg from
    tangential at design => ~4.4 deg NEGATIVE incidence ("designed with large negative
    incidence", CR section 9.3.2 / Fig 191). Reproduce with the module docstring's
    recipe above.
    """
    op = _solve_design_point()
    incidence = _vane_le_incidence_deg(op)

    POST_R9_INCIDENCE_DEG = (
        -4.361734120788512
    )  # ⭐ RE-PINNED BY THE SKIN-FRICTION (W-bar) FIX
    # (docs/centrifugal/45-jansen-skinfriction-fix.md). SAME MECHANISM as the E-R11 re-pin
    # below: the correction acts INSIDE THE IMPELLER, upstream of the vaneless space, so the
    # impeller-exit velocity triangle the diffuser inherits shifts and the swirl angle
    # reaching the vane LE shifts with it. But the magnitude is NOT tiny this time:
    # -3.906891 -> -4.361734 deg (-0.455 deg, 65x the E-R11 move).
    #
    # ⚠️ IT MOVED *TOWARD* NASA'S MEASURED -4.4 deg -- it now sits 0.04 deg off a number this
    # test has been 0.5 deg away from for its whole life. THIS WAS NOT SOUGHT AND IS NOT
    # EVIDENCE FOR THE FIX: nothing in the W-bar correction was chosen with the vane-LE
    # incidence in view (the fix was compelled by a normalisation failure in Kovar's Eq. (33),
    # and it makes Eckardt's PR WORSE). It is recorded because it is an INDEPENDENT
    # measurement the model was not fitted to, and it is the one place the fix landed on a
    # NASA number it had no business landing on. Reported, not banked.
    # (E-R11: -3.913939152888119 -> -3.906891. Pre-E-R3: -3.8441599944219433.)
    assert incidence == pytest.approx(POST_R9_INCIDENCE_DEG, abs=1e-6), (
        f"vane-LE incidence moved to {incidence:.3f} deg -- this slice must not touch "
        "incidence (it only changes the diffusion-factor/loading-loss vane count), and "
        "NASA's -4.4 deg is an independently constrained measurement, not a free knob"
    )
