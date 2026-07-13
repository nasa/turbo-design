"""Slice S4 (docs/centrifugal/17-tdd-plan.md) -- exit metal blockage + the mixing-A2 fix.

TWO CHANGES THAT MUST SHIP TOGETHER, ONE COMMIT. Shipping the first without the second
creates a silent decoupling -- docs/centrifugal/12-model-as-implemented.md S7.4's exact
prediction.

CHANGE 1 -- GEOMETRIC TE METAL BLOCKAGE. Both impeller blade rows (15 main + 15
splitter) reach the trailing edge, and both have a nonzero TANGENTIAL thickness there --
solid metal occupying the exit passage, exactly like the vaned diffuser's already-
accepted TE-thickness blockage (``tests/fixtures/hecc_stage.py``, ``blockage=0.031``).
Derived from NASA Appendix C by the SAME method already accepted for the LE thickness
t1 (``extract_hecc_le_thickness.py``):

    B_geom = (Z_main*t_main_TE + Z_splitter*t_splitter_TE) / (2*pi*r2)

``extract_hecc_te_blockage.py`` derives this at 98/99/100% chord (the
station-choice lever, exactly like t1's) and ADOPTS 100% (the true exit plane -- there
is no metal downstream of r2). B_geom ~ 0.0158, well inside the pre-registered
0.010-0.025 band.

CHANGE 2 -- THE MANDATORY COMPANION. ``ImpellerMixingAungier.delta_h`` computed
``A2 = 2*pi*r2*b2`` -- the UNBLOCKED geometric area -- hardcoded, while
``ImpellerLossState`` did not carry the impeller's configured blockage at all. The
moment ``Impeller.blockage != 0``, the velocity triangle (``Stage.solve``'s
``te_geom.with_blockage(...)``) and the mixing loss's own area silently disagree about
the SAME station. Fixed by adding ``ImpellerLossState.blockage2``, populated from
``Impeller.blockage`` at BOTH ``Stage.solve`` construction sites, and using
``A2 = 2*pi*r2*b2*(1-blockage2)`` inside the mixing loss.

WHY AERODYNAMIC BLOCKAGE REMAINS FORBIDDEN (data/coefficients.md; NOT tested here,
recorded for context): Oh, Yoon & Chung (1997) -- the loss set this module implements
-- carries NO blockage anywhere; exit non-uniformity is carried entirely by the mixing
loss, so an aerodynamic blockage on top would double-count it. Eckardt TM-75232 shows
the wake is momentum-RICH (not momentum-neutral), so a one-zone aerodynamic blockage
over-corrects work by construction on a backswept wheel. This slice adds METAL ONLY.

REFUTED IF: a zero-blockage run is not bit-identical to today (blockage2 default must
be 0.0 and must reproduce the old hardcoded ``2*pi*r2*b2`` exactly), or the mixing loss
does not respond to ``blockage2`` (the wiring is wrong), or B lands outside
0.010-0.025 from the coordinates (the extraction or station choice is wrong).
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests" / "fixtures"))

from turbodesign.centrifugal import Air, InletState  # noqa: E402
from turbodesign.centrifugal.losses import ImpellerLossState, ImpellerMixingAungier  # noqa: E402

# ---------------------------------------------------------------------------------
# A real, self-consistent HECC design-point impeller-exit operating state -- the SAME
# fixture as tests/unit/test_diffusion_factor_provenance.py (captured from a converged
# Stage.solve on the current working tree; see that file's docstring for the
# re-derivation recipe). Reused here, not re-solved, so these are unit tests of the
# mixing-loss closure in isolation.
# ---------------------------------------------------------------------------------
R1H = 0.037959334799999994
R1S = 0.10705720191138392
R2 = 0.21581
B2 = 0.015539618656433121
U1 = 183.26606553953252
VM1 = 139.06718837987157
T1 = 278.5282672219489
RHO1 = 1.12542622137953
U2 = 492.4220384078722
CM2 = 90.65796412570579
VT2 = 380.54268794749714
W2 = 143.99949832847747
T2 = 398.47020317955753
RHO2 = 2.5465574233145536
Z_EFF = 25.404764740917212
BACKSWEEP_DEG = -36.0
TIP_CLEARANCE = 0.000305
MDOT = 4.9269
OMEGA = 21789.0 * math.pi / 30.0
L_BLADE = 0.237875  # Impeller.l_blade_m, NASA-measured 3-D camber arc (slice S1)


def _hecc_state(**overrides) -> ImpellerLossState:
    base = dict(
        fluid=Air(),
        mdot=MDOT,
        omega=OMEGA,
        Z=Z_EFF,
        backsweep_deg=BACKSWEEP_DEG,
        r1h=R1H,
        r1s=R1S,
        U1=U1,
        Vm1=VM1,
        T1=T1,
        rho1=RHO1,
        r2=R2,
        b2=B2,
        tip_clearance=TIP_CLEARANCE,
        U2=U2,
        Cm2=CM2,
        Vt2=VT2,
        W2=W2,
        T2=T2,
        rho2=RHO2,
        L_blade=L_BLADE,
    )
    base.update(overrides)
    return ImpellerLossState(**base)


# =====================================================================================
# CHANGE 1 -- the derivation itself
# =====================================================================================


def test_te_metal_blockage_from_appendix_c():
    """B_geom, derived from NASA Appendix C by the t1 method, must land at
    0.016 +- 0.002, and an independently-coded cross-check (a DIFFERENT numerical
    method -- local quadratic interpolation instead of the adopted linear
    interpolation, both targeting the SAME common-chord station) must agree to within
    ~1% -- the same bar ``extract_hecc_le_thickness.py``'s t1 met.

    Imports the derivation, does not duplicate it (docs/centrifugal/17-tdd-plan.md's
    own rule: never copy configuration/derivations into a test).
    """
    from extract_hecc_te_blockage import (
        IN_TO_M,
        N_MAIN,
        N_SPLITTER,
        R2_M,
        _independent_quadratic_check,
        _te_thickness_per_section,
    )

    main_rows = _te_thickness_per_section("IMPELLER MAIN BLADE")
    split_rows = _te_thickness_per_section("IMPELLER SPLITTER BLADE")

    t_main = sum(r[1.00] for r in main_rows) / len(main_rows)
    t_split = sum(r[1.00] for r in split_rows) / len(split_rows)
    r2_in = R2_M / IN_TO_M
    B_adopted = (N_MAIN * t_main + N_SPLITTER * t_split) / (2.0 * math.pi * r2_in)

    assert B_adopted == pytest.approx(0.016, abs=0.002), (
        f"TE metal blockage {B_adopted:.5f} must land at 0.016 +- 0.002 -- if it does "
        "not, the extraction or the 100%-chord station choice is wrong, and the "
        "number must be reported as found, not adjusted to reach 0.016."
    )

    t_main_indep = _independent_quadratic_check(main_rows, 1.00)
    t_split_indep = _independent_quadratic_check(split_rows, 1.00)
    B_indep = (N_MAIN * t_main_indep + N_SPLITTER * t_split_indep) / (
        2.0 * math.pi * r2_in
    )
    rel_diff = abs(B_indep - B_adopted) / B_adopted
    assert rel_diff <= 0.01, (
        f"independent (quadratic-interpolation) read {B_indep:.5f} disagrees with the "
        f"adopted (linear-interpolation) read {B_adopted:.5f} by {100 * rel_diff:.2f}% "
        "-- more than the ~1% t1 standard."
    )


def test_te_blockage_station_choice_lever_is_visible():
    """98%/99%/100% chord must NOT agree -- the near-point TE closure (mirroring t1's
    near-point LE) makes 'the TE thickness' a station choice with a large lever, and
    that lever must be visible, not hidden behind the one adopted number."""
    from extract_hecc_te_blockage import (
        IN_TO_M,
        N_MAIN,
        N_SPLITTER,
        R2_M,
        _te_thickness_per_section,
    )

    main_rows = _te_thickness_per_section("IMPELLER MAIN BLADE")
    split_rows = _te_thickness_per_section("IMPELLER SPLITTER BLADE")
    r2_in = R2_M / IN_TO_M

    def B_at(frac: float) -> float:
        t_main = sum(r[frac] for r in main_rows) / len(main_rows)
        t_split = sum(r[frac] for r in split_rows) / len(split_rows)
        return (N_MAIN * t_main + N_SPLITTER * t_split) / (2.0 * math.pi * r2_in)

    B_100, B_99, B_98 = B_at(1.00), B_at(0.99), B_at(0.98)
    assert B_100 < B_99 < B_98, (
        "thickness must grow monotonically moving upstream from the TE"
    )
    assert B_98 > 3.0 * B_100, (
        "the 98%-chord station must be at least 3x the adopted 100%-chord blockage -- "
        "the lever this slice's docstring reports (measured: ~5-6x)"
    )


# =====================================================================================
# CHANGE 2 -- the mandatory companion: the mixing loss must see the SAME blockage
# =====================================================================================


def test_mixing_loss_tracks_the_state_blockage():
    """A2 must scale with (1 - blockage2). Before this slice's fix, A2 was hardcoded
    to 2*pi*r2*b2 and could not respond to ANY value of blockage2 -- this test MUST
    FAIL against the pre-fix code (either a TypeError, because ImpellerLossState had
    no ``blockage2`` field at all, or -- had a field been bolted on without touching
    ``delta_h`` -- silently return the SAME loss regardless of its value)."""
    model = ImpellerMixingAungier()

    dh_unblocked = model.delta_h(_hecc_state(blockage2=0.0))
    dh_blocked = model.delta_h(_hecc_state(blockage2=0.05))

    assert dh_blocked != pytest.approx(dh_unblocked, rel=1e-9), (
        "ImpellerMixingAungier.delta_h did not respond to blockage2 -- A2 is still "
        "hardcoded to the unblocked geometric area."
    )

    # Blocking the exit area RAISES the true meridional velocity the loss's W_out term
    # sees for the SAME Cm2 (mass continuity through a smaller area was already solved
    # for elsewhere; here we are only checking A2's own sensitivity within delta_h),
    # which must move the mixing loss -- direction not asserted here (it depends on
    # where W_out sits relative to W_sep), only that it MOVES.
    dh_more_blocked = model.delta_h(_hecc_state(blockage2=0.10))
    assert dh_more_blocked != pytest.approx(dh_blocked, rel=1e-9)


def _mixing_reference(state, blockage2: float) -> tuple[float, float, float]:
    """Hand-coded Aungier (1995) Eq. (26) -- the reference the model must reproduce.

    Returns (W_sep, W_out, dh). ``Wu2 = U2 - Vt2`` is the EXIT relative tangential, per
    E-R11/C1: this replaced the old ``Wu1xi = U1`` (Kovar et al. 2021 Eq. (45)'s addition,
    which the primary does not print -- see ImpellerMixingAungier's class docstring).
    """
    from turbodesign.centrifugal.losses import _camberline_length

    r2, b2 = state.r2, state.b2
    Z = state.Z if state.Z_exit is None else state.Z_exit
    D2 = 2.0 * r2
    Cu2, W2_ = state.Vt2, state.W2
    W1xi = math.hypot(state.Vm1, state.U1)
    L_tilde = _camberline_length(state)

    dW = 2.0 * math.pi * D2 * Cu2 / (Z * L_tilde)
    Deq = (W1xi + W2_ + dW) / 2.0 / W2_
    W_sep = W2_ if Deq <= 2.0 else W2_ * Deq / 2.0

    A2 = 2.0 * math.pi * r2 * b2 * (1.0 - blockage2)
    W_out = math.hypot(state.Cm2 * A2 / (math.pi * D2 * b2), state.U2 - Cu2)
    return W_sep, W_out, max(0.5 * (W_sep - W_out) ** 2, 0.0)


def test_zero_blockage_reduces_to_the_unblocked_aungier_formula():
    """REWRITTEN, E-R11/C1 (docs/centrifugal/42-r11-result.md).

    The OLD premise -- "blockage2 = 0.0 must reproduce the pre-slice hardcoded-A2 formula
    bit-for-bit, with W_out = hypot(Cm2*A2/(pi*D2*b2), U1)" -- IS GONE. That formula used
    the INLET blade speed U1 as the mixed-out relative tangential; Aungier (1995) Eq. (26),
    now HELD, prints a bare ``W_u`` inside an all-station-2 equation, so it is the EXIT
    relative tangential, Wu2 = U2 - Vt2. Asserting bit-identity against a formula the
    primary never printed would be pinning the defect.

    What slice S4 actually promised, and what is still true, is that ``blockage2`` is a
    pure ADDITION: at the field default (0.0) the A2 factor is inert, and delta_h reduces
    EXACTLY to Aungier's unblocked form. That is the invariant tested here.

    And the invariant now has TEETH it did not have before. At zero blockage,
    A2/(pi*D2*b2) = 1, so W_out = hypot(Cm2, U2-Vt2) = W2 IDENTICALLY -- the same W2 that
    is W_sep whenever D_eq <= 2. The unblocked mixing loss is therefore EXACTLY ZERO for
    any attached impeller: the term survives only through the D_eq > 2 stall branch. This
    is the C1 annihilation, asserted at its root rather than merely observed downstream.
    """
    model = ImpellerMixingAungier()

    state_default = (
        _hecc_state()
    )  # blockage2 not passed -- exercises the dataclass default
    state_explicit_zero = _hecc_state(blockage2=0.0)

    assert state_default.blockage2 == 0.0, "the field default must be exactly 0.0"
    dh_default = model.delta_h(state_default)
    assert dh_default == model.delta_h(state_explicit_zero), (
        "the default must be bit-for-bit the explicit zero"
    )

    W_sep, W_out, dh_ref = _mixing_reference(state_default, 0.0)

    # THE C1 IDENTITY: at zero blockage the mixed-out relative velocity IS W2.
    assert W_out == pytest.approx(state_default.W2, rel=1e-12)
    assert W_sep == pytest.approx(state_default.W2, rel=1e-12), (
        "this fixture must sit on the ATTACHED (D_eq <= 2) branch, or the identity above "
        "is not what is being tested"
    )
    assert dh_default == pytest.approx(dh_ref, rel=1e-12, abs=1e-9), (
        "blockage2=0.0 must reproduce Aungier Eq. (26) with A2 = 2*pi*r2*b2, bit-for-bit"
    )
    assert dh_default == 0.0, (
        "the unblocked mixing loss is EXACTLY zero on the attached branch -- E-R11/C1"
    )

    # ...and the reference is NOT vacuously comparing 0.0 to 0.0: the SAME hand-coded
    # Eq. (26), at a nonzero blockage, must reproduce the model on a nonzero number.
    state_blocked = _hecc_state(blockage2=0.05)
    _, _, dh_ref_blocked = _mixing_reference(state_blocked, 0.05)
    assert dh_ref_blocked > 0.0
    assert model.delta_h(state_blocked) == pytest.approx(dh_ref_blocked, rel=1e-12)


def test_both_loss_state_construction_sites_carry_blockage():
    """Stage.solve constructs ImpellerLossState at TWO sites (the internal-loss trial
    state inside the coupled continuity solve, and the parasitic-block state after
    convergence). Guard directly against the exact bug class this project keeps
    finding (Z_le, then Z_exit/has_splitters, now blockage2): populate ONE site and
    silently leave the other at the field default.

    Patches ``ImpellerLossState.__init__`` to capture EVERY state constructed during
    one ``Stage.solve`` call (the same technique
    ``tests/unit/test_diffusion_factor_provenance.py``'s docstring documents for
    re-deriving its own fixture) and asserts every single one carries the nonzero
    configured blockage -- not just the last (converged) one, and not just one of the
    two call sites.
    """
    import warnings

    import turbodesign.centrifugal.losses as losses_mod
    import turbodesign.centrifugal.solver as solver_mod
    from hecc_stage import (
        BACKSWEEP,
        IMPELLER,
        MDOT_DESIGN,
        P01,
        RPM,
        T01,
        diffusion_system,
    )
    from turbodesign.centrifugal import (
        Impeller,
        MeridionalPath,
        OhLossSet,
        Stage,
        WiesnerSlip,
    )

    captured: list[ImpellerLossState] = []
    original_init = losses_mod.ImpellerLossState.__init__

    def _capturing_init(self, *args, **kwargs):
        original_init(self, *args, **kwargs)
        captured.append(self)

    losses_mod.ImpellerLossState.__init__ = _capturing_init
    solver_mod.ImpellerLossState = losses_mod.ImpellerLossState
    try:
        impeller_kwargs = dict(IMPELLER)
        impeller_kwargs["blockage"] = 0.0158  # nonzero -- the derived TE metal value
        path = MeridionalPath.from_csv(
            REPO / "data" / "hecc" / "flowpath_hub.csv",
            REPO / "data" / "hecc" / "flowpath_shroud.csv",
        )
        stage = Stage(
            path=path,
            impeller=Impeller(backsweep_deg=BACKSWEEP, **impeller_kwargs),
            slip=WiesnerSlip(),
            losses=OhLossSet(),
            fluid=Air(),
            components=diffusion_system(),
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            stage.solve(mdot=MDOT_DESIGN, rpm=RPM, inlet=InletState(P0=P01, T0=T01))
    finally:
        losses_mod.ImpellerLossState.__init__ = original_init
        solver_mod.ImpellerLossState = losses_mod.ImpellerLossState

    assert len(captured) >= 2, (
        "expected at least one internal-loss trial state and one parasitic state"
    )
    for i, state in enumerate(captured):
        assert state.blockage2 == pytest.approx(0.0158), (
            f"ImpellerLossState #{i} (of {len(captured)}) was constructed with "
            f"blockage2={state.blockage2!r}, not the configured Impeller.blockage "
            "(0.0158) -- one of the TWO construction sites in Stage.solve is not "
            "populating blockage2."
        )
