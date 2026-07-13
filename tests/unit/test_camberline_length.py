"""Slice S1 (docs/centrifugal/17-tdd-plan.md) -- the camberline-length bug.

THE BUG. ``_camberline_length`` (and, before this slice, a deliberate inline duplicate
inside ``ImpellerSkinFrictionJansen``) computed its dominant term as
``phi = Vm1/U1`` = **0.759** for HECC -- the LE flow coefficient. The source (Gambini &
Vellini, *Turbomachinery: Fundamentals, Selection and Preliminary Design*, Springer
(2021), Ch. 6 -- paywalled, equation number UNVERIFIED) uses the GLOBAL flow
coefficient instead, ``phi_t = ((r1s/r2)**2 - (r1h/r2)**2) * (Vm1/U2)`` = **0.0608** for
HECC -- **12.5x smaller**. The shipped code's ``L_z = 0.580 m`` exceeds the impeller's
own diameter (``D2 = 0.432 m``): the model believed the impeller was longer than it is
wide.

THE DECISION (data/coefficients.md, docs/centrifugal/16-code-update-plan.md S1):
replace the closure with MEASURED geometry (NASA Appendix C gives it), falling back to
the CORRECTED closure only when no measurement is supplied. Two length quantities,
one bug: the through-blade length (``Impeller.l_blade_m``, feeds skin friction +
mixing) and the meridional length (``Impeller.l_main_m``, feeds leakage) are DIFFERENT
quantities that a single unsourced closure used to serve simultaneously.

REFUTED IF (docs/centrifugal/17-tdd-plan.md slice S1 / 15-experiments-prereg.md E1):
stage PR does not rise by >= 1.0 point when L~ is set to the measured blade length, or
Delta_h_sf does not fall by >= 900 J/kg. Either result would mean the loss is not
scaling as ``Delta_h_sf = 2*Cf*(L~/d_h)*W_bar**2`` requires, and the bug is not where
this slice's diagnosis says it is.
"""

from __future__ import annotations

import ast
import dataclasses
import inspect
import math
import sys
import warnings
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests" / "fixtures"))

from turbodesign.centrifugal import Air, InletState  # noqa: E402
from turbodesign.centrifugal.losses import (  # noqa: E402
    ImpellerLeakageAungier,
    ImpellerLossState,
    _axial_length_Lz,
    _camberline_length,
    _global_flow_coefficient,
)

# ---------------------------------------------------------------------------------
# A real, self-consistent HECC design-point operating state: LE/TE velocity triangle
# and geometry captured from a converged Stage.solve (tests/fixtures/hecc_stage.py's own
# machine, backsweep -36.0 deg NASA Fig 18, WiesnerSlip, OhLossSet, design mdot/rpm).
# Hardcoded here (rather than re-solving in every test) so these tests exercise the
# CLOSURE in isolation, independent of whatever the full stage solve's loss
# configuration happens to be at the time -- the numbers themselves are real, not
# invented; re-derive with:
#   uv run python -c "import sys; sys.path.insert(0,'tests/fixtures'); from hecc_stage \
#       import build, BACKSWEEP, MDOT_DESIGN, RPM, P01, T01; from \
#       turbodesign.centrifugal import InletState; import warnings; \
#       warnings.simplefilter('ignore'); op = build(BACKSWEEP).solve(mdot=MDOT_DESIGN, \
#       rpm=RPM, inlet=InletState(P0=P01, T0=T01)); print(op.stations.impeller_le.U, ...)"
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
CM2 = 91.81813852158302
VT2 = 379.6997719089894
W2 = 145.38528098199544
T2 = 398.27068854902996
RHO2 = 2.5465574233145536
Z_EFF = 25.404764740917212
BACKSWEEP_DEG = -36.0
TIP_CLEARANCE = 0.000305
MDOT = 4.9269
OMEGA = 21789.0 * math.pi / 30.0

D2 = 2.0 * R2  # 0.43162 m -- the impeller's own exit diameter


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
    )
    base.update(overrides)
    return ImpellerLossState(**base)


# ------------------------------------------------------------- L_z: the geometric absurdity


def test_Lz_cannot_exceed_the_impeller_diameter():
    """The falsifiable claim at the heart of this slice: the impeller cannot be longer
    (axially) than it is wide. Shipped: L_z = 0.580 m > D2 = 0.432 m. FAILS TODAY.
    """
    state = _hecc_state()
    L_z = _axial_length_Lz(state)
    assert L_z <= D2, (
        f"L_z = {L_z:.4f} m exceeds the impeller's own exit diameter D2 = {D2:.4f} m -- "
        "the model believes this impeller is longer (axially) than it is wide"
    )


def test_camberline_uses_the_global_flow_coefficient():
    """phi_t (the GLOBAL flow coefficient) must replace phi = Vm1/U1 (the LE one).

    phi_t = ((r1s/r2)**2 - (r1h/r2)**2) * (Vm1/U2) = 0.0608 for HECC (vs the shipped
    Vm1/U1 = 0.759 -- 12.5x too large). The corrected closure gives L_z ~ 0.104 m.
    """
    state = _hecc_state()

    phi_t = _global_flow_coefficient(state)
    assert phi_t == pytest.approx(0.0608, rel=0.03)

    phi_shipped = (
        state.Vm1 / state.U1
    )  # the BUG's flow coefficient, for comparison only
    assert phi_shipped == pytest.approx(0.759, rel=0.02)
    assert phi_shipped / phi_t == pytest.approx(12.5, rel=0.05), (
        "the whole point of this slice: the shipped flow coefficient is ~12.5x the "
        "correct (global) one"
    )

    L_z = _axial_length_Lz(state)
    assert L_z == pytest.approx(0.104, rel=0.05)


def test_the_old_flow_coefficient_form_is_unreachable():
    """No code path may return the shipped phi = Vm1/U1 form any more.

    ``_axial_length_Lz`` (the only function that computes L_z) must equal what you get
    by feeding it the GLOBAL flow coefficient, never the LE one -- checked by
    reconstructing both candidate L_z values from the closure's own published formula
    and asserting the function's return matches only the corrected candidate.
    """
    state = _hecc_state()
    r1h, r2 = state.r1h, state.r2

    def Lz_from_phi(phi: float) -> float:
        return 2.0 * r2 * (0.014 + 0.023 * r2 / r1h + 1.58 * phi)

    L_z_shipped_bug = Lz_from_phi(state.Vm1 / state.U1)
    L_z_corrected = Lz_from_phi(_global_flow_coefficient(state))

    L_z_actual = _axial_length_Lz(state)
    assert L_z_actual == pytest.approx(L_z_corrected, rel=1e-9)
    assert L_z_actual != pytest.approx(L_z_shipped_bug, rel=1e-3)


# ------------------------------------------------------------- l_blade_m: measured, not proxied


def test_blade_length_is_the_measured_arc_not_the_LE_angle_projection():
    """Impeller.l_blade_m (HECC) = 0.237875 m -- the MEASURED 3-D camberline arc.

    ds = sqrt(dm^2 + (r*dtheta)^2), integrated along each of NASA Appendix C's 11
    main-blade sections (Tables C.3-C.13) and span-averaged. It needs NO blade angle at
    all: it is the geometry.

    THIS TEST EXISTS TO REJECT A SPECIFIC WRONG ANSWER. The tempting derivation is
    l_main_m / mean(cos(beta1b)) = 0.292395 m, which divides the WHOLE blade's meridional
    length by the cosine of the INDUCER LEADING-EDGE angle (45.46 deg) -- a one-station
    quantity applied to the whole blade. HECC's blade is S-shaped (NASA Fig 18: a
    mid-chord beta minimum near 14-20 deg, where cos(beta) ~ 0.95, not 0.70), so it
    inflates the length by 23% and the skin-friction loss with it.

    That is this project's signature failure mode -- A QUANTITY VALID AT ONE STATION USED
    AT ANOTHER (docs/PHYSICS-RULES.md) -- and it was very nearly adopted here, because it
    happened to agree with a rough L_m/<cos beta> figure quoted in the slice-S1 plan.
    AGREEING WITH YOUR OWN ESTIMATE IS NOT VALIDATION.
    """
    from extract_hecc_blade_angles import _lines, blade_length_estimates, find_sections

    lines = _lines()
    est = blade_length_estimates(lines, find_sections(lines, "IMPELLER MAIN BLADE"))

    measured_arc = est["l_blade_m"]
    le_projection = est["l_blade_m_le_angle_projection"]

    assert measured_arc == pytest.approx(0.237875, rel=0.01)

    # The rejected proxy must stay visibly different -- if these ever converge, the arc
    # integration has silently degenerated into the projection.
    assert le_projection / measured_arc == pytest.approx(1.23, rel=0.03), (
        f"the LE-angle projection ({le_projection:.6f} m) should exceed the measured arc "
        f"({measured_arc:.6f} m) by ~23%. If it does not, one of the two is no longer "
        "computing what its name says."
    )

    from hecc_stage import IMPELLER  # tests/fixtures/hecc_stage.py's wired-in value

    assert "l_blade_m" in IMPELLER, (
        "tests/fixtures/hecc_stage.py's IMPELLER dict must supply l_blade_m -- without it "
        "skin friction and mixing fall back to the UNVERIFIED closure, not NASA's "
        "measured geometry"
    )
    assert IMPELLER["l_blade_m"] == pytest.approx(measured_arc, rel=0.01)
    assert IMPELLER["l_blade_m"] != pytest.approx(le_projection, rel=0.05), (
        "l_blade_m has regressed to the LE-angle projection"
    )


# ------------------------------------------------------------- the pre-registered result


def test_skin_friction_at_hecc_design_point():
    """Delta_h_sf: 2580 J/kg (shipped) -> ~1160 J/kg, at the full HECC design point.

    2580 is what the buggy closure gave (L~ = 0.522 m, an impeller longer than it is
    wide). 1160 is what NASA's measured camberline arc (0.237875 m) gives.

    REFUTED IF this does not land inside 1160 +- 150 J/kg: the loss would not be scaling
    with L~ the way Jansen's own form requires, and the bug is not where this slice's
    diagnosis says it is.
    """
    from hecc_stage import BACKSWEEP, MDOT_DESIGN, P01, RPM, T01, build

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        op = build(BACKSWEEP).solve(
            mdot=MDOT_DESIGN, rpm=RPM, inlet=InletState(P0=P01, T0=T01)
        )

    dh_sf = op.losses.internal["ImpellerSkinFrictionJansen"]
    # E-R3/D1 (docs/centrifugal/26-prereg-r3.md): d_h previously carried a spurious factor
    # 2.0 that Kovar et al. 2021 Eq. (29) does not have. Removing it RAISES d_h and so
    # LOWERS this loss: 1179.3 -> 913.6 J/kg. The band below was built around the buggy
    # value (~1160); it is re-centred on the sourced one. The 2580 in the message is the
    # pre-S1 shipped value, kept as the historical anchor.
    #
    # E-R11/C3 (docs/centrifugal/41-prereg-adoption.md, 42-r11-result.md): what
    # ``op.losses.internal`` now REPORTS is the f_c-SCALED loss -- Aungier's head-loss
    # correction (Eqs. 31-33) multiplies the INTERNAL sum, f_c = 1.2883 on HECC. So this
    # pin moves 915 -> 1175 J/kg. NOTE it is NOT a clean 913.6*1.2883 = 1177: the raw
    # (unscaled) loss ALSO moved, 913.6 -> 912.1 J/kg, because C3 feeds back through the
    # solve (more internal loss -> lower rho2 -> higher Cm2 -> a different W_bar). The
    # camberline arc L~ = 0.237875 m -- the quantity THIS file is about -- is untouched.
    #
    # ⭐ THE W-bar FIX (docs/centrifugal/45-jansen-skinfriction-fix.md) RE-CENTRES THIS PIN
    # 1175 -> 4427 J/kg, and the ratio (3.65x) is the whole story. W_bar was Kovar et al. 2021
    # Eq. (33), (2*W2 + W1s - W1h)/4 -- weights summing to 2 over a denominator of 4, so it
    # returned W/2 for a constant field and 103.6 m/s at HECC, BELOW the minimum (148.6 m/s) of
    # the velocities it averaged. It is now Oh, Yoon & Chung (1997) Table 6's printed five-term
    # form, W_bar = (V_1t + V_2 + W_1t + 2*W_1h + 3*W_2)/8, which normalises: 197.9 m/s.
    # L~ is UNTOUCHED, which is why this file still pins the number -- the band moved because
    # the VELOCITY did, not because the length did.
    assert dh_sf == pytest.approx(4427.0, abs=150.0), (
        f"Delta_h_sf = {dh_sf:.1f} J/kg -- expected ~4427 +- 150 (= f_c * ~3436; pre-W-bar-fix: "
        "1175; pre-E-R11: 915; pre-E-R3: 1179; pre-S1: 2580). If this misses, the loss is not "
        "scaling as Delta_h_sf = f_c * 2*Cf*(L~/d_h)*W_bar**2 requires, and the bug is not "
        "where this slice's diagnosis says it is."
    )


# ------------------------------------------------------------- leakage: a DIFFERENT quantity


def test_leakage_uses_the_meridional_length():
    """Leakage consumes L_meridional (Impeller.l_main_m), NOT L_blade -- two DIFFERENT
    quantities (docs/centrifugal/12-model-as-implemented.md S7.2; data/coefficients.md
    closes the old "Lm := L_tilde" proxy row).

    Two things are checked, for a specific reason stated here rather than hidden: this
    module's own Aungier-leakage algebra makes ``Delta_h_lk`` EXACTLY invariant to the
    numeric value of Lm (dP_cl ~ 1/Lm, Ucl ~ Lm**-0.5, mdot_cl ~ Lm**+0.5, and their
    product in Delta_h_lk cancels Lm identically -- verified both analytically and
    numerically: Delta_h_lk is bit-identical from Lm = 0.05 m to Lm = 5.0 m). So the
    WIRING (which field leakage reads) cannot be proven by its effect on the output --
    only by source inspection, which is what the second assertion does.
    """
    base = _hecc_state()

    # 1. Leakage must be insensitive to L_blade -- proving it does not read that field.
    state_small = dataclasses.replace(base, L_blade=0.05, L_meridional=0.209875)
    state_large = dataclasses.replace(base, L_blade=9.0, L_meridional=0.209875)
    dh_small = ImpellerLeakageAungier().delta_h(state_small)
    dh_large = ImpellerLeakageAungier().delta_h(state_large)
    assert dh_small == pytest.approx(dh_large, rel=1e-9), (
        "ImpellerLeakageAungier's output must not depend on L_blade at all"
    )

    # 2. The implementation must actually read state.L_meridional (not merely tolerate
    #    it being absent) -- checked by source inspection, since Delta_h_lk's own
    #    algebra is provably invariant to Lm's numeric value (see docstring above), so
    #    no behavioural test on the RETURNED NUMBER can distinguish the source field.
    src = inspect.getsource(ImpellerLeakageAungier.delta_h)
    assert "state.L_meridional" in src, (
        "ImpellerLeakageAungier.delta_h must read state.L_meridional for its Lm -- the "
        "MERIDIONAL length, a different quantity from the through-blade L_blade that "
        "skin friction and mixing use"
    )

    # 3. And the SHARED closure skin friction/mixing call must not itself read
    #    L_meridional -- the two consumers stay on separate fields. Walk the AST (not
    #    a raw substring check) so the docstring's prose -- which necessarily discusses
    #    both field names -- cannot trip the assertion.
    tree = ast.parse(inspect.getsource(_camberline_length))
    accessed = {
        node.attr
        for node in ast.walk(tree)
        if isinstance(node, ast.Attribute) and node.attr == "L_meridional"
    }
    assert not accessed, "_camberline_length must not read state.L_meridional"


# ------------------------------------------------------------- R1: dimensional defect


@pytest.mark.xfail(
    strict=True,
    reason=(
        "R1 (docs/centrifugal/20-review-fixes.md): ImpellerLeakageAungier.delta_h is "
        "dimensionally m/s, not J/kg (docs/PHYSICS-RULES.md rule 1) -- confirmed today "
        "at ~1.914812 J/kg, i.e. really 1.91 m/s. DO NOT make this pass by guessing the "
        "missing factor (e.g. multiplying by Ucl again) -- that would be exactly the "
        "magic-number failure this project exists to prevent. This xfail stays red "
        "until the formula is repaired against Aungier (2000)'s own primary source."
    ),
)
def test_leakage_delta_h_scales_as_a_specific_energy_not_a_velocity():
    """BROKEN, R1: Delta_h_lk = mdot_cl*Ucl/(2*mdot) is dimensionally
    [kg/s * m/s] / [kg/s] = m/s, not J/kg = m^2/s^2.

    Dimensional check via a scaling argument (no guessed factor involved): scale
    state.mdot and state.Vt2 (Cu2) UNIFORMLY by k, holding densities/geometry fixed.

        dP_cl ~ mdot * Cu2                    -> scales as k**2
        Ucl   ~ sqrt(dP_cl)                   -> scales as k**1
        mdot_cl ~ Ucl                          -> scales as k**1

    A genuine specific-energy (J/kg) quantity built the same way from Ucl must then
    scale as Ucl**2, i.e. k**2 (mdot_cl*Ucl**2/mdot ~ k*k**2/k = k**2). The shipped
    form instead returns mdot_cl*Ucl/mdot ~ k*k/k = k**1 -- linear in k, the scaling
    signature of a velocity, not a specific energy. Verified numerically today: k=2
    gives a 2.0x change in delta_h, not the 4.0x a dimensionally sound term requires.

    STRICT XFAIL: this assertion passes only once delta_h is dimensionally repaired.
    Until then it keeps the defect visible rather than letting it silently return a
    plausible-looking 1.9 J/kg.
    """
    k = 2.0
    base = _hecc_state()
    scaled = dataclasses.replace(base, mdot=base.mdot * k, Vt2=base.Vt2 * k)

    dh_base = ImpellerLeakageAungier().delta_h(base)
    dh_scaled = ImpellerLeakageAungier().delta_h(scaled)

    assert dh_scaled == pytest.approx(dh_base * k**2, rel=1e-6), (
        f"a dimensionally-correct J/kg leakage loss must scale as k**2 = {k**2} under "
        f"this uniform velocity scaling; got a factor of {dh_scaled / dh_base:.3f} "
        "instead (k**1 is the shipped, broken behaviour)"
    )
