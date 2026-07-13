"""Tests for the R2 loss-channel candidate sets (``turbodesign.centrifugal.candidates``).

``docs/centrifugal/25-prereg-r2.md`` S0/S4: ``OhLossSet`` is the FROZEN CONTROL and must stay
byte-identical -- this file's first test is a golden-value guard for exactly that, mirroring
``tests/unit/test_ledger_truth.py``'s own HECC design-point pin (imported values match:
psi 0.8192631815566565, PR 4.955847012745787). The remaining tests exercise the two NEW
candidate sets this experiment adds, without touching ``losses.py`` in any way.
"""

from __future__ import annotations

import math
import sys
import warnings
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "my_scripts"))

from hecc_stage import BACKSWEEP, MDOT_DESIGN, P01, RPM, T01, build  # noqa: E402

from turbodesign.centrifugal import (  # noqa: E402
    Air,
    Impeller,
    InletState,
    MeridionalPath,
    OhLossSet,
    Stage,
    WiesnerSlip,
)
from turbodesign.centrifugal.candidates import (  # noqa: E402
    ImpellerMixingBounded,
    ImpellerRecirculationCoppage,
    MixingBoundedLossSet,
    OhCoppageLossSet,
)
from turbodesign.centrifugal.losses import ImpellerLossState  # noqa: E402

# ------------------------------------------------------------- the frozen control, pinned

# Bit-identical to tests/unit/test_ledger_truth.py::PSI_AFTER_S5 / STAGE_PR_AFTER_S5 -- a
# SEPARATE pin, not a duplicate for its own sake: this file's whole purpose is to guard that
# nothing this experiment adds has, even accidentally, disturbed OhLossSet.
#
# E-R11 (docs/centrifugal/41-prereg-adoption.md; results 42-r11-result.md) re-pins BOTH.
# "Frozen" here means "not contaminated by THIS experiment's candidate sets" -- it does NOT
# mean "immune to an adopted, pre-registered correction to the shared loss library". The
# three adopted corrections (C1 mixing station -> exit; C2 the missing *U2 in leakage;
# C3 f_c, Aungier's head-loss correction, default-ON) all live INSIDE OhLossSet, so the
# control legitimately moves with them. Kept bit-identical to
# tests/unit/test_ledger_truth.py::PSI_AFTER_S5 / STAGE_PR_AFTER_S5, as before.
#
# ⭐ RE-PINNED AGAIN by the skin-friction (W-bar) fix -- Kovar et al. 2021 Eq. (33)'s velocity
# average did not normalise (weights summed to 2 over a denominator of 4) and was replaced by
# Oh, Yoon & Chung (1997) Table 6's five-term form. That correction ALSO lives inside OhLossSet,
# so this control legitimately moves with it, exactly as it did for C1/C2/C3.
# docs/centrifugal/45-jansen-skinfriction-fix.md.
PSI_FROZEN = 0.8127325809383643  # W-bar fix; was 0.82479529307983035 (E-R11)
STAGE_PR_FROZEN = 4.780524485318113  # W-bar fix; was 4.9424308497058096 (+5.50% -> +2.05% vs NASA)


def test_oh_loss_set_is_unchanged():
    """OhLossSet MUST remain byte-identical -- this experiment adds new objects, it does
    not edit the frozen one. If this ever fails, the control has been contaminated.
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # ImpellerRecirculationOh's own alpha2 tripwire
        op = build(BACKSWEEP).solve(mdot=MDOT_DESIGN, rpm=RPM, inlet=InletState(P0=P01, T0=T01))

    assert op.psi == pytest.approx(PSI_FROZEN, rel=1e-12)
    assert op.stage_PR == pytest.approx(STAGE_PR_FROZEN, rel=1e-12)


# ------------------------------------------------------------- a real operating state,
# for exercising the candidate models directly (not through a full Stage.solve)


def _hecc_te_state(**overrides) -> ImpellerLossState:
    """A converged HECC design-point ImpellerLossState -- captured the same way
    tests/unit/test_diffusion_factor_provenance.py captures its own fixture (see that
    file's module docstring for the exact re-derivation recipe), reused here so this file
    does not need its own copy of the whole Stage-solve machinery just to build one state.
    """
    base = dict(
        fluid=Air(),
        mdot=MDOT_DESIGN,
        omega=RPM * math.pi / 30.0,
        Z=25.404764740917212,  # Impeller.Z_eff, HECC 15+15
        backsweep_deg=BACKSWEEP,
        r1h=0.037959334799999994,
        r1s=0.10705720191138392,
        U1=183.26606553953252,
        Vm1=139.06718837987157,
        T1=278.5282672219489,
        rho1=1.12542622137953,
        r2=0.21581,
        b2=0.015539618656433121,
        tip_clearance=0.000305,
        U2=492.4220384078722,
        Cm2=90.65796412570579,
        Vt2=380.54268794749714,
        W2=143.99949832847747,
        T2=398.47020317955753,
        rho2=2.5465574233145536,
        Z_exit=30.0,
        has_splitters=True,
    )
    base.update(overrides)
    return ImpellerLossState(**base)


def test_recirculation_coppage_matches_the_transcribed_formula():
    """dh_rc = 0.02 * sqrt(tan(alpha2)) * Df^2 * U2^2 (Kovář et al. 2021 Eq. 46 / Galvas Eq. B73),
    checked directly against the closed form -- not merely "runs without error".
    """
    from turbodesign.centrifugal.losses import _diffusion_factor

    state = _hecc_te_state()
    model = ImpellerRecirculationCoppage()
    alpha2 = math.atan2(state.Vt2, state.Cm2)
    Df = _diffusion_factor(state)
    expected = 0.02 * math.sqrt(math.tan(alpha2)) * Df**2 * state.U2**2
    assert model.delta_h(state) == pytest.approx(expected, rel=1e-12)


@pytest.mark.parametrize(
    "cm2, vt2",
    [
        (90.0, -50.0),  # alpha2 < 0 (reversed swirl)
        (90.0, 0.0),  # alpha2 == 0.0 exactly -- the boundary is NOT admitted (0 < alpha2 required)
        (-1e-6, 90.0),  # alpha2 just over 90 deg
        (-10.0, 90.0),  # alpha2 ~= 96.3 deg
        (-50.0, 30.0),  # alpha2 ~= 149 deg
    ],
)
def test_recirculation_coppage_domain_guard_fires_loudly(cm2, vt2):
    """sqrt(tan(alpha2)) is undefined outside (0, 90) deg -- OUR OWN guard (neither Kovář et al.
    nor Galvas states one). Must RAISE, never silently clamp or return a complex-derived
    NaN/garbage value. (Cm2, Vt2) pairs are used directly, rather than round-tripping
    through degrees->atan2, so the boundary cases aren't at the mercy of floating-point
    rounding landing on the wrong side of the guard.
    """
    state = _hecc_te_state(Cm2=cm2, Vt2=vt2)
    with pytest.raises(ValueError, match="alpha2"):
        ImpellerRecirculationCoppage().delta_h(state)


def test_recirculation_coppage_is_parasitic_and_mixing_bounded_is_internal():
    """Rule 6: recirculation adds work with no pressure benefit (PARASITIC); mixing
    destroys pressure with no work term (INTERNAL). Getting this backwards inverts the
    physics -- assert it directly, not just "the set solves".
    """
    assert ImpellerRecirculationCoppage.kind == "parasitic"
    assert ImpellerMixingBounded.kind == "internal"


def test_oh_coppage_loss_set_differs_from_oh_loss_set_in_exactly_one_model():
    """OhCoppageLossSet swaps ONLY the recirculation model -- every other class must be
    the SAME class (by type), not merely "close".
    """
    oh_types = [type(m) for m in OhLossSet().models]
    cop_types = [type(m) for m in OhCoppageLossSet().models]

    assert len(oh_types) == len(cop_types)
    differences = [(a, b) for a, b in zip(oh_types, cop_types) if a is not b]
    assert len(differences) == 1, f"expected exactly one differing model, got {differences}"
    old_cls, new_cls = differences[0]
    assert old_cls.__name__ == "ImpellerRecirculationOh"
    assert new_cls.__name__ == "ImpellerRecirculationCoppage"


def test_mixing_bounded_loss_set_differs_from_oh_loss_set_in_exactly_one_model():
    """MixingBoundedLossSet swaps ONLY the mixing model -- recirculation stays Oh's own
    (frozen, untouched): this set isolates the INTERNAL channel exactly as
    OhCoppageLossSet isolates the PARASITIC channel (docs/centrifugal/25-prereg-r2.md).
    """
    oh_types = [type(m) for m in OhLossSet().models]
    mb_types = [type(m) for m in MixingBoundedLossSet().models]

    assert len(oh_types) == len(mb_types)
    differences = [(a, b) for a, b in zip(oh_types, mb_types) if a is not b]
    assert len(differences) == 1, f"expected exactly one differing model, got {differences}"
    old_cls, new_cls = differences[0]
    assert old_cls.__name__ == "ImpellerMixingAungier"
    assert new_cls.__name__ == "ImpellerMixingBounded"


def test_mixing_bounded_removes_the_1_over_z_collapse():
    """The whole point of the probe: two states identical except for blade count (Z_exit)
    must give the SAME mixing loss under ImpellerMixingBounded (both use z_reference), but
    DIFFERENT losses under the frozen ImpellerMixingAungier (which reads Z_exit directly).
    """
    from turbodesign.centrifugal.losses import ImpellerMixingAungier

    state_20 = _hecc_te_state(Z_exit=20.0, has_splitters=False)
    state_30 = _hecc_te_state(Z_exit=30.0, has_splitters=True)

    bounded_20 = ImpellerMixingBounded().delta_h(state_20)
    bounded_30 = ImpellerMixingBounded().delta_h(state_30)
    assert bounded_20 == pytest.approx(bounded_30, rel=1e-12)

    frozen_20 = ImpellerMixingAungier().delta_h(state_20)
    frozen_30 = ImpellerMixingAungier().delta_h(state_30)
    assert frozen_20 != pytest.approx(frozen_30, rel=1e-6)


def test_mixing_bounded_reuses_the_frozen_coefficients():
    """mixing_k and deq_threshold are UNCHANGED from the frozen ImpellerMixingAungier --
    only Z is fixed. Guards against a silent coefficient drift creeping into the probe.
    """
    from turbodesign.centrifugal.losses import ImpellerMixingAungier

    bounded = ImpellerMixingBounded()
    frozen = ImpellerMixingAungier()
    assert bounded.mixing_k == frozen.mixing_k
    assert bounded.deq_threshold == frozen.deq_threshold


def test_candidate_sets_solve_end_to_end_on_hecc():
    """Sanity: both candidate sets must actually solve a real Stage (not just their bare
    delta_h) without raising, at the HECC design point.
    """
    path = MeridionalPath.from_csv(
        REPO / "data" / "hecc" / "flowpath_hub.csv",
        REPO / "data" / "hecc" / "flowpath_shroud.csv",
    )
    from hecc_stage import IMPELLER, diffusion_system

    for losses in (OhCoppageLossSet(), MixingBoundedLossSet()):
        stage = Stage(
            path=path,
            impeller=Impeller(backsweep_deg=BACKSWEEP, **IMPELLER),
            slip=WiesnerSlip(),
            losses=losses,
            fluid=Air(),
            components=diffusion_system(),
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            op = stage.solve(mdot=MDOT_DESIGN, rpm=RPM, inlet=InletState(P0=P01, T0=T01))
        assert op.stage_PR is not None and op.stage_PR > 1.0
