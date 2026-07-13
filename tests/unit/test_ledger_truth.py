"""The design-point guard: the HECC number, pinned, so no slice moves it by accident.

WHY THIS FILE IMPORTS THE DRIVER INSTEAD OF COPYING IT
-----------------------------------------------------
The first version of this guard COPIED ``tests/fixtures/hecc_stage.py``'s ``IMPELLER`` dict,
with the comment *"kept byte-identical so this test builds the exact same machine the
driver does"*. That is a promise the code cannot keep, and it broke within the hour: slice
S1 corrected ``l_blade_m`` in the driver (0.292394 -> 0.237875 m), the copy kept the old
value, and this test went on **passing at rel=1e-12** while the machine it claimed to be
guarding had changed underneath it.

A regression guard that duplicates the configuration it guards **cannot fail for the
reason it exists**. It certifies a machine that no longer exists. That is precisely this
project's signature defect -- a term that cannot fire cannot be caught being wrong
(docs/PHYSICS-RULES.md) -- reproduced inside the very test written to prevent it.

So: import the driver. One machine, one definition. If the driver's geometry changes, this
test FAILS, loudly, and someone has to decide whether the change was intended and re-pin
the number in the same commit. That failure is the whole point of the file.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import pytest

from turbodesign.centrifugal import InletState

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "tests" / "fixtures"))

from hecc_stage import BACKSWEEP, MDOT_DESIGN, P01, RPM, T01, build  # noqa: E402

# Design-point truth, captured on the working tree after slice S1's camberline-length fix,
# with l_blade_m = 0.237875 m -- NASA's MEASURED 3-D camberline arc (Appendix C), not the
# leading-edge-angle projection (0.292395 m) that S1 first adopted and that this project's
# cross-station rule rejects.
#
# Slice S1 moved stage_PR from 4.9222 (+5.07%) to 4.9930 (+6.58%). THE MODEL GOT WORSE, and
# that is the finding: the buggy closure (L~ = 0.522 m -- an impeller longer than it is
# wide) was inflating Jansen skin friction (2580 -> 1160 J/kg) and MASKING ~1.5 points of
# the true pressure-ratio gap. The missing internal loss is ~6-7 kJ/kg, not ~4.75.
#
# Slice S2 (docs/centrifugal/17-tdd-plan.md; data/coefficients.md S3.1) moves stage_PR AGAIN,
# from 4.9930 (+6.58%) to 5.0040 (+6.82%) -- WORSE AGAIN, and pre-registered as such
# (docs/centrifugal/15-experiments-prereg.md E4). ``_diffusion_factor``'s K_BL/Z were an
# unsourced hybrid: Coppage/Galvas's D_f + the NO-SPLITTER constant 0.75, applied
# unconditionally, + Aungier's FRACTIONAL Z_eff (25.4048) borrowed from a DIFFERENT loss
# framework (wetted area / hydraulic diameter). Galvas himself (NASA TN D-7487 (1973), Eq.
# B59 + FORTRAN ``CONST1=0.75; IF(SPLT.EQ.1) CONST1=0.6``) prescribes K_BL = 0.6 for a
# SPLITTERED impeller and an INTEGER blade count at the exit (Z_exit = n_blades +
# n_splitters = 30 for HECC's 15+15) -- never a fractional count. The sourced
# (Z_exit=30, K_BL=0.6) pair lowers D_f (0.546 -> 0.529), which lowers the internal
# blade-loading loss (3618.4 -> 3397.7 J/kg, raising PR) and lowers the parasitic
# recirculation loss (12411.2 -> 11775.2 J/kg, lowering psi) -- psi and PR move in
# OPPOSITE directions, the signature docs/PHYSICS-RULES.md rule 6 predicts for one D_f
# driving one internal and one parasitic term.
#
# Slice S3 (docs/centrifugal/17-tdd-plan.md; tests/unit/test_inducer_throat.py) moves
# stage_PR AGAIN, from 5.0040 (+6.82%) to 4.9914 (+6.55%) -- BETTER this time, and NOT
# pre-registered: the plan predicted the design point would be UNCHANGED (Delta h_edf =
# 0 +- 50 J/kg, "the throat is matched to 0.3%"). ``ImpellerEntranceDiffusionAungier``
# was un-phantomed (Meroni et al. 2018 Eqs. 3-5; Wth solved from ``A_th = 0.020525 m2``
# instead of the old ``Wth := W1xi`` stand-in) and fires a real 248.3 J/kg at design
# under this codebase's established RMS convention.
#
# R11 (docs/centrifugal/20-review-fixes.md): that null is CONVENTION-SENSITIVE, not
# cleanly REFUTED -- an earlier version of this comment (and of docs/centrifugal/
# 13-literature-review.md Sec 1.1, commit 4fdb518) each independently explained the
# same fact two different ways ("bad brief velocities" vs "a convention mismatch").
# They are ONE fact, said once: data/coefficients.md Sec 3.1 is the canonical, fully
# worked account (raw 9 J/kg -> net ~0 under a mean-radius U1 convention, vs raw 309 ->
# net 248.3 under the RMS convention this module keeps -- the plan's own 200 J/kg bar
# sits INSIDE that band). RMS was kept because it is what every sibling loss inside the
# SAME ``delta_h`` already reads for "W1", not because it crosses the bar. The ROBUST,
# convention-independent part is the off-design growth (665.7 / 1155.2 / 1708.9 J/kg at
# 90/80/70% flow, tests/unit/test_inducer_throat.py) -- that is what makes the term
# load-bearing, not the design-point figure.
#
# Slice S4 (docs/centrifugal/17-tdd-plan.md; tests/unit/test_exit_blockage.py) moves
# stage_PR AGAIN, from 4.9914 (+6.55%) to 4.9798 (+6.30%) -- BETTER again, and this time
# TOWARD the measurement by CONSTRUCTION, not by tuning: two changes shipped in one
# commit, because shipping the first without the second would silently decouple them
# (docs/centrifugal/12-model-as-implemented.md S7.4's exact prediction). (1) A real,
# GEOMETRIC trailing-edge metal blockage (B = 0.015837, both blade rows' own tangential
# thickness at the true exit plane, derived from NASA Appendix C by the same method
# already accepted for the LE thickness t1 -- extract_hecc_te_blockage.py).
# (2) ``ImpellerMixingAungier`` had been computing its own exit area as the UNBLOCKED
# ``2*pi*r2*b2``, hardcoded, while the velocity triangle itself was already using the
# BLOCKED area (``Impeller.blockage`` existed and fed ``te_geom`` but ``ImpellerLossState``
# never carried it) -- a silent decoupling the instant blockage != 0.
# ``ImpellerLossState.blockage2``, populated from ``Impeller.blockage`` at BOTH
# ``Stage.solve`` construction sites, closes it.
#
# The pre-registered prediction (docs/centrifugal/17-tdd-plan.md S4) was psi ~ -0.27 pt,
# stage PR ~ -0.36 pt (+6.55% -> ~+6.2%). What was MEASURED: PR moved -0.25 pt
# (+6.55% -> +6.30%, matching the prediction closely) but psi moved -0.85 pt
# (+1.96% -> +1.11%) -- about 3x the predicted magnitude, though in the SAME direction
# (both toward the measurement). REPORTED, not adjusted to match the prereg.
#
# ---- THE MECHANISM. The FIRST TWO explanations of this number were WRONG. ----
#
# The NUMBER above has never been in dispute. The PHYSICS written down next to it was
# invented, twice:
#
#   (1) Commit 4cb610e's message, data/coefficients.md S3.1a, and THIS comment all said:
#       "the mixing loss's larger A2-driven delta_h feeds the coupled continuity solve's
#       Cm2." WRONG. A commit message cannot be edited -- data/coefficients.md S3.1b and
#       this comment are the CORRECTION OF RECORD.
#   (2) docs/centrifugal/20-review-fixes.md S R2 -- this project's OWN review document --
#       caught (1) and then proposed a SECOND wrong mechanism: "the internal-loss/
#       continuity feedback in te_residual amplifies dCm2 by ~2.5x through every internal
#       loss's state response." ALSO WRONG. Measured amplification: 1.05x. That section is
#       now marked REFUTED in its own document.
#
# MEASURED DECOMPOSITION (each step re-solved at the HECC design point; the mixing loss is
# un-blocked by forcing ImpellerLossState.blockage2 = 0.0 while leaving Impeller.blockage,
# and hence the velocity-triangle area, untouched):
#
#   S3 baseline (blockage = 0 everywhere)                          psi err +1.958%
#   blockage in the VELOCITY-TRIANGLE AREA ALONE (mixing unblocked)      +1.090%  (-0.867 pt)
#   the mixing-A2 fix ON TOP                                             +1.108%  (+0.018 pt)
#
# The A2 term the old text blamed is +0.018 pt -- OPPOSITE in sign, ~2% of the effect.
#
# WHERE THE -0.867 pt ACTUALLY LIVES. psi = work_actual/U2**2 = (work_euler +
# parasitic_total)/U2**2, exactly. Split at the design point:
#
#   Euler channel (U2*Vt2):        -0.275 pt   <- the prereg's own first-order chain
#                                                 (dCm2 = Cm2*B -> dVt2 = dCm2*tan36 -> dpsi).
#                                                 THE PREREG PREDICTED -0.27. IT WAS RIGHT.
#   PARASITIC channel:             -0.592 pt   <- THE ENTIRE OVERSHOOT. ImpellerRecirculationOh's
#                                                 sinh(3.5*alpha2**3): EXPONENTIAL in the exit
#                                                 swirl angle, responding to the shifted Cm2/Vt2.
#
# And there is NO dCm2 amplification: dCm2 = 1.511 m/s vs a naive Cm2*B = 1.436 m/s (1.05x).
# With the internal-loss/continuity coupling off entirely (losses=None -> _internal_losses_at
# returns {} and te_residual is pure continuity), dCm2 = 1.485 m/s (1.02x) and dpsi = -0.271 pt
# -- the prereg's chain, reproduced with no loss feedback at all. The ~3.2x total is real; it
# simply does not live in Cm2.
#
# WHY IT MATTERS BEYOND BOOKKEEPING: blockage reaches psi mainly through a PARASITIC channel,
# and parasitic loss adds WORK, NOT PRESSURE (docs/PHYSICS-RULES.md rule 6) -- so it is nearly
# invisible to PR. That independently explains falsification #4 (psi needs B ~ 0.09, PR needs
# B ~ 0.18: the two respond through DIFFERENT PHYSICS, so no single B satisfies both) and agrees
# with the parasitic probe (deleting 90% of parasitic work moves PR by 0.2%). See
# data/coefficients.md S3.1b.
#
# The prereg's own qualifying condition still holds ("REFUTED IF a zero-blockage run is not
# bit-identical... or the mixing loss does not respond to blockage") -- both do. What failed was
# not the experiment. It was the story we told about it, twice: verify against the source, not
# against our own summary -- and the summary was ours.
#
# ``blockage=0.0`` reproduces the S3 numbers above BIT-FOR-BIT
# (test_exit_blockage.py::test_zero_blockage_is_bit_identical, and verified directly
# against a full Stage.solve with Impeller.blockage forced to 0.0).
#
# Slice S5 (docs/centrifugal/17-tdd-plan.md; tests/unit/test_diffuser_vane_count.py) moves
# stage_PR AGAIN, from 4.9798 (+6.30%) to 4.9502 (+5.67%) -- BETTER, TOWARD the
# measurement, and this time by construction of a GEOMETRY fix, not a coefficient choice.
# ``VanedDiffuser._cascade`` built ONE cascade with n_vanes = n_vanes + n_splitters = 40
# and used it for EVERYTHING (deviation, incidence, friction, solidity, diffusion factor).
# But HECC's diffuser SPLITTER VANE leading edge sits at r = 0.2488 m, inside the
# r3 = 0.23133 -> r4 = 0.28459 m span -- only 20 (not 40) vanes exist over the first third
# of the passage, exactly where circulation per vane is largest. Solidity sits in
# Lieblein's Df DENOMINATOR, so crediting 40 vanes there UNDERSTATED the loading loss.
#
# THE FIX is station-dependent, not a blanket vane-count change:
#   - deviation (Carter's rule, a TE quantity) and friction (wetted-passage hydraulic
#     diameter) KEEP n_vanes + n_splitters = 40 -- both blade rows DO reach the TE. LEGAL,
#     UNCHANGED.
#   - the diffusion factor / loading loss (an LE -> TE quantity) uses
#     ``VanedDiffuser.Z_eff_loading``, the SAME Aungier length-weighted formula already
#     used for the impeller's Z_eff (Yang, Liu & Zhao 2023, Machines 11(1):118, Eq. 1),
#     with the vane CHORDS from NASA Appendix C Tables C.25-C.28
#     (extract_hecc_vane_angles.py): L_main = 0.0532632 m (2.097 in),
#     L_splitter = 0.0358396 m (1.411 in) -> Z_eff,vd = 33.4575.
#
# Measured: Df rose 1.078 -> 1.136 (sigma fell 1.314 -> 1.099, as Df's denominator
# predicts); the diffuser's own P0 loss (vaneless-space exit -> stage exit, the metric
# that reproduces NASA/CR-2014-218114/REV1 Table A.5's "total pressure loss from the
# diffuser to exit") rose from 6.320% to 6.876% -- BELOW NASA's MEASURED CAP of 7.17%
# (headroom 0.86 pt; this slice used 0.556 pt of it and left ~0.29 pt of the cap
# unused). Stage PR fell 0.63 percentage points (+6.30% -> +5.67%), inside the
# pre-registered -0.3 to -1.2 pt band. psi is UNTOUCHED (bit-identical: this slice
# changes only the stationary diffuser, downstream of the impeller) and the vane-LE
# incidence is UNTOUCHED (bit-identical: -3.854 deg before and after -- it is set by the
# unchanged 40-vane cascade's beta_le_deg and the unchanged vaneless-space exit state).
#
# R9 (docs/centrifugal/20-review-fixes.md) moves stage_PR AGAIN, from 4.9502 (+5.67%) to
# 4.9540 (+5.75%) -- WORSE, and PRE-REGISTERED as such: ``ImpellerClearanceJansen`` and
# ``ImpellerMixingAungier`` both consumed ``Z = Z_eff = 25.4048`` (Aungier's fractional,
# wetted-length splitter treatment) at the impeller EXIT, where 30 physical blades (15
# main + 15 splitters) actually exist -- the SEVENTH instance of this project's
# cross-station category error (data/coefficients.md's own S3.1 "OPEN, recorded, NOT
# changed by slice S2" row, from the day slice S2 fixed ``_diffusion_factor``'s Z but
# deliberately left these two alone). Both models are INTERNAL: more (correctly-counted)
# blades means less loss per blade, so fixing the category error can only RAISE delivered
# pressure and widen the gap to NASA -- the same direction every honestly-sourced fix in
# this file has moved things.
#
# PRE-REGISTERED (before running): Jansen clearance scales ~sqrt(1/Z), 947 -> ~871 J/kg
# (-8%); Aungier mixing's dW = 2*pi*D2*Cu2/(Z*L_tilde) scales as 1/Z, x(25.4048/30) =
# x0.847; combined, stage PR +5.67% -> ~+5.9%.
#
# MEASURED: ImpellerClearanceJansen 947.0726 -> 871.6625 J/kg (-7.96%) -- matches the
# prereg closely. ImpellerMixingAungier 1722.5463 -> 1725.3591 J/kg (+0.16%) -- the
# prereg's linear-in-Z reasoning for dW does NOT survive contact with the rest of the
# formula: with D_eq = W_max/W2 > 2 at this operating point, Delta h_mix = 0.5*(W2*D_eq/2
# - W_out)^2 depends on dW only through W_max = (W1xi+W2+dW)/2, a small perturbation on a
# sum dominated by W1xi+W2 -- so an 18% change in dW (dW itself DOES scale ~1/Z) barely
# moves the squared difference. The two internal losses were pre-registered to move in
# the SAME direction; only one did so by the predicted amount. Net: stage PR rose from
# +5.67% to +5.75%, NOT to the pre-registered ~+5.9% -- the mixing loss's near-zero
# response is why the combined move undershot the prediction. Reported, not adjusted to
# match the prereg.
#
# R10 (docs/centrifugal/20-review-fixes.md; data/coefficients.md S3.4a) moves stage_PR
# AGAIN, from 4.954036 (+5.749%) to 4.955847 (+5.788%) -- WORSE, and VERIFIED (not
# assumed) before this commit: the diffuser's skin-friction hydraulic diameter d_h used
# the unweighted n_vanes + n_splitters = 40 (cas.n_vanes), the OPPOSITE convention from
# the impeller's own friction d_h, which already uses the fractional, length-weighted
# Impeller.Z_eff = 25.4048 for the identical quantity class one component upstream. d_h
# is a WETTED-PERIMETER quantity -- it must count the passages that actually exist over
# the chord, length-weighted for the splitter's partial existence, exactly what
# VanedDiffuser.Z_eff_loading (33.4575) already computes and was, until this commit,
# wired only into the loading loss. Both prior choices (diffuser: 40; impeller: already
# length-weighted) happened to be the loss-INCREASING option -- the asymmetry that
# hides. Applying the same convention to the diffuser's d_h (cas.n_vanes ->
# cas_load.n_vanes) is the SIXTH time honest sourcing has moved this model further from
# NASA, not closer.
#
# MEASURED: stage PR +0.0018107 (4.954036341841492 -> 4.955847012745787); psi
# BIT-IDENTICAL (the diffuser sits downstream of where psi is measured); the diffuser's
# own P0-loss metric FALLS 6.8774% -> 6.8434% (it now uses LESS of NASA's 7.168%
# headroom: 0.291 pt -> 0.325 pt), so the R9-class hard cap does NOT trip.
#
# Re-pin these in the SAME COMMIT as any slice that legitimately moves them, and say in the
# commit message why they moved. Never widen the tolerance.
#
# E-R11 (docs/centrifugal/41-prereg-adoption.md, results in 42-r11-result.md) moves all
# three AGAIN, via three primary-sourced corrections adopted together:
#   C1  ImpellerMixingAungier's tangential moves from the INLET (U1) to the EXIT relative
#       tangential (U2 - Vt2), per Aungier (1995) Eq. (26) -- which ANNIHILATES the term
#       (HECC dh_mix 1735 -> 0.54 J/kg). Less internal loss => MORE pressure.
#   C2  ImpellerLeakageAungier gains the missing *U2 factor (Aungier Eq. 10 / Oh 1997
#       Table 6). PARASITIC, so it adds work, not pressure: stage PR moves 0.020 pp.
#   C3  f_c, Aungier's head-loss correction (Eqs. 31-33), NEW and default-ON, SCALES THE
#       INTERNAL LOSS SUM by 1.2883 at this point. More internal loss => LESS pressure.
# Net: PR 4.932766 (+5.30% vs NASA) -> 4.942431 (+5.50%). C1 alone would have given
# +7.23%; C3 claws most of that back. psi rises because the parasitic leakage term (C2)
# now carries real work.
#
# ⭐ RE-PINNED BY THE SKIN-FRICTION (W-bar) FIX -- docs/centrifugal/45-jansen-skinfriction-fix.md.
# ImpellerSkinFrictionJansen's velocity average was Kovar et al. 2021 Eq. (33),
# W_bar = (2*W2 + W1s - W1h)/4, whose weights sum to 2 over a denominator of 4: it returned W/2
# for a constant velocity field, and 103.6 m/s at HECC -- BELOW the minimum (148.6 m/s) of the
# velocities it averaged. Replaced with Oh, Yoon & Chung (1997) Table 6's printed five-term form,
# which normalises. Delta_h_sf 1213.8 -> 4426.7 J/kg; stage PR error +5.50% -> +2.05%.
# These pins are UPDATED, not relaxed: the tolerance is still rel=1e-12.
PSI_AFTER_S5 = 0.8127325809383643  # W-bar fix; was 0.82479529307983035 (E-R11)
STAGE_PR_AFTER_S5 = (
    4.780524485318113  # W-bar fix; was 4.9424308497058096 (+5.50% -> +2.05% vs NASA)
)
STAGE_ETA_POLY_REALGAS_AFTER_S5 = (
    0.853991995752914  # W-bar fix; was 0.86204978793140929
)


def test_hecc_design_point_is_bit_identical_after_s5():
    with warnings.catch_warnings():
        # Same suppression as tests/fixtures/hecc_stage.py -- this fixture is deep in
        # ImpellerRecirculationOh's exponential-warning region by design (see
        # data/coefficients.md S4); that is not what this test is checking.
        warnings.simplefilter("ignore")
        op = build(BACKSWEEP).solve(
            mdot=MDOT_DESIGN, rpm=RPM, inlet=InletState(P0=P01, T0=T01)
        )

    assert op.psi == pytest.approx(PSI_AFTER_S5, rel=1e-12)
    assert op.stage_PR == pytest.approx(STAGE_PR_AFTER_S5, rel=1e-12)
    assert op.stage_eta_poly_realgas == pytest.approx(
        STAGE_ETA_POLY_REALGAS_AFTER_S5, rel=1e-12
    )
