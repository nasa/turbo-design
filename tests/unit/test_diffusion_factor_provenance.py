"""Slice S2 (docs/centrifugal/17-tdd-plan.md) -- Galvas's contradiction.

THE FINDING. ``_diffusion_factor`` (``turbodesign/centrifugal/losses.py``) computes

    D_f = 1 - W2/W1s + (K_BL*dh_euler/U2**2) /
          [ (W1s/W2) * ( (Z/pi)*(1-D1s/D2) + 2*D1s/D2 ) ]

using ``K_BL = 0.75`` (unconditionally) and ``Z = state.Z``, i.e. ``Impeller.Z_eff`` --
Aungier's FRACTIONAL, splitter-aware effective blade count (25.4048 for HECC's 15+15
impeller).

Galvas, NASA TN D-7487 (1973), Eq. (B59) + p. 5-6, verbatim (verified directly against
``research/papers/galvas_1973_NASA-TN-D-7487_offdesign-centrifugal.pdf`` -- pdftotext,
this session):

    "a value of 0.75 is used for K_BL for conventional impellers and a value of 0.6 is
     used for impellers with splitters ... A parametric study of calculated diffusion
     factors with a variation in the number of blades indicated that changing the
     constant to 0.6 would compensate for the changing solidity near the exit."

and his FORTRAN listing (this session's pdftotext, ~line 1625-1630 of the extracted
text): ``CONST1=0.75`` ... ``IF(SPLT.EQ.1) CONST1=0.6`` -- read directly INTO ``DF=``.
And his input-list definition (same PDF, "Z  number of impeller blades at exit, Z3"):
Z is a BLADE COUNT, an integer quantity. Nowhere in Galvas does a fractional blade
count exist -- he explicitly considered varying Z for splitters and REJECTED it in
favour of the K_BL constant.

Our shipped expression is a hybrid that appears in NO source: Coppage/Galvas's D_f +
the NO-SPLITTER constant 0.75 + Aungier's fractional Z_eff from a DIFFERENT loss
framework (wetted area / hydraulic diameter, ``Impeller.Z_eff``'s own docstring).

THE FIX (this slice): adopt the one (Z, K_BL) pair Galvas's own code executes for a
splittered machine -- Z3 = 30 (an INTEGER; both blade rows reach HECC's TE) and
K_BL = 0.6 (the splitter branch). ``Z_eff`` (25.4048) is UNCHANGED and remains correct
for the loss models that are genuinely about wetted length/hydraulic diameter
(``ImpellerSkinFrictionJansen``'s ``d_h``) and for the slip model's ``Z``
(``data/coefficients.md`` S5, knob #6 -- retained, not turned). Only ``_diffusion_factor``
gets the integer exit count and the splitter-aware K_BL.

REFUTED IF (data/coefficients.md, docs/centrifugal/15-experiments-prereg.md E4): PR does
not rise when K_BL falls from 0.75 to 0.6, or ψ does not fall. The two MUST move in
OPPOSITE directions -- D_f drives one INTERNAL loss (blade loading, lowers PR when D_f
falls) and one PARASITIC loss (recirculation, lowers ψ when D_f falls); if they move
together, the internal/parasitic split (docs/PHYSICS-RULES.md rule 6) is broken.
"""

from __future__ import annotations

import math

import pytest

from turbodesign.centrifugal.losses import ImpellerLossState, _diffusion_factor
from turbodesign.centrifugal.state import Air

# ---------------------------------------------------------------------------------
# A real, self-consistent HECC design-point impeller-exit operating state -- captured
# from a converged Stage.solve on the CURRENT working tree (tests/fixtures/hecc_stage.py's
# machine: backsweep -36.0 deg, NASA Fig 18; WiesnerSlip; OhLossSet; design mdot/rpm),
# at the parasitic-block construction site (turbodesign/centrifugal/solver.py's
# ``state2``). Re-derive with:
#
#   uv run python -c "
#   import warnings, sys; sys.path.insert(0, 'tests/fixtures')
#   from hecc_stage import BACKSWEEP, MDOT_DESIGN, P01, RPM, T01, build
#   from turbodesign.centrifugal import InletState
#   import turbodesign.centrifugal.losses as L, turbodesign.centrifugal.solver as S
#   states = []
#   old = L.ImpellerLossState.__init__
#   def new(self, *a, **kw):
#       old(self, *a, **kw); states.append(self)
#   L.ImpellerLossState.__init__ = new; S.ImpellerLossState = L.ImpellerLossState
#   with warnings.catch_warnings():
#       warnings.simplefilter('ignore')
#       build(BACKSWEEP).solve(mdot=MDOT_DESIGN, rpm=RPM, inlet=InletState(P0=P01, T0=T01))
#   s = states[-1]
#   print(s.r1h, s.r1s, s.r2, s.b2, s.U1, s.Vm1, s.U2, s.Cm2, s.Vt2, s.W2, s.T2, s.Z)
#   "
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
Z_EFF = 25.404764740917212  # Impeller.Z_eff -- FRACTIONAL, splitter-aware, wetted-area sense
BACKSWEEP_DEG = -36.0
TIP_CLEARANCE = 0.000305
MDOT = 4.9269
OMEGA = 21789.0 * math.pi / 30.0


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


def _hand_computed_Df(state: ImpellerLossState, Z: float, K_BL: float) -> float:
    """The formula, evaluated by hand at an EXPLICIT (Z, K_BL) -- the reference
    against which ``_diffusion_factor``'s own field-driven selection is checked."""
    W2_, W1s = state.W2, state.W1s
    U2_, Vt2 = state.U2, state.Vt2
    D2, D1s = 2.0 * state.r2, 2.0 * state.r1s
    dh_euler = U2_ * Vt2
    denom = (W1s / W2_) * ((Z / math.pi) * (1.0 - D1s / D2) + 2.0 * D1s / D2)
    return 1.0 - W2_ / W1s + (K_BL * dh_euler / U2_**2) / denom


def test_splittered_impeller_uses_the_galvas_splitter_constant():
    """K_BL = 0.6 for a splittered impeller (Galvas Eq. B59 + FORTRAN
    ``IF(SPLT.EQ.1) CONST1=0.6``), NOT the shipped unconditional 0.75.

    Holds Z fixed at the integer exit count (30) and flips only ``has_splitters`` --
    isolating the K_BL branch from the Z branch (the next test isolates the reverse).
    """
    state_no_splitters = _hecc_state(Z_exit=30.0, has_splitters=False)
    state_splitters = _hecc_state(Z_exit=30.0, has_splitters=True)

    Df_no_splitters = _diffusion_factor(state_no_splitters)
    Df_splitters = _diffusion_factor(state_splitters)

    assert Df_no_splitters == pytest.approx(
        _hand_computed_Df(state_no_splitters, Z=30.0, K_BL=0.75)
    )
    assert Df_splitters == pytest.approx(
        _hand_computed_Df(state_splitters, Z=30.0, K_BL=0.6)
    )
    assert Df_splitters < Df_no_splitters, (
        "K_BL = 0.6 (splitter branch) must give a LOWER D_f than K_BL = 0.75 -- Galvas's "
        "own stated purpose ('compensate for the changing solidity near the exit')"
    )


def test_blade_count_in_Df_is_an_integer_exit_count():
    """D_f must use Z_exit = 30 (n_blades + n_splitters, an INTEGER), NOT the fractional
    Z_eff = 25.4048 that ``ImpellerLossState.Z`` carries for the slip model and the
    wetted-length loss models.

    Galvas's own input-list defines Z as "number of impeller blades at exit" -- a count.
    He never uses a splitter-weighted fraction inside D_f.
    """
    state = _hecc_state(Z_exit=30.0, has_splitters=True)

    Df = _diffusion_factor(state)
    Df_with_fractional_Z = _hand_computed_Df(state, Z=Z_EFF, K_BL=0.6)
    Df_with_integer_Z = _hand_computed_Df(state, Z=30.0, K_BL=0.6)

    assert Df == pytest.approx(Df_with_integer_Z)
    assert Df != pytest.approx(Df_with_fractional_Z, rel=1e-6), (
        "D_f must not be indifferent to feeding it the fractional Z_eff (25.4048) vs "
        "the integer exit count (30) -- if it were, state.Z_exit is not being read"
    )


def test_unsplittered_impeller_is_bit_identical():
    """An UNSPLITTERED impeller (Eckardt rotors O/A: n_splitters = 0) must be BIT-
    IDENTICAL after this change: K_BL = 0.75 (the conventional branch) and
    Z = n_blades (which already equals Z_eff when there are no splitters --
    ``Impeller.Z_eff`` returns ``float(n_blades)`` for ``n_splitters <= 0``).

    This is the load-bearing test: Kovář et al. 2021 validated this loss set against
    Eckardt with K_BL = 0.75. If an unsplittered machine's D_f moved, the whole
    multi-machine arbitration (docs/centrifugal/14-multi-machine-validation.md) would be
    comparing against a silently different physics model than the one Kovář et al. validated.
    """
    n_blades = (
        20.0  # representative unsplittered count (Eckardt-scale), not HECC's 15+15
    )
    state_old_api = _hecc_state(
        Z=n_blades
    )  # no Z_exit/has_splitters at all -- the OLD call site
    state_new_api = _hecc_state(Z=n_blades, Z_exit=n_blades, has_splitters=False)

    Df_old = _diffusion_factor(state_old_api)
    Df_new = _diffusion_factor(state_new_api)

    assert Df_new == pytest.approx(Df_old, rel=1e-12), (
        "an unsplittered impeller's D_f must not move: K_BL stays 0.75 and Z stays the "
        "(here, equal) blade count whether or not the caller populates the new fields"
    )
    assert Df_new == pytest.approx(
        _hand_computed_Df(state_new_api, Z=n_blades, K_BL=0.75)
    )


def test_Df_at_hecc_design_point():
    """The pre-registered number (data/coefficients.md, 15-experiments-prereg.md E4):
    D_f falls from ~0.546 (shipped: Z_eff=25.4048, K_BL=0.75) to ~0.524 (sourced:
    Z_exit=30, K_BL=0.6) at the HECC design point.
    """
    state = _hecc_state(Z_exit=30.0, has_splitters=True)
    Df = _diffusion_factor(state)
    assert Df == pytest.approx(0.524, abs=0.01)

    Df_shipped = _hecc_state()  # no new fields -- exercises the untouched fallback
    assert _diffusion_factor(Df_shipped) == pytest.approx(0.5463, rel=1e-3), (
        "sanity check: the SHIPPED (pre-slice) behaviour must still be reachable and "
        "must reproduce the documented 'before' value"
    )
