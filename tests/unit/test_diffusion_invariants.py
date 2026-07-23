"""The invariants a stationary component CANNOT violate.

These are cheap, they are exact, and they are the tests that catch the two error classes
that actually killed this codebase before: FRAME confusion and SILENT DEGENERACY.

1. T0 IS CONSERVED. A stationary row has U = 0, so it does no work. There is no physics
   by which a vaneless space, a vaned diffuser or an EGV can change the total temperature.
   If one does, it is a bug -- not a modelling choice. (This is the stationary counterpart
   of the rothalpy check on the rotor.)

2. LOSS-FREE MUST REDUCE TO THE ANALYTIC ANSWER. With Cf = 0 the vaneless space must
   reproduce the free vortex r*Vt = const EXACTLY, and leave P0 untouched. A component
   that cannot reproduce its own analytic limit is not trustworthy anywhere else.

3. ENTROPY MUST NOT FALL. Second law. It caught a manufactured-entropy bug in slice 1
   when nothing else did.

4. DEVIATION IS NOT SLIP. Carter's rule must be what sets the exit angle. docs/PHYSICS-RULES.md rule
   7: applying a slip model to a vaned diffuser is a physics error that returns a
   plausible number -- omega = 0 kills the Coriolis term the slip derivation rests on.
"""

from __future__ import annotations


import pytest

from turbodesign.centrifugal.diffusion import (
    ExitGuideVanes,
    VanedDiffuser,
    VanelessSpace,
    carters_deviation,
    state_from_totals,
)
from turbodesign.centrifugal.state import Air

FLUID = Air()
MDOT = 4.9269

# A representative HECC impeller-exit state: highly swirled, transonic in the ABSOLUTE
# frame but comfortably subsonic meridionally (docs/PHYSICS-RULES.md rule 5).
R2, B2 = 0.21581, 0.015467


def impeller_exit() -> "object":
    return state_from_totals(
        P0=520000.0, T0=475.0, Vm=110.0, Vt=440.0, r=R2, b=B2, fluid=FLUID, s=0.0
    )


# ------------------------------------------------------------------ T0 conservation


@pytest.mark.parametrize(
    "component",
    [
        VanelessSpace(r3=0.231331, b3=0.014199),
        VanedDiffuser(
            r3=0.231331,
            r4=0.284455,
            b=0.014199,
            n_vanes=20,
            n_splitters=20,
            beta_le_deg=74.0,
            beta_te_deg=50.0,
        ),
        ExitGuideVanes(
            r_mean=0.3065,
            b=0.020,
            n_vanes=60,
            chord=0.061925,
            beta_le_deg=45.0,
            beta_te_deg=20.0,
        ),
    ],
    ids=["vaneless", "vaned_diffuser", "egv"],
)
def test_a_stationary_component_cannot_change_T0(component):
    """U = 0 -> no work -> T0 conserved. EXACTLY. There is no physics that permits
    otherwise, so this is a bit-level invariant, not a tolerance."""
    inlet = impeller_exit()
    out = component.solve(inlet, MDOT, FLUID)
    assert out.T0 == pytest.approx(inlet.T0, rel=1e-12), (
        f"{component.name} changed T0 by {out.T0 - inlet.T0:+.6f} K. A stationary row "
        f"does NO WORK -- this is a frame error, not a loss."
    )


@pytest.mark.parametrize(
    "component",
    [
        VanelessSpace(r3=0.231331, b3=0.014199),
        VanedDiffuser(
            r3=0.231331,
            r4=0.284455,
            b=0.014199,
            n_vanes=20,
            n_splitters=20,
            beta_le_deg=74.0,
            beta_te_deg=50.0,
        ),
    ],
    ids=["vaneless", "vaned_diffuser"],
)
def test_a_stationary_component_cannot_raise_P0(component):
    """Losses destroy P0. Nothing here can create it."""
    inlet = impeller_exit()
    out = component.solve(inlet, MDOT, FLUID)
    assert out.P0 <= inlet.P0 * (1 + 1e-12)


@pytest.mark.parametrize(
    "component",
    [
        VanelessSpace(r3=0.231331, b3=0.014199),
        VanedDiffuser(
            r3=0.231331,
            r4=0.284455,
            b=0.014199,
            n_vanes=20,
            n_splitters=20,
            beta_le_deg=74.0,
            beta_te_deg=50.0,
        ),
    ],
    ids=["vaneless", "vaned_diffuser"],
)
def test_entropy_does_not_fall(component):
    """Second law. This is the check that caught slice 1 manufacturing entropy from an
    over-determined gas model when nothing else did."""
    inlet = impeller_exit()
    out = component.solve(inlet, MDOT, FLUID)
    assert out.s >= inlet.s - 1e-9


# ------------------------------------------------------------------ the analytic limit


def test_frictionless_vaneless_space_is_a_free_vortex():
    """Cf = 0 -> r*Vt = const, EXACTLY, and P0 untouched.

    A component that cannot reproduce its own analytic limit cannot be trusted with
    friction on. This is the vaneless space's equivalent of the loss-free impeller check.
    """
    inlet = impeller_exit()
    out = VanelessSpace(r3=0.231331, b3=B2, Cf=0.0).solve(inlet, MDOT, FLUID)

    assert out.r * out.Vt == pytest.approx(inlet.r * inlet.Vt, rel=1e-9), (
        "with no friction, angular momentum MUST be conserved"
    )
    assert out.P0 == pytest.approx(inlet.P0, rel=1e-9), "no friction -> no P0 loss"
    assert out.T0 == pytest.approx(inlet.T0, rel=1e-12)


def test_friction_reduces_swirl_and_total_pressure():
    """With friction on, the free-vortex result must be strictly degraded -- and in the
    right DIRECTION. Friction removes angular momentum; it cannot add it."""
    inlet = impeller_exit()
    ideal = VanelessSpace(r3=0.231331, b3=B2, Cf=0.0).solve(inlet, MDOT, FLUID)
    real = VanelessSpace(r3=0.231331, b3=B2, Cf=0.005).solve(inlet, MDOT, FLUID)

    assert real.r * real.Vt < ideal.r * ideal.Vt, "friction must REMOVE angular momentum"
    assert real.P0 < ideal.P0, "friction must destroy total pressure"
    assert real.T0 == pytest.approx(ideal.T0, rel=1e-12), "but it cannot change T0"


def test_the_vaneless_space_diffuses():
    """Its job: raise static pressure by trading velocity for it."""
    inlet = impeller_exit()
    out = VanelessSpace(r3=0.231331, b3=0.014199).solve(inlet, MDOT, FLUID)
    assert out.P > inlet.P
    assert out.V < inlet.V


# ------------------------------------------------------------------ deviation, not slip


def test_carters_rule_scales_with_camber_and_solidity():
    """delta = 0.26 * theta / sqrt(solidity). Zero camber -> zero deviation; more
    solidity -> tighter guidance -> less deviation."""
    assert carters_deviation(0.0, 1.5) == pytest.approx(0.0)
    assert carters_deviation(24.0, 1.5) > carters_deviation(12.0, 1.5)
    assert carters_deviation(24.0, 3.0) < carters_deviation(24.0, 1.5)


def test_the_flow_does_not_leave_at_the_vane_metal_angle():
    """Deviation is real: the exit flow angle must EXCEED the TE metal angle (the flow is
    under-turned). A model that returns the metal angle has silently dropped deviation --
    which is exactly the class of error that made the impeller's incidence loss zero."""
    d = VanedDiffuser(
        r3=0.231331,
        r4=0.284455,
        b=0.014199,
        n_vanes=20,
        n_splitters=20,
        beta_le_deg=74.0,
        beta_te_deg=50.0,
    )
    out = d.solve(impeller_exit(), MDOT, FLUID)
    assert out.alpha_deg > d.beta_te_deg, "flow must be UNDER-turned relative to the metal angle"


def test_a_vaned_diffuser_takes_no_slip_model():
    """docs/PHYSICS-RULES.md rule 7, enforced structurally: there is no slip parameter to pass.

    Slip rests on the Coriolis term; at omega = 0 it vanishes and the derivation collapses.
    A vaned diffuser is a CASCADE. This test exists so that a future 'helpful' refactor
    that adds a slip= kwarg fails loudly.
    """
    import dataclasses

    fields = {f.name for f in dataclasses.fields(VanedDiffuser)}
    assert "slip" not in fields and not any("slip" in f.lower() for f in fields), (
        "a slip model must NEVER be applicable to a vaned diffuser (docs/PHYSICS-RULES.md rule 7)"
    )


# ------------------------------------------------------------------ the pinch


def test_the_pinch_is_not_silently_dropped():
    """HECC pinches from b2 = 0.609 in to b3 = 0.559 in (-8.2%). A smaller exit area at
    the same mdot must show up as a HIGHER meridional velocity than the unpinched case."""
    inlet = impeller_exit()
    unpinched = VanelessSpace(r3=0.231331, b3=B2, Cf=0.0).solve(inlet, MDOT, FLUID)
    pinched = VanelessSpace(r3=0.231331, b3=0.014199, Cf=0.0).solve(inlet, MDOT, FLUID)
    assert pinched.b < unpinched.b
    assert pinched.Vm > unpinched.Vm, "a pinch must accelerate the meridional flow"
