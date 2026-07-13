# tests/test_velocity_triangle_oracle.py
"""Triangle closure, checked against the triangle rather than against the library.

This module does not import turbodesign.

compressor_math.py:261 is  Vr = W*sin(phi).  It should be  Vm*sin(phi).  W carries the
tangential component (W**2 = Vm**2 + Wt**2), so building Vr from W and then forming
V = sqrt(Vr**2 + Vt**2 + Vx**2) at :267 counts the tangential part twice.

Ten of the eleven sites in the library that build Vr from phi already use Vm (:85, :177,
:365, inlet.py:175, radeq.py:61, and five in turbine_math.py). One uses W.

The repository's own radial-equilibrium-derivation.md gives the same closure. Section 11
("Closure -- V_r and V_T in terms of V_M"):

    V_r = V_M sin(phi),    V_T = V_M tan(alpha),    V^2 = V_M^2 (1 + tan^2(alpha))

and the Section 13 code-mapping table repeats it literally:

    | V_r=V_M\\sin\\phi | `Vr = Vm*np.sin(phi)` | closure, Section 11 |
"""

import numpy as np
import pytest


def triangle_from_Vm(W, beta2, phi, U):
    """The closure as derived: Vr = Vm*sin(phi)."""
    Vm = W * np.cos(beta2)
    Wt = W * np.sin(beta2)
    Vr = Vm * np.sin(phi)
    Vx = Vm * np.cos(phi)
    Vt = Wt + U
    return Vm, Vr, Vx, Vt, Wt


def triangle_from_W(W, beta2, phi, U):
    """The closure as coded at compressor_math.py:261: Vr = W*sin(phi)."""
    Vm = W * np.cos(beta2)
    Wt = W * np.sin(beta2)
    Vr = W * np.sin(phi)
    Vx = Vm * np.cos(phi)
    Vt = Wt + U
    return Vm, Vr, Vx, Vt, Wt


@pytest.mark.parametrize("phi_deg", [0.0, 15.0, 45.0, 90.0])
def test_meridional_closure_holds_when_Vr_is_built_from_Vm(phi_deg):
    """Vm**2 == Vx**2 + Vr**2. An identity, not a correlation."""
    Vm, Vr, Vx, _, _ = triangle_from_Vm(
        250.0, np.radians(-37.5), np.radians(phi_deg), 300.0
    )
    assert Vm**2 == pytest.approx(Vx**2 + Vr**2, rel=1e-12)


def test_the_bug_is_invisible_at_phi_zero():
    """phi == 0 is the definition of an axial station, and the error is (W - Vm)*sin(phi).
    That is why this has never shown up in an axial case.
    """
    a = triangle_from_Vm(250.0, np.radians(-37.5), 0.0, 300.0)
    b = triangle_from_W(250.0, np.radians(-37.5), 0.0, 300.0)
    assert a == pytest.approx(b, rel=1e-15)


@pytest.mark.parametrize("phi_deg", [15.0, 45.0, 90.0])
def test_the_bug_breaks_meridional_closure_off_axis(phi_deg):
    """Off-axis, Vm**2 != Vx**2 + Vr**2: the triangle does not close."""
    Vm, Vr, Vx, _, _ = triangle_from_W(
        250.0, np.radians(-37.5), np.radians(phi_deg), 300.0
    )
    assert Vm**2 != pytest.approx(Vx**2 + Vr**2, rel=1e-6)


def test_the_error_is_W_over_Vm_at_ninety_degrees():
    """At phi = 90 deg the overstatement is W/Vm = 1/cos(beta2). At 37.5 deg of exit swirl
    that is 1.26, so Vr comes out 26% too large. The error grows smoothly from zero at
    phi = 0, so any station with a sloped endwall carries some of it.
    """
    beta2 = np.radians(-37.5)
    _, Vr_ok, _, _, _ = triangle_from_Vm(250.0, beta2, np.pi / 2, 300.0)
    _, Vr_bug, _, _, _ = triangle_from_W(250.0, beta2, np.pi / 2, 300.0)
    assert Vr_bug / Vr_ok == pytest.approx(1.0 / np.cos(beta2), rel=1e-12)
    assert Vr_bug / Vr_ok == pytest.approx(1.26, abs=0.01)
