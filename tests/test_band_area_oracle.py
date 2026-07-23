# tests/test_band_area_oracle.py
"""The streamtube band area, checked against geometry rather than against the library.

This module does not import turbodesign. It is an independent oracle: a code that agrees
with itself is not validated.

The band swept by a straight meridional segment revolved about the machine axis is the
lateral surface of a conical frustum. By Pappus's centroid theorem,

    A = 2*pi*r_centroid*L = pi*(r1 + r2)*sqrt(dx**2 + dr**2)

which reduces to pi*(r2**2 - r1**2) at dx = 0 (an annulus) and to 2*pi*r*b at dr = 0 (a
cylinder). One formula, valid at any meridional angle.
"""

import numpy as np
import pytest
from scipy.integrate import quad


def buggy_band_area(x1, r1, x2, r2):
    """Verbatim reproduction of flow_math.py:40 (upstream 4018a6c).

    S is a length in metres, so S/2 * dx**2 is a cubic metre, and it is added to r1*dx,
    a square metre. tests/test_current_area_formula.py checks that this reproduction still
    matches the library.
    """
    dx = x2 - x1
    S = r2 - r1  # metres -- not a slope
    C = np.sqrt(1.0 + ((r2 - r1) / dx) ** 2)
    return 2 * np.pi * C * (S / 2 * dx**2 + r1 * dx)


def buggy_branching_area(x1, r1, x2, r2):
    """Verbatim reproduction of the branch guard at flow_math.py:32 (upstream 4018a6c):
    below the 1e-5 threshold on dx, the annulus formula pi*(r2**2 - r1**2); at or above
    it, buggy_band_area. Reproduced branch and all, so the test below exercises the actual
    discontinuity rather than a formula that resembles it.
    """
    dx = x2 - x1
    if np.abs(dx) < 1e-5:
        return np.pi * (r2**2 - r1**2)
    return buggy_band_area(x1, r1, x2, r2)


def pappus_band_area(x1, r1, x2, r2):
    return np.pi * (r1 + r2) * np.hypot(x2 - x1, r2 - r1)


def slope_variant_band_area(x1, r1, x2, r2):
    """The obvious repair to buggy_band_area: redefine S as the slope dr/dx rather than a
    length, leaving C and everything else unchanged. It fixes the dimensional error but
    not the sign, which is what test_area_is_never_negative_for_a_reversed_cut measures.
    """
    dx = x2 - x1
    S = (r2 - r1) / dx  # a slope now, not metres
    C = np.sqrt(1.0 + ((r2 - r1) / dx) ** 2)
    return 2 * np.pi * C * (S / 2 * dx**2 + r1 * dx)


def quadrature_band_area(x1, r1, x2, r2):
    """A second implementation of the same surface-of-revolution integral, evaluated
    numerically rather than through Pappus's closed form. For a straight-line generatrix --
    the only case this file exercises -- that closed form is this integral; the two are
    independent implementations of one theorem, not two theorems (see
    test_pappus_agrees_with_quadrature_on_a_cone).
    """
    dx, dr = x2 - x1, r2 - r1
    speed = np.hypot(dx, dr)
    val, _ = quad(
        lambda t: 2 * np.pi * (r1 + t * dr) * speed,
        0.0,
        1.0,
        epsabs=1e-14,
        epsrel=1e-14,
    )
    return float(val)


def test_pappus_agrees_with_quadrature_on_a_cone():
    """Two independent implementations, agreeing to machine precision -- not two
    independent theorems. For a straight-line generatrix, Pappus's centroid theorem is the
    closed form of the exact integral quadrature_band_area evaluates numerically; they are
    provably the same quantity, computed two ways. The cross-check catches a typo in either
    implementation (a wrong centroid, a missing factor of 2, wrong bounds, a transposed
    r1/r2). It would not catch an error shared by both -- a wrong definition of the surface
    itself -- because for this shape there is only one integral to get right.
    """
    a = pappus_band_area(0.0, 0.25, 0.10, 0.35)
    b = quadrature_band_area(0.0, 0.25, 0.10, 0.35)
    assert a == pytest.approx(b, rel=1e-12)


def test_a_cylinder_does_not_discriminate():
    """With dr = 0 the buggy formula is exact, so a cylinder passes against it. Any check
    of the area formula that uses one establishes nothing.
    """
    args = (0.0, 0.2, 0.05, 0.2)
    assert buggy_band_area(*args) == pytest.approx(pappus_band_area(*args), rel=1e-12)
    assert pappus_band_area(*args) == pytest.approx(2 * np.pi * 0.2 * 0.05, rel=1e-12)


def test_a_cone_discriminates():
    """The discriminating case: a sloped endwall, 0.25 -> 0.35 in radius over 0.10 axial.
    The buggy formula under-reports the area by more than 10%.
    """
    args = (0.0, 0.25, 0.10, 0.35)
    err = (buggy_band_area(*args) - pappus_band_area(*args)) / pappus_band_area(*args)
    assert err < -0.10, (
        f"expected the bug to under-report by more than 10%, got {err:.2%}"
    )


def test_pappus_reduces_to_the_annulus_when_dx_is_zero():
    a = pappus_band_area(0.0, 0.25, 0.0, 0.35)
    assert a == pytest.approx(np.pi * (0.35**2 - 0.25**2), rel=1e-12)


def test_pappus_is_continuous_across_the_branch_guard():
    """flow_math.py:32 guards on abs(dx) < 1e-5 -- on dx, not on the meridional angle. A
    2-micron change in a coordinate (9e-6 -> 1.1e-5) crosses the guard and switches formula:
    the annulus pi*(r2**2 - r1**2) below the threshold, buggy_band_area at or above it. That
    switch moves the result by about -16.7%. Pappus, built from hypot and multiplication, is
    continuous, and moves by about 2e-9 relative.
    """
    x1, r1, r2 = 0.0, 0.25, 0.35
    below = buggy_branching_area(x1, r1, 9e-6, r2)
    above = buggy_branching_area(x1, r1, 1.1e-5, r2)
    cliff = (above - below) / below
    assert cliff < -0.15, (
        f"expected the branch guard to produce a double-digit-percent jump for a "
        f"2-micron change in dx, got {cliff:.2%}"
    )

    p_below = pappus_band_area(x1, r1, 9e-6, r2)
    p_above = pappus_band_area(x1, r1, 1.1e-5, r2)
    p_change = (p_above - p_below) / p_below
    assert abs(p_change) < 1e-6, (
        f"Pappus should be continuous across the guard; got a relative change of "
        f"{p_change:.2e}"
    )


def test_area_is_never_negative_for_a_reversed_cut():
    """Pappus cannot go negative: hypot is unsigned, so a reversed cut (decreasing x) gives
    the same area as the forward cut. The slope-variant repair of buggy_band_area is not so
    lucky: C*dx = sign(dx)*L, so a cut ordered with decreasing x flips its sign. Measured
    below, and it is the reason to delete the branch rather than patch S.
    """
    forward = pappus_band_area(0.0, 0.25, 0.10, 0.35)
    backward = pappus_band_area(0.10, 0.35, 0.0, 0.25)
    assert forward > 0 and backward > 0
    assert forward == pytest.approx(backward, rel=1e-12)

    slope_forward = slope_variant_band_area(0.0, 0.25, 0.10, 0.35)
    slope_backward = slope_variant_band_area(0.10, 0.35, 0.0, 0.25)
    # The dimensional half of the repair is right: the forward magnitude matches Pappus.
    assert slope_forward == pytest.approx(forward, rel=1e-12)
    # The sign half is not: the reversed cut returns a negative area, which Pappus -- and
    # an area -- cannot.
    assert slope_backward < 0, (
        f"expected the slope-variant repair to still yield a negative area on a "
        f"reversed cut, got {slope_backward}"
    )
    assert slope_backward == pytest.approx(-forward, rel=1e-12)
