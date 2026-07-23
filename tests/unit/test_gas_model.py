"""The gas model must be thermodynamically self-consistent.

An ideal gas has TWO independent constants. Supplying cp, gamma AND R over-determines
it, and if they do not satisfy cp = gamma*R/(gamma-1), the model manufactures entropy
from nothing -- a loss-free impeller then generates ds > 0, which is physically
impossible and invisible to every test that does not check entropy explicitly.

This actually happened: the Slice 1 spec asked for cp=1005, gamma=1.4, R=287
(inconsistent: 1.4*287/0.4 = 1004.5), and a LOSSLESS impeller generated
ds = +0.185 J/(kg K). These tests exist so it cannot happen again.
"""

import math

import pytest

from turbodesign.centrifugal import Air


def test_gamma_is_derived_not_an_input():
    """Air must not accept gamma as an independent constant."""
    with pytest.raises(TypeError):
        Air(cp=1005.0, R=287.0, gamma=1.4)  # type: ignore[call-arg]


def test_the_ideal_gas_identity_holds():
    """cp == gamma*R/(gamma-1), by construction rather than by luck."""
    air = Air()
    assert air.cp == pytest.approx(air.gamma * air.R / (air.gamma - 1.0), rel=1e-12)


def test_the_isentropic_exponent_is_unambiguous():
    """cp/R and gamma/(gamma-1) are the SAME number.

    They are equal only when the gas is consistent. When it is not, the two forms
    disagree and the answer depends on which one the author happened to type.
    """
    air = Air()
    assert air.cp / air.R == pytest.approx(air.gamma / (air.gamma - 1.0), rel=1e-12)


def test_cv_is_consistent():
    air = Air()
    assert air.cv == pytest.approx(air.cp - air.R, rel=1e-12)
    assert air.gamma == pytest.approx(air.cp / air.cv, rel=1e-12)


def test_speed_of_sound():
    """Needed for the Slice 10 choke criterion (on the MERIDIONAL Mach number)."""
    air = Air()
    assert air.speed_of_sound(288.15) == pytest.approx(
        math.sqrt(air.gamma * air.R * 288.15), rel=1e-12
    )
    assert air.speed_of_sound(288.15) == pytest.approx(340.3, abs=1.0)


@pytest.mark.parametrize("cp,R", [(1005.0, 287.0), (1004.7, 287.05), (1150.0, 260.0)])
def test_identity_holds_for_any_valid_gas(cp, R):
    """Not just the default -- the invariant is structural."""
    air = Air(cp=cp, R=R)
    assert air.cp == pytest.approx(air.gamma * air.R / (air.gamma - 1.0), rel=1e-12)
    assert air.cp / air.R == pytest.approx(air.gamma / (air.gamma - 1.0), rel=1e-12)
