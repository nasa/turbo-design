# tests/test_relative_frame_oracle.py
"""The ideal relative-frame total pressure, checked against thermodynamics.

This module does not import turbodesign for its arithmetic: it is an independent oracle,
and a code that agrees with itself is not validated. The one exception is the last test,
which reads a converged solve through the `example_spool` fixture in order to compare the
library's own numbers against the closure derived here.

compressor_math.py:246 sets the ideal exit relative total pressure to the inlet value,
P0R2_ideal = P0R1. That holds only when U2 == U1. Three lines below, :249 computes T0R2
from rothalpy using the local U = omega*r. The two lines cannot both be right. The ideal
exit relative total pressure is isentropic in the rotating frame between T0R1 and T0R2:

    P0R2_ideal = P0R1 * (T0R2/T0R1)**(gamma/(gamma-1))

A note on the closure. The exponent form above assumes constant cp. This library's gas
model is Cantera, i.e. cp varies with temperature, and the rigorous isentropic closure for
a variable-cp ideal gas is

    ln(P0R2/P0R1)|_s = (1/R) * integral_{T0R1}^{T0R2} cp(T)/T dT

which collapses to the constant-gamma form only when cp is constant. The two are not
interchangeable by assumption, so the difference is measured here rather than asserted
away (test_the_two_closures_agree_to_within_a_tenth_of_the_delta). Over the temperature
span the examples actually cover it is 1.6% of the pressure-ratio delta, which is small
enough to quote the constant-gamma number -- provided that choice is stated. Re-run this
file if the gas model or the temperature range changes.
"""

import numpy as np
import pytest
from scipy.integrate import quad

R_AIR = 287.05
CP_AIR = 1004.0  # J/kg-K, near 290 K


def ideal_P0R_ratio_constant_gamma(T0R1: float, T0R2: float, gamma: float) -> float:
    return (T0R2 / T0R1) ** (gamma / (gamma - 1.0))


def ideal_P0R_ratio_variable_cp(
    T0R1: float, T0R2: float, cp_of_T, R: float = R_AIR
) -> float:
    integral, _ = quad(lambda T: cp_of_T(T) / T, T0R1, T0R2, epsabs=1e-12, epsrel=1e-12)
    return float(np.exp(integral / R))


def T0R2_from_rothalpy(T0R1: float, U1: float, U2: float, cp: float = CP_AIR) -> float:
    """Rothalpy I = h + W**2/2 - U**2/2 is conserved across a rotor. With h0_rel = cp*T0R,

        cp*T0R1 - U1**2/2 = cp*T0R2 - U2**2/2   ->   T0R2 = T0R1 + (U2**2 - U1**2)/(2*cp)

    The -U**2/2 term is what carries the radius change; it is identically zero only when
    U2 == U1.
    """
    return T0R1 + (U2**2 - U1**2) / (2.0 * cp)


def test_variable_cp_reduces_to_constant_gamma_when_cp_is_constant():
    """The two closures must agree exactly in the limit. If they do not, the oracle is wrong."""
    gamma = 1.4
    cp = gamma * R_AIR / (gamma - 1.0)
    a = ideal_P0R_ratio_constant_gamma(288.0, 320.0, gamma)
    b = ideal_P0R_ratio_variable_cp(288.0, 320.0, lambda T: cp)
    assert b == pytest.approx(a, rel=1e-12)


def test_the_ideal_ratio_is_one_when_the_blade_speed_is_unchanged():
    """One direction of the equivalence, and the only case in which :246 is right.

    U2 == U1 makes the rothalpy term vanish, so T0R2 == T0R1 and the ideal ratio is 1
    under both closures -- which is what compressor_math.py:246 hard-codes.
    """
    U = 291.0
    T0R2 = T0R2_from_rothalpy(288.0, U, U)
    assert T0R2 == 288.0
    assert ideal_P0R_ratio_constant_gamma(288.0, T0R2, 1.4) == 1.0
    assert ideal_P0R_ratio_variable_cp(
        288.0, T0R2, lambda T: 1004.0 + 0.04 * T
    ) == pytest.approx(1.0, rel=1e-15)


def test_the_ideal_ratio_departs_from_one_when_the_blade_speed_changes():
    """The other direction, and the defect. A rotor whose radius contracts by 10% along a
    streamline drops U from 290 to 261 m/s. Rothalpy then puts T0R2 8 K below T0R1, and the
    ideal relative total pressure ratio is 0.907 -- not the 1.0 that :246 assigns it.
    """
    T0R1, U1, U2 = 288.0, 290.0, 261.0
    T0R2 = T0R2_from_rothalpy(T0R1, U1, U2)
    assert T0R2 == pytest.approx(280.04, abs=0.01)

    ratio = ideal_P0R_ratio_constant_gamma(T0R1, T0R2, 1.4)
    assert ratio == pytest.approx(0.9066, abs=1e-4)
    assert abs(ratio - 1.0) > 0.05, (
        f"a 10% radius change should move the ideal ratio by several percent, got {ratio}"
    )


def test_the_two_closures_agree_to_within_a_tenth_of_the_delta():
    """Bound the disagreement between the constant-gamma and variable-cp closures; do not
    assume it is negligible.

    The range used is the one the library actually spans across the rotor of Mattingly
    example 9.1 at the hub streamline: T0R goes from 287.147118 K to 288.256793 K, i.e.
    1.11 K. cp(T) for air is represented by a linear fit whose slope is a bounding
    sensitivity, not a fitted real-gas value; the NASA-Glenn/McBride polynomial is not
    reproduced here. Over that span the two closures give 1.013529 and 1.013739. The
    quantity that matters is not the ratio but the delta from unity -- 0.013529 -- and the
    closures differ by 1.6% of it. That is small enough to quote the constant-gamma
    number, so long as the closure used is stated rather than left implicit.
    """
    cp0, dcp_dT = (
        1004.0,
        0.04,
    )  # J/kg-K, J/kg-K**2 -- bounding slope, not a fitted value
    T1, T2 = 287.147118, 288.256793  # example 9.1, hub streamline, across the rotor
    const = ideal_P0R_ratio_constant_gamma(T1, T2, 1.402560)
    varcp = ideal_P0R_ratio_variable_cp(T1, T2, lambda T: cp0 + dcp_dT * T)
    spread = abs(varcp - const) / (const - 1.0)  # relative to the delta, not to 1.0
    assert spread < 0.10, (
        f"the two closures disagree by {spread:.1%} of the delta ({const=}, {varcp=}). "
        f"A published delta must state which closure produced it."
    )


def test_the_meanline_is_not_a_no_op_under_a_variable_cp_gas(example_spool):
    """At the meanline of example 9.1 the radius is unchanged across the rotor, so the
    argument above says T0R2/T0R1 should be exactly 1 and :246 should be exactly right
    there. The library returns 0.99939 instead, because its gas is variable-cp: cp is not
    the same on the two rows, so cp*T0R1 - U**2/2 = cp*T0R2 - U**2/2 no longer forces
    T0R1 == T0R2.

    That is 6e-4 in T0R and about 2e-3 in the ideal pressure ratio. It is small, but it is
    not zero, and a test that pinned the meanline at 1 to within 1e-9 would be pinning the
    wrong number. The value is read from the solve rather than hardcoded, so it tracks the
    library.
    """
    loc = example_spool("mattingly-axial-compressor/example9.1.py")
    inlet, rotor = loc["spool"].blade_rows[0], loc["rotor"]

    mean = 1  # hub, mean, shroud
    assert rotor.r[mean] == pytest.approx(inlet.r[mean], rel=1e-12), (
        "example 9.1 is supposed to hold the mean radius constant across the rotor"
    )

    tau = rotor.T0R[mean] / inlet.T0R[mean]
    delta = ideal_P0R_ratio_constant_gamma(1.0, tau, 1.402560) - 1.0
    assert abs(delta) > 1e-3, (
        f"T0R is unchanged across the meanline to within {abs(tau - 1.0):.1e}; if this "
        f"ever becomes exact, the gas model has changed and this test's premise with it"
    )
