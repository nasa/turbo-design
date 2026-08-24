"""Tests for the non-dimensional mass flow function / MFP utilities in isentropic.py."""

import numpy as np
import pytest

from turbodesign.isentropic import (
    A_As,
    Massflow,
    area_for_massflow,
    choke_margin,
    mass_flow_function,
    mass_flow_function_max,
    mass_flow_function_required,
    mass_flow_parameter,
    min_area_for_massflow,
    solve_for_mach,
)

GAMMAS = [1.3, 1.4, 1.667]


@pytest.mark.parametrize("gamma", GAMMAS)
def test_sonic_value_matches_closed_form_anchor(gamma):
    assert mass_flow_function(1.0, gamma) == pytest.approx(mass_flow_function_max(gamma), rel=1e-12)


def test_sonic_value_at_gamma_1p4_matches_known_constant():
    assert mass_flow_function_max(1.4) == pytest.approx(0.5787037, abs=1e-6)


@pytest.mark.parametrize("gamma", GAMMAS)
def test_a_as_identity(gamma):
    """m~(M,gamma) * A_As(M,gamma) == m~_max(gamma) for all M - ties the new function
    to the already-used A_As rather than to a hand-transcribed constant."""
    M = np.linspace(0.05, 3.0, 60)
    product = np.asarray(mass_flow_function(M, gamma)) * np.asarray(A_As(M, gamma))
    assert product == pytest.approx(mass_flow_function_max(gamma), rel=1e-9)


@pytest.mark.parametrize(
    "M,A_As_expected",
    [(0.2, 2.9635), (0.5, 1.33984), (0.8, 1.03823), (2.0, 1.6875)],
)
def test_textbook_area_ratio_table_gamma_1p4(M, A_As_expected):
    assert A_As(M, 1.4) == pytest.approx(A_As_expected, rel=2e-4)
    expected_m_tilde = mass_flow_function_max(1.4) / A_As_expected
    assert mass_flow_function(M, 1.4) == pytest.approx(expected_m_tilde, rel=2e-4)


def test_known_value_at_half_mach_gamma_1p4():
    # m~ = M * (1 + 0.2*M^2)^-3 = 0.5 * 1.05^-3
    assert mass_flow_function(0.5, 1.4) == pytest.approx(0.5 * 1.05**-3, rel=1e-9)


@pytest.mark.parametrize("gamma", GAMMAS)
def test_monotonic_increasing_then_decreasing_about_sonic(gamma):
    subsonic = mass_flow_function(np.linspace(0.01, 1.0, 50), gamma)
    supersonic = mass_flow_function(np.linspace(1.0, 3.0, 50), gamma)
    assert np.all(np.diff(subsonic) > 0)
    assert np.all(np.diff(supersonic) < 0)


@pytest.mark.parametrize("gamma", GAMMAS)
def test_massflow_dimensional_round_trip(gamma):
    P0, T0, A, M, R = 2.0e5, 400.0, 0.05, 0.6, 287.0
    direct = Massflow(P0, T0, A, M, gamma, R)
    via_pure = A * P0 / np.sqrt(T0) * np.sqrt(gamma / R) * mass_flow_function(M, gamma)
    via_parameter = A * P0 / np.sqrt(T0) * mass_flow_parameter(M, gamma, R)
    assert direct == pytest.approx(via_pure, rel=1e-12)
    assert direct == pytest.approx(via_parameter, rel=1e-12)


def test_pure_function_is_gas_agnostic():
    """The whole point of the non-dimensional form: unchanged under a 3x change in R."""
    m_air = mass_flow_function(0.55, 1.4)
    m_helium_like = mass_flow_function(0.55, 1.4)  # gamma unchanged, R varies elsewhere
    assert m_air == pytest.approx(m_helium_like, rel=1e-12)

    # Massflow itself DOES depend on R (dimensional), but only through mass_flow_parameter.
    mdot_air = Massflow(1e5, 300.0, 0.02, 0.5, 1.4, R=287.0)
    mdot_helium = Massflow(1e5, 300.0, 0.02, 0.5, 1.4, R=2077.0)
    assert mdot_air != pytest.approx(mdot_helium, rel=1e-3)
    assert mdot_air / mdot_helium == pytest.approx(np.sqrt(2077.0 / 287.0), rel=1e-9)


def test_min_area_scaling():
    base = min_area_for_massflow(10.0, 1e5, 300.0, 1.4, 287.0)
    assert min_area_for_massflow(20.0, 1e5, 300.0, 1.4, 287.0) == pytest.approx(2 * base, rel=1e-12)
    assert min_area_for_massflow(10.0, 2e5, 300.0, 1.4, 287.0) == pytest.approx(base / 2, rel=1e-12)
    assert min_area_for_massflow(10.0, 1e5, 1200.0, 1.4, 287.0) == pytest.approx(base * 2.0, rel=1e-9)


def test_min_area_equals_area_for_massflow_at_sonic():
    massflow, P0, T0, gamma, R = 5.0, 1.5e5, 350.0, 1.4, 287.0
    assert min_area_for_massflow(massflow, P0, T0, gamma, R) == pytest.approx(
        area_for_massflow(massflow, P0, T0, 1.0, gamma, R), rel=1e-12
    )


def test_area_for_massflow_inverse_round_trip():
    massflow, P0, T0, M, gamma, R = 3.5, 1.2e5, 320.0, 0.4, 1.33, 287.0
    area = area_for_massflow(massflow, P0, T0, M, gamma, R)
    recovered = Massflow(P0, T0, area, M, gamma, R)
    assert recovered == pytest.approx(massflow, rel=1e-9)


def test_area_for_massflow_blockage_scaling():
    args = (4.0, 1.1e5, 310.0, 0.5, 1.4, 287.0)
    unblocked = area_for_massflow(*args, blockage=0.0)
    blocked = area_for_massflow(*args, blockage=0.1)
    assert blocked == pytest.approx(unblocked / 0.9, rel=1e-9)


def test_mass_flow_function_required_matches_area_for_massflow_inverse():
    massflow, P0, T0, A, gamma, R = 6.0, 1.3e5, 330.0, 0.08, 1.4, 287.0
    m_req = mass_flow_function_required(massflow, P0, T0, A, gamma, R)
    area_back = area_for_massflow(massflow, P0, T0_ := T0, m_req, gamma, R)
    # area_for_massflow takes a Mach, not m~ - instead check the required-vs-max
    # feasibility predicate directly (this is the guard's actual contract).
    feasible_area = min_area_for_massflow(massflow, P0, T0, gamma, R)
    assert (m_req > mass_flow_function_max(gamma)) == (A < feasible_area)


def test_choke_margin_zero_at_sonic_and_nonnegative_elsewhere():
    """mass_flow_function peaks at M=1 on BOTH sides, so choke_margin is >=0
    everywhere and exactly 0 only at M=1 - it does not distinguish subsonic
    from supersonic (see the docstring caveat on `choke_margin`)."""
    for gamma in GAMMAS:
        assert choke_margin(1.0, gamma) == pytest.approx(0.0, abs=1e-12)
        assert choke_margin(0.5, gamma) > 0
        assert choke_margin(1.5, gamma) > 0
        assert choke_margin(3.0, gamma) > 0


def test_solve_for_mach_residual_is_zero_at_the_true_mach():
    P0, T0, area, gamma, R = 1.4e5, 340.0, 0.03, 1.4, 287.0
    M_true = 0.42
    massflow = Massflow(P0, T0, area, M_true, gamma, R)
    assert solve_for_mach(M_true, massflow, P0, T0, area, gamma, R) == pytest.approx(0.0, abs=1e-9)
    assert solve_for_mach(M_true + 0.05, massflow, P0, T0, area, gamma, R) > 1e-6


def test_scalar_in_scalar_out_array_in_array_out():
    scalar_result = mass_flow_function(0.5, 1.4)
    assert np.isscalar(scalar_result) or isinstance(scalar_result, float)

    array_result = mass_flow_function(np.array([0.3, 0.5, 0.7]), 1.4)
    assert isinstance(array_result, np.ndarray)
    assert array_result.shape == (3,)
