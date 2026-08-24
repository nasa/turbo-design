from typing import Union
import numpy as np
import numpy.typing as npt
import math

ArrayLike = Union[float, npt.NDArray[np.float64]]


def _maybe_return_scalar(result: npt.NDArray[np.float64], *inputs: object) -> ArrayLike:
    """Return a float when all driving inputs were scalar, else an ndarray."""
    if all(np.isscalar(inp) for inp in inputs):
        return float(np.asarray(result))
    return np.asarray(result, dtype=float)


def IsenP(M:ArrayLike,gamma:float) -> ArrayLike:
    """Computes the ratio P0/Ps

    Args:
        M (np.ndarray): Mach Number
        gamma (float): specific heat ratio

    Returns:
        float: P0/P ratio 
    """
    M_arr = np.asarray(M, dtype=float)
    result = np.power((1+(gamma-1)/2.0 * M_arr*M_arr),gamma/(gamma-1))
    return _maybe_return_scalar(result, M)


def FindMachP0P(P0_P:ArrayLike,gamma:float) -> ArrayLike:
    """Finds the mach number given a P0/P ratio

    Args:
        P0_P (np.ndarray): ratio of total to static pressure
        gamma (float): specific heat ratio

    Returns:
        float: [description]
    """
    n = (gamma-1)/gamma
    P0_P_arr = np.asarray(P0_P, dtype=float)
    c = 2.0/(gamma-1) * (np.power(P0_P_arr,n) - 1.0)
    M = np.sqrt(c)
    return _maybe_return_scalar(M, P0_P)
    


def IsenT(M:ArrayLike,gamma:float) -> ArrayLike:
    """Computes T0/Ts

    Args:
        M (np.ndarray): _description_
        gamma (float): _description_

    Returns:
        float: Ratio of T0/Ts
    """
    M_arr = np.asarray(M, dtype=float)
    result = (1.0+(gamma-1.0)/2.0 *M_arr*M_arr)
    return _maybe_return_scalar(result, M)


def A_As(M:ArrayLike,gamma:float) -> ArrayLike:
    """Computes the ratio of Area to Throat Area give a given mach number and gamma 

    Args:
        M (np.ndarray): Mach Number
        gamma (float): Specific Heat Ratio 

    Returns:
        float: Area to throat area ratio 
    """
    a = (gamma+1.0)/(2.0*(gamma-1.0))
    temp1 = np.power((gamma+1.0)/2.0,-a)
    M_arr = np.asarray(M, dtype=float)
    temp2 = np.power((1+(gamma-1)/2*M_arr*M_arr),a)/M_arr
    result = temp1*temp2
    return _maybe_return_scalar(result, M)


def mass_flow_function(M:ArrayLike,gamma:float) -> ArrayLike:
    """Non-dimensional mass flow function (Mattingly), a pure function of Mach
    number and gamma only - gas- and scale-agnostic.

        m~ = M * (1 + (gamma-1)/2 * M^2) ^ (-(gamma+1)/(2*(gamma-1)))

    Note: this is the "pure" non-dimensional form (no sqrt(gamma) factor).
    `mass_flow_parameter` layers the dimensional sqrt(gamma/R) scaling on top
    of this to recover Mattingly's mdot*sqrt(Tt)/(A*Pt) quantity.

    Identity: mass_flow_function(M,gamma) * A_As(M,gamma) == mass_flow_function_max(gamma)
    so choke_margin(M,gamma) == 1 - 1/A_As(M,gamma) - the fraction of annulus
    area in excess of the sonic throat area. This ratio is quadratically flat
    near M=1 (margin is only ~0.04 at M=0.8, ~0.01 at M=0.9) so it should be
    reported alongside M, not in place of it, near the choke point.

    Args:
        M (np.ndarray): Mach Number
        gamma (float): specific heat ratio

    Returns:
        float: Non-dimensional mass flow function m~
    """
    M_arr = np.asarray(M, dtype=float)
    gamma_val = float(gamma)
    expo = -(gamma_val + 1.0) / (2.0 * (gamma_val - 1.0))
    result = M_arr * np.power(1.0 + (gamma_val - 1.0) / 2.0 * M_arr * M_arr, expo)
    return _maybe_return_scalar(result, M)


def mass_flow_function_max(gamma:float) -> float:
    """Sonic (M=1) value of `mass_flow_function`, in closed form.

    Args:
        gamma (float): specific heat ratio

    Returns:
        float: m~_max = ((gamma+1)/2) ^ (-(gamma+1)/(2*(gamma-1)))
    """
    gamma_val = float(gamma)
    expo = -(gamma_val + 1.0) / (2.0 * (gamma_val - 1.0))
    return float(np.power((gamma_val + 1.0) / 2.0, expo))


def mass_flow_parameter(M:ArrayLike,gamma:float,R:float=287.0) -> ArrayLike:
    """Mattingly's dimensional-per-sqrt(R) mass flow parameter:

        MFP = sqrt(gamma/R) * mass_flow_function(M,gamma)

    such that mdot = A*P0/sqrt(T0) * MFP(M,gamma,R). See `Massflow`.

    Args:
        M (np.ndarray): Mach Number
        gamma (float): specific heat ratio
        R (float): Ideal gas constant [J/(kg*K)]. Defaults to 287 (air).

    Returns:
        float: Mass flow parameter
    """
    gamma_val = float(gamma)
    R_val = float(R)
    m_tilde = np.asarray(mass_flow_function(M, gamma_val), dtype=float)
    result = np.sqrt(gamma_val / R_val) * m_tilde
    return _maybe_return_scalar(result, M)


def mass_flow_function_required(massflow:ArrayLike,P0:ArrayLike,T0:ArrayLike,A:ArrayLike,gamma:float,R:float,blockage:float=0.0) -> ArrayLike:
    """Non-dimensional mass flow function required to pass `massflow` through
    area `A` at the given total conditions - the inverse of `Massflow`.

    Args:
        massflow (float): Target massflow [kg/s]
        P0 (float): Total pressure [Pa]
        T0 (float): Total temperature [K]
        A (float): Flow area [m^2]
        gamma (float): specific heat ratio
        R (float): Ideal gas constant [J/(kg*K)]
        blockage (float): Fractional area blockage (0 to 1). Defaults to 0.

    Returns:
        float: Required non-dimensional mass flow function m~_req. Compare
        against `mass_flow_function_max(gamma)` to test feasibility (M<=1).
    """
    massflow_arr = np.asarray(massflow, dtype=float)
    P0_arr = np.asarray(P0, dtype=float)
    T0_arr = np.asarray(T0, dtype=float)
    A_arr = np.asarray(A, dtype=float)
    gamma_val = float(gamma)
    R_val = float(R)
    A_eff = A_arr * (1.0 - blockage)
    result = massflow_arr * np.sqrt(T0_arr) / (A_eff * P0_arr * np.sqrt(gamma_val / R_val))
    return _maybe_return_scalar(result, massflow, P0, T0, A)


def choke_margin(M:ArrayLike,gamma:float) -> ArrayLike:
    """Fraction of flow capacity remaining before choking (M=1).

        choke_margin = 1 - mass_flow_function(M,gamma) / mass_flow_function_max(gamma)

    `mass_flow_function` peaks at M=1 on BOTH sides (it decreases moving away
    from M=1 whether subsonic or supersonic), so this margin is >= 0 for any
    M and exactly 0 only at M=1 - it does not distinguish subsonic from
    supersonic. In this codebase every Mach this is evaluated at comes from a
    bounded subsonic solve ([0, 1]), so in practice it reads as "how far
    below choked capacity," but the function itself is direction-agnostic.

    Args:
        M (np.ndarray): Mach Number
        gamma (float): specific heat ratio

    Returns:
        float: Choke margin
    """
    gamma_val = float(gamma)
    m_tilde = np.asarray(mass_flow_function(M, gamma_val), dtype=float)
    m_tilde_max = mass_flow_function_max(gamma_val)
    result = 1.0 - m_tilde / m_tilde_max
    return _maybe_return_scalar(result, M)


def min_area_for_massflow(massflow:ArrayLike,P0:ArrayLike,T0:ArrayLike,gamma:float,R:float,blockage:float=0.0) -> ArrayLike:
    """Minimum (sonic, M=1) annulus area needed to pass `massflow` at the
    given total conditions. Equivalent to `area_for_massflow(..., M=1.0, ...)`.

    Args:
        massflow (float): Target massflow [kg/s]
        P0 (float): Total pressure [Pa]
        T0 (float): Total temperature [K]
        gamma (float): specific heat ratio
        R (float): Ideal gas constant [J/(kg*K)]
        blockage (float): Fractional area blockage (0 to 1). Defaults to 0.

    Returns:
        float: Minimum flow area [m^2]
    """
    massflow_arr = np.asarray(massflow, dtype=float)
    P0_arr = np.asarray(P0, dtype=float)
    T0_arr = np.asarray(T0, dtype=float)
    gamma_val = float(gamma)
    R_val = float(R)
    m_tilde_max = mass_flow_function_max(gamma_val)
    result = massflow_arr * np.sqrt(T0_arr) / (P0_arr * np.sqrt(gamma_val / R_val) * m_tilde_max * (1.0 - blockage))
    return _maybe_return_scalar(result, massflow, P0, T0)


def area_for_massflow(massflow:ArrayLike,P0:ArrayLike,T0:ArrayLike,M:ArrayLike,gamma:float,R:float,blockage:float=0.0) -> ArrayLike:
    """Annulus area required to pass `massflow` at the given total conditions
    and target Mach number. Gas-agnostic and scale-agnostic sizing relation
    (the workhorse for scaling a component's annulus from a non-dimensional
    design point to a new massflow/gas/scale).

    Args:
        massflow (float): Target massflow [kg/s]
        P0 (float): Total pressure [Pa]
        T0 (float): Total temperature [K]
        M (float): Target Mach number at this station
        gamma (float): specific heat ratio
        R (float): Ideal gas constant [J/(kg*K)]
        blockage (float): Fractional area blockage (0 to 1). Defaults to 0.

    Returns:
        float: Required flow area [m^2]
    """
    massflow_arr = np.asarray(massflow, dtype=float)
    P0_arr = np.asarray(P0, dtype=float)
    T0_arr = np.asarray(T0, dtype=float)
    gamma_val = float(gamma)
    R_val = float(R)
    m_tilde = np.asarray(mass_flow_function(M, gamma_val), dtype=float)
    result = massflow_arr * np.sqrt(T0_arr) / (P0_arr * np.sqrt(gamma_val / R_val) * m_tilde * (1.0 - blockage))
    return _maybe_return_scalar(result, massflow, P0, T0, M)


def Massflow(P0:ArrayLike,T0:ArrayLike,A:ArrayLike,M:ArrayLike,gamma:float,R:float=287) -> ArrayLike:
    """Massflow rate calculation

    Args:
        P0 (float): Inlet Total Pressure (Pa)
        T0 (float): Inlet Total Temperature (K)
        A (float): Area (m^2)
        M (float): Mach Number
        gamma (float): Ratio of specific heats
        R (float): Ideal Gas Constant. Defaults to 287 J/(KgK).

    Returns:
        float: Massflow rate [kg/s]
    """
    P0_arr = np.asarray(P0, dtype=float)
    T0_arr = np.asarray(T0, dtype=float)
    A_arr = np.asarray(A, dtype=float)
    gamma_val = float(gamma)
    R_val = float(R)
    mdot = A_arr * P0_arr / np.sqrt(T0_arr) * np.asarray(mass_flow_parameter(M, gamma_val, R_val), dtype=float)
    return _maybe_return_scalar(mdot, P0, T0, A, M)


def solve_for_mach(M: float, massflow: float, P0: float, T0: float, area: float, gamma: float, R: float) -> float:
    """Residual between desired and estimated massflow for a guessed Mach number.

    Args:
        M (float): Mach number guess (dimensionless).
        massflow (float): Target massflow [kg/s].
        P0 (float): Total pressure [Pa].
        T0 (float): Total temperature [K].
        area (float): Flow area [m^2].
        gamma (float): Specific heat ratio Cp/Cv [-].
        R (float): Gas constant [J/(kg·K)].

    Returns:
        float: Absolute massflow residual [kg/s].
    """
    estimate = Massflow(P0, T0, area, M, gamma, R)
    residual = np.abs(massflow - estimate)
    return residual
