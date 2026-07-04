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
        float: Nusselt Number
    """
    P0_arr = np.asarray(P0, dtype=float)
    T0_arr = np.asarray(T0, dtype=float)
    A_arr = np.asarray(A, dtype=float)
    M_arr = np.asarray(M, dtype=float)
    gamma_val = float(gamma)
    R_val = float(R)
    mdot = A_arr * P0_arr/np.sqrt(T0_arr) * np.sqrt(gamma_val/R_val) * M_arr \
        *np.power(1.0+(gamma_val-1.0)/2.0 * M_arr*M_arr, -(gamma_val+1.0)/(2.0*(gamma_val-1.0)))
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
    expo = -(gamma + 1.0) / (2.0 * (gamma - 1.0))
    
    estimate = area* P0/np.sqrt(T0)*np.sqrt(gamma / R)*M*np.power(1.0 + (gamma - 1.0) / 2.0 * M * M, expo)
    residual = np.abs(massflow - estimate)
    return residual
