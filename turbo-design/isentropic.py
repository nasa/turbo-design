import numpy as np
import math


def IsenP(M:np.ndarray,gamma:float) -> float:
    """Computes the ratio P0/Ps

    Args:
        M (np.ndarray): Mach Number
        gamma (float): specific heat ratio

    Returns:
        float: P0/P ratio 
    """
    return np.power((1+(gamma-1)/2.0 * M*M),gamma/(gamma-1))


def FindMachP0P(P0_P:np.ndarray,gamma:float) -> float:
    """Finds the mach number given a P0/P ratio

    Args:
        P0_P (np.ndarray): ratio of total to static pressure
        gamma (float): specific heat ratio

    Returns:
        float: [description]
    """
    n = (gamma-1)/gamma
    c = 2.0/(gamma-1) * (np.power(P0_P,n) - 1.0)

    M = np.sqrt(c)
    return M # Subsonic and supersonic solution
    


def IsenT(M:np.ndarray,gamma:float) -> float:
    """Computes T0/Ts

    Args:
        M (np.ndarray): _description_
        gamma (float): _description_

    Returns:
        float: Ratio of T0/Ts
    """
    return (1.0+(gamma-1.0)/2.0 *M*M)


def A_As(M:np.ndarray,gamma:float) -> float:
    """Computes the ratio of Area to Throat Area give a given mach number and gamma 

    Args:
        M (np.ndarray): Mach Number
        gamma (float): Specific Heat Ratio 

    Returns:
        float: Area to throat area ratio 
    """
    a = (gamma+1.0)/(2.0*(gamma-1.0))
    temp1 = np.power((gamma+1.0)/2.0,a)
    temp2 = np.power((1+(gamma-1)/2*M*M),-a)/M
    return temp1*temp2


def Massflow(P0:float,T0:float,A:float,M:float,gamma:float,R:float=287):
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
    mdot = A * P0/np.sqrt(T0) * np.sqrt(gamma/R) * M \
        *np.power(1.0+(gamma-1.0)/2.0 * M*M, -(gamma+1.0)/(2.0*(gamma-1.0)))
    
    return mdot