from __future__ import annotations

from typing import Optional

import numpy as np
import numpy.typing as npt
from scipy.optimize import minimize_scalar

from pyturbo.helper import convert_to_ndarray

from .bladerow import BladeRow
from .enums import LossType, RowType
from .isentropic import IsenP, IsenT, solve_for_mach
from .turbine_math import T0_coolant_weighted_average, compute_massflow

__all__ = ["stator_calc", "rotor_calc"]


def stator_calc(
    row: BladeRow,
    upstream: BladeRow,
    downstream: Optional[BladeRow] = None,
    calculate_vm: bool = True,
    static_defined: bool = False,
) -> None:
    """Solve compressor stator exit conditions."""

    def stator_calculation(Yp: npt.NDArray) -> npt.NDArray:
        row.Yp = Yp

        T0_coolant_local = 0.0
        if row.coolant is not None:
            T0_coolant_local = T0_coolant_weighted_average(row)
        row.T0 = upstream.T0 - T0_coolant_local

        if static_defined:
            row.P0 = upstream.P0 - row.Yp * (upstream.P0 - row.P)
        else:  # Compressor
            row.P0 = row.P0_is - row.Yp * (upstream.P0 - upstream.P)

        residual_error = np.zeros(len(row.area)-1)

        if downstream is not None:
            row.P0_P = float((row.P0 / downstream.P).mean())
            row.rp = ((row.P - downstream.P) / (upstream.P0 - downstream.P)).mean()

        deviation_func = getattr(row, "deviation_function", None)
        deviation = deviation_func(row, upstream) if callable(deviation_func) else 0.0
        row.deviation[:] = deviation
        
        if calculate_vm:
            # Get the static pressure from the massflow distribution.
            M = np.zeros(len(row.area))
            # Initial massflow fraction
            streamline_massflow = np.diff(row.massflow)
            
            for j in range(1, len(row.area)):
                res = minimize_scalar(solve_for_mach,bounds=[0.01, 1.0],
                    args=(
                        streamline_massflow[j - 1],
                        row.P0[j],
                        row.T0[j],
                        row.area[j],
                        row.gamma,
                        row.R,
                    ),
                )
                M[j] = res.x
                residual_error[j-1] = res.fun
            M[0] = 1.0 / (len(M) - 1) * M[1:].sum()  # Preserve the same average.
            row.M = M
            P0_P = IsenP(M, row.gamma)
            row.P = row.P0 / P0_P

            T0_T = IsenT(row.M, row.gamma)
            row.T = row.T0 / T0_T
            row.V = row.M * np.sqrt(row.gamma * row.R * row.T)
            row.Vm = row.V * np.cos(row.alpha2)
            row.Vx = row.Vm * np.cos(row.phi)
            row.Vr = row.Vm * np.sin(row.phi)
            row.Vt = row.Vm * np.tan(row.alpha2 + row.deviation)
        else:  # We know Vm, P0, T0, P
            row.Vx = row.Vm * np.cos(row.phi)
            row.Vr = row.Vm * np.sin(row.phi)
            row.Vt = row.Vm * np.tan(row.alpha2 + row.deviation)
            row.V = np.sqrt(row.Vx**2 + row.Vr**2 + row.Vt**2)
            row.T = row.P / (row.R * row.rho)  # We know P, this is a guess
            row.M = row.V / np.sqrt(row.gamma * row.R * row.T)

        if upstream.row_type == RowType.Rotor:
            row.alpha1 = upstream.alpha2 + upstream.deviation
            
        row.beta1 = upstream.beta2
        row.rho = row.P / (row.R * row.T)
        row.U = row.omega * row.r
        row.Wt = row.Vt - row.U
        row.P0_stator_inlet = upstream.P0
        row.entropy_rise = 0.5 * (row.Cp + upstream.Cp) * np.log(row.T / upstream.T) - row.R * np.log(row.P / upstream.P)
        return row.entropy_rise, residual_error

    loss_type = getattr(row.loss_function, "loss_type", None)

    if loss_type == LossType.Pressure:
        row.Yp = row.loss_function(row, upstream)
        stator_calculation(row.Yp)
    elif loss_type == LossType.Entropy:
        desired_entropy_rise = convert_to_ndarray(row.loss_function(row, upstream))
        if len(desired_entropy_rise) == 1:  # bulk value
            fun = lambda s: np.abs(desired_entropy_rise - stator_calculation(s))
            x = minimize_scalar(fun, bounds=[0.01, 0.4])
            row.Yp[:] = x
        elif len(desired_entropy_rise) > 1:  # entropy rise array
            fun = lambda s, j: np.abs(desired_entropy_rise[j] - stator_calculation(s)[j])
            for j in range(len(desired_entropy_rise)):
                x = minimize_scalar(fun, bounds=[0.01, 0.4], args=(j))
                row.Yp[j] = x
            stator_calculation(row.Yp)
    else:  # LossType.Enthalpy
        row.Yp[:] = 0
        stator_calculation(row.Yp)


def _log_rotor_failure(reason: str,row:BladeRow) -> None:
    def _fmt(val):
        try:
            return np.array2string(np.asarray(val), precision=5)
        except Exception:
            return str(val)

    print(f"[RotorCalc] Failure detected: {reason}")
    print(f"    row.T0R: {_fmt(row.T0R)}")
    print(f"    row.T: {_fmt(row.T)}")
    print(f"    row.W: {_fmt(row.W)}")
    print(f"    row.M: {_fmt(getattr(row, 'M', np.nan))}")
    print(f"    row.M_rel: {_fmt(getattr(row, 'M_rel', np.nan))}")
    print(f"    row.Yp: {_fmt(getattr(row, 'Yp', np.nan))}")
    if np.any(row.T >= row.T0R):
        print("    Note: T should be less than T0R.")
    yp_val = getattr(row, "Yp", None)
    if yp_val is not None and np.any(yp_val > 0.3):
        print(
            "    Note: row.Yp exceeded 0.3 which may indicate an issue with the design or loss model."
        )
        
def rotor_streamline_calculation(Yp: npt.NDArray,row:BladeRow,upstream:BladeRow,calculate_vm:bool=True):
    """Calculate the streamline quantities for the rotor

    Args:
        Yp (npt.NDArray): _description_

    Raises:
        ValueError: _description_

    Returns:
        npt.NDArray: _description_
    """
    row.Yp = Yp
    
    row.P0R = row.P0R_is - row.Yp * (upstream.P0R - upstream.P)
    deviation_func = getattr(row, "deviation_function", None)
    deviation = deviation_func(row, upstream) if callable(deviation_func) else 0.0
    
    # For compressors we use P0_ratio
    row.P0_is = upstream.P0 * row.P0_ratio
    row.P0 = row.P0_is - row.Yp * (upstream.P0 - upstream.P)
    
    M = np.linspace(0,1,10)
    P0_P = IsenP(M,row.gamma) # Lets vary the mach number to adjust the static pressure
    P = row.P0 / P0_P
    
    expo = -(row.gamma + 1.0) / (2.0 * (row.gamma - 1.0))
    
    estimate = row.area* row.P0/np.sqrt(T0)*np.sqrt(row.gamma / row.R)*M*np.power(1.0 + (row.gamma - 1.0) / 2.0 * M * M, expo)
    
    
    # row.deviation = deviation
    # residual_error = np.zeros(len(row.area)-1)
    # if calculate_vm:
    #     streamline_massflow = np.diff(row.massflow)
    #     M_rel = np.zeros(len(row.area))
    #     for j in range(1, len(row.area)):
    #         res = minimize_scalar(
    #             solve_for_mach,
    #             bounds=[0.01, 1.0],
    #             args=(
    #                 streamline_massflow[j - 1],
    #                 row.P0[j],
    #                 row.T0[j],
    #                 row.area[j],
    #                 row.gamma,
    #                 row.R,
    #             ),
    #         )
    #         M_rel[j] = res.x
    #         residual_error[j-1] = res.fun
    #     M_rel[0] = 1.0 / (len(M_rel) - 1) * M_rel[1:].sum()
    #     row.M_rel = M_rel
    #     P0_P = IsenP(M_rel, row.gamma)
    #     row.P = row.P0 / P0_P

    #     P0R_P = row.P0R / row.P
    #     T0R_T = P0R_P ** ((row.gamma - 1.0) / row.gamma)
    #     row.T = row.T0R / T0R_T
    #     row.W = np.sqrt(2.0 * row.Cp * (row.T0R - row.T))
    #     nan_in_velocity = np.isnan(np.sum(row.W))
    #     temp_issue = np.any(row.T >= row.T0R)
    #     high_loss = np.any(getattr(row, "Yp", 0) > 0.3)
    #     if nan_in_velocity:
    #         reason = "nan detected in relative velocity"
    #         if temp_issue:
    #             reason += "; T >= T0R shouldn't happen because of T-s diagram"
    #         if high_loss:
    #             reason += "; Yp > 0.3 This could be a problem with the loss model;"
    #         _log_rotor_failure(reason)
    #         raise ValueError("nan detected")
    #     row.Vr = row.W * np.sin(row.phi)
    #     row.Vm = row.W * np.cos(row.beta2 + row.deviation)
    #     row.Wt = row.W * np.sin(row.beta2)
    #     row.Vx = row.Vm * np.cos(row.phi)
    #     row.Vt = row.Wt + row.U
    #     row.V = np.sqrt(row.Vr**2 + row.Vt**2 + row.Vx**2)
    #     row.M = row.V / np.sqrt(row.gamma * row.R * row.T)
    #     row.Vm = np.sqrt(row.Vx**2 + row.Vr**2)
    #     row.T0 = row.T + row.V**2 / (2.0 * row.Cp)
    #     row.alpha2 = np.arctan2(row.Vt, row.Vm)
    # else:  # We know Vm, P0, T0
    #     row.Vr = row.Vm * np.sin(row.phi)
    #     row.Vx = row.Vm * np.cos(row.phi)
    #     row.W = np.sqrt(2.0 * row.Cp * (row.T0R - row.T))
    #     row.Wt = row.W * np.sin(row.beta2 + row.deviation)
    #     row.U = row.omega * row.r
    #     row.Vt = row.Wt + row.U
    #     row.alpha2 = np.arctan2(row.Vt, row.Vm)
    #     row.V = np.sqrt(row.Vm**2 * (1.0 + np.tan(row.alpha2) ** 2))
    #     row.M = row.V / np.sqrt(row.gamma * row.R * row.T)

    # T0_T = 1.0 + (row.gamma - 1.0) / 2.0 * row.M**2
    # row.P0 = row.P * T0_T ** (row.gamma / (row.gamma - 1.0))
    # row.P0_P = (row.P0_stator_inlet / row.P).mean()

    # row.M_rel = row.W / np.sqrt(row.gamma * row.R * row.T)
    # row.T0 = row.T + row.V**2 / (2.0 * row.Cp)
    # row.entropy_rise = 0.5 * (row.Cp + upstream.Cp) * np.log(
    #     row.T / upstream.T
    # ) - row.R * np.log(row.P / upstream.P)
    # return row.entropy_rise, residual_error

def rotor_calc(
    row: BladeRow,
    upstream: BladeRow,
    calculate_vm: bool = True,
) -> None:
    """Calculates quantities given beta2 

    Args:
        row (BladeRow): Rotor Row
        upstream (BladeRow): Stator Row or Rotor Row
    """
    def _log_rotor_failure(reason:str):
        def _fmt(val):
            try:
                return np.array2string(np.asarray(val), precision=5)
            except Exception:
                return str(val)

        print(f"[RotorCalc] Failure detected: {reason}")
        print(f"    row.T0R: {_fmt(row.T0R)}")
        print(f"    row.T: {_fmt(row.T)}")
        print(f"    row.W: {_fmt(row.W)}")
        print(f"    row.M: {_fmt(getattr(row,'M', np.nan))}")
        print(f"    row.M_rel: {_fmt(getattr(row,'M_rel', np.nan))}")
        print(f"    row.Yp: {_fmt(getattr(row,'Yp', np.nan))}")
        if np.any(row.T >= row.T0R):
            print("    Note: T should be less than T0R.")
        yp_val = getattr(row,'Yp', None)
        if yp_val is not None and np.any(yp_val > 0.3):
            print("    Note: row.Yp exceeded 0.3 which may indicate an issue with the design or loss model.")
    
    row.P0_is = upstream.P0*row.P0_ratio
    row.P0 = row.P0_is - row.Yp * (upstream.P0 - upstream.P)
    
    def calculate_vm_func(M):
        P0_P = IsenP(M,row.gamma) # Lets vary the mach number to adjust the static pressure
        row.P = row.P0 / P0_P

        upstream_radius = upstream.r
        row.U = row.omega*row.r
        # Upstream Relative Frame Calculations 
        upstream.U = upstream.rpm*np.pi/30 * upstream_radius # rad/s 
        upstream.Wt = upstream.Vt - upstream.U
        upstream.W = np.sqrt(upstream.Vx**2 + upstream.Wt**2 + upstream.Vr**2)
        upstream.beta2 = np.arctan2(upstream.Wt,upstream.Vm)
        upstream.T0R = upstream.T+upstream.W**2/(2*upstream.Cp)
        upstream.P0R = upstream.P * (upstream.T0R/upstream.T)**((upstream.gamma)/(upstream.gamma-1))      
        upstream.M_rel = upstream.W/np.sqrt(upstream.gamma*upstream.R*upstream.T)
        
        upstream_rothalpy = upstream.T0R*upstream.Cp - 0.5*upstream.U**2 # H01R - 1/2 U1^2 
        if np.any(upstream_rothalpy < 0):
            print('U is too high, reduce RPM or radius')
        # Rotor Exit Calculations
        row.beta1 = upstream.beta2
        #row.Yp # Evaluated earlier 
        row.P0R = upstream.P0R - row.Yp*(upstream.P0R-row.P)
        
        # Total Relative Temperature stays constant through the rotor. Adjust for change in radius from rotor inlet to exit
        row.T0R = upstream.T0R # (upstream_rothalpy + 0.5*row.U**2)/row.Cp # - T0_coolant_weighted_average(row) 
        P0R_P = row.P0R / row.P
        T0R_T = P0R_P**((row.gamma-1)/row.gamma)
        row.T = (row.T0R/T0R_T)     # Exit static temperature
        row.W = np.sqrt(2*row.Cp*(row.T0R-row.T)) #! nan popups here a lot for radial machines 
        nan_in_velocity = np.isnan(np.sum(row.W))
        temp_issue = np.any(row.T >= row.T0R)
        high_loss = np.any(getattr(row,'Yp',0) > 0.3)
        if nan_in_velocity:
            # Need to adjust T
            reason = "nan detected in relative velocity"
            if temp_issue:
                reason += "; T >= T0R shouldn't happen because of T-s diagram'"
            if high_loss:
                reason += "; Yp > 0.3 This could be a problem with the loss model;"
            _log_rotor_failure(reason)
            raise ValueError(f'nan detected')
        row.Vr = row.W*np.sin(row.phi)
        row.Vm = row.W*np.cos(row.beta2)
        row.Wt = row.W*np.sin(row.beta2)
        row.Vx = row.Vm*np.cos(row.phi)
        row.Vt = row.Wt + row.U 
        row.V = np.sqrt(row.Vr**2+row.Vt**2+row.Vx**2)
        row.M = row.V/np.sqrt(row.gamma*row.R*row.T)
        row.Vm = np.sqrt(row.Vx**2+row.Vr**2)
        row.T0 = row.T + row.V**2/(2*row.Cp)
        row.P0 = row.P*(row.T0/row.T)**(row.gamma/(row.gamma-1))
        row.alpha2 = np.arctan2(row.Vt,row.Vm)
        
        row.M_rel = row.W/np.sqrt(row.gamma*row.R*row.T)
        row.T0 = row.T+row.V**2/(2*row.Cp)
        compute_massflow(row)
        return np.abs(upstream.total_massflow - row.total_massflow)
    
    if calculate_vm:
        
    
    else: # We know Vm, P0, T0
        row.Vr = row.Vm*np.sin(row.phi)
        row.Vx = row.Vm*np.cos(row.phi)
        
        row.W = np.sqrt(2*row.Cp*(row.T0R-row.T))
        row.Wt = row.W*np.sin(row.beta2)
        row.U = row.omega * row.r 
        row.Vt = row.Wt+row.U
        
        row.alpha2 = np.arctan2(row.Vt,row.Vm)
        row.V = np.sqrt(row.Vm**2*(1+np.tan(row.alpha2)**2))
        
        row.M = row.V/np.sqrt(row.gamma*row.R*row.T)
        T0_T = (1+(row.gamma-1)/2 * row.M**2)
        row.P0 = row.P * T0_T**(row.gamma/(row.gamma-1))
    
    row.M_rel = row.W/np.sqrt(row.gamma*row.R*row.T)
    row.T0 = row.T+row.V**2/(2*row.Cp)
    compute_massflow(row)