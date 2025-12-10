from __future__ import annotations

from typing import Optional

import numpy as np
import numpy.typing as npt
from scipy.optimize import minimize_scalar

from pyturbo.helper import convert_to_ndarray

from .bladerow import BladeRow
from .enums import LossType, RowType
from .isentropic import IsenP, IsenT, solve_for_mach
from .turbine_math import T0_coolant_weighted_average
from .flow_math import compute_massflow, compute_streamline_areas

__all__ = ["stator_calc", "rotor_calc", "polytropic_efficiency"]


def polytropic_efficiency(pi: float, tau: float, gamma: float) -> float:
    """Compute polytropic efficiency from pressure/temperature ratios.

    Args:
        pi: Total-pressure ratio (pt,out / pt,in).
        tau: Total-temperature ratio (Tt,out / Tt,in).
        gamma: Ratio of specific heats.

    Returns:
        Polytropic efficiency consistent with the provided ratios.
    """
    pi_safe = max(pi, 1e-8)
    tau_safe = max(tau, 1e-8)
    ln_tau = np.log(tau_safe)
    if abs(ln_tau) < 1e-12:
        return 1.0
    return ((gamma - 1.0) / gamma) * np.log(pi_safe) / ln_tau


def stator_calc(
    row: BladeRow,
    upstream: BladeRow,
    downstream: Optional[BladeRow] = None,
    calculate_vm: bool = True,
    static_defined: bool = False,
) -> None:
    """Solve compressor stator exit conditions.

    Args:
        row: Current stator blade row being solved.
        upstream: Upstream blade row providing inlet total conditions.
        downstream: Optional downstream row used for reaction/P0_P bookkeeping.
        calculate_vm: If True, iterates Mach to satisfy massflow; if False, assumes Vm known.
        static_defined: Treat P as prescribed (turbine-like) when True; otherwise compressor mode.
    """

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
    elif loss_type == LossType.Polytropic:
        target_eta_poly = float(row.loss_function(row, upstream))

        def _objective(y: float) -> float:
            yp_array = np.full_like(row.Yp, y)
            stator_calculation(yp_array)
            pi = float((row.P0 / upstream.P0).mean())
            tau = float((row.T0 / upstream.T0).mean())
            eta_poly = polytropic_efficiency(pi, tau, row.gamma)
            return abs(eta_poly - target_eta_poly)

        res = minimize_scalar(_objective, bounds=[0.0, 0.6], method="bounded")
        row.Yp = np.full_like(row.Yp, res.x)
        stator_calculation(row.Yp)
    else:  # LossType.Enthalpy
        row.Yp[:] = 0
        stator_calculation(row.Yp)


def rotor_calc(
    row: BladeRow,
    upstream: BladeRow,
    calculate_vm: bool = True,
    static_defined: bool = False,
) -> None:
    """Solve compressor rotor exit conditions.

    Args:
        row: Rotor blade row being solved.
        upstream: Upstream blade row providing inlet relative/absolute conditions.
        calculate_vm: If True, iterates Mach to satisfy massflow; if False, assumes Vm known.
        static_defined: Treat P as prescribed (turbine-like) when True; otherwise compressor mode.
    """
    row.P0_is = upstream.P0*row.P0_ratio
    row.P0 = row.P0_is - row.Yp * (upstream.P0 - upstream.P)
    
    # Upstream relative frame
    upstream.U = upstream.rpm * np.pi / 30 * upstream.r
    upstream.Wt = upstream.Vt - upstream.U
    upstream.W = np.sqrt(upstream.Vx ** 2 + upstream.Wt ** 2 + upstream.Vr ** 2)
    upstream.beta2 = np.arctan2(upstream.Wt, upstream.Vm)
    upstream.T0R = upstream.T + upstream.W ** 2 / (2 * upstream.Cp)
    upstream.P0R = upstream.P * (upstream.T0R / upstream.T) ** (upstream.gamma / (upstream.gamma - 1))
    upstream.M_rel = upstream.W / np.sqrt(upstream.gamma * upstream.R * upstream.T)
    upstream_rothalpy = upstream.T0R * upstream.Cp - 0.5 * upstream.U ** 2
    
    if np.any(upstream_rothalpy < 0):
        print('U is too high, reduce RPM or radius')
        
    def calculate_vm_func(M_rel: float, apply: bool = False):
        """Compute Vm of rotor locally guessing relative mach number at rotor exit. This calculation is done to balance the massflow

        Args:
            M (float): Relative mach number at rotor exit 
            apply (bool, optional): Apply calculations. Defaults to False.
        """
        # Use local scratch copies to avoid polluting row during optimizer iterations
        P0R_local = upstream.P0R - row.Yp * (upstream.P0R - upstream.P)
        T0R_local = upstream.T0R
        
        P_local = P0R_local / IsenP(M_rel, row.gamma)
        U_local = row.omega * row.r
        
        P0R_P = P0R_local / P_local
        T0R_T = P0R_P ** ((row.gamma - 1) / row.gamma)
        T_local = T0R_local / T0R_T
        W_local = np.sqrt(2 * row.Cp * (T0R_local - T_local))

        if np.isnan(W_local).any() or np.any(T_local >= T0R_local):
            return np.inf

        Vr_local = W_local * np.sin(row.phi)
        Vm_local = W_local * np.cos(row.beta2)
        Wt_local = W_local * np.sin(row.beta2)
        Vx_local = Vm_local * np.cos(row.phi)
        Vt_local = Wt_local + U_local
        V_local = np.sqrt(Vr_local ** 2 + Vt_local ** 2 + Vx_local ** 2)
        M_local = V_local / np.sqrt(row.gamma * row.R * T_local)
        T0_local = T_local + V_local ** 2 / (2 * row.Cp)
        P0_local = P_local * (T0_local / T_local) ** (row.gamma / (row.gamma - 1))

        # compute massflow using locals
        rho_local = P_local / (row.R * T_local)
        total_area, streamline_area = compute_streamline_areas(row)
        massflow_local = np.zeros_like(row.massflow)
        for j in range(1, len(row.percent_hub_shroud)):
            Vm_seg = 0.5 * (Vm_local[j] + Vm_local[j - 1])
            rho_seg = 0.5 * (rho_local[j] + rho_local[j - 1])
            massflow_local[j] = Vm_seg * rho_seg * streamline_area[j] * (1 - row.blockage) + massflow_local[j - 1]
        total_massflow_local = massflow_local[-1]

        if apply:
            row.P = P_local
            row.T = T_local
            row.W = W_local
            row.Vr = Vr_local
            row.Vm = Vm_local
            row.Wt = Wt_local
            row.Vx = Vx_local
            row.Vt = Vt_local
            row.V = V_local
            row.M = M_local
            row.T0 = T0_local
            row.P0 = P0_local
            row.alpha2 = np.arctan2(row.Vt, row.Vm)
            row.M_rel = W_local / np.sqrt(row.gamma * row.R * T_local)
            row.total_massflow = total_massflow_local
            row.total_massflow_no_coolant = total_massflow_local
            row.massflow = massflow_local
            row.total_area = total_area
            row.area = streamline_area
        return np.abs(upstream.total_massflow - total_massflow_local)
    
    if calculate_vm:
        res = minimize_scalar(calculate_vm_func, bounds=[0.01,1])
        # apply best solution to row
        calculate_vm_func(res.x, apply=True)
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
