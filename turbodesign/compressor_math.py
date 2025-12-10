from __future__ import annotations

from typing import Optional

import numpy as np
import numpy.typing as npt
from scipy.optimize import minimize_scalar, minimize

from pyturbo.helper import convert_to_ndarray

from .bladerow import BladeRow, compute_gas_constants
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


def stator_calc(row: BladeRow, upstream: BladeRow, calculate_vm: bool = True) -> None:
    """Solve compressor stator exit conditions by matching exit massflow."""


    def calculate_vm_func(M_guess: float, apply: bool = False) -> float:
        """Solve stator for a guessed Mach; returns massflow residual."""
        T0_coolant_local = T0_coolant_weighted_average(row) if row.coolant is not None else 0.0
        T0_local = upstream.T0 - T0_coolant_local
        P0_local = row.P0_is - row.Yp * (upstream.P0 - upstream.P)

        deviation_func = getattr(row, "deviation_function", None)
        deviation = deviation_func(row, upstream) if callable(deviation_func) else 0.0

        M_local = np.full_like(row.area, M_guess, dtype=float)
        P0_P = IsenP(M_local, row.gamma)
        P_local = P0_local / P0_P
        T0_T = IsenT(M_local, row.gamma)
        T_local = T0_local / T0_T
        V_local = M_local * np.sqrt(row.gamma * row.R * T_local)
        Vm_local = V_local * np.cos(row.alpha2)
        Vx_local = Vm_local * np.cos(row.phi)
        Vr_local = Vm_local * np.sin(row.phi)
        Vt_local = Vm_local * np.tan(row.alpha2 + deviation)

        rho_local = P_local / (row.R * T_local)
        U_local = row.omega * row.r
        Wt_local = Vt_local - U_local
        alpha1_local = upstream.alpha2 + upstream.deviation if upstream.row_type == RowType.Rotor else row.alpha1
        entropy_rise_local = 0.5 * (row.Cp + upstream.Cp) * np.log(T_local / upstream.T) - row.R * np.log(P_local / upstream.P)

        # massflow integration (include blockage and optional coolant)
        total_area, streamline_area = compute_streamline_areas(row)
        massflow_local = np.zeros_like(row.massflow)
        massflow_fraction = np.linspace(0, 1, len(row.percent_hub_shroud))
        for j in range(1, len(row.percent_hub_shroud)):
            Vm_seg = 0.5 * (Vm_local[j] + Vm_local[j - 1])
            rho_seg = 0.5 * (rho_local[j] + rho_local[j - 1])
            massflow_local[j] = Vm_seg * rho_seg * streamline_area[j] * (1 - row.blockage) + massflow_local[j - 1]
        total_massflow_no_coolant = massflow_local[-1]
        if row.coolant is not None:
            massflow_local += massflow_fraction * row.coolant.massflow_percentage * total_massflow_no_coolant
        total_massflow_local = massflow_local[-1]

        if apply:
            row.T0 = T0_local
            row.P0 = P0_local
            row.M = M_local
            row.P = P_local
            row.T = T_local
            row.V = V_local
            row.Vm = Vm_local
            row.Vx = Vx_local
            row.Vr = Vr_local
            row.Vt = Vt_local
            row.alpha1 = alpha1_local
            row.beta1 = upstream.beta2
            row.deviation[:] = deviation
            row.rho = rho_local
            row.U = U_local
            row.Wt = Wt_local
            row.P0_stator_inlet = upstream.P0
            row.entropy_rise = entropy_rise_local
            row.total_area = total_area
            row.area = streamline_area
            row.massflow = massflow_local
            row.total_massflow_no_coolant = total_massflow_no_coolant
            row.total_massflow = total_massflow_local
        target_massflow = getattr(upstream, "total_massflow", total_massflow_local)
        return abs(target_massflow - total_massflow_local)

    if calculate_vm:
        res = minimize_scalar(calculate_vm_func,bounds=[0.01,1])
        calculate_vm_func(res.x, apply=True)
    else: # We know Vm, P0, T0, P 
        row.Vx = row.Vm*np.cos(row.phi)
        row.Vr = row.Vm*np.sin(row.phi)
        row.Vt = row.Vm*np.tan(row.alpha2)
        row.V = np.sqrt(row.Vx**2 + row.Vr**2 + row.Vt**2)
        row.T = row.P/(row.R*row.rho)   # We know P, this is a guess
        row.M = row.V/np.sqrt(row.gamma*row.R*row.T)

def rotor_calc(
    row: BladeRow,
    upstream: BladeRow,
    calculate_vm: bool = True,
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

        # compute massflow using locals (include blockage and optional coolant)
        rho_local = P_local / (row.R * T_local)
        total_area, streamline_area = compute_streamline_areas(row)
        massflow_local = np.zeros_like(row.massflow)
        massflow_fraction = np.linspace(0, 1, len(row.percent_hub_shroud))
        for j in range(1, len(row.percent_hub_shroud)):
            Vm_seg = 0.5 * (Vm_local[j] + Vm_local[j - 1])
            rho_seg = 0.5 * (rho_local[j] + rho_local[j - 1])
            massflow_local[j] = Vm_seg * rho_seg * streamline_area[j] * (1 - row.blockage) + massflow_local[j - 1]
        total_massflow_no_coolant = massflow_local[-1]
        if row.coolant is not None:
            massflow_local += massflow_fraction * row.coolant.massflow_percentage * total_massflow_no_coolant
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
            row.P0R = P0R_local
            row.T0R = T0R_local
            row.alpha2 = np.arctan2(row.Vt, row.Vm)
            row.M_rel = W_local / np.sqrt(row.gamma * row.R * T_local)
            row.total_massflow = total_massflow_local
            row.total_massflow_no_coolant = total_massflow_no_coolant
            row.massflow = massflow_local
            row.total_area = total_area
            row.area = streamline_area
            row.entropy_rise = 0.5 * (row.Cp + upstream.Cp) * np.log(T_local / upstream.T) - row.R * np.log(P_local / upstream.P)
    
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
    compute_gas_constants(row)
