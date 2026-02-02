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
from .outlet import OutletType

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
    loss_fn = getattr(row, "loss_function", None)
    loss_type = getattr(loss_fn, "loss_type", LossType.Pressure)
    target_eta_poly = None
    target_entropy = None

    if loss_type == LossType.Pressure and callable(loss_fn):
        row.Yp = loss_fn(row, upstream)  # type: ignore[arg-type]
    else:
        row.Yp[:] = 0
        if loss_type == LossType.Polytropic:
            if np.any(row.eta_poly):
                target_eta_poly = float(np.mean(row.eta_poly))
            elif callable(loss_fn):
                target_eta_poly = float(loss_fn(row, upstream))  # type: ignore[arg-type]
        elif loss_type == LossType.Entropy:
            if np.any(row.entropy_rise):
                target_entropy = float(np.mean(row.entropy_rise))
            elif callable(loss_fn):
                target_entropy = float(loss_fn(row, upstream))  # type: ignore[arg-type]

    def calculate_vm_func(M_guess: float, apply: bool = False) -> float:
        """Solve stator for a guessed Mach; returns massflow residual."""
        T0_coolant_local = T0_coolant_weighted_average(row) if row.coolant is not None else 0.0
        T0_local = upstream.T0 - T0_coolant_local
        
        P0_local = row.P0
        if row.row_type == RowType.IGV:
            row.P0 = row.P0_is - row.Yp * (upstream.P0 - upstream.P)
        else:
            row.P0_is = P0_local + row.Yp * (upstream.P0 - upstream.P)

        deviation_func = getattr(row, "deviation_function", None)
        deviation = deviation_func(row, upstream) if callable(deviation_func) else 0.0
        deviation_rad = np.radians(deviation)

        M_local = np.full_like(row.area, M_guess, dtype=float)
        P0_P = IsenP(M_local, row.gamma)
        P_local = P0_local / P0_P
        T0_T = IsenT(M_local, row.gamma)
        T_local = T0_local / T0_T
        V_local = M_local * np.sqrt(row.gamma * row.R * T_local)
        Vm_local = V_local * np.cos(row.alpha2)
        Vx_local = Vm_local * np.cos(row.phi)
        Vr_local = Vm_local * np.sin(row.phi)
        Vt_local = Vm_local * np.tan(row.alpha2 + deviation_rad)

        rho_local = P_local / (row.R * T_local)
        U_local = row.omega * row.r
        Wt_local = Vt_local - U_local
        alpha1_local = upstream.alpha2 + upstream.deviation if upstream.row_type == RowType.Rotor else row.alpha1
        entropy_rise_local = 0.5 * (row.Cp + upstream.Cp) * np.log(T_local / upstream.T) - row.R * np.log(P_local / upstream.P)

        # massflow integration (include blockage and optional coolant)
        total_area, streamline_area = compute_streamline_areas(row)
        n_streams = len(row.percent_hub_shroud)
        massflow_local = np.zeros(n_streams, dtype=float)
        massflow_fraction = np.array([1.0]) if n_streams <= 1 else np.linspace(0, 1, n_streams)
        if n_streams <= 1:
            massflow_local[0] = Vm_local[0] * rho_local[0] * streamline_area[0] * (1 - row.blockage)
            total_massflow_no_coolant = massflow_local[0]
            if row.coolant is not None:
                massflow_local += massflow_fraction * row.coolant.massflow_percentage * total_massflow_no_coolant
            total_massflow_local = massflow_local[-1]
        else:
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
            row.deviation[:] = deviation_rad
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
            # pi_local: stage total-pressure ratio (pt_out/pt_in); tau_local: total-temperature ratio (Tt_out/Tt_in)
            pi_local = float(np.mean(row.P0) / np.mean(upstream.P0)) if np.all(row.P0) else 1.0
            tau_local = float(np.mean(row.T0) / np.mean(upstream.T0)) if np.all(row.T0) else 1.0
            row.eta_poly = polytropic_efficiency(pi_local, tau_local, row.gamma)
            tau_is = (row.P0_is / upstream.P0) ** ((row.gamma - 1.0) / row.gamma)
            row.T0_is = upstream.T0 * tau_is
        target_massflow = getattr(upstream, "total_massflow", total_massflow_local)
        return abs(target_massflow - total_massflow_local)

    def solve_massflow_for_current_loss() -> None:
        res = minimize_scalar(calculate_vm_func, bounds=[0.01, 1], method="bounded")
        calculate_vm_func(res.x, apply=True)

    if calculate_vm:
        if loss_type == LossType.Polytropic and target_eta_poly is not None:
            def obj(y: float) -> float:
                row.Yp[:] = y
                res_local = minimize_scalar(calculate_vm_func, bounds=[0.01, 1], method="bounded")
                calculate_vm_func(res_local.x, apply=True)
                return abs(float(row.eta_poly) - target_eta_poly)

            res_y = minimize_scalar(obj, bounds=[0.0, 0.95], method="bounded")
            row.Yp[:] = res_y.x
            solve_massflow_for_current_loss()
        elif loss_type == LossType.Entropy and target_entropy is not None:
            def obj_entropy(y: float) -> float:
                row.Yp[:] = y
                res_local = minimize_scalar(calculate_vm_func, bounds=[0.01, 1], method="bounded")
                calculate_vm_func(res_local.x, apply=True)
                return abs(float(np.mean(row.entropy_rise)) - target_entropy)

            res_y = minimize_scalar(obj_entropy, bounds=[0.0, 0.95], method="bounded")
            row.Yp[:] = res_y.x
            solve_massflow_for_current_loss()
        else:
            solve_massflow_for_current_loss()
    else:  # We know Vm, P0, T0, P
        row.Vx = row.Vm * np.cos(row.phi)
        row.Vr = row.Vm * np.sin(row.phi)
        row.Vt = row.Vm * np.tan(row.alpha2)
        row.V = np.sqrt(row.Vx ** 2 + row.Vr ** 2 + row.Vt ** 2)
        row.T = row.P / (row.R * row.rho)   # We know P, this is a guess
        row.M = row.V / np.sqrt(row.gamma * row.R * row.T)

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
    """
    loss_fn = getattr(row, "loss_function", None)
    loss_type = getattr(loss_fn, "loss_type", LossType.Pressure)
    target_eta_poly = None
    target_entropy = None

    if loss_type == LossType.Pressure and callable(loss_fn):
        row.Yp = loss_fn(row, upstream)  # type: ignore[arg-type]
    else:
        row.Yp[:] = 0
        if loss_type == LossType.Polytropic:
            if np.any(row.eta_poly):
                target_eta_poly = float(np.mean(row.eta_poly))
            elif callable(loss_fn):
                target_eta_poly = float(loss_fn(row, upstream))  # type: ignore[arg-type]
        elif loss_type == LossType.Entropy:
            if np.any(row.entropy_rise):
                target_entropy = float(np.mean(row.entropy_rise))
            elif callable(loss_fn):
                target_entropy = float(loss_fn(row, upstream))  # type: ignore[arg-type]

    # Use the frozen target (if available) so diagnostic code can overwrite row.P0_ratio
    # without changing the initial guess used by this solver.
    P0_ratio_target = getattr(row, "P0_ratio_target", 0.0) or row.P0_ratio
    row.P0 = upstream.P0 * P0_ratio_target
    row.P0_is = row.P0 + row.Yp * (upstream.P0 - upstream.P)
    
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
        deviation_func = getattr(row, "deviation_function", None)
        deviation_val = deviation_func(row, upstream) if callable(deviation_func) else 0.0
        deviation_rad = np.radians(deviation_val)

        P0R_local = upstream.P0R - row.Yp * (upstream.P0R - upstream.P)
        # Use rothalpy conservation: I = Cp*T0R - U^2/2 = const across rotor
        U_local = row.omega * row.r
        T0R_local = (upstream_rothalpy + 0.5 * U_local ** 2) / row.Cp

        P_local = P0R_local / IsenP(M_rel, row.gamma)

        P0R_P = P0R_local / P_local
        T0R_T = P0R_P ** ((row.gamma - 1) / row.gamma)
        T_local = T0R_local / T0R_T
        W_local = np.sqrt(2 * row.Cp * (T0R_local - T_local))

        if np.isnan(W_local).any() or np.any(T_local >= T0R_local):
            return np.inf

        Vr_local = W_local * np.sin(row.phi)
        beta2_eff = row.beta2 + deviation_rad
        Vm_local = W_local * np.cos(beta2_eff)
        Wt_local = W_local * np.sin(beta2_eff)
        Vx_local = Vm_local * np.cos(row.phi)
        Vt_local = Wt_local + U_local
        V_local = np.sqrt(Vr_local ** 2 + Vt_local ** 2 + Vx_local ** 2)
        M_local = V_local / np.sqrt(row.gamma * row.R * T_local)
        T0_local = T_local + V_local ** 2 / (2 * row.Cp)
        P0_local = P_local * (T0_local / T_local) ** (row.gamma / (row.gamma - 1))

        # compute massflow using locals (include blockage and optional coolant)
        rho_local = P_local / (row.R * T_local)
        total_area, streamline_area = compute_streamline_areas(row)
        n_streams = len(row.percent_hub_shroud)
        massflow_local = np.zeros(n_streams, dtype=float)
        massflow_fraction = np.array([1.0]) if n_streams <= 1 else np.linspace(0, 1, n_streams)
        if n_streams <= 1:
            massflow_local[0] = Vm_local[0] * rho_local[0] * streamline_area[0] * (1 - row.blockage)
            total_massflow_no_coolant = massflow_local[0]
            if row.coolant is not None:
                massflow_local += massflow_fraction * row.coolant.massflow_percentage * total_massflow_no_coolant
            total_massflow_local = massflow_local[-1]
        else:
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
            row.U = U_local
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
            row.deviation[:] = deviation_rad
            # pi_local: stage total-pressure ratio (pt_out/pt_in); tau_local: total-temperature ratio (Tt_out/Tt_in)
            pi_local = float(np.mean(row.P0) / np.mean(upstream.P0)) if np.all(row.P0) else 1.0
            tau_local = float(np.mean(row.T0) / np.mean(upstream.T0)) if np.all(row.T0) else 1.0
            row.eta_poly = polytropic_efficiency(pi_local, tau_local, row.gamma)
            tau_is = (row.P0_is / upstream.P0) ** ((row.gamma - 1.0) / row.gamma)
            row.T0_is = upstream.T0 * tau_is
    
        return np.abs(upstream.total_massflow - total_massflow_local)
    
    def solve_massflow_for_current_loss() -> None:
        res = minimize_scalar(calculate_vm_func, bounds=[0.01, 1], method="bounded")
        calculate_vm_func(res.x, apply=True)

    if calculate_vm:
        if loss_type == LossType.Polytropic and target_eta_poly is not None:
            def obj(y: float) -> float:
                row.Yp[:] = y
                res_local = minimize_scalar(calculate_vm_func, bounds=[0.01, 1], method="bounded")
                calculate_vm_func(res_local.x, apply=True)
                return abs(float(row.eta_poly) - target_eta_poly)

            res_y = minimize_scalar(obj, bounds=[0.0, 0.95], method="bounded")
            row.Yp[:] = res_y.x
            solve_massflow_for_current_loss()
        elif loss_type == LossType.Entropy and target_entropy is not None:
            def obj_entropy(y: float) -> float:
                row.Yp[:] = y
                res_local = minimize_scalar(calculate_vm_func, bounds=[0.01, 1], method="bounded")
                calculate_vm_func(res_local.x, apply=True)
                return abs(float(np.mean(row.entropy_rise)) - target_entropy)

            res_y = minimize_scalar(obj_entropy, bounds=[0.0, 0.95], method="bounded")
            row.Yp[:] = res_y.x
            solve_massflow_for_current_loss()
        else:
            solve_massflow_for_current_loss()
    else: # We know Vm from radeq, beta2 from blade angle, T0R from rothalpy
        deviation_func = getattr(row, "deviation_function", None)
        deviation_val = deviation_func(row, upstream) if callable(deviation_func) else 0.0
        deviation_rad = np.radians(deviation_val)
        beta2_eff = row.beta2 + deviation_rad

        row.U = row.omega * row.r
        row.T0R = (upstream_rothalpy + 0.5 * row.U ** 2) / row.Cp
        row.P0R = upstream.P0R - row.Yp * (upstream.P0R - upstream.P)

        row.Vr = row.Vm * np.sin(row.phi)
        row.Vx = row.Vm * np.cos(row.phi)

        # Compute W from velocity triangle (geometric closure)
        row.W = row.Vm / np.cos(beta2_eff)
        row.Wt = row.W * np.sin(beta2_eff)
        row.Vt = row.Wt + row.U

        row.alpha2 = np.arctan2(row.Vt, row.Vm)
        row.V = np.sqrt(row.Vm ** 2 * (1 + np.tan(row.alpha2) ** 2))

        # Update T from energy conservation in rotating frame
        row.T = row.T0R - row.W ** 2 / (2 * row.Cp)

        row.M = row.V / np.sqrt(row.gamma * row.R * row.T)

    # Compute T0 first, then derive P0 from P0R to keep the velocity triangle consistent.
    # P0/P0R = (T0/T0R)^(gamma/(gamma-1)) always holds since both reference the same static state.
    row.M_rel = row.W / np.sqrt(row.gamma * row.R * row.T)
    row.T0 = row.T + row.V ** 2 / (2 * row.Cp)
    row.P0 = row.P0R * (row.T0 / row.T0R) ** (row.gamma / (row.gamma - 1))
    compute_gas_constants(row)
