"""OTAC compressor loss model translations (best effort).

This module ports the OTAC ``*.int`` compressor loss models into the current
Python API. Geometry/flow mappings assume ``upstream`` is ``FL_IR`` and ``row``
is ``FL_OR``. Many correlations require design parameters (e.g., blade counts,
clearances); these are exposed as constructor arguments with pragmatic defaults.
"""

from __future__ import annotations

import warnings
import numpy as np
import numpy.typing as npt

from ...bladerow import BladeRow
from ...enums import LossType
from ..losstype import LossBaseClass


def _mean(val, default: float = 0.0) -> float:
    try:
        arr = np.asarray(val)
        return float(np.mean(arr)) if arr.size else default
    except Exception:
        return default


def _mag(val) -> float:
    arr = np.asarray(val)
    return float(np.linalg.norm(arr))


def _span(row: BladeRow) -> float:
    r = np.asarray(row.r)
    return float(np.max(r) - np.min(r))


def _hub(row: BladeRow) -> float:
    return float(np.min(np.asarray(row.r)))


def _tip(row: BladeRow) -> float:
    return float(np.max(np.asarray(row.r)))


# ---------------------------------------------------------------------------
# High-level placeholders for axial/turbine/NASA correlations
# ---------------------------------------------------------------------------


class AxialCompressorAungier(LossBaseClass):
    """Aungier axial-compressor pressure-loss (omega)."""

    def __init__(self):
        super().__init__(LossType.Pressure)

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        Vm1 = _mean(getattr(upstream, "Vm", upstream.V))
        Vm2 = _mean(getattr(row, "Vm", row.V))
        Vt1 = _mean(getattr(upstream, "Vt", upstream.V))
        Vt2 = _mean(getattr(row, "Vt", row.V))
        M2 = _mean(getattr(row, "M", 0.0))

        Vm1 = max(Vm1, 1e-6)
        turning = abs(np.arctan2(Vt2, Vm2) - np.arctan2(Vt1, Vm1))

        Df = abs(Vt2 - Vt1) / Vm1 + max(0.0, 1 - Vm2 / Vm1)
        Deq = Df  # stand-in for equivalent diffusion factor

        omega_prof = 0.04 + 0.6 * Df**2 + 0.2 * max(0.0, Deq - 0.5)
        omega_sec = 0.01 * (turning**2)
        tip_clearance = float(getattr(row, "tip_clearance", 0.0))
        omega_tip = 0.02 * tip_clearance
        omega_shock = 0.0
        if M2 > 0.9:
            omega_shock = 0.02 * (M2 - 0.9) ** 2

        omega = omega_prof + omega_sec + omega_tip + omega_shock
        omega = np.maximum(omega, 0.0)
        return np.full_like(row.r, omega, dtype=float)


class AxialCompressorEntropy(LossBaseClass):
    """Koch-Smith axial-compressor loss (entropy-based)."""

    def __init__(self):
        super().__init__(LossType.Entropy)

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        Vm1 = _mean(getattr(upstream, "Vm", upstream.V))
        Vm2 = _mean(getattr(row, "Vm", row.V))
        Vt1 = _mean(getattr(upstream, "Vt", upstream.V))
        Vt2 = _mean(getattr(row, "Vt", row.V))
        M2 = _mean(getattr(row, "M", 0.0))

        Vm1 = max(Vm1, 1e-6)
        turning = abs(np.arctan2(Vt2, Vm2) - np.arctan2(Vt1, Vm1))
        Df = abs(Vt2 - Vt1) / Vm1 + max(0.0, 1 - Vm2 / Vm1)

        ds_over_cp = 0.03 + 0.5 * Df**2 + 0.1 * (turning**2)
        if M2 > 1.0:
            ds_over_cp += 0.05 * (M2 - 1.0) ** 2

        return np.full_like(row.r, ds_over_cp, dtype=float)


class AxialCompressorWrightMiller(LossBaseClass):
    """Wright-Miller axial compressor loss (pressure coefficient)."""

    def __init__(self):
        super().__init__(LossType.Pressure)

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        beta1 = _mean(np.degrees(getattr(row, "beta1", upstream.beta2)))
        beta2 = _mean(np.degrees(getattr(row, "beta2", row.beta2)))
        turning = abs(beta2 - beta1)
        V_ratio = _mean(getattr(row, "V", row.V)) / max(_mean(getattr(upstream, "V", upstream.V)), 1e-6)
        omega = abs(np.tan(np.radians(beta1)) - np.tan(np.radians(beta2)))
        omega *= (0.5 + 0.5 * V_ratio) * (1 + turning / 90.0)
        omega = np.maximum(omega, 0.0)
        return np.full_like(row.r, omega, dtype=float)


class AxialTurbineAinleyMathiesonOTAC(LossBaseClass):
    """Placeholder for OTAC axial turbine Ainley-Mathieson."""

    def __init__(self):
        super().__init__(LossType.Pressure)

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        warnings.warn("AxialTurbineAinleyMathiesonOTAC not yet translated; returning zeros.", stacklevel=2)
        return np.zeros_like(row.r)


class AxialTurbineKackerOkapuuOTAC(LossBaseClass):
    """Placeholder for OTAC axial turbine Kacker-Okapuu."""

    def __init__(self):
        super().__init__(LossType.Pressure)

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        warnings.warn("AxialTurbineKackerOkapuuOTAC not yet translated; returning zeros.", stacklevel=2)
        return np.zeros_like(row.r)


class NASA23B20(LossBaseClass):
    """NASA23B/20 compressor loss correlation (simplified profile + incidence)."""

    def __init__(self):
        super().__init__(LossType.Pressure)

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        beta1 = _mean(np.degrees(upstream.beta2))
        beta2 = _mean(np.degrees(row.beta2))
        turning = abs(beta2 - beta1)
        Vm1 = _mean(getattr(upstream, "Vm", upstream.V))
        Vm2 = _mean(getattr(row, "Vm", row.V))
        Vm1 = max(Vm1, 1e-6)
        DF = abs(_mean(getattr(row, "Vt", row.V)) - _mean(getattr(upstream, "Vt", upstream.V))) / Vm1 + max(
            0.0, 1 - Vm2 / Vm1
        )
        omega_min = 0.06 + 0.4 * DF**2
        i_ml = 1.5
        incidence = abs(beta1 - i_ml)
        omega = omega_min * (1 + (incidence / 15) ** 2) * (1 + turning / 120)
        return np.full_like(row.r, omega, dtype=float)


class NASA74A(LossBaseClass):
    """NASA 74A 5-stage axial compressor loss (simplified)."""

    def __init__(self):
        super().__init__(LossType.Pressure)

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        Vm1 = _mean(getattr(upstream, "Vm", upstream.V))
        Vm2 = _mean(getattr(row, "Vm", row.V))
        Vm1 = max(Vm1, 1e-6)
        Vt1 = _mean(getattr(upstream, "Vt", upstream.V))
        Vt2 = _mean(getattr(row, "Vt", row.V))
        DF = abs(Vt2 - Vt1) / Vm1 + max(0.0, 1 - Vm2 / Vm1)
        if len(np.asarray(row.r)) > 1:
            pct_rad = np.linspace(0, 1, len(row.r))
        else:
            pct_rad = np.array([0.5])
        base = 0.02 + 0.08 * DF**2
        radial_factor = 1 + 0.3 * pct_rad
        omega = base * radial_factor
        return np.array(omega, dtype=float)


class RadialInput(LossBaseClass):
    """Radial input helper (no additional loss)."""

    def __init__(self):
        super().__init__(LossType.Pressure)

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        return np.zeros_like(row.r)


class DiffuserVanelessStanitz(LossBaseClass):
    """Stanitz vaneless diffuser loss (pressure)."""

    def __init__(self):
        super().__init__(LossType.Enthalpy)

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        r_in = _tip(upstream)
        r_out = _tip(row)
        length = max(r_out - r_in, 1e-6)
        Dh = 2 * _span(row) if _span(row) > 0 else 1.0
        Re = abs(_mean(row.V)) * Dh * _mean(getattr(row, "rho", 1.0)) / max(getattr(row, "mu", 1.0), 1e-6)
        Cf = 0.026 / (Re ** 0.2) if Re > 0 else 0.0
        loss = Cf * length / Dh * (_mean(row.V) ** 2) / 2
        loss = np.maximum(loss, 0.0)
        return np.full_like(row.r, loss, dtype=float)


# ---------------------------------------------------------------------------
# Impeller correlations (translated)
# ---------------------------------------------------------------------------


class ImpellerBladeLoadingAungier(LossBaseClass):
    """Aungier blade loading loss."""

    def __init__(self, number_of_blades: int = 12, splitter_le: float = 0.0, loading_coefficient: float = 1.0):
        super().__init__(LossType.Enthalpy)
        self.number_of_blades = number_of_blades
        self.splitter_le = splitter_le
        self.loading_coefficient = loading_coefficient

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        mean_blades = self.number_of_blades / 2 + (self.number_of_blades / 2) * (1 - self.splitter_le)
        radius_exit = _tip(row)
        radius_inlet = _hub(upstream)
        lb = float(getattr(row, "chord", radius_exit - radius_inlet))

        dVrel = 2 * np.pi * (2 * radius_exit) * _mean(row.U) * self.loading_coefficient / (mean_blades * lb)
        dh_blade = dVrel**2 / 48.0

        kbar = (_mean(row.phi) - _mean(upstream.phi)) / max(lb, 1e-6)
        bbar = ((radius_exit - radius_inlet) + _span(row)) / 2.0
        Vrelbar = (_mag(getattr(upstream, "W", upstream.V)) + _mag(getattr(row, "W", row.V))) / 2.0
        dh_hub_shroud = (kbar * bbar * Vrelbar) ** 2 / 12.0

        loss = dh_blade + dh_hub_shroud
        loss = np.maximum(loss, 0.0)
        return np.full_like(row.r, loss, dtype=float)


class ImpellerBladeLoadingCoppage(LossBaseClass):
    """Coppage blade loading loss."""

    def __init__(
        self,
        number_of_blades: int = 12,
        splitter_le: float = 0.0,
        loading_coefficient: float = 1.0,
        surge_vrel: float = 1.0,
    ):
        super().__init__(LossType.Enthalpy)
        self.number_of_blades = number_of_blades
        self.splitter_le = splitter_le
        self.loading_coefficient = loading_coefficient
        self.surge_vrel = surge_vrel

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        Vrel_out = _mag(getattr(row, "W", row.V))
        Kbl = 0.75 if self.splitter_le == 0 else 0.6
        r_tip_in = _tip(upstream)
        r_exit = _tip(row)
        denom = (self.surge_vrel / max(Vrel_out, 1e-6)) * (
            (self.number_of_blades / np.pi) * (1 - r_tip_in / max(r_exit, 1e-6)) + 2 * r_tip_in / max(r_exit, 1e-6)
        )
        Df = 1 - Vrel_out / max(self.surge_vrel, 1e-6) + Kbl * self.loading_coefficient / max(denom, 1e-6)
        dh = 0.05 * Df**2 * (_mean(row.U) ** 2)
        dh = np.maximum(dh, 0.0)
        return np.full_like(row.r, dh, dtype=float)


class ImpellerClearanceJansen(LossBaseClass):
    """Jansen clearance loss."""

    def __init__(
        self,
        number_of_blades: int = 12,
        tip_clearance_axial: float = 0.0,
        exit_blade_height: float | None = None,
        loss_modifier: float = 1.0,
    ):
        super().__init__(LossType.Enthalpy)
        self.number_of_blades = number_of_blades
        self.tip_clearance_axial = tip_clearance_axial
        self.exit_blade_height = exit_blade_height
        self.loss_modifier = loss_modifier

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        h_exit = self.exit_blade_height if self.exit_blade_height is not None else _span(row)
        Vtheta = _mean(getattr(row, "Vt", row.V))
        Vm = _mean(getattr(row, "Vm", row.V))
        r_tip_in = _tip(upstream)
        r_hub_in = _hub(upstream)
        r_exit = _tip(row)
        ratio = abs(r_exit - r_tip_in)
        rho_ratio = _mean(getattr(row, "rho", 1.0)) / max(_mean(getattr(upstream, "rho", 1.0)), 1e-6)
        inner = (4 * np.pi) / (h_exit * self.number_of_blades) * (
            (r_tip_in**2 - r_hub_in**2) / (ratio * (1 + rho_ratio))
        ) * Vtheta * Vm
        dh = 0.6 * (self.tip_clearance_axial / max(h_exit, 1e-6)) * Vtheta * np.sqrt(max(inner, 0.0))
        dh *= self.loss_modifier
        dh = np.maximum(dh, 0.0)
        return np.full_like(row.r, dh, dtype=float)


class ImpellerDiscFrictionDaily(LossBaseClass):
    """Daily & Nece disc friction loss."""

    def __init__(self, bf_gap: float | None = None, loss_modifier: float = 1.0):
        super().__init__(LossType.Enthalpy)
        self.bf_gap = bf_gap
        self.loss_modifier = loss_modifier

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        radius_exit = _tip(row)
        span = _span(row)
        gap = self.bf_gap if self.bf_gap is not None else float(getattr(row, "tip_clearance", 0.0) * span)
        gap = gap if gap > 0 else 1e-5

        rho_avg = _mean([getattr(upstream, "rho", 0.0), getattr(row, "rho", 0.0)], 1.0)
        mu = float(getattr(row, "mu", 0.0) or getattr(upstream, "mu", 0.0) or 1.0)

        U = np.asarray(getattr(row, "U", row.omega * row.r))
        W_in = _mag(getattr(upstream, "W", upstream.V))
        if W_in == 0:
            W_in = _mag(getattr(upstream, "V", 1e-6))

        Re = np.abs(U) * radius_exit * rho_avg / mu
        Re = np.maximum(Re, 1e-6)

        f_df = np.where(
            Re < 3e5,
            3.7 * (gap / radius_exit) ** 0.1 / np.sqrt(Re),
            0.102 * (gap / radius_exit) ** 0.1 / (Re ** 0.2),
        )

        dh_disc = f_df * rho_avg * (radius_exit**2) * (U**3) / (4 * W_in)
        loss = self.loss_modifier * dh_disc

        cp = _mean([row.Cp, upstream.Cp], row.Cp)
        ht_delta = cp * (np.asarray(row.T0) - np.asarray(upstream.T0))
        cap = 0.25 * ht_delta
        loss = np.minimum(loss, cap)
        loss = np.maximum(loss, 0.0)

        return np.array(loss, dtype=float)


class ImpellerIncidenceAungier(LossBaseClass):
    """Aungier incidence loss."""

    def __init__(self, blade_inlet_angle_deg: float | None = None, loss_modifier: float = 1.0):
        super().__init__(LossType.Enthalpy)
        self.blade_inlet_angle_deg = blade_inlet_angle_deg
        self.loss_modifier = loss_modifier

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        r_tip = _tip(upstream)
        r_hub = _hub(upstream)
        omega = float(getattr(upstream, "omega", 0.0))
        V = _mag(upstream.V)
        Vm = _mean(getattr(upstream, "Vm", upstream.V))
        alpha = _mean(np.degrees(upstream.alpha1))

        U_tip = omega * r_tip
        U_hub = omega * r_hub
        Vtheta_tip_rel = U_tip - _mean(getattr(upstream, "Vt", 0.0))
        Vtheta_hub_rel = U_hub - _mean(getattr(upstream, "Vt", 0.0))
        Vrel_tip = np.hypot(V * np.cos(np.radians(alpha)), Vtheta_tip_rel)
        Vrel_hub = np.hypot(V * np.cos(np.radians(alpha)), Vtheta_hub_rel)

        blade_beta = self.blade_inlet_angle_deg
        if blade_beta is None:
            beta_metal = getattr(upstream, "beta1_metal", None)
            blade_beta = _mean(beta_metal if beta_metal is not None else np.degrees(upstream.beta1))

        target = abs(Vm / max(np.cos(np.radians(blade_beta)), 1e-6))

        dh_hub = 0.4 * (Vrel_hub - target) ** 2
        dh_tip = 0.4 * (Vrel_tip - target) ** 2
        dh_mean = 0.4 * (_mag(getattr(upstream, "W", upstream.V)) - target) ** 2
        dh_inc = (dh_hub + dh_tip + 10 * dh_mean) / 12.0
        dh_inc = np.maximum(dh_inc * self.loss_modifier, 0.0)
        return np.full_like(row.r, dh_inc, dtype=float)


class ImpellerIncidenceConrad(LossBaseClass):
    """Conrad incidence loss."""

    def __init__(
        self,
        leading_edge_thickness: float | None = None,
        number_of_blades: int = 12,
        blade_inlet_angle_deg: float | None = None,
        radius_tip_inlet: float | None = None,
        radius_hub_inlet: float | None = None,
        f_incidence: float = 1.0,
        loss_modifier_inc: float = 1.0,
    ):
        super().__init__(LossType.Enthalpy)
        self.leading_edge_thickness = leading_edge_thickness
        self.number_of_blades = number_of_blades
        self.blade_inlet_angle_deg = blade_inlet_angle_deg
        self.radius_tip_inlet = radius_tip_inlet
        self.radius_hub_inlet = radius_hub_inlet
        self.f_incidence = f_incidence
        self.loss_modifier_inc = loss_modifier_inc

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        r_tip = self.radius_tip_inlet if self.radius_tip_inlet is not None else _tip(upstream)
        r_hub = self.radius_hub_inlet if self.radius_hub_inlet is not None else _hub(upstream)
        le_thickness = self.leading_edge_thickness
        if le_thickness is None:
            chord = float(getattr(row, "chord", 0.0) or 1.0)
            le_thickness = 0.02 * chord

        blade_inlet_angle = self.blade_inlet_angle_deg
        if blade_inlet_angle is None:
            beta_metal = getattr(upstream, "beta1_metal", None)
            blade_inlet_angle = _mean(beta_metal if beta_metal is not None else np.degrees(upstream.beta1))

        area = float(np.pi * (r_tip**2 - r_hub**2))
        beta_opt_rad = np.arctan(
            area
            / (area - le_thickness * (r_tip - r_hub) * self.number_of_blades / 2.0)
            * np.tan(np.radians(blade_inlet_angle))
        )
        beta_opt_deg = np.degrees(beta_opt_rad)

        beta_in_deg = _mean(np.degrees(upstream.beta1))
        Vrel = _mag(getattr(upstream, "W", upstream.V))
        if Vrel == 0:
            Vrel = _mag(upstream.V)

        Wui = Vrel * np.sin(np.radians(abs(beta_in_deg - beta_opt_deg)))
        dh_incidence = self.f_incidence * 0.5 * Wui**2

        loss = self.loss_modifier_inc * dh_incidence
        loss = np.maximum(loss, 0.0)
        return np.full_like(row.r, loss, dtype=float)


class ImpellerLeakageAungier(LossBaseClass):
    """Aungier tip-leakage loss."""

    def __init__(
        self,
        number_of_blades: int = 12,
        splitter_le: float = 0.0,
        seal_clearance: float | None = None,
        loading_coefficient: float = 1.0,
        loss_modifier: float = 1.0,
    ):
        super().__init__(LossType.Enthalpy)
        self.number_of_blades = number_of_blades
        self.splitter_le = splitter_le
        self.seal_clearance = seal_clearance
        self.loading_coefficient = loading_coefficient
        self.loss_modifier = loss_modifier

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        mean_blades = self.number_of_blades / 2 + (self.number_of_blades / 2) * (1 - self.splitter_le)
        r_exit = _tip(row)
        r_in = _hub(upstream)
        r_tip_in = _tip(upstream)
        bwidth = _span(row)
        rbar = (r_in + r_exit) / 2.0
        bbar = ((r_tip_in - r_in) + bwidth) / 2.0
        lb = float(getattr(row, "chord", r_exit - r_in))

        Vtheta_out = _mean(getattr(row, "Vt", row.V))
        Vtheta_in = _mean(getattr(upstream, "Vt", upstream.V))
        rho_out = _mean(getattr(row, "rho", 1.0))
        deltaP = rho_out * abs(r_exit * Vtheta_out - r_in * Vtheta_in) / max(lb, 1e-6)
        seal_gap = self.seal_clearance if self.seal_clearance is not None else float(getattr(row, "tip_clearance", 0.0) * bwidth)

        Ucl = 0.816 * np.sqrt(max(2 * deltaP / max(rho_out, 1e-6), 0.0))
        mdot_cl = rho_out * mean_blades * seal_gap * lb * Ucl
        W_out = _mag(getattr(row, "W", row.V))
        U_out = _mean(row.U)
        dh_leakage = (mdot_cl * Ucl) / max(2 * W_out * max(U_out, 1e-6), 1e-6) * U_out**2

        loss = self.loss_modifier * dh_leakage
        loss = np.maximum(loss, 0.0)
        return np.full_like(row.r, loss, dtype=float)


class ImpellerMixingAungier(LossBaseClass):
    """Aungier & Dean mixing loss."""

    def __init__(
        self,
        number_of_blades: int = 12,
        splitter_le: float = 0.0,
        loading_coefficient: float = 1.0,
        area_exit_factor: float = 0.9,
        loss_modifier: float = 1.0,
    ):
        super().__init__(LossType.Enthalpy)
        self.number_of_blades = number_of_blades
        self.splitter_le = splitter_le
        self.loading_coefficient = loading_coefficient
        self.area_exit_factor = area_exit_factor
        self.loss_modifier = loss_modifier

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        mean_blades = self.number_of_blades / 2 + (self.number_of_blades / 2) * (1 - self.splitter_le)
        r_exit = _tip(row)
        lb = float(getattr(row, "chord", r_exit - _hub(upstream)))
        dVrel = 2 * np.pi * (2 * r_exit) * _mean(row.U) * self.loading_coefficient / (mean_blades * lb)
        Vrel_max = (_mag(getattr(upstream, "W", upstream.V)) + _mag(getattr(row, "W", row.V)) + dVrel) / 2
        Vrel_out = _mag(getattr(row, "W", row.V))
        Deq = Vrel_max / max(Vrel_out, 1e-6)

        if Deq <= 2:
            Vrel_sep = Vrel_out
        else:
            Vrel_sep = Vrel_out * Deq / 2

        Vm = _mean(getattr(row, "Vm", row.V))
        Vtheta_rel = _mean(getattr(row, "Wt", row.Vt if hasattr(row, "Vt") else row.V))
        area_exit = float(np.pi * (r_exit**2 - _hub(row) ** 2))
        Vrel_out_eff = np.hypot(Vm * area_exit * self.area_exit_factor / max(np.pi * _span(row) * r_exit * 2, 1e-6), Vtheta_rel)

        dh_mix = 0.5 * (Vrel_sep - Vrel_out_eff) ** 2
        cp = _mean([row.Cp, upstream.Cp], row.Cp)
        cap = 0.3 * cp * (_mean(row.T0) - _mean(upstream.T0))
        dh_mix = np.clip(dh_mix, 0.0, cap)

        loss = self.loss_modifier * dh_mix
        return np.full_like(row.r, loss, dtype=float)


class ImpellerMixingJohnston(LossBaseClass):
    """Johnston & Dean mixing loss."""

    def __init__(self, bstar: float = 1.0, loss_modifier: float = 1.0):
        super().__init__(LossType.Enthalpy)
        self.bstar = bstar
        self.loss_modifier = loss_modifier

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        omega = float(getattr(upstream, "omega", 0.0))
        r_tip = _tip(upstream)
        Vtheta = _mean(getattr(upstream, "Vt", upstream.V))
        V = _mag(upstream.V)
        alpha = _mean(np.degrees(upstream.alpha1))
        Vtheta_tip_rel = omega * r_tip - Vtheta
        Vrel_tip = np.hypot(V * np.cos(np.radians(alpha)), Vtheta_tip_rel)
        Vrel_out = _mag(getattr(row, "W", row.V))
        wake = 1 - (1 / 0.45) * (Vrel_out / max(Vrel_tip, 1e-6))
        wake = np.clip(wake, 0.0, 0.99)

        Vtheta_out = _mean(getattr(row, "Vt", row.V))
        Vm_out = _mean(getattr(row, "Vm", row.V))
        coeff = 1.0 / (1 + (Vtheta_out / max(Vm_out, 1e-6)) ** 2)
        dh = coeff * ((1 - wake - self.bstar) / max(1 - wake, 1e-6)) ** 2 * (_mag(row.V) ** 2 / 2)
        dh = np.maximum(dh, 0.0) * self.loss_modifier
        return np.full_like(row.r, dh, dtype=float)


class ImpellerPrescribed(LossBaseClass):
    """Conrad incidence loss (prescribed percentage of enthalpy rise)."""

    def __init__(self, loss_pct: float = 0.0):
        super().__init__(LossType.Enthalpy)
        self.loss_pct = loss_pct

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        cp = _mean([row.Cp, upstream.Cp], row.Cp)
        ht_in = cp * np.asarray(upstream.T0)
        ht_out = cp * np.asarray(row.T0)
        dh = ht_out - ht_in
        loss = self.loss_pct * dh
        loss = np.maximum(loss, 0)  # clamp negative due to numerical noise
        return np.array(loss, dtype=float)


class ImpellerRecirculationAungier(LossBaseClass):
    """Aungier recirculation loss."""

    def __init__(
        self,
        number_of_blades: int = 12,
        splitter_le: float = 0.0,
        loading_coefficient: float = 1.0,
        lb: float | None = None,
        loss_modifier: float = 1.0,
    ):
        super().__init__(LossType.Enthalpy)
        self.number_of_blades = number_of_blades
        self.splitter_le = splitter_le
        self.loading_coefficient = loading_coefficient
        self.lb = lb
        self.loss_modifier = loss_modifier

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        mean_blades = self.number_of_blades / 2 + (self.number_of_blades / 2) * (1 - self.splitter_le)
        r_exit = _tip(row)
        lb = self.lb if self.lb is not None else float(getattr(row, "chord", r_exit - _hub(upstream)))
        dVrel = 2 * np.pi * (2 * r_exit) * _mean(row.U) * self.loading_coefficient / (mean_blades * lb)
        Vrel_max = (_mag(getattr(upstream, "W", upstream.V)) + _mag(getattr(row, "W", row.V)) + dVrel) / 2
        Deq = Vrel_max / max(_mag(getattr(row, "W", row.V)), 1e-6)

        Vtheta_rel = _mean(getattr(row, "Wt", row.Vt if hasattr(row, "Vt") else row.V))
        Vm = _mean(getattr(row, "Vm", row.V))
        beta = np.degrees(np.arctan2(Vtheta_rel, max(Vm, 1e-6)))
        dh = (Deq / 2 - 1) * (abs(Vtheta_rel) / max(Vm, 1e-6) - 2 * (1 / max(np.tan(np.radians(beta)), 1e-6))) * (
            _mean(row.U) ** 2
        )
        cp = _mean([row.Cp, upstream.Cp], row.Cp)
        cap = 0.25 * cp * (_mean(row.T0) - _mean(upstream.T0))
        dh = np.clip(dh, 0.0, cap) * self.loss_modifier
        return np.full_like(row.r, dh, dtype=float)


class ImpellerRecirculationOh(LossBaseClass):
    """Oh recirculation loss."""

    def __init__(
        self,
        splitter_le: float = 0.0,
        loading_coefficient: float = 1.0,
        number_of_blades: int = 12,
        radius_tip_inlet: float | None = None,
        radius_exit: float | None = None,
        surge_vrel: float = 1.0,
        loss_modifier: float = 1.0,
    ):
        super().__init__(LossType.Enthalpy)
        self.splitter_le = splitter_le
        self.loading_coefficient = loading_coefficient
        self.number_of_blades = number_of_blades
        self.radius_tip_inlet = radius_tip_inlet
        self.radius_exit = radius_exit
        self.surge_vrel = surge_vrel
        self.loss_modifier = loss_modifier

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        r_tip_in = self.radius_tip_inlet if self.radius_tip_inlet is not None else _tip(upstream)
        r_exit = self.radius_exit if self.radius_exit is not None else _tip(row)
        omega = float(getattr(upstream, "omega", 0.0))
        Vtheta = _mean(getattr(upstream, "Vt", upstream.V))
        alpha = _mean(np.degrees(upstream.alpha1))
        V = _mag(upstream.V)

        U_tip = omega * r_tip_in
        Vtheta_tip_rel = U_tip - Vtheta
        Vrel_tip = np.hypot(V * np.cos(np.radians(alpha)), Vtheta_tip_rel)

        Kbl = 0.75 if self.splitter_le == 0 else 0.6
        Df = 1 - _mag(getattr(row, "W", row.V)) / max(self.surge_vrel, 1e-6) + Kbl * self.loading_coefficient / max(
            (self.surge_vrel / max(_mag(getattr(row, "W", row.V)), 1e-6))
            * ((self.number_of_blades / np.pi) * (1 - r_tip_in / max(r_exit, 1e-6)) + 2 * r_tip_in / max(r_exit, 1e-6)),
            1e-6,
        )

        dh = 8e-5 * np.sinh(3.5 * (np.radians(_mean(row.alpha2)) ** 3)) * Df**2 * (_mean(row.U) ** 2)
        cp = _mean([row.Cp, upstream.Cp], row.Cp)
        cap = 0.5 * cp * (_mean(row.T0) - _mean(upstream.T0))
        dh = np.clip(dh * self.loss_modifier, 0.0, cap)
        return np.full_like(row.r, dh, dtype=float)


class ImpellerSkinFrictionCoppage(LossBaseClass):
    """Coppage skin friction loss (simplified)."""

    def __init__(
        self,
        number_of_blades: int = 12,
        roughness: float = 1e-5,
        loss_modifier: float = 1.0,
        splitters: bool = False,
    ):
        super().__init__(LossType.Enthalpy)
        self.number_of_blades = number_of_blades
        self.roughness = roughness
        self.loss_modifier = loss_modifier
        self.splitters = splitters

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        r_exit = _tip(row)
        bwidth = _span(row)
        lambda_ratio = _hub(upstream) / max(_tip(upstream), 1e-6)
        Dh = 2 * bwidth if bwidth > 0 else 1.0
        U = _mean(row.U)
        mu = float(getattr(row, "mu", 0.0) or getattr(upstream, "mu", 0.0) or 1.0)
        rho = _mean(getattr(row, "rho", 1.0))
        Re = abs(U) * Dh * rho / mu
        Re = max(Re, 1e3)
        Cf = 0.26 / (Re ** 0.25)
        if self.roughness > 0:
            Cf += 0.11 * (self.roughness / max(Dh, 1e-6)) ** 0.25
        Wbar = 0.125 * (
            _mag(upstream.V)
            + _mag(row.V)
            + _mag(getattr(upstream, "W", upstream.V))
            + 2 * _mag(getattr(row, "W", row.V))
            + 3 * _mag(getattr(row, "W", row.V))
        )
        Ksf = 7.0 if self.splitters else 5.6
        dh = Ksf * Cf * (Wbar**2) / 2
        dh = np.maximum(dh * self.loss_modifier, 0.0)
        return np.full_like(row.r, dh, dtype=float)


class ImpellerSkinFrictionJansen(LossBaseClass):
    """Jansen skin friction loss (simplified Casey/Colebrook form)."""

    def __init__(self, roughness: float = 1e-5, loss_modifier: float = 1.0):
        super().__init__(LossType.Enthalpy)
        self.roughness = roughness
        self.loss_modifier = loss_modifier

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        bwidth = _span(row)
        Dh = 2 * bwidth if bwidth > 0 else 1.0
        U = _mean(row.U)
        mu = float(getattr(row, "mu", 0.0) or getattr(upstream, "mu", 0.0) or 1.0)
        rho = _mean(getattr(row, "rho", 1.0))
        Re = abs(U) * Dh * rho / mu
        Re = max(Re, 1e3)
        # Swamee-Jain approximation for Colebrook
        Cf = 0.25 / (np.log10(self.roughness / (3.7 * Dh) + 5.74 / (Re**0.9)) ** 2)
        Wbar = _mag(getattr(row, "W", row.V))
        dh = Cf * (Wbar**2) / 2
        dh = np.maximum(dh * self.loss_modifier, 0.0)
        return np.full_like(row.r, dh, dtype=float)


class ImpellerVarious(LossBaseClass):
    """Aggregate loss using multiple sub-correlations."""

    def __init__(self):
        super().__init__(LossType.Enthalpy)
        # Compose key submodels with default parameters
        self.blade_loading = ImpellerBladeLoadingCoppage()
        self.clearance = ImpellerClearanceJansen()
        self.mixing = ImpellerMixingJohnston()
        self.disc_friction = ImpellerDiscFrictionDaily()
        self.leakage = ImpellerLeakageAungier()
        self.recirculation = ImpellerRecirculationOh()
        self.incidence = ImpellerIncidenceConrad()
        self.skin_friction = ImpellerSkinFrictionJansen()

    def __call__(self, row: BladeRow, upstream: BladeRow) -> npt.NDArray:
        components = [
            self.blade_loading(row, upstream),
            self.clearance(row, upstream),
            self.mixing(row, upstream),
            self.disc_friction(row, upstream),
            self.leakage(row, upstream),
            self.recirculation(row, upstream),
            self.incidence(row, upstream),
            self.skin_friction(row, upstream),
        ]
        total = np.zeros_like(row.r, dtype=float)
        for comp in components:
            total = total + np.asarray(comp, dtype=float)
        return total


__all__ = [
    "AxialCompressorAungier",
    "AxialCompressorEntropy",
    "AxialCompressorWrightMiller",
    "AxialTurbineAinleyMathiesonOTAC",
    "AxialTurbineKackerOkapuuOTAC",
    "DiffuserVanelessStanitz",
    "ImpellerBladeLoadingAungier",
    "ImpellerBladeLoadingCoppage",
    "ImpellerClearanceJansen",
    "ImpellerDiscFrictionDaily",
    "ImpellerIncidenceAungier",
    "ImpellerIncidenceConrad",
    "ImpellerLeakageAungier",
    "ImpellerMixingAungier",
    "ImpellerMixingJohnston",
    "ImpellerPrescribed",
    "ImpellerRecirculationAungier",
    "ImpellerRecirculationOh",
    "ImpellerSkinFrictionCoppage",
    "ImpellerSkinFrictionJansen",
    "ImpellerVarious",
    "NASA23B20",
    "NASA74A",
    "RadialInput",
]
