"""Shaft-coupled matching between a compressor and a turbine.

Terminology note: in this codebase "Spool" (`CompressorSpool`, `TurbineSpool`)
means one component's blade-row stack - a compressor OR a turbine section -
not a physical rotating shaft. `ShaftMatch` is named for the shaft that
couples two of them: given an already-sized `CompressorSpool`, it sizes and
solves the matching `TurbineSpool` on the same shaft, enforcing shaft power
balance and mass-flow continuity (with simple fuel-air-ratio / bleed-fraction
scalars).

Currently only the compressor-known direction is implemented (the classic
gas-generator matching problem: compressor -> combustor -> turbine). Sizing a
compressor from a known turbine is not yet supported.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field
from typing import Callable, List, Optional, Tuple

import numpy as np
from cantera.composite import Solution
from scipy.optimize import minimize_scalar

from .bladerow import BladeRow
from .compressor_spool import CompressorSpool
from .inlet import Inlet
from .isentropic import area_for_massflow, Massflow
from .outlet import Outlet
from .passage import Passage
from .turbine_spool import TurbineSpool

__all__ = ["ComponentSizingSpec", "SizingEstimate", "MatchResult", "ShaftMatch", "turbine_massflow"]


def turbine_massflow(compressor_massflow: float, fuel_air_ratio: float = 0.0, bleed_fraction: float = 0.0) -> float:
    """Turbine inlet massflow implied by compressor massflow plus fuel/bleed.

    `bleed_fraction` is spool-level extraction upstream of the turbine and is
    NOT the same as the per-row `Coolant` mechanism on `BladeRow` - combining
    the two (to avoid double-counting) is out of scope here.

    Args:
        compressor_massflow: Compressor inlet massflow [kg/s].
        fuel_air_ratio: Fuel-to-air ratio added in the combustor (>= 0).
        bleed_fraction: Fraction of compressor flow bled off before the
            turbine (0 <= bleed_fraction < 1).

    Returns:
        Turbine inlet massflow [kg/s].
    """
    return compressor_massflow * (1.0 + fuel_air_ratio - bleed_fraction)


@dataclass
class ComponentSizingSpec:
    """Design targets and meanline knobs for the component being sized."""

    inlet_T0: float                                   # e.g. combustor exit T4 [K] - cannot be derived
    inlet_P0: Optional[float] = None                  # default: known-side exit P0 * (1 - combustor_pressure_loss)
    combustor_pressure_loss: float = 0.04

    exit_static_pressure: Optional[float] = None      # exactly one of these two is required
    expansion_pressure_ratio: Optional[float] = None  # P0_exit / P0_inlet

    num_stages: int = 1
    eta_total_guess: float = 0.89
    stage_loading: float = 1.8                        # Delta_h0 / U^2 per stage
    ngv_choked: bool = True                            # size the first-stage NGV throat at M=1 (standard practice)
    inlet_mach: float = 0.15                           # only used when ngv_choked is False
    mean_mach: float = 0.6
    exit_mach: float = 0.45
    mean_radius: Optional[float] = None                # default: derived from stage_loading and U
    blockage: float = 0.0

    fluid: Optional[Solution] = None                  # hot-side gas; default: known component's fluid


@dataclass
class SizingEstimate:
    """Result of `ShaftMatch.size()` - a 0-D/1-D sizing pass, no solver involved."""

    massflow: float
    inlet_P0: float
    inlet_T0: float
    exit_P0: float
    exit_T0: float
    exit_P: float
    required_power: float
    rpm: float
    mean_radius: float
    blade_speed: float
    gamma: float
    R: float
    station_area: np.ndarray   # [inlet, mean, exit] annulus areas [m^2]
    station_M: np.ndarray      # target Mach at each station - station_M[0] == 1.0 when ngv_choked
    warnings: List[str] = field(default_factory=list)


@dataclass
class MatchResult:
    """Result of `ShaftMatch.match()`."""

    converged: bool
    iterations: int
    massflow_compressor: float
    massflow_turbine: float
    power_compressor: float
    power_turbine: float
    massflow_residual: float
    power_residual: float
    area_scale: float
    compressor: CompressorSpool
    turbine: TurbineSpool
    sizing: SizingEstimate


def _hot_gas_properties(fluid: Solution, T0: float, P0: float) -> Tuple[float, float, float]:
    """Snapshot-and-restore Cp/gamma/R lookup from a shared Cantera Solution.

    Uses total conditions as a proxy for static conditions, adequate for a
    first-cut sizing estimate. Mirrors the save/restore discipline in
    `bladerow.compute_gas_constants` - `fluid` is typically a shared object,
    so mutating it in place without restoring would leak state.

    Args:
        fluid: Cantera Solution to query.
        T0: Temperature to evaluate properties at [K].
        P0: Pressure to evaluate properties at [Pa].

    Returns:
        Tuple of (Cp [J/(kg*K)], gamma [-], R [J/(kg*K)]).
    """
    saved_state = fluid.state
    try:
        fluid.TP = T0, P0
        Cp = float(fluid.cp)
        Cv = float(fluid.cv)
    finally:
        fluid.state = saved_state
    return Cp, Cp / Cv, Cp - Cv


def _clone_row(row: BladeRow, fluid: Optional[Solution]) -> BladeRow:
    """Deep-copy a blade row without deep-copying its shared Cantera fluid.

    Every `ShaftMatch.build()` call must construct fresh rows rather than
    mutate and re-solve existing ones: `interpolate_streamline_quantities`
    and `compute_power` both mutate rows in place, and reusing a solved row
    across outer iterations would leak state between area-scale guesses.
    """
    ref_fluid = getattr(row, "fluid", None)
    memo = {id(ref_fluid): ref_fluid} if ref_fluid is not None else {}
    cloned = copy.deepcopy(row, memo)
    cloned.fluid = fluid if fluid is not None else ref_fluid
    return cloned


class ShaftMatch:
    """Sizes and solves a `TurbineSpool` to match an already-solved `CompressorSpool`
    on a common shaft.

    Two phases:
        1. `size()` - cheap 0-D/1-D algebra (power balance, mass continuity,
           non-dimensional-MFP-based area sizing). No solver involved.
        2. `build()` / `match()` - constructs a full `TurbineSpool` by scaling
           a user-supplied reference turbine's annulus (angles, loss models,
           and stage layout are reused verbatim; only radii are scaled), then
           iterates the annulus area scale until both the mass-flow and
           shaft-power residuals converge.

    Args:
        known: Already-solved `CompressorSpool` defining the shaft's power
            demand and inlet massflow. Sizing a compressor from a known
            turbine is not yet supported.
        spec: Design targets/knobs for the turbine being sized.
        reference: An existing, already-constructed `TurbineSpool` (typically
            from another example/design) whose `Passage` and blade rows are
            scaled to match `spec`/`known`. Its angles, loss models, chords,
            and stage count are reused as-is.
        fuel_air_ratio: Combustor fuel-to-air ratio (mass basis), >= 0.
        bleed_fraction: Fraction of compressor flow bled upstream of the
            turbine, in [0, 1).
        mechanical_efficiency: Shaft mechanical efficiency, in (0, 1].
        rpm: Shaft speed [rad/s... no, rev/min, matching `CompressorSpool.rpm`].
            Defaults to `known.rpm` (single shaft).
    """

    def __init__(
        self,
        known: CompressorSpool,
        spec: ComponentSizingSpec,
        reference: TurbineSpool,
        fuel_air_ratio: float = 0.0,
        bleed_fraction: float = 0.0,
        mechanical_efficiency: float = 1.0,
        rpm: Optional[float] = None,
    ) -> None:
        if not isinstance(known, CompressorSpool):
            raise NotImplementedError(
                "ShaftMatch currently only supports sizing a turbine from an "
                "already-sized CompressorSpool (compressor -> combustor -> "
                "turbine). Sizing a compressor from a known TurbineSpool is "
                "not yet implemented."
            )
        if fuel_air_ratio < 0:
            raise ValueError("fuel_air_ratio must be >= 0.")
        if not (0.0 <= bleed_fraction < 1.0):
            raise ValueError("bleed_fraction must be in [0, 1).")
        if not (0.0 < mechanical_efficiency <= 1.0):
            raise ValueError("mechanical_efficiency must be in (0, 1].")
        if (spec.exit_static_pressure is None) == (spec.expansion_pressure_ratio is None):
            raise ValueError(
                "ComponentSizingSpec must set exactly one of "
                "exit_static_pressure or expansion_pressure_ratio."
            )

        self.known = known
        self.spec = spec
        self.reference = reference
        self.fuel_air_ratio = float(fuel_air_ratio)
        self.bleed_fraction = float(bleed_fraction)
        self.mechanical_efficiency = float(mechanical_efficiency)
        self.rpm = float(rpm) if rpm is not None else float(known.rpm)
        if self.rpm <= 0:
            raise ValueError("rpm must be positive (set explicitly or via known.rpm).")

    def size(self) -> SizingEstimate:
        """Phase 1: pure-algebra sizing pass. No solver, no geometry construction.

        Returns:
            SizingEstimate with the turbine's target massflow, station total
            conditions, required power, blade speed/radius, and station
            annulus areas (via the non-dimensional mass flow function).
        """
        compressor = self.known
        spec = self.spec
        rows = compressor._all_rows()

        mdot_c = float(rows[1].total_massflow_no_coolant)
        exit_row = rows[-2]
        P0_3 = float(np.mean(exit_row.P0))
        T0_3 = float(np.mean(exit_row.T0))
        P_c = compressor.total_power()
        if P_c <= 0:
            raise ValueError(
                f"CompressorSpool.total_power() = {P_c:.3f} W is not positive; "
                "cannot size a matching turbine. Has `known.solve()` been called?"
            )

        mdot_t = turbine_massflow(mdot_c, self.fuel_air_ratio, self.bleed_fraction)
        P_req = P_c / self.mechanical_efficiency

        T0_4 = float(spec.inlet_T0)
        P0_4 = float(spec.inlet_P0) if spec.inlet_P0 is not None else P0_3 * (1.0 - spec.combustor_pressure_loss)

        fluid = spec.fluid if spec.fluid is not None else compressor.fluid
        Cp_hot, gamma_hot, R_hot = _hot_gas_properties(fluid, T0_4, P0_4)

        warnings: List[str] = []
        dT0 = P_req / (mdot_t * Cp_hot)
        if dT0 / T0_4 > 0.5:
            warnings.append(
                f"Required temperature drop dT0={dT0:.1f} K is >50% of inlet "
                f"T0={T0_4:.1f} K; single-stage sizing assumptions may not hold."
            )
        T0_5 = T0_4 - dT0

        if spec.expansion_pressure_ratio is not None:
            pr = float(spec.expansion_pressure_ratio)
            T0_5_ideal = T0_4 * pr ** ((gamma_hot - 1.0) / gamma_hot)
            implied_eta = dT0 / max(T0_4 - T0_5_ideal, 1e-9)
            if not (0.7 <= implied_eta <= 0.95):
                warnings.append(
                    f"Implied turbine efficiency {implied_eta:.3f} from "
                    "expansion_pressure_ratio is outside the typical range [0.70, 0.95]."
                )
            P0_5 = P0_4 * pr
        else:
            eta = spec.eta_total_guess
            T0_5_ideal = T0_4 - dT0 / max(eta, 1e-6)
            P0_5 = P0_4 * (T0_5_ideal / T0_4) ** (gamma_hot / (gamma_hot - 1.0))

        num_stages = max(spec.num_stages, 1)
        dh0_stage = Cp_hot * dT0 / num_stages
        if spec.mean_radius is not None:
            r_mean = float(spec.mean_radius)
            U = r_mean * self.rpm * np.pi / 30.0
            implied_loading = dh0_stage / max(U**2, 1e-9)
            if not (1.0 <= implied_loading <= 2.5):
                warnings.append(
                    f"Implied stage loading {implied_loading:.2f} at fixed "
                    f"mean_radius={r_mean:.4f} m is outside the typical range [1.0, 2.5]."
                )
        else:
            U = float(np.sqrt(dh0_stage / max(spec.stage_loading, 1e-6)))
            if U > 500.0:
                warnings.append(f"Blade speed U={U:.1f} m/s exceeds a common material limit (~500 m/s).")
            r_mean = U / (self.rpm * np.pi / 30.0)

        if spec.exit_static_pressure is not None:
            exit_P = float(spec.exit_static_pressure)
        else:
            exit_P = P0_5 / (1.0 + (gamma_hot - 1.0) / 2.0 * spec.exit_mach**2) ** (gamma_hot / (gamma_hot - 1.0))

        # The first-stage NGV throat is conventionally choked (M=1) - it is the engine's
        # metering orifice, and this is the classic "choked-MFP" sizing relation:
        # A_4 = mdot*sqrt(R*T04) / (P04 * m~_max). `inlet_mach` is only used as an override
        # when ngv_choked=False (e.g. sizing a non-metering station instead).
        inlet_M = 1.0 if spec.ngv_choked else spec.inlet_mach
        station_T0 = np.array([T0_4, 0.5 * (T0_4 + T0_5), T0_5])
        station_P0 = np.array([P0_4, np.sqrt(P0_4 * P0_5), P0_5])
        station_M = np.array([inlet_M, spec.mean_mach, spec.exit_mach])
        station_area = area_for_massflow(mdot_t, station_P0, station_T0, station_M, gamma_hot, R_hot, spec.blockage)
        station_area = np.atleast_1d(np.asarray(station_area, dtype=float))

        # Self-consistency: feeding each station's area back through Massflow at its
        # target Mach must reproduce mdot_t (this is what the unit tests check).
        check = Massflow(station_P0, station_T0, station_area * (1.0 - spec.blockage), station_M, gamma_hot, R_hot)
        if not np.allclose(check, mdot_t, rtol=1e-6):
            raise AssertionError(f"Sizing self-consistency check failed: {check} != {mdot_t}")

        return SizingEstimate(
            massflow=mdot_t,
            inlet_P0=P0_4, inlet_T0=T0_4,
            exit_P0=P0_5, exit_T0=T0_5, exit_P=exit_P,
            required_power=P_req,
            rpm=self.rpm, mean_radius=r_mean, blade_speed=U,
            gamma=gamma_hot, R=R_hot,
            station_area=station_area, station_M=station_M,
            warnings=warnings,
        )

    def build(self, sizing: SizingEstimate, area_scale: float = 1.0, exit_P: Optional[float] = None) -> TurbineSpool:
        """Phase 2a: construct a fresh `TurbineSpool` scaled from `self.reference`.

        Reuses the reference's angles, loss models, chords, and stage layout
        verbatim; only the annulus radii are scaled (about each station's
        local mean radius, holding the mean radius fixed so blade speed is
        preserved) to hit `sizing`'s target areas.

        Args:
            sizing: Result of `size()`.
            area_scale: Multiplier on the reference annulus height at every
                passage control point. 1.0 reproduces the reference geometry.
            exit_P: Exit static pressure override [Pa]. Defaults to
                `sizing.exit_P`. `match()` treats this as its second free
                knob (area_scale mainly sets flow capacity/massflow; exit_P
                mainly sets the expansion ratio/power) since, for a
                fixed-angle turbine in pressure-balance mode, `spool.massflow`
                is only a solver seed - the achieved massflow and power are
                both outputs of the geometry and boundary pressures, not
                independently dictated.

        Returns:
            A new, unsolved `TurbineSpool` (call `.solve()` or use `match()`).
        """
        ref = self.reference
        ref_passage = ref.passage
        s = max(float(area_scale), 1e-6)

        r_hub = np.asarray(ref_passage.rhub_pts, dtype=float)
        r_shroud = np.asarray(ref_passage.rshroud_pts, dtype=float)
        r_mean = 0.5 * (r_hub + r_shroud)
        new_r_hub = r_mean - (r_mean - r_hub) * s
        new_r_shroud = r_mean + (r_shroud - r_mean) * s

        passage = Passage(
            ref_passage.xhub_pts.tolist(), new_r_hub.tolist(),
            ref_passage.xshroud_pts.tolist(), new_r_shroud.tolist(),
            passageType=ref_passage.passageType,
            zero_phi=getattr(ref_passage, "zero_phi", False),
        )

        fluid = self.spec.fluid if self.spec.fluid is not None else self.known.fluid
        rows = [_clone_row(row, fluid) for row in ref.rows]

        inlet = Inlet(hub_location=ref.inlet.hub_location)
        inlet.alpha2 = np.array(ref.inlet.alpha2, dtype=float, copy=True)
        inlet.init_total(sizing.inlet_P0, sizing.inlet_T0, M=self.spec.inlet_mach)

        outlet = Outlet(num_streamlines=ref.num_streamlines, location=ref.outlet.location)
        outlet.init_static(sizing.exit_P if exit_P is None else float(exit_P), percent_radii=[0.5])

        turbine = TurbineSpool(
            passage, sizing.massflow, inlet, outlet, rows,
            num_streamlines=ref.num_streamlines, fluid=fluid, rpm=self.rpm,
        )
        # solve_radeq.adjust_streamlines is called with 2 args by TurbineSpool
        # but requires 3 (TypeError) - keep this off on every built spool.
        turbine.adjust_streamlines = False
        return turbine

    def match(
        self,
        tol_rel: float = 1e-3,
        area_bounds: Tuple[float, float] = (0.05, 50.0),
        exit_pressure_frac_bounds: Tuple[float, float] = (0.2, 0.95),
        max_iter: int = 15,
    ) -> MatchResult:
        """Phase 2b: solve for the (area scale, exit pressure) pair that closes
        both the mass-flow and shaft-power residuals.

        For a fixed-angle turbine in pressure-balance mode, `spool.massflow`
        is only a solver seed, not an enforced constraint: the achieved
        massflow and power are both outputs of the annulus geometry and the
        boundary pressures. `build()` therefore exposes two knobs instead of
        one: the annulus area scale (which mainly sets flow capacity, i.e.
        massflow) and the exit static pressure (which mainly sets the
        expansion ratio, i.e. power). These are solved as nested 1-D bounded
        root-finds - an outer search over area scale targeting the massflow
        residual, with an inner search over exit-pressure fraction (of inlet
        P0) targeting the power residual at each area scale - using the same
        `minimize_scalar(..., method="bounded")` pattern as
        `CompressorSpool.solve_massflow_for_pressure_ratio`.

        Args:
            tol_rel: Relative tolerance on both residuals for `converged`.
                Search precision (the optimizers' `xatol`) is independent of
                this - loosening `tol_rel` reports success sooner, it does
                not make the underlying search coarser.
            area_bounds: (lower, upper) bounds for the area scale search.
            exit_pressure_frac_bounds: (lower, upper) bounds for exit static
                pressure as a fraction of inlet total pressure.
            max_iter: Maximum bounded-solver iterations, applied to both the
                inner and outer search.

        Returns:
            MatchResult with both residuals reported honestly - a converged
            massflow residual does not imply the power residual also closed
            (or vice versa); `converged` requires both.
        """
        sizing = self.size()
        mdot_required = sizing.massflow
        power_required = sizing.required_power
        P0_inlet = sizing.inlet_P0
        search_xatol = min(tol_rel, 1e-4)  # search precision independent of the reported tolerance

        last: dict = {}

        def inner_objective(exit_frac: float, area_scale: float) -> float:
            exit_P = exit_frac * P0_inlet
            turbine = self.build(sizing, area_scale=area_scale, exit_P=exit_P)
            turbine.solve()
            mdot_achieved = float(turbine._all_rows()[1].total_massflow_no_coolant)
            power_achieved = float(turbine.total_power())
            last.update(turbine=turbine, area_scale=area_scale, exit_P=exit_P, mdot=mdot_achieved, power=power_achieved)
            return abs(power_achieved - power_required) / power_required

        def outer_objective(area_scale: float) -> float:
            minimize_scalar(
                inner_objective, bounds=exit_pressure_frac_bounds, method="bounded",
                args=(area_scale,), options={"xatol": search_xatol, "maxiter": max_iter},
            )
            resid_mdot = (last["mdot"] - mdot_required) / mdot_required
            resid_power = (last["power"] - power_required) / power_required
            # Massflow alone is nearly flat over a wide range of area_scale (the inner
            # exit-pressure search can usually bring achieved massflow close regardless),
            # so minimizing only resid_mdot leaves the outer search free to settle
            # anywhere along that flat region - including points with a poor power
            # match. Combine both residuals so the outer search actually favors area
            # scales where the inner solve's power match is also good.
            return float(np.hypot(resid_mdot, resid_power))

        res_outer = minimize_scalar(
            outer_objective, bounds=area_bounds, method="bounded",
            options={"xatol": search_xatol, "maxiter": max_iter},
        )
        area_scale = float(res_outer.x)
        if abs(last.get("area_scale", np.nan) - area_scale) > 1e-9:
            outer_objective(area_scale)  # ensure `last` reflects the returned optimum

        turbine = last["turbine"]
        mdot_achieved = last["mdot"]
        power_achieved = last["power"]
        resid_mdot = (mdot_achieved - mdot_required) / mdot_required
        resid_power = (power_achieved - power_required) / power_required
        converged = abs(resid_mdot) <= tol_rel and abs(resid_power) <= tol_rel

        return MatchResult(
            converged=converged,
            iterations=int(res_outer.nfev),
            massflow_compressor=float(self.known._all_rows()[1].total_massflow_no_coolant),
            massflow_turbine=mdot_achieved,
            power_compressor=float(self.known.total_power()),
            power_turbine=power_achieved,
            massflow_residual=float(resid_mdot),
            power_residual=float(resid_power),
            area_scale=area_scale,
            compressor=self.known,
            turbine=turbine,
            sizing=sizing,
        )

    def export(self, prefix: str, result: MatchResult) -> None:
        """Export the compressor, turbine, and match summary to JSON.

        Args:
            prefix: Filename prefix; writes `{prefix}_compressor.json`,
                `{prefix}_turbine.json`, and `{prefix}_shaft_match.json`.
            result: Output of `match()`.
        """
        import json

        result.compressor.export_properties(f"{prefix}_compressor.json")
        result.turbine.export_properties(f"{prefix}_turbine.json")

        class NumpyEncoder(json.JSONEncoder):
            def default(self, obj):  # type: ignore[override]
                if isinstance(obj, np.ndarray):
                    return obj.tolist()
                return super().default(obj)

        summary = {
            "converged": result.converged,
            "iterations": result.iterations,
            "massflow_compressor": result.massflow_compressor,
            "massflow_turbine": result.massflow_turbine,
            "power_compressor": result.power_compressor,
            "power_turbine": result.power_turbine,
            "massflow_residual": result.massflow_residual,
            "power_residual": result.power_residual,
            "area_scale": result.area_scale,
            "sizing": {
                "massflow": result.sizing.massflow,
                "inlet_P0": result.sizing.inlet_P0,
                "inlet_T0": result.sizing.inlet_T0,
                "exit_P0": result.sizing.exit_P0,
                "exit_T0": result.sizing.exit_T0,
                "exit_P": result.sizing.exit_P,
                "required_power": result.sizing.required_power,
                "rpm": result.sizing.rpm,
                "mean_radius": result.sizing.mean_radius,
                "blade_speed": result.sizing.blade_speed,
                "station_area": result.sizing.station_area,
                "station_M": result.sizing.station_M,
                "warnings": result.sizing.warnings,
            },
        }
        with open(f"{prefix}_shaft_match.json", "w") as f:
            json.dump(summary, f, indent=4, cls=NumpyEncoder)
