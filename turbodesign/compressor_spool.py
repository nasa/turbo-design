# type: ignore[arg-type, reportUnknownArgumentType]
from __future__ import annotations
from turtle import down, up
from typing import Dict, List, Union, Optional, Tuple
import json

import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt

from cantera.composite import Solution
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar

# --- Project-local imports
from .bladerow import BladeRow, interpolate_streamline_quantities
from .enums import RowType, LossType
from .loss.turbine import TD2
from .passage import Passage
from .inlet import Inlet
from .outlet import Outlet, OutletType
from .compressor_math import rotor_calc, stator_calc, polytropic_efficiency
from .flow_math import compute_massflow, compute_streamline_areas, compute_power
from .turbine_math import (
    inlet_calc,
    compute_gas_constants,
    compute_reynolds,
)
from .solve_radeq import adjust_streamlines, radeq
from pyturbo.helper import convert_to_ndarray

# Default fraction of the stator-to-stator pressure rise attributed to a rotor
DEFAULT_ROTOR_PRESSURE_FRACTION = 0.5

class CompressorSpool:
    """Used to design compressors 

    This class (formerly named *Spool*) encapsulates both the generic geometry/plotting
    utilities from the original base spool and the turbine-solving logic that lived
    in the turbine-specific spool implementation.

    Notes on differences vs. the two-class design:
    - `field(default_factory=...)` was previously used on a non-dataclass attribute
      (`t_streamline`). Here it's handled in `__init__` to avoid a silent bug.
    - `fluid` defaults to `Solution('air.yaml')` if not provided.
    - All turbine-specific methods (initialize/solve/massflow balancing/etc.) are
      preserved here. If you ever add a *CompressorSpool* in the future, consider
      splitting turbine/compressor behaviors behind a strategy/solver object.
    """

    # Class-level defaults (avoid mutable defaults here!)
    rows: List[BladeRow]
    massflow: float
    rpm: float

    # Types/attributes documented for linters; values set in __init__
    passage: Passage
    t_streamline: npt.NDArray
    num_streamlines: int

    _fluid: Solution
    _adjust_streamlines: bool

    def __init__(
        self,
        passage: Passage,
        massflow: float,
        inlet: Inlet,
        outlet: Outlet,
        rows: List[BladeRow],
        num_streamlines: int = 3,
        fluid: Optional[Solution] = None,
        rpm: float = -1,
        rotor_pressure_fraction: float = DEFAULT_ROTOR_PRESSURE_FRACTION,
    ) -> None:
        """Initialize a compressor spool

        Args:
            passage: Passage defining hub and shroud
            massflow: massflow at spool inlet
            inlet: Inlet object
            outlet: Outlet object
            rows: List of blade rows between inlet and outlet
            num_streamlines: number of streamlines used through the meridional passage
            fluid: cantera gas solution; defaults to air.yaml if None
            rpm: RPM for the entire spool. Individual rows can override later.
            rotor_pressure_fraction: Fraction of total pressure rise in rotors (0.0 to 1.0)
        """
        self.passage = passage
        self.massflow = massflow
        self.inlet = inlet
        self.outlet = outlet
        self.rows = rows
        self.num_streamlines = num_streamlines
        self._fluid = fluid if fluid is not None else Solution("air.yaml")
        self.rpm = rpm
        self.rotor_pressure_fraction = float(np.clip(rotor_pressure_fraction, 0.0, 1.0))

        # Previously this used dataclasses.field on a non-dataclass; do it explicitly
        self.t_streamline = np.zeros((10,), dtype=float)
        self._adjust_streamlines = True
        self.convergence_history: List[Dict] = []

        # Assign IDs, RPMs, and axial chords where appropriate
        for i, br in enumerate(self._all_rows()):
            br.id = i
            if not isinstance(br, (Outlet)):
                br.rpm = rpm
                br.axial_chord = br.hub_location * self.passage.hub_length
            # Freeze any configured P0_ratio targets for later use (diagnostics may overwrite P0_ratio).
            if getattr(br, "P0_ratio_target", 0.0) == 0 and getattr(br, "P0_ratio", 0.0) != 0:
                br.P0_ratio_target = br.P0_ratio
            if isinstance(br, BladeRow) and br.row_type == RowType.Rotor:
                setattr(br, "rotor_pressure_fraction", getattr(br, "rotor_pressure_fraction", self.rotor_pressure_fraction))

        # Propagate initial fluid to rows
        for br in self._all_rows():
            br.fluid = self._fluid

    def _all_rows(self) -> List[BladeRow]:
        """Convenience to iterate inlet + interior rows + outlet."""
        return [self.inlet, *self.rows, self.outlet]

    @property
    def blade_rows(self) -> List[BladeRow]:
        """Backwards-compatible combined row list."""
        return self._all_rows()

    def set_rotor_pressure_fraction(self, value: float) -> None:
        """Update default pressure split fraction for all rotor rows."""
        self.rotor_pressure_fraction = float(np.clip(value, 0.0, 1.0))
        for row in self.rows:
            if row.row_type == RowType.Rotor:
                setattr(row, "rotor_pressure_fraction", self.rotor_pressure_fraction)

    # ------------------------------
    # Properties
    # ------------------------------
    @property
    def fluid(self) -> Optional[Solution]:
        return self._fluid

    @fluid.setter
    def fluid(self, newFluid: Solution) -> None:
        """Change the gas used in the spool and cascade to rows."""
        self._fluid = newFluid
        for br in self._all_rows():
            br.fluid = self._fluid

    @property
    def adjust_streamlines(self) -> bool:
        return self._adjust_streamlines

    @adjust_streamlines.setter
    def adjust_streamlines(self, val: bool) -> None:
        self._adjust_streamlines = val

    # ------------------------------
    # Row utilities
    # ------------------------------
    def set_blade_row_rpm(self, index: int, rpm: float) -> None:
        self.rows[index].rpm = rpm

    def set_blade_row_type(self, blade_row_index: int, rowType: RowType) -> None:
        self.rows[blade_row_index].row_type = rowType

    def set_blade_row_exit_angles(
        self,
        radius: Dict[int, List[float]],
        beta: Dict[int, List[float]],
        IsSupersonic: bool = False,
    ) -> None:
        """Set intended exit flow angles for rows (useful when geometry is fixed)."""
        for k, v in radius.items():
            self.rows[k].radii_geom = v
        for k, v in beta.items():
            self.rows[k].beta_geom = v
            self.rows[k].beta_fixed = True
        for br in self._all_rows():
            br.solution_type = "supersonic" if IsSupersonic else "subsonic"

    # ------------------------------
    # Streamline setup/geometry
    # ------------------------------
    def initialize_streamlines(self) -> None:
        """Initialize streamline storage per row and compute curvature."""
        for row in self._all_rows():
            row.phi = np.zeros((self.num_streamlines,))
            row.rm = np.zeros((self.num_streamlines,))
            row.r = np.zeros((self.num_streamlines,))
            row.m = np.zeros((self.num_streamlines,))

            t_radial = np.array([0.5]) if self.num_streamlines == 1 else np.linspace(0, 1, self.num_streamlines)
            self.calculate_streamline_curvature(row, t_radial)

            if self.num_streamlines == 1:
                area = self.passage.get_area(row.hub_location)
                row.total_area = area
                row.area = np.array([area])

            # Ensure a loss model exists on blade rows
            if not isinstance(row, (Inlet, Outlet)) and row.loss_function is None:
                row.loss_function = TD2()

        # With radii known, couple blade geometry (pitch/chord/stagger) if specified
        for row in self._all_rows():
            if isinstance(row, BladeRow) and row.row_type not in (RowType.Inlet, RowType.Outlet):
                try:
                    row.synchronize_blade_geometry()
                except Exception:
                    pass

    def calculate_streamline_curvature(
        self, row: BladeRow, t_radial: Union[List[float], npt.NDArray]
    ) -> None:
        """Interpolate passage curvature metrics onto a blade row.

        Args:
            row: BladeRow to populate with phi, rm, r, and m along streamlines.
            t_radial: Parametric hub-to-shroud locations (0–1) at which to sample curvature.
        """
        for i, tr in enumerate(t_radial):
            t_s, x_s, r_s = self.passage.get_streamline(tr)
            phi, rm, r = self.passage.streamline_curvature(x_s, r_s)
            row.phi[i] = float(interp1d(t_s, phi)(row.hub_location))
            row.rm[i] = float(interp1d(t_s, rm)(row.hub_location))
            row.r[i] = float(interp1d(t_s, r)(row.hub_location))
            row.m[i] = float(
                interp1d(t_s, self.passage.get_m(tr, resolution=len(t_s)))(row.hub_location)
            )
        chord = np.asarray(row.chord, dtype=float)
        mean_chord = float(np.mean(chord)) if chord.size else 0.0
        if row.num_blades and mean_chord != 0:
            mean_r = float(np.mean(row.r))
            pitch = 2 * np.pi * mean_r / row.num_blades
            row.pitch_to_chord = pitch / mean_chord

    # ------------------------------
    # initialization/solve
    # ------------------------------
    def initialize(self) -> None:
        """Initialize massflow and thermodynamic state through rows (compressor).

        Sets inlet totals, interpolates geometry, propagates gas properties, and
        runs per-row calcs to seed the solver.
        """
        rows = self._all_rows()

        # Inlet
        W0 = self.massflow
        inlet: Inlet = self.inlet
        if self.fluid:
            inlet.__initialize_fluid__(self.fluid)  # type: ignore[arg-type]
        else:
            inlet.__initialize_fluid__(  # type: ignore[call-arg]
                R=rows[1].R,
                gamma=rows[1].gamma,
                Cp=rows[1].Cp,
            )

        inlet.total_massflow = W0
        inlet.total_massflow_no_coolant = W0
        inlet.massflow = np.array([W0]) if self.num_streamlines == 1 else np.linspace(0, 1, self.num_streamlines) * W0
        
        inlet.__interpolate_quantities__(self.num_streamlines)  # type: ignore[attr-defined]
        inlet.__initialize_velocity__(self.passage, self.num_streamlines)  # type: ignore[attr-defined]
        interpolate_streamline_quantities(inlet, self.passage, self.num_streamlines)

        compute_gas_constants(inlet, self.fluid)
        inlet_calc(inlet)

        for row in rows:
            interpolate_streamline_quantities(row, self.passage, self.num_streamlines)

        # Pass T0, P0 to downstream rows
        for i in range(1, len(rows) - 1):
            upstream = rows[i - 1]
            downstream = rows[i + 1] if i + 1 < len(rows) else None

            row = rows[i]
            if row.coolant is not None:
                T0c = row.coolant.T0
                P0c = row.coolant.P0
                W0c = row.coolant.massflow_percentage * self.massflow
                Cpc = row.coolant.Cp
            else:
                T0c = 100
                P0c = 0
                W0c = 0
                Cpc = 0

            T0 = upstream.T0
            P0 = upstream.P0
            Cp = upstream.Cp

            T0 = (W0 * Cp * T0 + W0c * Cpc * T0c) / (Cpc * W0c + Cp * W0)
            P0 = (W0 * Cp * P0 + W0c * Cpc * P0c) / (Cpc * W0c + Cp * W0)
            Cp = (W0 * Cp + W0c * Cpc) / (W0c + W0) if (W0c + W0) != 0 else Cp

            if row.row_type == RowType.Stator:
                T0 = upstream.T0
            else:
                T0 = upstream.T0 - row.power / (Cp * (W0 + W0c))

            W0 += W0c
            row.T0 = T0
            row.P0 = P0
            row.Cp = Cp
            row.total_massflow = W0
            row.massflow = np.array([row.total_massflow]) if self.num_streamlines == 1 else np.linspace(0, 1, self.num_streamlines) * row.total_massflow

            # Pass gas constants
            row.rho = upstream.rho
            row.gamma = upstream.gamma
            row.R = upstream.R

            total_area, streamline_area = compute_streamline_areas(row)
            row.total_area = total_area
            row.area = streamline_area
            if row.row_type == RowType.Stator or row.row_type == RowType.IGV:
                if row.row_type == RowType.IGV:
                    row.P0_is = upstream.P0
                stator_calc(row, upstream, calculate_vm=True)  # type: ignore[arg-type]
            elif row.row_type == RowType.Rotor:
                # Align rotor ideal P0 target with downstream stator if provided (stage-level target)
                if downstream and downstream.row_type == RowType.Stator:
                    downstream.P0 = upstream.P0*downstream.P0_ratio
                    downstream.Yp = downstream.loss_function(row, upstream)
                    downstream.P0_is = downstream.P0 + downstream.Yp * (upstream.P0-upstream.P)
                    row.P0_ratio = downstream.P0_ratio
                else:
                    row.P0 = row.P0_ratio * upstream.P0
                rotor_calc(row, upstream,calculate_vm=True)
                compute_power(row, upstream, is_compressor=True)

    def solve(self) -> None:
        """Run streamline initialization and solve the compressor flow field.

        The solution method is determined by the outlet configuration:
        - If outlet.outlet_type is massflow_static_pressure: use angle matching
        - Otherwise: use pressure balance
        """
        self.initialize_streamlines()
        self.initialize()

        if self.outlet.outlet_type == OutletType.massflow_static_pressure:
            print("Using angle matching mode: blade exit angles will be adjusted to match specified massflow")
            self._angle_match()
        else:
            print("Using pressure balance mode: blade exit angles are fixed, total pressures will be adjusted")
            self.balance_pressure()

    def solve_angle_match(self) -> None:
        """Explicit angle-matching solve by temporarily setting outlet type."""
        prev_type = self.outlet.outlet_type
        prev_massflow = getattr(self.outlet, 'total_massflow', None)
        try:
            if prev_massflow is None:
                self.outlet.total_massflow = self.massflow
            self.outlet.outlet_type = OutletType.massflow_static_pressure
            self.solve()
        finally:
            self.outlet.outlet_type = prev_type
            if prev_massflow is None and hasattr(self.outlet, 'total_massflow'):
                delattr(self.outlet, 'total_massflow')

    def solve_balance_pressure(self) -> None:
        """Explicit pressure-balance solve by temporarily setting outlet type."""
        prev_type = self.outlet.outlet_type
        try:
            self.outlet.outlet_type = OutletType.total_pressure
            self.solve()
        finally:
            self.outlet.outlet_type = prev_type

    def overall_pressure_ratio(self) -> float:
        """Compute overall total pressure ratio (inlet to last internal row)."""
        rows = self._all_rows()
        if len(rows) < 2:
            return 1.0
        return float(np.mean(np.mean(rows[-2].P0 / self.inlet.P0) ))

    def overall_polytropic_efficiency(self) -> float:
        """Compute overall polytropic efficiency from inlet to last internal row."""
        rows = self._all_rows()
        if len(rows) < 2:
            return 0.0
        pi = float(np.mean(rows[-2].P0) / np.mean(self.inlet.P0))
        tau = float(np.mean(rows[-2].T0)/np.mean(self.inlet.T0))
        gamma = float(np.mean(self.inlet.gamma)) if hasattr(self.inlet, "gamma") else 1.4
        if tau <= 0 or abs(np.log(tau)) < 1e-12 or pi <= 1.0:
            return 0.0
        return ((gamma - 1.0) / gamma) * np.log(pi) / np.log(tau)

    def solve_massflow_for_pressure_ratio(self, target_pr: float, bounds: tuple[float, float], meanline: bool = False) -> tuple[float, float]:
        """Solve inlet massflow to hit a target overall total-pressure ratio.

        Args:
            target_pr: desired overall P0 ratio (inlet / last internal row).
            bounds: (lower, upper) bounds for massflow during search.
            meanline: if True, force a single streamline and disable streamline adjustment.

        Returns:
            Tuple of (converged massflow, achieved pressure ratio).
        """
        if meanline:
            self.num_streamlines = 1
            self._adjust_streamlines = False

        lower, upper = bounds
        if lower <= 0 or upper <= 0 or lower >= upper:
            raise ValueError("Massflow bounds must be positive and (lower < upper).")

        def objective(mdot: float) -> float:
            self.massflow = mdot
            self.solve()
            achieved = self.overall_pressure_ratio()
            return (achieved - target_pr) ** 2

        res = minimize_scalar(objective, bounds=bounds, method="bounded")
        self.massflow = float(res.x)
        self.solve()
        achieved = self.overall_pressure_ratio()
        return self.massflow, achieved
    
    def balance_pressure(self) -> None:
        """Balance Pressure assumes we know:
            1. The blade angles
            2. Total Pressure Ratio 
            3. Massflow
            
            We find the static pressures in between the blade rows such that massflow is balanced.
            Implemented by marching rows (compressor mode) without guessing pressure ratios.
        """
        rows = self._all_rows()

        print("Looping to converge massflow (compressor)")
        loop_iter = 0
        max_iter = 10
        prev_err = 1e9
        self.convergence_history = []  # Reset convergence history
        while loop_iter < max_iter:
            for i in range(1, len(rows) - 1):
                row = rows[i]
                upstream = rows[i - 1]
                downstream = rows[i + 1] if i + 1 < len(rows) else None

                if row.row_type == RowType.Inlet:
                    row.Yp = 0
                    continue

                if row.row_type == RowType.Rotor:
                    rotor_calc(row, upstream, calculate_vm=True)
                    if self.num_streamlines > 1:
                        row = radeq(row, upstream, downstream)
                        compute_gas_constants(row, self.fluid)
                        rotor_calc(row, upstream, calculate_vm=False)
                elif row.row_type == RowType.Stator or row.row_type == RowType.IGV:
                    if row.row_type == RowType.IGV:
                        row.P0_is = upstream.P0
                    stator_calc(row, upstream, calculate_vm=True)
                    if self.num_streamlines > 1:
                        row = radeq(row, upstream, downstream)
                        compute_gas_constants(row, self.fluid)
                        stator_calc(row, upstream, calculate_vm=False)

                compute_gas_constants(row, self.fluid)
                compute_power(row, upstream, is_compressor=True)

            target = rows[1].total_massflow_no_coolant
            self.inlet.massflow = np.array([target]) if self.num_streamlines == 1 else np.linspace(0, 1, self.num_streamlines) * target
            self.inlet.total_massflow_no_coolant = rows[1].total_massflow_no_coolant
            self.inlet.total_massflow = rows[1].total_massflow_no_coolant
            self.inlet.calculated_massflow = self.inlet.total_massflow_no_coolant
            inlet_calc(self.inlet)

            # if self.adjust_streamlines:
            #     adjust_streamlines(rows[:-1], self.passage, np.linspace(0, 1, self.num_streamlines))

            self.outlet.transfer_quantities(rows[-2])
            self.outlet.P = self.outlet.get_static_pressure(self.outlet.percent_hub_shroud)

            err = self._massflow_std(rows[1:-1])
            loop_iter += 1
            print(f"Loop {loop_iter} massflow convergence error:{err}")

            # Store convergence history
            self.convergence_history.append({
                'iteration': loop_iter,
                'massflow_std': float(err),
                'massflow_change': float(abs(err - prev_err)),
                'relative_change': float(abs((err - prev_err) / max(err, 1e-6))),
                'massflow': float(rows[1].total_massflow_no_coolant)
            })

            denom = max(err, 1e-6)
            if abs((err - prev_err) / denom) <= 0.05:
                break
            prev_err = err

        compute_reynolds(rows, self.passage)

    @staticmethod
    def _massflow_std(blade_rows: List[BladeRow]) -> float:
        """Compute standard deviation of massflow across rows for diagnostics."""
        totals = []
        for row in blade_rows:
            if hasattr(row, "total_massflow_no_coolant"):
                totals.append(row.total_massflow_no_coolant)
            elif len(getattr(row, "massflow", [])) > 0:
                totals.append(row.massflow[-1])
        return float(np.std(totals)) if totals else 0.0
    # ------------------------------
    # Massflow / angle matching
    # ------------------------------
    def _angle_match(self) -> None:
        """Match massflow between streamtubes by tweaking exit angles."""
        blade_rows = self._all_rows()
        self.convergence_history = []  # Reset convergence history
        prev_err = 1e9

        for iter_num in range(3):
            for i, row in enumerate(blade_rows):
                # Only adjust blade rows; skip inlet/outlet and other utility rows
                if row.row_type not in (RowType.Rotor, RowType.Stator):
                    continue

                upstream = blade_rows[i - 1] if i > 0 else blade_rows[i]
                downstream = blade_rows[i + 1] if i < len(blade_rows) - 1 else None

                if row.row_type == RowType.Stator:
                    bounds = [0, 80]
                elif row.row_type == RowType.Rotor:
                    bounds = [-80, 0]
                else:
                    bounds = [0, 0]

                for j in range(1, self.num_streamlines):
                    res = minimize_scalar(
                        match_massflow_objective,
                        bounds=bounds,
                        args=(j, row, upstream, downstream, self.fluid),
                        options={'xatol': 1e-3},
                        method="bounded",
                    )
                    if row.row_type == RowType.Rotor:
                        row.beta2[j] = np.radians(res.x)
                        row.beta2[0] = 1 / (len(row.beta2) - 1) * row.beta2[1:].sum()
                    elif row.row_type == RowType.Stator:
                        row.alpha2[j] = np.radians(res.x)
                        row.alpha2[0] = 1 / (len(row.alpha2) - 1) * row.alpha2[1:].sum()
                compute_gas_constants(upstream, self.fluid)
                compute_gas_constants(row, self.fluid)
                compute_massflow(row)
                compute_power(row, upstream, is_compressor=True)

            # Track convergence history
            err = self._massflow_std(blade_rows[1:-1])
            self.convergence_history.append({
                'iteration': iter_num + 1,
                'massflow_std': float(err),
                'massflow_change': float(abs(err - prev_err)),
                'relative_change': float(abs((err - prev_err) / max(err, 1e-6))),
                'massflow': float(blade_rows[1].total_massflow_no_coolant)
            })
            prev_err = err
            print(f"Angle match iteration {iter_num + 1}, massflow std: {err:.6f}")


    # ------------------------------
    # Export / Plotting
    # ------------------------------
    def export_properties(self, filename: str = "compressor_spool.json") -> None:
        """Export compressor spool properties and blade row data to JSON file.

        Exports comprehensive compressor design data including blade row properties,
        streamline coordinates, efficiency metrics, pressure ratios, stage loading,
        and power calculations for each stage. Useful for post-processing and result
        archiving.

        Args:
            filename: Output JSON file path (default: "compressor_spool.json")

        Returns:
            None. Writes JSON file to specified path.

        Example:
            >>> spool.export_properties("r35_compressor_results.json")
        """
        blade_rows = self._all_rows()
        blade_rows_out = []
        degree_of_reaction = []
        total_total_efficiency = []
        total_static_efficiency = []
        stage_loading = []
        euler_power = []
        enthalpy_power = []
        x_streamline = np.zeros((self.num_streamlines, len(blade_rows)))
        r_streamline = np.zeros((self.num_streamlines, len(blade_rows)))
        massflow = []

        for indx, row in enumerate(blade_rows):
            blade_rows_out.append(row.to_dict())
            if row.row_type == RowType.Rotor:
                degree_of_reaction.append(
                    (
                        (blade_rows[indx - 1].P - row.P)
                        / (blade_rows[indx - 2].P - row.P)
                    ).mean()
                )
                total_total_efficiency.append(row.eta_total)
                total_static_efficiency.append(row.eta_static)
                stage_loading.append(row.stage_loading)
                euler_power.append(row.euler_power)
                enthalpy_power.append(row.power)
            if row.row_type not in (RowType.Inlet, RowType.Outlet):
                massflow.append(row.massflow[-1])

            for j, p in enumerate(row.percent_hub_shroud):
                t, x, r = self.passage.get_streamline(p)
                x_streamline[j, indx] = float(interp1d(t, x)(row.percent_hub))
                r_streamline[j, indx] = float(interp1d(t, r)(row.percent_hub))

        Pratio_Total_Total = np.mean(self.inlet.P0 / blade_rows[-2].P0)
        Pratio_Total_Static = np.mean(self.inlet.P0 / blade_rows[-2].P)
        flow_fn_massflow = float(np.mean(massflow)) if massflow else 0.0
        FlowFunction = flow_fn_massflow * np.sqrt(self.inlet.T0.mean()) * float(np.mean(self.inlet.P0)) / 1000
        CorrectedSpeed = self.rpm * np.pi / 30 / np.sqrt(self.inlet.T0.mean())
        EnergyFunction = (
            (self.inlet.T0 - blade_rows[-2].T0)
            * 0.5
            * (self.inlet.Cp + blade_rows[-2].Cp)
            / self.inlet.T0
        )
        EnergyFunction = np.mean(EnergyFunction)

        # English-unit conversions
        massflow_kg_s = float(np.mean(massflow)) if massflow else 0.0
        massflow_lbm_s = massflow_kg_s / 0.45359237
        euler_power_hp = [p / 745.7 for p in euler_power]
        enthalpy_power_hp = [p / 745.7 for p in enthalpy_power]

        data = {
            "blade_rows": blade_rows_out,
            "massflow": massflow_kg_s,
            "massflow_lbm_s": massflow_lbm_s,
            "rpm": self.rpm,
            "r_streamline": r_streamline.tolist(),
            "x_streamline": x_streamline.tolist(),
            "rhub": self.passage.rhub_pts.tolist(),
            "rshroud": self.passage.rshroud_pts.tolist(),
            "xhub": self.passage.xhub_pts.tolist(),
            "xshroud": self.passage.xshroud_pts.tolist(),
            "num_streamlines": self.num_streamlines,
            "euler_power": euler_power,
            "euler_power_hp": euler_power_hp,
            "enthalpy_power": enthalpy_power,
            "enthalpy_power_hp": enthalpy_power_hp,
            "total-total_efficiency": total_total_efficiency,
            "total-static_efficiency": total_static_efficiency,
            "stage_loading": stage_loading,
            "degree_of_reaction": degree_of_reaction,
            "Pratio_Total_Total": float(Pratio_Total_Total),
            "Pratio_Total_Static": float(Pratio_Total_Static),
            "FlowFunction": float(FlowFunction),
            "CorrectedSpeed": float(CorrectedSpeed),
            "EnergyFunction": float(EnergyFunction),
            "eta_polytropic_overall": float(self.overall_polytropic_efficiency()),
            "units": {
                "massflow": {"metric": "kg/s", "english": "lbm/s"},
                "rpm": {"metric": "rpm", "english": "rpm"},
                "euler_power": {"metric": "W", "english": "hp"},
                "enthalpy_power": {"metric": "W", "english": "hp"},
                "Pratio_Total_Total": {"metric": "—", "english": "—"},
                "Pratio_Total_Static": {"metric": "—", "english": "—"},
                "FlowFunction": {"metric": "kg/s·K^0.5·Pa", "english": "lbm/s·R^0.5·psf"},
                "CorrectedSpeed": {"metric": "rad/s·K^-0.5", "english": "rad/s·R^-0.5"},
                "EnergyFunction": {"metric": "—", "english": "—"},
            },
        }

        class NumpyEncoder(json.JSONEncoder):
            def default(self, obj):  # type: ignore[override]
                if isinstance(obj, np.ndarray):
                    return obj.tolist()
                return super().default(obj)

        with open(filename, "w") as f:
            json.dump(data, f, indent=4, cls=NumpyEncoder)

    def plot(self) -> None:
        """Plot hub/shroud and streamlines."""
        blade_rows = self._all_rows()
        plt.figure(num=1, clear=True, dpi=150, figsize=(15, 10))
        plt.plot(
            self.passage.xhub_pts,
            self.passage.rhub_pts,
            label="hub",
            linestyle="solid",
            linewidth=2,
            color="black",
        )
        plt.plot(
            self.passage.xshroud_pts,
            self.passage.rshroud_pts,
            label="shroud",
            linestyle="solid",
            linewidth=2,
            color="black",
        )

        hub_length = np.sum(
            np.sqrt(np.diff(self.passage.xhub_pts) ** 2 + np.diff(self.passage.rhub_pts) ** 2)
        )
        x_streamline = np.zeros((self.num_streamlines, len(self.blade_rows)))
        r_streamline = np.zeros((self.num_streamlines, len(self.blade_rows)))
        for i in range(len(blade_rows)):
            x_streamline[:, i] = blade_rows[i].x
            r_streamline[:, i] = blade_rows[i].r

        for i in range(1, len(blade_rows) - 1):
            plt.plot(x_streamline[:, i], r_streamline[:, i], "--b", linewidth=1.5)

        for i, row in enumerate(blade_rows):
            plt.plot(row.x, row.r, linestyle="dashed", linewidth=1.5, color="blue", alpha=0.4)
            plt.plot(x_streamline[:, i], r_streamline[:, i], "or")

            if i == 0:
                pass
            else:
                upstream = blade_rows[i - 1]
                if upstream.row_type == RowType.Inlet:
                    cut_line1, _, _ = self.passage.get_cutting_line(
                        (row.hub_location * hub_length + (0.5 * row.blade_to_blade_gap * row.axial_chord) - row.axial_chord)
                        / hub_length
                    )
                else:
                    cut_line1, _, _ = self.passage.get_cutting_line(
                        (upstream.hub_location * hub_length) / hub_length
                    )
                cut_line2, _, _ = self.passage.get_cutting_line(
                    (row.hub_location * hub_length - (0.5 * row.blade_to_blade_gap * row.axial_chord)) / hub_length
                )

            if row.row_type == RowType.Stator:
                x1, r1 = cut_line1.get_point(np.linspace(0, 1, 10))
                plt.plot(x1, r1, "m")
                x2, r2 = cut_line2.get_point(np.linspace(0, 1, 10))
                plt.plot(x2, r2, "m")
                x_text = (x1 + x2) / 2
                r_text = (r1 + r2) / 2
                plt.text(x_text.mean(), r_text.mean(), "Stator", fontdict={"fontsize": "xx-large"})
            elif row.row_type == RowType.Rotor:
                x1, r1 = cut_line1.get_point(np.linspace(0, 1, 10))
                plt.plot(x1, r1, color="brown")
                x2, r2 = cut_line2.get_point(np.linspace(0, 1, 10))
                plt.plot(x2, r2, color="brown")
                x_text = (x1 + x2) / 2
                r_text = (r1 + r2) / 2
                plt.text(x_text.mean(), r_text.mean(), "Rotor", fontdict={"fontsize": "xx-large"})

        plt.axis("scaled")
        plt.savefig("Meridional.png", transparent=False, dpi=150)
        plt.show()

    def plot_velocity_triangles(self) -> None:
        """Plot velocity triangles for each blade row (turbines).
        """
        blade_rows = self._all_rows()
        prop = dict(arrowstyle="-|>,head_width=0.4,head_length=0.8", shrinkA=0, shrinkB=0)

        for j in range(self.num_streamlines):
            x_start = 0.0
            y_max = 0.0
            y_min = 0.0
            plt.figure(num=1, clear=True)
            for i in range(1, len(blade_rows) - 1):
                row = blade_rows[i]
                x_end = x_start + row.Vm.mean()
                dx = x_end - x_start

                Vt = row.Vt[j]
                Wt = row.Wt[j]
                U = row.U[j]

                y_max = max(y_max, Vt, Wt)
                y_min = min(y_min, Vt, Wt)

                # V
                plt.annotate("", xy=(x_end, Vt), xytext=(x_start, 0), arrowprops=prop)
                plt.text((x_start + x_end) / 2, Vt / 2 * 1.1, "V", fontdict={"fontsize": "xx-large"})

                # W
                plt.annotate("", xy=(x_end, Wt), xytext=(x_start, 0), arrowprops=prop)
                plt.text((x_start + x_end) / 2, Wt / 2 * 1.1, "W", fontdict={"fontsize": "xx-large"})

                if abs(Vt) > abs(Wt):
                    plt.annotate("", xy=(x_end, Wt), xytext=(x_end, 0), arrowprops=prop)  # Wt
                    plt.text(x_end + dx * 0.1, Wt / 2, "Wt", fontdict={"fontsize": "xx-large"})

                    plt.annotate("", xy=(x_end, U + Wt), xytext=(x_end, Wt), arrowprops=prop)  # U
                    plt.text(x_end + dx * 0.1, (Wt + U) / 2, "U", fontdict={"fontsize": "xx-large"})
                else:
                    plt.annotate("", xy=(x_end, Vt), xytext=(x_end, 0), arrowprops=prop)  # Vt
                    plt.text(x_end + dx * 0.1, Vt / 2, "Vt", fontdict={"fontsize": "xx-large"})

                    plt.annotate("", xy=(x_end, Wt + U), xytext=(x_end, Wt), arrowprops=prop)  # U
                    plt.text(x_end + dx * 0.1, Wt + U / 2, "U", fontdict={"fontsize": "xx-large"})

                y = y_min if -np.sign(Vt) > 0 else y_max
                plt.text((x_start + x_end) / 2, -np.sign(Vt) * y * 0.95, row.row_type.name, fontdict={"fontsize": "xx-large"})
                x_start += row.Vm[j]
                plt.axis([0, x_end + dx, y_min, y_max])
            plt.ylabel("Tangental Velocity [m/s]")
            plt.xlabel("Vm [m/s]")
            plt.title(f"Velocity Triangles for Streamline {j}")
            plt.savefig(f"streamline_{j:04d}.png", transparent=False, dpi=150)

    def save_convergence_history(self, filename: str = "convergence_history.jsonl") -> None:
        """Save convergence history to JSONL file.

        Writes the convergence history collected during solve() to a JSON Lines file,
        where each line is a JSON object representing one iteration.

        Args:
            filename: Output JSONL file path (default: "convergence_history.jsonl")

        Returns:
            None. Writes JSONL file to specified path.

        Example:
            >>> spool.solve()
            >>> spool.save_convergence_history("compressor_convergence.jsonl")
        """
        import json
        from pathlib import Path

        output_path = Path(filename)
        with open(output_path, 'w') as f:
            for entry in self.convergence_history:
                f.write(json.dumps(entry) + '\n')
        print(f"Convergence history saved to {output_path}")

    def plot_convergence(self, save_to_file: Optional[Union[bool, str]] = None) -> None:
        """Plot convergence history showing massflow error vs iteration.

        Displays a semi-log plot of the massflow standard deviation error across
        iterations. If convergence history is empty, warns user.

        Args:
            save_to_file: If True, saves to "convergence.png". If string, saves to that filename.
                         If None/False, displays plot without saving.

        Returns:
            None. Either displays plot or saves to file.

        Example:
            >>> spool.solve()
            >>> spool.plot_convergence()  # Display plot
            >>> spool.plot_convergence(save_to_file=True)  # Save to convergence.png
            >>> spool.plot_convergence(save_to_file="my_convergence.png")  # Save to custom file
        """
        if not self.convergence_history:
            print("Warning: No convergence history available. Run solve() first.")
            return

        iterations = [entry['iteration'] for entry in self.convergence_history]
        massflow_std = [entry['massflow_std'] for entry in self.convergence_history]
        relative_change = [entry['relative_change'] for entry in self.convergence_history]

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

        # Plot massflow std deviation
        ax1.semilogy(iterations, massflow_std, 'o-', linewidth=2, markersize=8)
        ax1.set_xlabel('Iteration', fontsize=12)
        ax1.set_ylabel('Massflow Std Dev [kg/s]', fontsize=12)
        ax1.set_title('Convergence History: Massflow Standard Deviation', fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3)

        # Plot relative change
        ax2.semilogy(iterations, relative_change, 's-', color='orange', linewidth=2, markersize=8)
        ax2.set_xlabel('Iteration', fontsize=12)
        ax2.set_ylabel('Relative Change', fontsize=12)
        ax2.set_title('Convergence History: Relative Change', fontsize=14, fontweight='bold')
        ax2.axhline(y=0.05, color='r', linestyle='--', label='Convergence Threshold (0.05)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()

        if save_to_file:
            filename = "convergence.png" if save_to_file is True else str(save_to_file)
            plt.savefig(filename, dpi=150, bbox_inches='tight')
            print(f"Convergence plot saved to {filename}")
        else:
            plt.show()


def outlet_pressure(percents: List[float], inletP0: float, outletP: float) -> npt.NDArray:
    """Linearly interpolate total pressure values along the spool."""
    percents_arr = convert_to_ndarray(percents)
    return inletP0 + (outletP - inletP0) * percents_arr


def match_massflow_objective(exit_angle: float, index: int, row: BladeRow, upstream: BladeRow, downstream: Optional[BladeRow] = None, fluid: Optional[Solution] = None) -> float:
    """Objective for adjusting exit angle to match a target massflow slice."""
    if row.row_type not in (RowType.Rotor, RowType.Stator):
        return 0.0

    lt = getattr(row, "loss_function", None)
    loss_type = getattr(lt, "loss_type", None)

    if loss_type == LossType.Pressure and callable(lt):
        row.Yp = lt(row, upstream)  # type: ignore[arg-type]

    if row.row_type == RowType.Rotor:
        row.beta2[index] = np.radians(exit_angle)
        rotor_calc(row, upstream)
    elif row.row_type == RowType.Stator:
        row.alpha2[index] = np.radians(exit_angle)
        stator_calc(row, upstream)

    if fluid is not None:
        compute_gas_constants(upstream, fluid)
        compute_gas_constants(row, fluid)

    compute_massflow(row)
    compute_power(row, upstream)

    # drive radial distribution of massflow linearly by index using upstream total as target
    target_total = None
    for candidate in ("total_massflow_no_coolant", "total_massflow"):
        val = getattr(upstream, candidate, None)
        if val is not None and val != 0:
            target_total = val
            break
    if target_total is None:
        target_total = row.total_massflow if getattr(row, "total_massflow", 0) != 0 else row.massflow[-1]

    target = target_total * index / max(len(row.massflow) - 1, 1)
    return float(np.abs(target - row.massflow[index]))
