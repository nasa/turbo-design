# type: ignore[arg-type, reportUnknownArgumentType]
from __future__ import annotations

from multiprocessing import Value
import stat
from turtle import down
from typing import Dict, List, Union, Optional
import json

import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt

from cantera.composite import Solution
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar, fmin_slsqp

# --- Project-local imports
from .bladerow import BladeRow, interpolate_streamline_quantities
from .enums import RowType, LossType
from .outlet import OutletType
from .loss.turbine import TD2
from .passage import Passage
from .inlet import Inlet
from .outlet import Outlet
from .turbine_math import (
    inlet_calc,
    rotor_calc,
    stator_calc,
    compute_power,
    compute_gas_constants,
    compute_reynolds,
)
from .flow_math import compute_massflow, compute_streamline_areas, compute_power
from .solve_radeq import adjust_streamlines, radeq
from pyturbo.helper import line2D, convert_to_ndarray


class TurbineSpool:
    """Used with turbines

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
    ) -> None:
        """Initialize a (turbine) spool

        Args:
            passage: Passage defining hub and shroud
            massflow: massflow at spool inlet
            inlet: Inlet object
            outlet: Outlet object
            rows: Blade rows between inlet and outlet (stators/rotors only)
            num_streamlines: number of streamlines used through the meridional passage
            fluid: cantera gas solution; defaults to air.yaml if None
            rpm: RPM for the entire spool. Individual rows can override later.
        """
        self.passage = passage
        self.massflow = massflow
        self.num_streamlines = num_streamlines
        self._fluid = fluid
        self.rpm = rpm

        self.inlet = inlet
        self.outlet = outlet
        if self.outlet.outlet_type != OutletType.static_pressure:
            assert "Outlet needs to be statically defined for turbine calculation"
        self.rows = rows
        self.t_streamline = np.zeros((10,), dtype=float)
        self._adjust_streamlines = True
        self.convergence_history: List[Dict] = []

        # Assign IDs, RPMs, and axial chords where appropriate
        for i, br in enumerate(self._all_rows()):
            br.id = i
            if not isinstance(br, (Inlet, Outlet)):
                br.rpm = rpm
                br.axial_chord = br.hub_location * self.passage.hub_length

        # Propagate initial fluid to rows
        if self._fluid is not None:
            for br in self._all_rows():
                br.fluid = self._fluid

    def _all_rows(self) -> List[BladeRow]:
        """Convenience to iterate inlet + interior rows + outlet."""
        return [self.inlet, *self.rows, self.outlet]

    @property
    def blade_rows(self) -> List[BladeRow]:
        """Backwards-compatible combined row list."""
        return self._all_rows()

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
        for br in self.blade_rows:
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
        self, row: BladeRow, t_hub_shroud: Union[List[float], npt.NDArray]
    ) -> None:
        """Calculates the streamline curvature

        Args:
            row (BladeRow):  current blade row 
            t_radial (Union[List[float], npt.NDArray]): percent along line from hub to shroud
        """
        for i, tr in enumerate(t_hub_shroud):
            t_s, x_s, r_s = self.passage.get_streamline(tr)
            phi, rm, r = self.passage.streamline_curvature(x_s, r_s)
            row.phi[i] = float(interp1d(t_s, phi)(row.hub_location))
            row.rm[i] = float(interp1d(t_s, rm)(row.hub_location))
            row.r[i] = float(interp1d(t_s, r)(row.hub_location))
            row.m[i] = float(
                interp1d(t_s, self.passage.get_m(tr, resolution=len(t_s)))(row.hub_location)
            )
        # Back-compute pitch_to_chord if blade count is specified and chord is nonzero
        chord = np.asarray(row.chord, dtype=float)
        mean_chord = float(np.mean(chord)) if chord.size else 0.0
        if row.num_blades and mean_chord != 0:
            mean_r = float(np.mean(row.r))
            pitch = 2 * np.pi * mean_r / row.num_blades
            row.pitch_to_chord = pitch / mean_chord

    def solve_for_static_pressure(self,upstream:BladeRow,row:BladeRow):
        """Solve for static pressure at blade row exit using isentropic flow relations.

        Uses massflow-area-Mach number relation to find static pressure from known
        total conditions. Attempts both subsonic and supersonic solutions and selects
        the subsonic solution.

        Args:
            upstream: Upstream blade row providing inlet conditions
            row: Current blade row where static pressure is being solved

        Returns:
            None. Updates row.M, row.T, and row.P in-place.
        """
        if row.row_type == RowType.Stator:
            b = row.total_area * row.P0 / np.sqrt(row.T0) * np.sqrt(row.gamma/row.R)
        else:
            b = row.total_area * row.P0R / np.sqrt(row.T0R) * np.sqrt(row.gamma/row.R)

        solve_for_M = upstream.total_massflow / b
        fun = lambda M : np.abs(solve_for_M - M*(1+(row.gamma-1)/2 * M**2) ** (-(row.gamma+1)/(2*(row.gamma-1))))
        M_subsonic = minimize_scalar(fun,0.1, bounds=[0,1])
        M_supersonic = minimize_scalar(fun,1.5, bounds=[1,5])
        row.M = M_subsonic
        if row.row_type == RowType.Stator:
            row.T = row.T0/IsenT(M_subsonic,row.gamma)
        else: 
            row.T = row.T0R/IsenT(M_subsonic,row.gamma)
        a = np.sqrt(row.T*row.gamma*row.R)
        row.P = row.total_massflow * row.R*row.T / (row.total_area * row.M * a) 
        # When total conditions are defined we calculate static pressure
        if row.row_type == RowType.Stator:
            row.P = upstream.P0 - (upstream.P0 - row.P0) / row.Yp 
        else:
            row.P = upstream.P0R - (upstream.P0R - row.P0R) / row.Yp 
    # ------------------------------
    # initialization/solve
    # ------------------------------
    def initialize(self) -> None:
        """Initialize massflow and thermodynamic state through rows (turbines)."""
        blade_rows = self._all_rows()
        Is_static_defined = (self.outlet.outlet_type == OutletType.static_pressure) or (self.outlet.outlet_type == OutletType.massflow_static_pressure)

        # Inlet
        W0 = self.massflow
        inlet = self.inlet
        if self.fluid:
            inlet.__initialize_fluid__(self.fluid)  # type: ignore[arg-type]
        elif inlet.gamma is not None:
            inlet.__initialize_fluid__(R=inlet.R, gamma=inlet.gamma, Cp=inlet.Cp)  # type: ignore[call-arg]
        elif blade_rows[1].gamma is not None:
            inlet.__initialize_fluid__(  # type: ignore[call-arg]
                R=blade_rows[1].R,
                gamma=blade_rows[1].gamma,
                Cp=blade_rows[1].Cp,
            )

        inlet.total_massflow = W0
        inlet.total_massflow_no_coolant = W0
        inlet.massflow = np.array([W0]) if self.num_streamlines == 1 else np.linspace(0, 1, self.num_streamlines) * W0

        inlet.__interpolate_quantities__(self.num_streamlines)  # type: ignore[attr-defined]
        inlet.__initialize_velocity__(self.passage, self.num_streamlines)  # type: ignore[attr-defined]
        interpolate_streamline_quantities(inlet, self.passage, self.num_streamlines)

        inlet_calc(inlet)

        for i,row in enumerate(blade_rows):
            interpolate_streamline_quantities(row, self.passage, self.num_streamlines)
        
        outlet = self.outlet
        for j in range(self.num_streamlines):
            percents = np.zeros(shape=(len(blade_rows) - 2)) + 0.3
            percents[-1] = 1
            if Is_static_defined:
                Ps_range = step_pressures(percents=percents, inletP0=inlet.P0[j], outletP=outlet.P[j])
                for i in range(1, len(blade_rows) - 1):
                    blade_rows[i].P[j] = Ps_range[i - 1]
            else:
                P0_range = step_pressures(percents=percents, inletP0=inlet.P0[j], outletP=outlet.P0[j])
                for i in range(1, len(blade_rows) - 1):
                    if blade_rows[i].row_type == RowType.Stator:
                        blade_rows[i].P0[j] = P0_range[i - 1]
                    else:
                        blade_rows[i].P0R[j] = P0_range[i - 1]
                
        # Pass T0, P0 to downstream rows
        for i in range(1, len(blade_rows) - 1):
            upstream = blade_rows[i - 1]
            downstream = blade_rows[i + 1] if i + 1 < len(blade_rows) else None

            row = blade_rows[i]
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
                
            # Adjust for Coolant
            T0 = (W0 * upstream.Cp * upstream.T0 + W0c * Cpc * T0c) / (Cpc * W0c + upstream.Cp * W0)
            # P0 = (W0 * upstream.Cp * upstream.P0 + W0c * Cpc * P0c) / (Cpc * W0c + upstream.Cp * W0)
            Cp = (W0 * upstream.Cp + W0c * Cpc) / (W0c + W0) if (W0c + W0) != 0 else upstream.Cp
            # Adjust for power 
            if row.row_type == RowType.Rotor:
                T0 = T0 - row.power / (Cp * (W0 + W0c))

            W0 += W0c
            row.T0 = T0
            # row.P0 = P0
            row.Cp = Cp
            row.total_massflow = W0
            row.massflow = np.array([row.total_massflow]) if self.num_streamlines == 1 else np.linspace(0, 1, self.num_streamlines) * row.total_massflow

            # Pass gas constants
            row.rho = upstream.rho
            row.gamma = upstream.gamma
            row.R = upstream.R
            
            if row.loss_function.loss_type == LossType.Pressure:
                row.Yp = row.loss_function(row, upstream)
            elif row.loss_function.loss_type == LossType.Enthalpy: 
                row.Yp = 0
                    
            if row.row_type == RowType.Stator:
                stator_calc(row, upstream, downstream,True,Is_static_defined)  # type: ignore[arg-type]
                compute_massflow(row)
            elif row.row_type == RowType.Rotor:
                rotor_calc(row, upstream,True,Is_static_defined)
                compute_massflow(row)
                compute_power(row, upstream)

    def solve(self) -> None:
        """Solve for exit angles/pressures to satisfy chosen massflow constraint."""
        self.initialize_streamlines()
        self.initialize()

        if self.outlet.outlet_type == OutletType.massflow_static_pressure:
            print("Using angle matching mode: blade exit angles will be adjusted to match specified massflow")
            self._angle_match()
        else:
            print("Using pressure balance mode: blade exit angles are fixed, static pressures will be adjusted")
            self._balance_pressure()


    def total_power(self) -> float:
        """Return total turbine power extracted (sum over rotor rows)."""
        total = 0.0
        for row in self._all_rows():
            if getattr(row, "row_type", None) == RowType.Rotor:
                total += float(getattr(row, "power", 0.0) or 0.0)
        return total

    def solve_massflow_for_power(self, target_power: float, massflow_guess: Optional[float] = None, tol_rel: float = 1e-3, max_iter: int = 8, relax: float = 0.7, bounds: tuple[float, float] = (1e-6, 1e9)) -> tuple[float, float]:
        """Power-driven closure: iterate inlet massflow to hit a target turbine power.

        This uses a simple algebraic update (no additional nested optimizer):
            mdot_next = mdot_current * (P_target / P_current)

        The inner flow solution still uses the existing pressure-balance method to
        maintain a consistent massflow between rows for the current guess.

        Args:
            target_power: Desired turbine power [W]. Use a positive value for power extracted.
            massflow_guess: Optional starting guess for inlet massflow [kg/s]. Defaults to `self.massflow`.
            tol_rel: Relative tolerance on power error.
            max_iter: Maximum outer iterations.
            relax: Under-relaxation factor (0–1) for massflow updates.
            bounds: (lower, upper) bounds for massflow during updates.

        Returns:
            Tuple of (achieved_massflow_kg_s, achieved_power_W).
        """
        target = float(target_power)
        if target <= 0:
            raise ValueError("target_power must be positive for turbine power-based solve.")

        lower, upper = bounds
        if lower <= 0 or upper <= 0 or lower >= upper:
            raise ValueError("Massflow bounds must be positive and (lower < upper).")

        mdot = float(self.massflow if massflow_guess is None else massflow_guess)
        mdot = float(np.clip(mdot, lower, upper))

        # Temporarily store original outlet type and ensure pressure balance mode
        prev_outlet_type = self.outlet.outlet_type
        prev_massflow = getattr(self.outlet, 'total_massflow', None)
        self.outlet.outlet_type = OutletType.static_pressure

        try:
            for _ in range(max_iter):
                # Important: prevent a previous computed `row.power` from being treated as an input
                # in `initialize()` when power is not a design target.
                for r in self.rows:
                    if r.row_type == RowType.Rotor:
                        r.power = 0.0
                        r.power_mean = 0.0

                self.massflow = mdot
                self.solve()

                achieved_power = self.total_power()
                achieved_mdot = float(getattr(self._all_rows()[1], "total_massflow_no_coolant", mdot) or mdot)

                if achieved_power <= 0 or not np.isfinite(achieved_power):
                    raise ValueError(f"Non-physical power encountered during solve (power={achieved_power}).")

                err_rel = abs(achieved_power - target) / target
                if err_rel <= tol_rel:
                    return achieved_mdot, achieved_power

                mdot_update = achieved_mdot * (target / achieved_power)
                mdot = float(np.clip(relax * mdot_update + (1.0 - relax) * achieved_mdot, lower, upper))

            return float(getattr(self._all_rows()[1], "total_massflow_no_coolant", self.massflow) or self.massflow), self.total_power()
        finally:
            self.outlet.outlet_type = prev_outlet_type
            if prev_massflow is not None:
                self.outlet.total_massflow = prev_massflow

    # ------------------------------
    # Massflow matching/balancing
    # ------------------------------
    def _angle_match(self) -> None:
        """Match massflow between streamtubes by tweaking exit angles."""
        rows = self._all_rows()
        massflow_target = np.linspace(0,rows[-1].total_massflow,self.num_streamlines)

        self.convergence_history = []  # Reset convergence history
        past_err = -100.0
        loop_iter = 0
        err = 1e-3

        print("Looping to converge massflow (angle matching)")
        while (np.abs((err - past_err) / err) > 0.05) and (loop_iter < 10):
            for i in range(1,len(rows)-1):
                upstream = rows[i - 1] if i > 0 else rows[i]
                downstream = rows[i + 1] if i < len(rows) - 1 else None

                # Use custom massflow target if defined, otherwise use default
                if rows[i].massflow_target is not None:
                    current_massflow_target = rows[i].massflow_target
                else:
                    current_massflow_target = massflow_target

                if rows[i].row_type == RowType.Stator:
                    bounds = [0, 80]
                elif rows[i].row_type == RowType.Rotor:
                    bounds = [-80, 0]
                else:
                    bounds = [0, 0]

                for j in range(1, self.num_streamlines):
                    res = minimize_scalar(
                        massflow_loss_function,
                        bounds=bounds,
                        args=(j, rows[i], upstream, current_massflow_target[j], downstream),
                        options={'xatol': 1e-4},
                        method="bounded",
                    )
                    if rows[i].row_type == RowType.Rotor:
                        rows[i].beta2[j] = np.radians(res.x)
                        rows[i].beta2[0] = 1 / (len(rows[i].beta2) - 1) * rows[i].beta2[1:].sum()
                    elif rows[i].row_type == RowType.Stator:
                        rows[i].alpha2[j] = np.radians(res.x)
                        rows[i].alpha2[0] = 1 / (len(rows[i].alpha2) - 1) * rows[i].alpha2[1:].sum()
                compute_gas_constants(upstream, self.fluid)
                compute_gas_constants(rows[i], self.fluid)

            # Adjust inlet to match massflow found at first blade row
            target = rows[1].total_massflow_no_coolant
            self.inlet.massflow = np.array([target]) if self.num_streamlines == 1 else (np.linspace(0, 1, self.num_streamlines) * target)
            self.inlet.total_massflow_no_coolant = rows[1].total_massflow_no_coolant
            self.inlet.total_massflow = rows[1].total_massflow_no_coolant
            self.inlet.calculated_massflow = self.inlet.total_massflow_no_coolant
            inlet_calc(self.inlet)

            if self.adjust_streamlines:
                adjust_streamlines(rows, self.passage)

            # Track convergence history
            past_err = err
            err = self.__massflow_std__(rows[1:-1])
            loop_iter += 1

            self.convergence_history.append({
                'iteration': loop_iter,
                'massflow_std': float(err),
                'massflow_change': float(abs(err - past_err)),
                'relative_change': float(abs((err - past_err) / max(err, 1e-6))),
                'massflow': float(rows[1].total_massflow_no_coolant)
            })
            print(f"Angle match iteration {loop_iter}, massflow std: {err:.6f}")

        compute_reynolds(rows, self.passage)

    @staticmethod
    def __massflow_std__(blade_rows: List[BladeRow]) -> float:
        """Calculate massflow standard deviation across blade rows.

        Computes the standard deviation of total massflow (without coolant) across
        all blade rows. Used as a convergence criterion for pressure balance and
        angle matching iterations. Warns if deviation exceeds 1.0 kg/s.

        Args:
            blade_rows: List of all blade rows (inlet, stators, rotors, outlet)

        Returns:
            float: Two times the standard deviation of massflow [kg/s]
        """
        total_massflow = []
        massflow_stage = []
        stage_ids = list({row.stage_id for row in blade_rows if row.stage_id >= 0})

        for row in blade_rows:
            total_massflow.append(row.total_massflow_no_coolant)
            sign = 1
            for s in stage_ids:
                for r in blade_rows:
                    if r.stage_id == s and r.row_type == RowType.Rotor:
                        massflow_stage.append(sign * r.total_massflow_no_coolant)
                        sign *= -1
            if len(stage_ids) % 2 == 1 and massflow_stage:
                massflow_stage.append(massflow_stage[-1] * sign)
        deviation = np.std(total_massflow) * 2
        if deviation > 1.0:
            print("high massflow deviation detected")
        return np.std(total_massflow) * 2

    def _balance_pressure(self) -> None:
        """Balance massflow between rows using radial equilibrium."""
        rows = self._all_rows()
        past_err = -100.0
        loop_iter = 0
        err = 1e-3
        self.convergence_history = []  # Reset convergence history
        
        def balance_loop(
            x0: List[float],
            rows: List[BladeRow],
            P0: List[float],
            P_or_P0: List[float],
        ) -> float:
            """Runs through the calclulation and outputs the standard deviation of massflow

            Args:
                x0 (List[float]): Array of percent breakdown (P0 to P) or (P0 to P0_exit)
                rows (List[BladeRow]): _description_
                P0 (npt.NDArray): _description_
                P_or_P0 (npt.NDArray): _description_

            Returns:
                float: _description_
            """
            nonlocal err, past_err, loop_iter
            static_defined = (self.outlet.outlet_type == OutletType.static_pressure)
            P_exit = P_or_P0
            for j in range(self.num_streamlines):
                Ps_guess = step_pressures(x0, P0[j], P_exit[j])
                for i in range(1, len(rows) - 2):
                    rows[i].P[j] = float(Ps_guess[i - 1])
            rows[-2].P[:] = P_exit[-1]
            
            # Loop through massflow calculation for all rows
            for i in range(1, len(rows) - 1):
                row = rows[i]
                upstream = rows[i - 1] if i > 0 else rows[i]
                downstream = rows[i + 1]

                if row.row_type == RowType.Inlet:
                    row.Yp = 0
                else:
                    if row.loss_function.loss_type == LossType.Pressure:  # type: ignore[union-attr]
                        row.Yp = row.loss_function(row, upstream)  # type: ignore[assignment]
                        for _ in range(2):
                            if row.row_type == RowType.Rotor:
                                rotor_calc(row, upstream, 
                                        calculate_vm=True,outlet_type=OutletType.static_pressure if static_defined else OutletType.total_pressure)
                                if self.num_streamlines > 1:
                                    row = radeq(row, upstream, downstream)
                                    compute_gas_constants(row, self.fluid)
                                    rotor_calc(row, upstream, 
                                            calculate_vm=False,outlet_type=OutletType.static_pressure if static_defined else OutletType.total_pressure)
                            elif row.row_type == RowType.Stator:
                                stator_calc(row, upstream, downstream, 
                                            calculate_vm=True,outlet_type=OutletType.static_pressure if static_defined else OutletType.total_pressure)
                                if self.num_streamlines > 1:
                                    row = radeq(row, upstream, downstream)
                                    compute_gas_constants(row, self.fluid)
                                    stator_calc(row, upstream, downstream, 
                                                calculate_vm=False,outlet_type=OutletType.static_pressure if static_defined else OutletType.total_pressure)
                            compute_gas_constants(row, self.fluid)
                            compute_massflow(row)
                            compute_power(row, upstream)

                    elif row.loss_function.loss_type == LossType.Enthalpy: 
                        if row.row_type == RowType.Rotor:
                            row.Yp = 0
                            rotor_calc(row,upstream,calculate_vm=True)
                            eta_total = float(row.loss_function(row,upstream))
                            
                            def find_yp(Yp,row,upstream):
                                row.Yp = Yp
                                rotor_calc(row,upstream,calculate_vm=True)
                                row = radeq(row,upstream)
                                compute_gas_constants(row,self.fluid)
                                rotor_calc(row,upstream,calculate_vm=False)
                                return abs(row.eta_total - eta_total)
                            
                            res = minimize_scalar(find_yp,bounds=[0,0.6],args=(row,upstream))
                            row.Yp = res.x
                        elif row.row_type == RowType.Stator:
                            row.Yp = 0
                            stator_calc(row,upstream,downstream,calculate_vm=True)
                            if self.num_streamlines > 1:
                                row = radeq(row,upstream) 
                                compute_gas_constants(row,self.fluid)
                                stator_calc(row,upstream,downstream,calculate_vm=False)
                        compute_gas_constants(row,self.fluid)
                        compute_massflow(row)
                        compute_power(row,upstream)
            print(x0)
            
            past_err = err
            err = self.__massflow_std__(rows[1:-1])
            loop_iter += 1
            
            # Store convergence history
            self.convergence_history.append({
                'iteration': loop_iter,
                'massflow_std': float(err),
                'massflow_change': float(abs(err - past_err)),
                'relative_change': float(abs((err - past_err) / max(err, 1e-6))),
                'massflow': float(rows[1].total_massflow_no_coolant)
            })
            
            return self.__massflow_std__(rows[1:-1])

        pressure_ratio_ranges: List[tuple] = []
        pressure_ratio_guess: List[float] = []
        for i in range(1, len(rows) - 2):
            bounds = tuple(float(v) for v in rows[i].inlet_to_outlet_pratio)
            pressure_ratio_ranges.append(bounds)
            pressure_ratio_guess.append(float(np.mean(bounds)))

        if self.outlet.outlet_type != OutletType.static_pressure:
            raise ValueError("For turbine calculations, please define outlet using init_static")
        
        print("Looping to converge massflow")
        while (np.abs((err - past_err) / err) > 0.05) and (loop_iter < 10):
            if len(pressure_ratio_ranges) == 1: # Single stage, use minimize scalar 
                x = minimize_scalar(
                    fun=balance_loop,
                    args=(rows, self.inlet.P0, self.outlet.P),
                    bounds=pressure_ratio_ranges[0],
                    tol=1e-4,
                    method="bounded")
                print(x)
            else:   # Multiple stages, use slsqp
                x = fmin_slsqp(
                    func=balance_loop,
                    args=(rows, self.inlet.P0, self.outlet.P),
                    bounds=pressure_ratio_ranges,
                    x0=pressure_ratio_guess,
                    epsilon=1e-4,
                    iter=200)
                pressure_ratio_guess = x.tolist()
                
            # Adjust inlet to match massflow found at first blade row
            target = rows[1].total_massflow_no_coolant
            self.inlet.massflow = np.array([target]) if self.num_streamlines == 1 else (np.linspace(0, 1, self.num_streamlines) * target)
            self.inlet.total_massflow_no_coolant = rows[1].total_massflow_no_coolant
            self.inlet.total_massflow = rows[1].total_massflow_no_coolant
            self.inlet.calculated_massflow = self.inlet.total_massflow_no_coolant
            inlet_calc(self.inlet)

            if self.adjust_streamlines:
                adjust_streamlines(rows[:-1], self.passage)

            self.outlet.transfer_quantities(rows[-2])  # outlet
            self.outlet.P = self.outlet.get_static_pressure(self.outlet.percent_hub_shroud)


        compute_reynolds(rows, self.passage)

    # ------------------------------
    # Export / Plotting
    # ------------------------------
    def export_properties(self, filename: str = "turbine_spool.json") -> None:
        """Export turbine spool properties and blade row data to JSON file.

        Exports comprehensive turbine design data including blade row properties,
        streamline coordinates, efficiency metrics, degree of reaction, stage loading,
        and power calculations for each stage. Useful for post-processing and result
        archiving.

        Args:
            filename: Output JSON file path (default: "turbine_spool.json")

        Returns:
            None. Writes JSON file to specified path.

        Example:
            >>> spool.export_properties("eee_hpt_results.json")
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
        # Use scalarized inlet conditions to avoid shape mismatches with per-row massflow
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

        with open(filename, "w", encoding="utf-8") as f:
            json.dump(data, f, indent=4, cls=NumpyEncoder, ensure_ascii=False)

    def plot(self) -> None:
        """Plot hub/shroud and streamlines with improved labels and formatting."""
        blade_rows = self._all_rows()
        fig, ax = plt.subplots(1, 1, figsize=(16, 8), dpi=150)

        # Plot hub and shroud with thicker lines
        ax.plot(
            self.passage.xhub_pts,
            self.passage.rhub_pts,
            label="Hub",
            linestyle="solid",
            linewidth=3,
            color="black",
            zorder=10
        )
        ax.plot(
            self.passage.xshroud_pts,
            self.passage.rshroud_pts,
            label="Shroud",
            linestyle="solid",
            linewidth=3,
            color="black",
            zorder=10
        )

        hub_length = np.sum(
            np.sqrt(np.diff(self.passage.xhub_pts) ** 2 + np.diff(self.passage.rhub_pts) ** 2)
        )

        # Prepare streamline data
        x_streamline = np.zeros((self.num_streamlines, len(blade_rows)))
        r_streamline = np.zeros((self.num_streamlines, len(blade_rows)))
        for i in range(len(blade_rows)):
            x_streamline[:, i] = blade_rows[i].x
            r_streamline[:, i] = blade_rows[i].r

        # Plot streamlines connecting blade rows
        for i in range(1, len(blade_rows) - 1):
            ax.plot(x_streamline[:, i], r_streamline[:, i],
                   linestyle="--", linewidth=1.2, color="gray", alpha=0.6, zorder=1)

        # Track label positions to avoid overlaps
        label_positions = []

        for i, row in enumerate(blade_rows):
            # Plot blade row exit locations
            ax.plot(row.x, row.r, linestyle="none", marker="o",
                   markersize=6, color="red", alpha=0.7, zorder=5)

            # Label inlet
            if row.row_type == RowType.Inlet:
                x_pos = row.x.mean()
                r_pos = row.r.mean()
                ax.axvline(x=x_pos, color='green', linestyle=':', linewidth=2, alpha=0.7, zorder=2)
                ax.text(x_pos, self.passage.rshroud_pts.max() * 1.05, 'INLET',
                       fontsize=12, fontweight='bold', ha='center', va='bottom',
                       bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgreen', alpha=0.7))
                label_positions.append((x_pos, 'INLET'))

            # Plot blade rows with proper labels
            elif row.row_type in [RowType.Stator, RowType.Rotor]:
                if i > 0:
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

                    # Plot blade leading and trailing edges
                    if row.row_type == RowType.Stator:
                        color = 'purple'
                        label = f'Stator {row.stage_id + 1}'
                    else:
                        color = 'brown'
                        label = f'Rotor {row.stage_id + 1}'

                    x1, r1 = cut_line1.get_point(np.linspace(0, 1, 10))
                    ax.plot(x1, r1, color=color, linewidth=2.5, alpha=0.8, zorder=3)
                    x2, r2 = cut_line2.get_point(np.linspace(0, 1, 10))
                    ax.plot(x2, r2, color=color, linewidth=2.5, alpha=0.8, zorder=3)

                    # Mark exit location with vertical line
                    x_exit = row.x.mean()
                    ax.axvline(x=x_exit, color=color, linestyle='--',
                             linewidth=1.5, alpha=0.5, zorder=2)

                    # Add exit label at top
                    ax.text(x_exit, self.passage.rshroud_pts.max() * 1.02, f'{label} Exit',
                           fontsize=10, ha='center', va='bottom', rotation=0,
                           color=color, fontweight='bold')

            # Label outlet
            elif row.row_type == RowType.Outlet:
                x_pos = row.x.mean()
                ax.axvline(x=x_pos, color='blue', linestyle=':', linewidth=2, alpha=0.7, zorder=2)
                ax.text(x_pos, self.passage.rshroud_pts.max() * 1.05, 'OUTLET',
                       fontsize=12, fontweight='bold', ha='center', va='bottom',
                       bbox=dict(boxstyle='round,pad=0.5', facecolor='lightblue', alpha=0.7))

        # Formatting
        ax.set_xlabel('Axial Distance [m]', fontsize=13, fontweight='bold')
        ax.set_ylabel('Radial Distance [m]', fontsize=13, fontweight='bold')
        ax.set_title(f'Meridional View - {self.num_streamlines} Streamlines',
                    fontsize=14, fontweight='bold', pad=40)
        ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
        ax.legend(loc='upper left', fontsize=11, framealpha=0.9)
        ax.set_aspect('equal', adjustable='box')

        plt.tight_layout()
        plt.savefig("Meridional.png", transparent=False, dpi=200, bbox_inches='tight')
        plt.show()

    def plot_velocity_triangles(self) -> None:
        """Plot velocity triangles for each blade row with improved styling and annotations."""
        blade_rows = self._all_rows()

        # Define arrow properties for different velocity types
        prop_V = dict(arrowstyle="-|>,head_width=0.5,head_length=1.0",
                     shrinkA=0, shrinkB=0, color='blue', lw=2.5)
        prop_W = dict(arrowstyle="-|>,head_width=0.5,head_length=1.0",
                     shrinkA=0, shrinkB=0, color='red', lw=2.5)
        prop_U = dict(arrowstyle="-|>,head_width=0.5,head_length=1.0",
                     shrinkA=0, shrinkB=0, color='green', lw=2.5)
        prop_component = dict(arrowstyle="-|>,head_width=0.4,head_length=0.8",
                             shrinkA=0, shrinkB=0, color='gray', lw=1.5, linestyle='--')

        for j in range(self.num_streamlines):
            x_start = 0.0
            y_max = 0.0
            y_min = 0.0

            fig, ax = plt.subplots(1, 1, figsize=(14, 8), dpi=150)

            for i in range(1, len(blade_rows) - 1):
                row = blade_rows[i]
                x_end = x_start + row.Vm[j]
                dx = x_end - x_start

                Vt = row.Vt[j]
                Wt = row.Wt[j]
                U = row.U[j]
                Vm = row.Vm[j]

                y_max = max(y_max, Vt, Wt, U + Wt, U + Vt)
                y_min = min(y_min, Vt, Wt, 0)

                # Draw absolute velocity V (blue)
                ax.annotate("", xy=(x_end, Vt), xytext=(x_start, 0), arrowprops=prop_V, zorder=5)
                v_mag = np.sqrt(Vm**2 + Vt**2)
                ax.text((x_start + x_end) / 2, Vt / 2 + np.sign(Vt) * 15,
                       f"V={v_mag:.1f}", fontsize=12, fontweight='bold',
                       ha='center', color='blue',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.7))

                # Draw relative velocity W (red)
                ax.annotate("", xy=(x_end, Wt), xytext=(x_start, 0), arrowprops=prop_W, zorder=5)
                w_mag = np.sqrt(Vm**2 + Wt**2)
                ax.text((x_start + x_end) / 2, Wt / 2 - np.sign(Wt) * 15,
                       f"W={w_mag:.1f}", fontsize=12, fontweight='bold',
                       ha='center', color='red',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='lightcoral', alpha=0.7))

                # Draw velocity components and U
                if abs(Vt) > abs(Wt):
                    # Draw Wt component
                    ax.annotate("", xy=(x_end, Wt), xytext=(x_end, 0), arrowprops=prop_component, zorder=3)
                    ax.text(x_end + dx * 0.08, Wt / 2, f"Wt={Wt:.1f}",
                           fontsize=10, ha='left', color='gray')

                    # Draw U (blade speed)
                    ax.annotate("", xy=(x_end, U + Wt), xytext=(x_end, Wt), arrowprops=prop_U, zorder=4)
                    ax.text(x_end + dx * 0.08, (Wt + U + Wt) / 2, f"U={U:.1f}",
                           fontsize=11, ha='left', fontweight='bold', color='green',
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgreen', alpha=0.7))
                else:
                    # Draw Vt component
                    ax.annotate("", xy=(x_end, Vt), xytext=(x_end, 0), arrowprops=prop_component, zorder=3)
                    ax.text(x_end + dx * 0.08, Vt / 2, f"Vt={Vt:.1f}",
                           fontsize=10, ha='left', color='gray')

                    # Draw U (blade speed)
                    ax.annotate("", xy=(x_end, Wt), xytext=(x_end, Vt), arrowprops=prop_U, zorder=4)
                    ax.text(x_end + dx * 0.08, (Vt + Wt) / 2, f"U={U:.1f}",
                           fontsize=11, ha='left', fontweight='bold', color='green',
                           bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgreen', alpha=0.7))

                # Draw Vm component (dashed horizontal)
                ax.plot([x_start, x_end], [0, 0], 'k--', linewidth=1.5, alpha=0.5, zorder=2)
                ax.text((x_start + x_end) / 2, -5, f"Vm={Vm:.1f}",
                       fontsize=10, ha='center', va='top', color='black')

                # Add blade row label
                label_y = y_min - (y_max - y_min) * 0.15 if Vt > 0 else y_max + (y_max - y_min) * 0.15
                stage_label = f"{row.row_type.name} {row.stage_id + 1}"
                ax.text((x_start + x_end) / 2, label_y, stage_label,
                       fontsize=13, ha='center', fontweight='bold',
                       bbox=dict(boxstyle='round,pad=0.5',
                               facecolor='lightyellow' if row.row_type == RowType.Stator else 'lightcoral',
                               edgecolor='black', linewidth=2))

                # Add separation line between blade rows
                if i < len(blade_rows) - 2:
                    ax.axvline(x=x_end, color='gray', linestyle=':', linewidth=1, alpha=0.5, zorder=1)

                x_start = x_end

            # Formatting
            margin = (y_max - y_min) * 0.2
            ax.set_ylim([y_min - margin, y_max + margin])
            ax.set_xlim([0, x_end * 1.1])

            ax.set_ylabel('Tangential Velocity [m/s]', fontsize=13, fontweight='bold')
            ax.set_xlabel('Meridional Velocity Vm [m/s]', fontsize=13, fontweight='bold')
            ax.set_title(f'Velocity Triangles - Streamline {j} (r={blade_rows[1].r[j]:.4f} m)',
                        fontsize=14, fontweight='bold', pad=20)

            ax.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
            ax.axhline(y=0, color='black', linewidth=1.5, zorder=2)

            # Add legend
            from matplotlib.patches import FancyArrow
            legend_elements = [
                plt.Line2D([0], [0], color='blue', linewidth=2.5, label='V (Absolute Velocity)'),
                plt.Line2D([0], [0], color='red', linewidth=2.5, label='W (Relative Velocity)'),
                plt.Line2D([0], [0], color='green', linewidth=2.5, label='U (Blade Speed)')
            ]
            ax.legend(handles=legend_elements, loc='upper right', fontsize=10, framealpha=0.9)

            plt.tight_layout()
            plt.savefig(f"streamline_{j:04d}.png", transparent=False, dpi=200, bbox_inches='tight')
            plt.close()

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
            >>> spool.save_convergence_history("turbine_convergence.jsonl")
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
        ax1.set_xlabel('Iteration', fontsize=16)
        ax1.set_ylabel('2× Massflow Std Dev [kg/s]', fontsize=16)
        ax1.set_title('Convergence History: Massflow Standard Deviation', fontsize=14, fontweight='bold')
        ax1.grid(True, alpha=0.3)

        # Plot relative change
        ax2.semilogy(iterations, relative_change, 's-', color='orange', linewidth=2, markersize=8)
        ax2.set_xlabel('Iteration', fontsize=16)
        ax2.set_ylabel(r'Massflow Residual $\left|\frac{err_{n-1} - err_n}{err_n}\right|$', fontsize=16)
        ax2.set_title('Convergence History: Relative Error Change', fontsize=14, fontweight='bold')
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()

        if save_to_file:
            filename = "convergence.png" if save_to_file is True else str(save_to_file)
            plt.savefig(filename, dpi=150, bbox_inches='tight')
            print(f"Convergence plot saved to {filename}")
        else:
            plt.show()


# ------------------------------
# Helper functions (kept module-level)
# ------------------------------
def massflow_loss_function(
    exit_angle: float,
    index: int,
    row: BladeRow,
    upstream: BladeRow,
    massflow_target:float,
    downstream: Optional[BladeRow] = None,
    fluid: Optional[Solution] = None
) -> float:
    if row.row_type == RowType.Inlet:
        row.Yp = 0
    else:
        if row.loss_function.loss_type == LossType.Pressure:  # type: ignore[union-attr]
            row.Yp = row.loss_function(row, upstream)  # type: ignore[assignment]
            if row.row_type == RowType.Rotor:
                row.beta2[index] = np.radians(exit_angle)
                rotor_calc(row, upstream)
            elif row.row_type == RowType.Stator:
                row.alpha2[index] = np.radians(exit_angle)
                stator_calc(row, upstream, downstream)
            compute_gas_constants(upstream, fluid)
            compute_gas_constants(row, fluid)
        elif row.loss_function.loss_type == LossType.Enthalpy:  # type: ignore[union-attr]
            if row.row_type == RowType.Rotor:
                row.Yp = 0
                row.beta2[index] = np.radians(exit_angle)
                rotor_calc(row, upstream)
                T0_drop = row.loss_function(row, upstream)  # type: ignore[arg-type]
                T0_target = row.T0.mean() - T0_drop

                def find_yp(Yp):
                    row.Yp = Yp
                    rotor_calc(row, upstream)
                    compute_gas_constants(upstream, fluid)
                    compute_gas_constants(row, fluid)
                    return abs(row.T0.mean() - T0_target)

                res = minimize_scalar(find_yp, bounds=[0, 0.6], method="bounded")
                row.Yp = res.x
            elif row.row_type == RowType.Stator:
                row.Yp = 0
                row.alpha2[index] = np.radians(exit_angle)
                stator_calc(row, upstream, downstream)
                compute_gas_constants(upstream, fluid)
                compute_gas_constants(row, fluid)

    compute_massflow(row)
    compute_power(row, upstream)

    if row.row_type != RowType.Inlet:
        T03_is = upstream.T0 * (row.P0 / upstream.P0) ** ((row.gamma - 1) / row.gamma)
        row.eta_total = (upstream.T0.mean() - row.T0.mean()) / (upstream.T0.mean() - T03_is.mean())

    return float(np.abs(massflow_target - row.massflow[index]))


def step_pressures(percents: List[float], inletP0: float, outletP: float) -> npt.NDArray:
    """Map a list of percents [0..1] to each row's outlet static pressure."""
    percents_arr = convert_to_ndarray(percents)
    Ps = np.zeros((len(percents_arr),))
    for i in range(len(percents_arr)):
        Ps[i] = float(interp1d((0, 1), (inletP0, outletP))(percents_arr[i]))
        inletP0 = Ps[i]
    return Ps
