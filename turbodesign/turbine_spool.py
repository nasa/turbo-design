# type: ignore[arg-type, reportUnknownArgumentType]
from __future__ import annotations

from multiprocessing import Value
import stat
from typing import Dict, List, Union, Optional
import json

import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt

from cantera.composite import Solution
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar, fmin_slsqp
from sympy import true

# --- Project-local imports
from .bladerow import BladeRow, interpolate_streamline_quantities
from .enums import RowType, MassflowConstraint, LossType, PassageType
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
from .flow_math import compute_massflow, compute_streamline_areas
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
    massflow_constraint: MassflowConstraint
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
        massflow_constraint: MassflowConstraint = MassflowConstraint.AngleMatch,
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
            massflow_constraint: AngleMatch (adjust turning) or PressureBalance (radial eq).
        """
        self.passage = passage
        self.massflow = massflow
        self.num_streamlines = num_streamlines
        self._fluid = fluid if fluid is not None else Solution("air.yaml")
        self.massflow_constraint = massflow_constraint
        self.rpm = rpm

        self.inlet = inlet
        self.outlet = outlet
        if not self.outlet.static_defined:
            assert "Outlet needs to be statically defined for turbine calculation"
        self.rows = rows
        self.t_streamline = np.zeros((10,), dtype=float)
        self._adjust_streamlines = True

        # Assign IDs, RPMs, and axial chords where appropriate
        for i, br in enumerate(self._all_rows()):
            br.id = i
            if not isinstance(br, (Inlet, Outlet)):
                br.rpm = rpm
                br.axial_chord = br.hub_location * self.passage.hub_length

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

            t_radial = np.linspace(0, 1, self.num_streamlines)
            self.calculate_streamline_curvature(row, t_radial)

            # Ensure a loss model exists on blade rows
            if not isinstance(row, (Inlet, Outlet)) and row.loss_function is None:
                row.loss_function = TD2()

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
        if row.num_blades and row.chord != 0:
            mean_r = float(row.r.mean())
            pitch = 2 * np.pi * mean_r / row.num_blades
            row.pitch_to_chord = pitch / row.chord

    def solve_for_static_pressure(self,upstream:BladeRow,row:BladeRow):
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
        Is_static_defined = self.outlet.static_defined

        # Inlet
        W0 = self.massflow
        inlet = self.inlet
        if self.fluid:
            inlet.__initialize_fluid__(self.fluid)  # type: ignore[arg-type]
        else:
            inlet.__initialize_fluid__(  # type: ignore[call-arg]
                R=blade_rows[1].R,
                gamma=blade_rows[1].gamma,
                Cp=blade_rows[1].Cp,
            )

        inlet.total_massflow = W0
        inlet.total_massflow_no_coolant = W0
        inlet.massflow = np.linspace(0, 1, self.num_streamlines) * W0

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
            row.massflow = np.linspace(0, 1, self.num_streamlines) * row.total_massflow

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

        if self.massflow_constraint == MassflowConstraint.AngleMatch:
            self._angle_match()
        elif self.massflow_constraint == MassflowConstraint.PressureBalance:
            self._balance_pressure()

    # ------------------------------
    # Massflow matching/balancing
    # ------------------------------
    def _angle_match(self) -> None:
        """Match massflow between streamtubes by tweaking exit angles."""
        blade_rows = self._all_rows()
        for _ in range(3):
            for i, row in enumerate(blade_rows):
                upstream = blade_rows[i - 1] if i > 0 else blade_rows[i]
                downstream = blade_rows[i + 1] if i < len(blade_rows) - 1 else None

                if row.row_type == RowType.Stator:
                    bounds = [0, 80]
                elif row.row_type == RowType.Rotor:
                    bounds = [-80, 0]
                else:
                    bounds = [0, 0]

                if row.row_type != RowType.Inlet:
                    for j in range(1, self.num_streamlines):
                        res = minimize_scalar(
                            massflow_loss_function,
                            bounds=bounds,
                            args=(j, row, upstream, downstream),
                            tol=1e-3,
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

            adjust_streamlines(blade_rows, self.passage)
        compute_reynolds(blade_rows, self.passage)

    @staticmethod
    def __massflow_std__(blade_rows: List[BladeRow]) -> float:
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
            static_defined = self.outlet.static_defined
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
                                        calculate_vm=True,static_defined=static_defined)
                                row = radeq(row, upstream, downstream)
                                compute_gas_constants(row, self.fluid)
                                rotor_calc(row, upstream, 
                                        calculate_vm=False,static_defined=static_defined)
                            elif row.row_type == RowType.Stator:
                                stator_calc(row, upstream, downstream, 
                                            calculate_vm=True,static_defined=static_defined)
                                row = radeq(row, upstream, downstream)
                                compute_gas_constants(row, self.fluid)
                                stator_calc(row, upstream, downstream, 
                                            calculate_vm=False,static_defined=static_defined)
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
                            row = radeq(row,upstream) 
                            compute_gas_constants(row,self.fluid)
                            stator_calc(row,upstream,downstream,calculate_vm=False)
                        compute_gas_constants(row,self.fluid)
                        compute_massflow(row)
                        compute_power(row,upstream)
            print(x0)
            return self.__massflow_std__(rows[1:-1])

        pressure_ratio_ranges: List[tuple] = []
        pressure_ratio_guess: List[float] = []
        for i in range(1, len(rows) - 2):
            bounds = tuple(float(v) for v in rows[i].inlet_to_outlet_pratio)
            pressure_ratio_ranges.append(bounds)
            pressure_ratio_guess.append(float(np.mean(bounds)))

        if not self.outlet.static_defined:
            raise ValueError("For turbine calculations, please define outlet using init_static")
        
        print("Looping to converge massflow")
        past_err = -100.0
        loop_iter = 0
        err = 1e-3
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
            self.inlet.massflow = (np.linspace(0, 1, self.num_streamlines) * rows[1].total_massflow_no_coolant)
            self.inlet.total_massflow_no_coolant = rows[1].total_massflow_no_coolant
            self.inlet.total_massflow = rows[1].total_massflow_no_coolant
            self.inlet.calculated_massflow = self.inlet.total_massflow_no_coolant
            inlet_calc(self.inlet)

            if self.adjust_streamlines:
                adjust_streamlines(rows[:-1], self.passage)

            self.outlet.transfer_quantities(rows[-2])  # outlet
            self.outlet.P = self.outlet.get_static_pressure(self.outlet.percent_hub_shroud)

            past_err = err
            err = self.__massflow_std__(rows)
            loop_iter += 1
            print(f"Loop {loop_iter} massflow convergenced error:{err}")

        compute_reynolds(rows, self.passage)

    # ------------------------------
    # Export / Plotting
    # ------------------------------
    def export_properties(self, filename: str = "turbine_spool.json") -> None:
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
        x_streamline = np.zeros((self.num_streamlines, len(blade_rows)))
        r_streamline = np.zeros((self.num_streamlines, len(blade_rows)))
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


# ------------------------------
# Helper functions (kept module-level)
# ------------------------------



def massflow_loss_function(
    exit_angle: float,
    index: int,
    row: BladeRow,
    upstream: BladeRow,
    downstream: Optional[BladeRow] = None,
    fluid: Optional[Solution] = None,
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
            upstream = compute_gas_constants(upstream, fluid)
            row = compute_gas_constants(row, fluid)
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
        T3_is = upstream.T0 * (1 / row.P0_P) ** ((row.gamma - 1) / row.gamma)
        a = np.sqrt(row.gamma * row.R * T3_is)
        T03_is = T3_is * (1 + (row.gamma - 1) / 2 * (row.V / a) ** 2)
        row.eta_total = (upstream.T0.mean() - row.T0.mean()) / (upstream.T0.mean() - T03_is.mean())

    # drive radial distribution of massflow linearly by index
    target = row.total_massflow * index / (len(row.massflow) - 1)
    return float(np.abs(target - row.massflow[index]))


def step_pressures(percents: List[float], inletP0: float, outletP: float) -> npt.NDArray:
    """Map a list of percents [0..1] to each row's outlet static pressure."""
    percents_arr = convert_to_ndarray(percents)
    Ps = np.zeros((len(percents_arr),))
    for i in range(len(percents_arr)):
        Ps[i] = float(interp1d((0, 1), (inletP0, outletP))(percents_arr[i]))
        inletP0 = Ps[i]
    return Ps
