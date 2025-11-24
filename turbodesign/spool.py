# type: ignore[arg-type, reportUnknownArgumentType]
from __future__ import annotations

from typing import Dict, List, Union, Optional
import json
import copy

import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt

from cantera.composite import Solution
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar, fmin_slsqp

# --- Project-local imports
from .bladerow import BladeRow, interpolate_streamline_radii
from .enums import RowType, MassflowConstraint, LossType, PassageType
from .loss.turbine import TD2
from .passage import Passage
from .inlet import Inlet
from .outlet import Outlet
from .td_math import (
    inlet_calc,
    rotor_calc,
    stator_calc,
    compute_massflow,
    compute_power,
    compute_gas_constants,
    compute_reynolds,
)
from .solve_radeq import adjust_streamlines, radeq
from pyturbo.helper import line2D, convert_to_ndarray


class Spool:
    """Used with either compressor or turbines

    This class encapsulates both the generic geometry/plotting utilities
    from the original base *Spool* and the turbine-solving logic that lived
    in *TurbineSpool*.

    Notes on differences vs. the two-class design:
    - `field(default_factory=...)` was previously used on a non-dataclass attribute
      (`t_streamline`). Here it's handled in `__init__` to avoid a silent bug.
    - `fluid` defaults to `Solution('air.yaml')` if not provided.
    - All turbine-specific methods (initialize/solve/massflow balancing/etc.) are
      preserved here. If you ever add a *CompressorSpool* in the future, consider
      splitting turbine/compressor behaviors behind a strategy/solver object.
    """

    # Class-level defaults (avoid mutable defaults here!)
    blade_rows: List[BladeRow]
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
        rows: List[BladeRow],
        num_streamlines: int = 3,
        fluid: Optional[Solution] = None,
        rpm: float = -1,
        massflow_constraint: MassflowConstraint = MassflowConstraint.MatchMassFlow,
    ) -> None:
        """Initialize a (turbine) spool

        Args:
            passage: Passage defining hub and shroud
            massflow: massflow at spool inlet
            rows: List of blade rows (include inlet at [0] and outlet at [-1])
            num_streamlines: number of streamlines used through the meridional passage
            fluid: cantera gas solution; defaults to air.yaml if None
            rpm: RPM for the entire spool. Individual rows can override later.
            massflow_constraint: MatchMassFlow (adjust turning) or BalanceMassFlow (radial eq).
        """
        self.passage = passage
        self.massflow = massflow
        self.blade_rows = rows
        self.num_streamlines = num_streamlines
        self._fluid = fluid if fluid is not None else Solution("air.yaml")
        self.massflow_constraint = massflow_constraint
        self.rpm = rpm

        # Previously this used dataclasses.field on a non-dataclass; do it explicitly
        self.t_streamline = np.zeros((10,), dtype=float)
        self._adjust_streamlines = True

        # Assign IDs, RPMs, and axial chords where appropriate
        for i, br in enumerate(self.blade_rows):
            br.id = i
            if not isinstance(br, (Inlet, Outlet)):
                br.rpm = rpm
                br.axial_chord = br.hub_location * self.passage.hub_length

        # Propagate initial fluid to rows
        for br in self.blade_rows:
            br.fluid = self._fluid

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
        self.blade_rows[index].rpm = rpm

    def set_blade_row_type(self, blade_row_index: int, rowType: RowType) -> None:
        self.blade_rows[blade_row_index].row_type = rowType

    def set_blade_row_exit_angles(
        self,
        radius: Dict[int, List[float]],
        beta: Dict[int, List[float]],
        IsSupersonic: bool = False,
    ) -> None:
        """Set intended exit flow angles for rows (useful when geometry is fixed)."""
        for k, v in radius.items():
            self.blade_rows[k].radii_geom = v
        for k, v in beta.items():
            self.blade_rows[k].beta_geom = v
            self.blade_rows[k].beta_fixed = True
        for br in self.blade_rows:
            br.solution_type = (
                SolutionType.supersonic if IsSupersonic else SolutionType.subsonic
            )

    # ------------------------------
    # Streamline setup/geometry
    # ------------------------------
    def initialize_streamlines(self) -> None:
        """Initialize streamline storage per row and compute curvature."""
        for row in self.blade_rows:
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
        self, row: BladeRow, t_radial: Union[List[float], npt.NDArray]
    ) -> None:
        for i, tr in enumerate(t_radial):
            t_s, x_s, r_s = self.passage.get_streamline(tr)
            phi, rm, r = self.passage.streamline_curvature(x_s, r_s)
            row.phi[i] = float(interp1d(t_s, phi)(row.hub_location))
            row.rm[i] = float(interp1d(t_s, rm)(row.hub_location))
            row.r[i] = float(interp1d(t_s, r)(row.hub_location))
            row.m[i] = float(
                interp1d(t_s, self.passage.get_m(tr, resolution=len(t_s)))(row.hub_location)
            )

    # ------------------------------
    # initialization/solve
    # ------------------------------
    def initialize(self) -> None:
        """Initialize massflow and thermodynamic state through rows (turbines)."""
        Is_static_defined = self.blade_rows[-1].static_defined

        # Inlet
        W0 = self.massflow
        inlet: Inlet = self.blade_rows[0]  # type: ignore[assignment]
        if self.fluid:
            inlet.__initialize_fluid__(self.fluid)  # type: ignore[arg-type]
        else:
            inlet.__initialize_fluid__(  # type: ignore[call-arg]
                R=self.blade_rows[1].R,
                gamma=self.blade_rows[1].gamma,
                Cp=self.blade_rows[1].Cp,
            )

        inlet.total_massflow = W0
        inlet.total_massflow_no_coolant = W0
        inlet.massflow = np.linspace(0, 1, self.num_streamlines) * W0

        inlet.__interpolate_quantities__(self.num_streamlines)  # type: ignore[attr-defined]
        inlet.__initialize_velocity__(self.passage, self.num_streamlines)  # type: ignore[attr-defined]
        interpolate_streamline_radii(inlet, self.passage, self.num_streamlines)

        compute_gas_constants(inlet, self.fluid)
        inlet_calc(inlet)

        for row in self.blade_rows:
            interpolate_streamline_radii(row, self.passage, self.num_streamlines)

        outlet: Outlet = self.blade_rows[-1]  # type: ignore[assignment]
        for j in range(self.num_streamlines):
            P0 = inlet.get_total_pressure(inlet.percent_hub_shroud[j])  # type: ignore[attr-defined]
            percents = np.zeros(shape=(len(self.blade_rows) - 2)) + 0.3
            percents[-1] = 1
            if Is_static_defined:
                Ps_range = outlet_pressure(percents=percents, inletP0=inlet.P0[j], outletP=outlet.P[j])
                for i in range(1, len(self.blade_rows) - 1):
                    self.blade_rows[i].P[j] = Ps_range[i - 1]
            else:
                P0_range = outlet_pressure(percents=percents, inletP0=inlet.P0[j], outletP=outlet.P[j])
                for i in range(1, len(self.blade_rows) - 1):
                    self.blade_rows[i].P0[j] = P0_range[i - 1]
                    

        # Pass T0, P0 to downstream rows
        for i in range(1, len(self.blade_rows) - 1):
            upstream = self.blade_rows[i - 1]
            downstream = self.blade_rows[i + 1] if i + 1 < len(self.blade_rows) else None

            row = self.blade_rows[i]
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
            row.massflow = np.linspace(0, 1, self.num_streamlines) * row.total_massflow

            # Pass gas constants
            row.rho = upstream.rho
            row.gamma = upstream.gamma
            row.R = upstream.R

            if row.row_type == RowType.Stator:
                stator_calc(row, upstream, downstream,Is_static_defined)  # type: ignore[arg-type]
                compute_massflow(row)
            elif row.row_type == RowType.Rotor:
                rotor_calc(row, upstream,Is_static_defined)
                compute_massflow(row)
                compute_power(row, upstream)

    def solve(self) -> None:
        """Solve for exit angles/pressures to satisfy chosen massflow constraint."""
        self.initialize_streamlines()
        self.initialize()

        if self.massflow_constraint == MassflowConstraint.MatchMassFlow:
            self.__match_massflow()
        elif self.massflow_constraint == MassflowConstraint.BalanceMassFlow:
            self.__balance_massflow()

    # ------------------------------
    # Massflow matching/balancing
    # ------------------------------
    def __match_massflow(self) -> None:
        """Match massflow between streamtubes by tweaking exit angles."""
        for _ in range(3):
            for i, row in enumerate(self.blade_rows):
                upstream = self.blade_rows[i - 1] if i > 0 else self.blade_rows[i]
                downstream = self.blade_rows[i + 1] if i < len(self.blade_rows) - 1 else None

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

            adjust_streamlines(self.blade_rows, self.passage)
        compute_reynolds(self.blade_rows, self.passage)

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

    def __balance_massflow(self) -> None:
        """Balance massflow between rows using radial equilibrium."""

        def balance_massflows(
            x0: List[float],
            blade_rows: List[BladeRow],
            P0: npt.NDArray,
            P: npt.NDArray,
            balance_mean_pressure: bool = True,
        ) -> float:
            if balance_mean_pressure:
                for j in range(self.num_streamlines):
                    Ps = outlet_pressure(x0, P0[j], P[j])
                    for i in range(1, len(blade_rows) - 2):
                        blade_rows[i].P[j] = float(Ps[i - 1])
                blade_rows[-2].P = P
            else:
                for i in range(1, len(blade_rows) - 1):
                    for j in range(self.num_streamlines):
                        blade_rows[i].P[j] = P[j] * x0[(i - 1) * self.num_streamlines + j]

            calculate_massflows(blade_rows, True, self.fluid)
            return self.__massflow_std__(blade_rows[1:-1])

        outlet_P = []
        outlet_P_guess = []
        for i in range(1, len(self.blade_rows) - 2):
            outlet_P.append(self.blade_rows[i].inlet_to_outlet_pratio)
            outlet_P_guess.append(np.mean(self.blade_rows[i].inlet_to_outlet_pratio))

        print("Looping to converge massflow")
        past_err = -100.0
        loop_iter = 0
        err = 1e-3
        while (np.abs((err - past_err) / err) > 0.05) and (loop_iter < 10):
            if len(outlet_P) == 1:
                res = minimize_scalar(
                    fun=balance_massflows,
                    args=(self.blade_rows, self.blade_rows[0].P0, self.blade_rows[-1].P),
                    bounds=outlet_P[0],
                    tol=1e-3,
                    method="bounded",
                )
                x = res.x
            else:
                x = fmin_slsqp(
                    func=balance_massflows,
                    args=(self.blade_rows, self.blade_rows[0].P0, self.blade_rows[-1].P),
                    bounds=outlet_P,
                    x0=outlet_P_guess,
                    epsilon=1e-3,
                    iter=100,
                )
                outlet_P_guess = x  # type: ignore[assignment]

            # Adjust inlet to match massflow found at first blade row
            self.blade_rows[0].massflow = (
                np.linspace(0, 1, self.num_streamlines)
                * self.blade_rows[1].total_massflow_no_coolant
            )
            self.blade_rows[0].total_massflow_no_coolant = self.blade_rows[1].total_massflow_no_coolant
            self.blade_rows[0].total_massflow = self.blade_rows[1].total_massflow_no_coolant
            self.blade_rows[0].calculated_massflow = self.blade_rows[0].total_massflow_no_coolant
            inlet_calc(self.blade_rows[0])

            if self.adjust_streamlines:
                adjust_streamlines(self.blade_rows[:-1], self.passage)

            self.blade_rows[-1].transfer_quantities(self.blade_rows[-2])  # outlet
            self.blade_rows[-1].P = self.blade_rows[-1].get_static_pressure(
                self.blade_rows[-1].percent_hub_shroud
            )

            past_err = err
            err = self.__massflow_std__(self.blade_rows)
            loop_iter += 1
            print(f"Loop {loop_iter} massflow convergenced error:{err}")

        compute_reynolds(self.blade_rows, self.passage)

    # ------------------------------
    # Export / Plotting
    # ------------------------------
    def export_properties(self, filename: str = "turbine_spool.json") -> None:
        blade_rows_out = []
        degree_of_reaction = []
        total_total_efficiency = []
        total_static_efficiency = []
        stage_loading = []
        euler_power = []
        enthalpy_power = []
        x_streamline = np.zeros((self.num_streamlines, len(self.blade_rows)))
        r_streamline = np.zeros((self.num_streamlines, len(self.blade_rows)))
        massflow = []

        for indx, row in enumerate(self.blade_rows):
            blade_rows_out.append(row.to_dict())
            if row.row_type == RowType.Rotor:
                degree_of_reaction.append(
                    (
                        (self.blade_rows[indx - 1].P - row.P)
                        / (self.blade_rows[indx - 2].P - row.P)
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

        Pratio_Total_Total = np.mean(self.blade_rows[0].P0 / self.blade_rows[-2].P0)
        Pratio_Total_Static = np.mean(self.blade_rows[0].P0 / self.blade_rows[-2].P)
        FlowFunction = (
            np.mean(massflow) * np.sqrt(self.blade_rows[0].T0) * self.blade_rows[0].P0 / 1000
        )
        CorrectedSpeed = self.rpm * np.pi / 30 / np.sqrt(self.blade_rows[0].T0.mean())
        EnergyFunction = (
            (self.blade_rows[0].T0 - self.blade_rows[-2].T0)
            * 0.5
            * (self.blade_rows[0].Cp + self.blade_rows[-2].Cp)
            / self.blade_rows[0].T0
        )
        EnergyFunction = np.mean(EnergyFunction)

        data = {
            "blade_rows": blade_rows_out,
            "massflow": float(np.mean(massflow)) if massflow else 0.0,
            "rpm": self.rpm,
            "r_streamline": r_streamline.tolist(),
            "x_streamline": x_streamline.tolist(),
            "rhub": self.passage.rhub_pts.tolist(),
            "rshroud": self.passage.rshroud_pts.tolist(),
            "xhub": self.passage.xhub_pts.tolist(),
            "xshroud": self.passage.xshroud_pts.tolist(),
            "num_streamlines": self.num_streamlines,
            "euler_power": euler_power,
            "enthalpy_power": enthalpy_power,
            "total-total_efficiency": total_total_efficiency,
            "total-static_efficiency": total_static_efficiency,
            "stage_loading": stage_loading,
            "degree_of_reaction": degree_of_reaction,
            "Pratio_Total_Total": float(Pratio_Total_Total),
            "Pratio_Total_Static": float(Pratio_Total_Static),
            "FlowFunction": float(FlowFunction),
            "CorrectedSpeed": float(CorrectedSpeed),
            "EnergyFunction": float(EnergyFunction),
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
        for i in range(len(self.blade_rows)):
            x_streamline[:, i] = self.blade_rows[i].x
            r_streamline[:, i] = self.blade_rows[i].r

        for i in range(1, len(self.blade_rows) - 1):
            plt.plot(x_streamline[:, i], r_streamline[:, i], "--b", linewidth=1.5)

        for i, row in enumerate(self.blade_rows):
            plt.plot(row.x, row.r, linestyle="dashed", linewidth=1.5, color="blue", alpha=0.4)
            plt.plot(x_streamline[:, i], r_streamline[:, i], "or")

            if i == 0:
                pass
            else:
                upstream = self.blade_rows[i - 1]
                if upstream.row_type == RowType.Inlet:
                    cut_line1, _, _ = self.passage.get_cutting_line(
                        (row.location * hub_length + (0.5 * row.blade_to_blade_gap * row.axial_chord) - row.axial_chord)
                        / hub_length
                    )
                else:
                    cut_line1, _, _ = self.passage.get_cutting_line(
                        (upstream.location * hub_length) / hub_length
                    )
                cut_line2, _, _ = self.passage.get_cutting_line(
                    (row.location * hub_length - (0.5 * row.blade_to_blade_gap * row.axial_chord)) / hub_length
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
        """Plot velocity triangles for each blade row (turbines)."""
        prop = dict(arrowstyle="-|>,head_width=0.4,head_length=0.8", shrinkA=0, shrinkB=0)

        for j in range(self.num_streamlines):
            x_start = 0.0
            y_max = 0.0
            y_min = 0.0
            plt.figure(num=1, clear=True)
            for i in range(1, len(self.blade_rows) - 1):
                row = self.blade_rows[i]
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

def calculate_massflows(
    blade_rows: List[BladeRow],
    calculate_vm: bool = False,
    fluid: Optional[Solution] = None,
) -> None:
    for i in range(1, len(blade_rows) - 1):
        row = blade_rows[i]
        upstream = blade_rows[i - 1] if i > 0 else blade_rows[i]
        downstream = blade_rows[i + 1]

        if row.row_type == RowType.Inlet:
            row.Yp = 0
        else:
            if row.loss_function.loss_type == LossType.Pressure:  # type: ignore[union-attr]
                row.Yp = row.loss_function(row, upstream)  # type: ignore[assignment]
                for _ in range(2):
                    if row.row_type == RowType.Rotor:
                        rotor_calc(row, upstream, calculate_vm=True)
                        row = radeq(row, upstream, downstream)
                        compute_gas_constants(row, fluid)
                        rotor_calc(row, upstream, calculate_vm=False)
                    elif row.row_type == RowType.Stator:
                        stator_calc(row, upstream, downstream, calculate_vm=True)
                        row = radeq(row, upstream, downstream)
                        compute_gas_constants(row, fluid)
                        stator_calc(row, upstream, downstream, calculate_vm=False)
                    compute_gas_constants(row, fluid)
                    compute_massflow(row)
                    compute_power(row, upstream)

            elif row.loss_function.loss_type == LossType.Enthalpy: 
                if row.row_type == RowType.Rotor:
                    row.Yp = 0
                    rotor_calc(row,upstream,calculate_vm=calculate_vm)
                    eta_total = float(row.loss_function(row,upstream))
                    def find_yp(Yp,row,upstream):
                        row.Yp = Yp
                        rotor_calc(row,upstream,calculate_vm=True)
                        row = radeq(row,upstream)
                        compute_gas_constants(row,fluid)
                        rotor_calc(row,upstream,calculate_vm=False)
                        return abs(row.eta_total - eta_total)
                    
                    res = minimize_scalar(find_yp,bounds=[0,0.6],args=(row,upstream))
                    row.Yp = res.x
                elif row.row_type == RowType.Stator:
                    row.Yp = 0
                    stator_calc(row,upstream,downstream,calculate_vm=True)
                    row = radeq(row,upstream) 
                    compute_gas_constants(row,fluid)
                    stator_calc(row,upstream,downstream,calculate_vm=False)
                compute_gas_constants(row,fluid)
                compute_massflow(row)
                compute_power(row,upstream)


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


def outlet_pressure(percents: List[float], inletP0: float, outletP: float) -> npt.NDArray:
    """Map a list of percents [0..1] to each row's outlet static pressure."""
    percents_arr = convert_to_ndarray(percents)
    Ps = np.zeros((len(percents_arr),))
    for i in range(len(percents_arr)):
        Ps[i] = float(interp1d((0, 1), (inletP0, outletP))(percents_arr[i]))
        inletP0 = Ps[i]
    return Ps
