# type: ignore[arg-type, reportUnknownArgumentType]
from __future__ import annotations

from typing import Dict, List, Union, Optional
import json

import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt

from cantera.composite import Solution
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar, fmin_slsqp

from turbodesign import loss
from turbodesign.loss import losstype

# --- Project-local imports
from .bladerow import BladeRow, interpolate_streamline_radii
from .enums import RowType, MassflowConstraint, LossType, PassageType
from .loss.turbine import TD2
from .passage import Passage
from .inlet import Inlet
from .outlet import Outlet
from .td_math import (
    inlet_calc,
    compute_massflow,
    compute_power,
    compute_gas_constants,
    compute_reynolds,
    T0_coolant_weighted_average
)
from .solve_radeq import adjust_streamlines, radeq
from pyturbo.helper import line2D, convert_to_ndarray



def stator_calc(row:BladeRow,upstream:BladeRow,downstream:Optional[BladeRow]=None,calculate_vm:bool=True, static_defined:bool=False):
    """Given P0, T0, P, alpha2 of stator calculate all other quantities

    Usage:
        Set row.P0 = upstream.P0 - any pressure loss
        row.T0 = upstream.T0 - any cooling
        row.P = row.rp*(row.P0 - rotor.P) + rotor.P 
        Set alpha2 
        
    Args:
        row (BladeRow): Stator Row
        upstream (BladeRow): Stator or Rotor Row 
        downstream (BladeRow): Stator or Rotor Row. Defaults to None
        calculate_vm (bool): True to calculate the meridional velocity. False, do not calculate this and let radeq calculate it
        static_defined (bool): True if static conditions defined at the outlet. False if total conditions defined at outlet
    """
 
    # Static Pressure is assumed
    row.T0 = upstream.T0 - T0_coolant_weighted_average(row)
    loss_type = getattr(row.loss_function, "loss_type", None)
    
    if loss_type == LossType.Pressure:
        row.Yp = row.loss_function(row,upstream)
        stator_calculation(row.Yp)
    elif loss_type == LossType.Entropy: 
        desired_entropy_rise = convert_to_ndarray(row.loss_function(row,upstream))
        if len(desired_entropy_rise) == 1: # If entropy rise is a bulk value
            fun = lambda s: np.abs(desired_entropy_rise - stator_calculation(s))
            x = minimize_scalar(fun,bounds=[0.01,0.4])
            row.Yp = x
        elif len(desired_entropy_rise) > 1: # In case entropy rise is an array 
            fun = lambda s,j: np.abs(desired_entropy_rise[j] - stator_calculation(s)[j])
            for j in range(len(desired_entropy_rise)):
                x = minimize_scalar(fun,bounds=[0.01,0.4],args=(j))
                row.Yp[j] = x
            stator_calculation(row.Yp)
    else: # loss_type == LossType.Enthalpy: 
        # Enthalpy loss is typically calculated for a stage so I bulk the loss after the rotor. This is a really terrible way of doing loss but this is here for those legacy loss models. 
        row.Yp = 0
        stator_calculation(row.Yp)


        
    def stator_calculation(Yp:npt.NDArray):
        row.Yp = Yp
        if static_defined:
            row.P0 = upstream.P0 - row.Yp*(upstream.P0-row.P) # When static conditions are defined, use it to calculate P0
        else:
            row.P0 = row.P0_is - row.Yp*(upstream.P0 - upstream.P)
            
        if downstream is not None:
            row.P0_P = float((row.P0/downstream.P).mean())
            row.rp = ((row.P-downstream.P)/(upstream.P0-downstream.P)).mean()
            
        if calculate_vm:
            row.M = ((row.P0/row.P)**((row.gamma-1)/row.gamma) - 1) * 2/(row.gamma-1)
            row.M = np.sqrt(row.M)
            T0_T = (1+(row.gamma-1)/2 * row.M**2)
            row.T = row.T0/T0_T
            row.V = row.M*np.sqrt(row.gamma*row.R*row.T)
            row.Vm = row.V*np.cos(row.alpha2)
            row.Vx = row.Vm*np.cos(row.phi)
            row.Vr = row.Vm*np.sin(row.phi)
            row.Vt = row.Vm*np.tan(row.alpha2)
        else: # We know Vm, P0, T0, P
            row.Vx = row.Vm*np.cos(row.phi)
            row.Vr = row.Vm*np.sin(row.phi)
            row.Vt = row.Vm*np.tan(row.alpha2)
            row.V = np.sqrt(row.Vx**2 + row.Vr**2 + row.Vt**2)
            row.T = row.P/(row.R*row.rho)   # We know P, this is a guess
            row.M = row.V/np.sqrt(row.gamma*row.R*row.T)
            
        if upstream.row_type == RowType.Rotor:
            row.alpha1 = upstream.alpha2 # Upstream rotor absolute frame flow angle
        row.beta1 = upstream.beta2
        row.rho = row.P/(row.R*row.T)
        row.U = row.omega*row.r
        row.Wt = row.Vt-row.U
        row.P0_stator_inlet = upstream.P0
        row.entropy_rise = 0.5*(row.Cp+upstream.Cp)*np.log(row.T/upstream.T) - row.R * np.log(row.P/upstream.P)    
        return row.entropy_rise

def solve_for_mach(M:float, row:BladeRow,area:float, massflow:float):
    expo = -(row.gamma+1)/(2*(row.gamma-1))
    return np.abs(massflow - area * row.P0/row.T0 * np.sqrt(row.gamma/row.R) * M * (1+(row.gamma-1)/2*M**2)**expo)

def rotor_calc(row:BladeRow,upstream:BladeRow,calculate_vm:bool=True,static_defined:bool=False):
    """Calculates quantities given beta2 

    Args:
        row (BladeRow): Rotor Row
        upstream (BladeRow): Stator Row or Rotor Row
        calculate_vm (bool): True to calculate the meridional velocity. False, do not calculate this and let radeq calculate it
        static_defined (bool): True if static conditions defined at the outlet. False if total conditions defined at outlet
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

    row.P0_stator_inlet = upstream.P0_stator_inlet
    ## P0_P is assumed 
    # row.P = row.P0_stator_inlet*1/row.P0_P
    
    upstream_radius = upstream.r
    # Upstream Relative Frame Calculations 
    upstream.U = upstream.rpm*np.pi/30 * upstream_radius # rad/s 
    upstream.Wt = upstream.Vt - upstream.U
    upstream.W = np.sqrt(upstream.Vx**2 + upstream.Wt**2 + upstream.Vr**2)
    upstream.beta2 = np.arctan2(upstream.Wt,upstream.Vm)
    upstream.T0R = upstream.T+upstream.W**2/(2*upstream.Cp)
    upstream.P0R = upstream.P * (upstream.T0R/upstream.T)**((upstream.gamma)/(upstream.gamma-1))      
    upstream.M_rel = upstream.W/np.sqrt(upstream.gamma*upstream.R*upstream.T)
    upstream_rothalpy = upstream.T0R*upstream.Cp - 0.5*upstream.U**2 # H01R - 1/2 U1^2 
    row.U = row.omega*row.r
    
    if np.any(upstream_rothalpy < 0):
        print('U is too high, reduce RPM or radius')
    
    # Rotor Exit Calculations
    row.beta1 = upstream.beta2
    
    row.T0R = upstream.T0R - T0_coolant_weighted_average(row)
    loss_type = getattr(row.loss_function, "loss_type", None)
    
    if loss_type == LossType.Pressure:
        row.Yp = row.loss_function(row,upstream)
        rotor_calculation(row.Yp)
    elif loss_type == LossType.Entropy: 
        desired_entropy_rise = convert_to_ndarray(row.loss_function(row,upstream))
        if len(desired_entropy_rise) == 1: # If entropy rise is a bulk value
            fun = lambda s: np.abs(desired_entropy_rise - rotor_calculation(s))
            x = minimize_scalar(fun,bounds=[0.01,0.4])
            row.Yp = x
        elif len(desired_entropy_rise) > 1: # In case entropy rise is an array 
            fun = lambda s,j: np.abs(desired_entropy_rise[j] - rotor_calculation(s)[j])
            for j in range(len(desired_entropy_rise)):
                x = minimize_scalar(fun,bounds=[0.01,0.4],args=(j))
                row.Yp[j] = x
            rotor_calculation(row.Yp)
    else: # loss_type == LossType.Enthalpy: 
        # Enthalpy loss is typically calculated for a stage so I bulk the loss after the rotor. This is a really terrible way of doing loss but this is here for those legacy loss models. 
        rotor_calculation(row.Yp)
    
    def rotor_calculation(Yp:npt.NDArray): 
        row.Yp = Yp
        # Total Relative Temperature stays constant through the rotor. Adjust for change in radius from rotor inlet to exit
        row.T0R = upstream.T0R # (upstream_rothalpy + 0.5*row.U**2)/row.Cp # - T0_coolant_weighted_average(row) 
        
        # ---- Compressor ----
        row.P0R = upstream.P0R - row.Yp*(upstream.P0R-upstream.P)
        M_rel = row.r*0
        for j in range(1,row.r):
            M_rel[j] = minimize_scalar(solve_for_mach,bounds=[0.01,1],args=(row,row.area*row.percent_hub_shroud,np.diff(upstream.massflow)))
        M_rel[0] = 1 / (len(M_rel) - 1) * M_rel[1:].sum() # This way average stays the same 
        # We basically need static pressure to do the rest of the calculations 
        row.P = row.P0R/(1+(row.gamma-1)/2*M_rel**2)**(row.gamma/(row.gamma-1)) 
        # ---- end compressor ----
        
        P0R_P = row.P0R / row.P
        T0R_T = P0R_P**((row.gamma-1)/row.gamma)
        row.T = (row.T0R/T0R_T)     # Exit static temperature
        if calculate_vm:    # Calculates the T0 at the exit
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
            row.alpha2 = np.arctan2(row.Vt,row.Vm)
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
        row.P0_P = (row.P0_stator_inlet/row.P).mean()

        row.M_rel = row.W/np.sqrt(row.gamma*row.R*row.T)
        row.T0 = row.T+row.V**2/(2*row.Cp)
        row.entropy_rise = 0.5*(row.Cp+upstream.Cp)*np.log(row.T/upstream.T) - row.R * np.log(row.P/upstream.P)    
        return row.entropy_rise    
    
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
            rows: List of blade rows between inlet and outlet
            num_streamlines: number of streamlines used through the meridional passage
            fluid: cantera gas solution; defaults to air.yaml if None
            rpm: RPM for the entire spool. Individual rows can override later.
            massflow_constraint: AngleMatch (adjust turning) or PressureBalance (radial eq).
        """
        self.passage = passage
        self.massflow = massflow
        self.inlet = inlet
        self.outlet = outlet
        self.rows = rows
        self.num_streamlines = num_streamlines
        self._fluid = fluid if fluid is not None else Solution("air.yaml")
        self.massflow_constraint = massflow_constraint
        self.rpm = rpm

        # Previously this used dataclasses.field on a non-dataclass; do it explicitly
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
        if row.num_blades and row.chord != 0:
            mean_r = float(row.r.mean())
            pitch = 2 * np.pi * mean_r / row.num_blades
            row.pitch_to_chord = pitch / row.chord

    # ------------------------------
    # initialization/solve
    # ------------------------------
    def initialize(self) -> None:
        """Initialize massflow and thermodynamic state through rows (turbines)."""
        blade_rows = self._all_rows()
        Is_static_defined = self.outlet.static_defined # This is set when you initialize the outlet 

        # Inlet
        W0 = self.massflow
        inlet: Inlet = self.inlet
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
        interpolate_streamline_radii(inlet, self.passage, self.num_streamlines)

        compute_gas_constants(inlet, self.fluid)
        inlet_calc(inlet)

        for row in blade_rows:
            interpolate_streamline_radii(row, self.passage, self.num_streamlines)

        outlet: Outlet = self.outlet
        for j in range(self.num_streamlines):
            P0 = inlet.get_total_pressure(inlet.percent_hub_shroud[j])  # type: ignore[attr-defined]
            percents = np.zeros(shape=(len(blade_rows) - 2)) + 0.3
            percents[-1] = 1
            if Is_static_defined:
                Ps_range = outlet_pressure(percents=percents, inletP0=inlet.P0[j], outletP=outlet.P[j])
                for i in range(1, len(blade_rows) - 1):
                    blade_rows[i].P[j] = Ps_range[i - 1]
            else:
                P0_range = outlet_pressure(percents=percents, inletP0=inlet.P0[j], outletP=outlet.P[j])
                for i in range(1, len(blade_rows) - 1):
                    blade_rows[i].P0[j] = P0_range[i - 1]
                    

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
        if self.massflow_constraint == MassflowConstraint.AngleMatch:
            pass
        elif self.massflow_constraint == MassflowConstraint.PressureBalance: # Balances the Total pressure
            self._pressure_balance()
    
    def _pressure_balance(self):
        blade_rows = self._all_rows()
        outlet_P = []
        outlet_P_guess = []
        for i in range(1, len(blade_rows) - 2):
            outlet_P.append(blade_rows[i].inlet_to_outlet_pratio) # inlet_to_outlet_pratio is a default property  
            outlet_P_guess.append(np.mean(blade_rows[i].inlet_to_outlet_pratio))
            
        
    # ------------------------------
    # Export / Plotting
    # ------------------------------
    def export_properties(self, filename: str = "compressor_spool.json") -> None:
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
            
            elif row.loss_function.loss_type == LossType.Enthalpy: #! Need to think about this since Yp is calculated inside stator_calc and rotor_calc
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
