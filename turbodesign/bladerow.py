from dataclasses import dataclass, field, Field
from typing import Any, Callable, List, Optional, Sequence, Tuple, Union
from .enums import RowType, PowerType
import numpy as np 
import numpy.typing as npt
from scipy.interpolate import interp1d
from .arrayfuncs import convert_to_ndarray
from cantera import Solution, composite
from .coolant import Coolant
from pyturbo.helper import line2D
from pyturbo.aero.airfoil2D import Airfoil2D
from .loss import LossBaseClass
from .deviation.deviation_base import DeviationBaseClass
from .passage import Passage
from .arrayfuncs import safe_interpolate
    

@dataclass(eq=False)
class BladeRow:
    """A single blade row (stator or rotor) in a turbomachine stage.

    Attributes
    ----------
    **Core Configuration**

    id : int
        Row identifier.
    stage_id : int
        Stage identifier.
    row_type : RowType
        Stator or Rotor.
    loss_function : LossBaseClass, optional
        Loss model applied to this row.
    deviation_function : DeviationBaseClass, optional
        Deviation model applied to this row.
    cutting_line : line2D, optional
        Line perpendicular to the streamline.
    rp : float
        Degree of reaction.
    hub_location : float
        Hub axial location.
    shroud_location : float
        Shroud axial location.

    **Fluid Properties**

    R : float
        Ideal gas constant, J/(kg K). Default 287.15.
    gamma : float
        Ratio of specific heats Cp/Cv. Default 1.33.
    Cp : float
        Specific heat at constant pressure, J/(kg K). Default 1019.
    Cv : float
        Specific heat at constant volume, J/(kg K).
    mu : float
        Dynamic viscosity, Pa s.

    **Mass Flow**

    total_massflow : float
        Total mass flow including upstream cooling, kg/s.
    massflow : ndarray
        Mass flow distribution per radial station.
    total_massflow_no_coolant : float
        Inlet mass flow without coolant, kg/s.
    massflow_target : ndarray, optional
        Custom mass flow distribution for angle matching, kg/s.

    **Streamline Geometry**

    percent_hub : float
        Where blade row is defined along the hub (0-1).
    percent_hub_shroud : ndarray
        Percent streamline length from hub to shroud.
    x : ndarray
        Axial coordinates.
    r : ndarray
        Radial coordinates.
    m : ndarray
        Meridional coordinates.
    total_area : float
        Total annular flow area.
    area : ndarray
        Flow area per streamline.

    **Row Efficiency**

    eta_total : float
        Total-to-total isentropic efficiency.
    eta_static : float
        Total-to-static isentropic efficiency.
    eta_poly : float
        Polytropic efficiency.
    stage_loading : float
        Stage loading coefficient (work per stage).

    **Flow Angles** *(radians)*

    alpha1, alpha2 : ndarray
        Absolute flow angles at inlet and exit.
    beta1, beta2 : ndarray
        Relative flow angles at inlet and exit.
    deviation : ndarray
        Flow deviation from metal angle.
    beta1_fixed, beta2_fixed : bool
        Whether inlet/exit geometry is already defined.

    **Velocities**

    Vm : ndarray
        Meridional velocity.
    Vx : ndarray
        Axial velocity.
    Vt : ndarray
        Tangential (swirl) velocity.
    Vr : ndarray
        Radial velocity.
    V : ndarray
        Absolute velocity magnitude.
    U : ndarray
        Blade peripheral velocity.
    W : ndarray
        Relative velocity magnitude.
    Wt : ndarray
        Relative tangential velocity.
    M : ndarray
        Absolute Mach number.
    M_rel : ndarray
        Relative Mach number.
    omega : float
        Angular velocity, rad/s.

    **Thermodynamic Quantities**

    P0 : ndarray
        Total pressure, Pa.
    T0 : ndarray
        Total temperature, K.
    P : ndarray
        Static pressure, Pa.
    T : ndarray
        Static temperature, K.
    rho : ndarray
        Density, kg/m^3.
    P0R : ndarray
        Relative total pressure, Pa.
    T0R : ndarray
        Relative total temperature, K.
    entropy_rise : ndarray
        Entropy rise across row.

    **Performance**

    power : float
        Power, W.
    P0_P : float
        Total-to-static pressure ratio.
    P0_ratio : float
        Total-to-total pressure ratio.
    flow_coefficient : float
        Flow coefficient (Vm/U).
    Reynolds : float
        Reynolds number.
    Yp : ndarray
        Pressure loss coefficient.

    **Blade Geometry** *(set via properties)*

    axial_chord : float
        Axial chord length. Set via property.
    aspect_ratio : float
        Height-to-chord ratio. Set via property.
    pitch_to_chord : float
        Pitch-to-chord ratio. Set via property.
    stagger : float
        Stagger angle, degrees. Set via property.
    num_blades : int
        Blade count. Set via property.
    tip_clearance : float
        Clearance as fraction of span. Set via property.
    te_pitch : float
        Trailing-edge-to-pitch ratio. Set via property.
    blade_to_blade_gap : float
        Inter-row gap as fraction of chord. Set via property.
    """

    id: int = 0
    stage_id: int = 0
    row_type: RowType = RowType.Stator
    loss_function: Optional[LossBaseClass] = None
    deviation_function: Optional[DeviationBaseClass] = None
    cutting_line: Optional[line2D] = None         # Line perpendicular to the streamline
    rp: float = 0.4              # Degree of Reaction
    hub_location: float = 0.0
    shroud_location: float = 0.0

    # Fluid
    R: float = 287.15           # Ideal Gas constant J/(Kg K)
    gamma: float = 1.33         # Ratio of Cp/Cv
    Cp: float = 1019            # Cp J/(Kg*K)
    Cv: float = 1019/1.14       # Cv J/(Kg*K)
    _coolant: Optional[Coolant] = None     # Coolant Fluid
    mu: float = 0

    total_massflow: float = 0    # Massflow spool + all upstream cooling flow [kg/s]
    massflow: npt.NDArray = field(default_factory=lambda: np.array([0]))  # Massflow per radii
    total_massflow_no_coolant: float = 0     # Inlet massflow
    massflow_target: Optional[npt.NDArray] = None  # Custom massflow distribution for angle matching [kg/s]
    # ----------------------------------

    # Streamline Properties 
    percent_hub: float = 0 # Where blade row is defined along the hub. 
    percent_hub_shroud: npt.NDArray = field(default_factory=lambda: np.array([0]))    # Percent streamline length from hub to shroud.
    x: npt.NDArray = field(default_factory=lambda: np.array([0]))       # x - coordinates (useful for computing axial chord)
    r: npt.NDArray = field(default_factory=lambda: np.array([0]))       # Radius - coordinates 
    m: npt.NDArray = field(default_factory=lambda: np.array([0]))       # meridional 
    total_area: float = 0
    area: npt.NDArray = field(default_factory=lambda: np.array([0]))
    # Calculated massflow is the massflow computed after radial eq solver
    calculated_massflow: float = 0

    # Row Efficiency (calculated or specified)
    eta_total: float = 0     # Total to Total
    eta_static: float = 0    # Total to static
    eta_poly: float = 0     # Polytropic efficiency (per row if applicable)
    stage_loading: float = 0 # stage loading how much work done per stage

    alpha1: npt.NDArray = field(default_factory=lambda: np.array([0]))               # Blade inlet absolute flow angle
    alpha2: npt.NDArray = field(default_factory=lambda: np.array([0]))               # Blade exit absolute flow angle
    
    beta1: npt.NDArray = field(default_factory=lambda: np.array([0]))                 # Blade inlet relative flow angle
    beta2: npt.NDArray = field(default_factory=lambda: np.array([0]))                 # Blade exit relative flow angle
    
    deviation: npt.NDArray = field(default_factory=lambda: np.array([0])) 
    _beta1_metal: npt.NDArray = field(default_factory=lambda: np.array([0]))           # blade inlet metal angle
    beta1_metal_radii: npt.NDArray = field(default_factory=lambda: np.array([0]))      # radii where metal angle is defined
    
    _beta2_metal: npt.NDArray = field(default_factory=lambda: np.array([0]))           # blade exit metal angle
    beta2_metal_radii: npt.NDArray = field(default_factory=lambda: np.array([0]))      # radii where metal angle is defined
    
    beta1_fixed: bool = False    # Geometry already defined. This affects the inlet flow angle
    beta2_fixed: bool = False    # Geometry already defined. This affects the exit flow angle

    # Velocities 
    Vm: npt.NDArray = field(default_factory=lambda: np.array([0]))               # Meridional velocity
    Vx: npt.NDArray = field(default_factory=lambda: np.array([0]))               # Axial Velocity
    Vt: npt.NDArray = field(default_factory=lambda: np.array([0]))               # Tangential Velocity
    Vr: npt.NDArray = field(default_factory=lambda: np.array([0]))                # Radial velocity 
    V: npt.NDArray = field(default_factory=lambda: np.array([0]))                # Absolute Velocity in 3D coordinate system
    V2: npt.NDArray = field(default_factory=lambda: np.array([0]))               # Absolute Velocity in Theta-Axial plane
    M: npt.NDArray = field(default_factory=lambda: np.array([0]))                # Mach Number
    M_rel: npt.NDArray = field(default_factory=lambda: np.array([0]))            # Relative Mach Number
    U: npt.NDArray = field(default_factory=lambda: np.array([0]))                # Peripheral velocity
    W: npt.NDArray = field(default_factory=lambda: np.array([0]))                # Relative Velocity in Theta-Axial plane
    Wt: npt.NDArray = field(default_factory=lambda: np.array([0]))               # Relative Tangential Velocity    

    _rpm: float = field(default=0, init=False, repr=False)
    omega: float = 0 # angular velocity rad/s
    
    P0_stator_inlet: npt.NDArray = field(default_factory=lambda: np.array([0]))              # Every quantity is an exit quantity, This is used  for efficiency calcs
    T0_stator_inlet: npt.NDArray = field(default_factory=lambda: np.array([0]))              # Every quantity is an exit quantity, This is used  for efficiency calcs
    P0: npt.NDArray = field(default_factory=lambda: np.array([0]))              # Total Quantities 
    P0_is: npt.NDArray = field(default_factory=lambda: np.array([0]))
    T0: npt.NDArray = field(default_factory=lambda: np.array([0]))          
    T0_is: npt.NDArray = field(default_factory=lambda: np.array([0])) 
    P0R: npt.NDArray = field(default_factory=lambda: np.array([0]))             # Relative Total Pressure (Pa)
    P0R_is: npt.NDArray = field(default_factory=lambda: np.array([0]))
    T0R: npt.NDArray = field(default_factory=lambda: np.array([0]))
    
    # Static Quantities
    P: npt.NDArray = field(default_factory=lambda: np.array([0]))
    T: npt.NDArray = field(default_factory=lambda: np.array([0]))
    T_is: npt.NDArray = field(default_factory=lambda: np.array([0]))
    rho: npt.NDArray = field(default_factory=lambda: np.array([0]))
    entropy_rise: npt.NDArray = field(default_factory=lambda: np.array([0]))
    
    # Related to streamline curvature
    phi: npt.NDArray = field(default_factory=lambda: np.array([0]))                      # Inclination angle x,r plane. AY td2.f
    rm: npt.NDArray = field(default_factory=lambda: np.array([0]))                      # Curvature
    incli_curve_radii: npt.NDArray = field(default_factory=lambda: np.array([0]))       # radius at which curvature was evaluated
    mprime: npt.NDArray = field(default_factory=lambda: np.array([0]))                   # Mprime distance
    
    Yp: npt.NDArray = field(default_factory=lambda: np.array([0]))                       # Pressure loss
    blockage: float = 0 
    flow_coefficient: float = 0     # Vm/U or similar nondimensional flow coefficient
    power: float = 0                 # Watts 
    power_mean: float = 0
    power_distribution: npt.NDArray = field(default_factory=lambda: np.array([0]))  # How power is divided by radius.
    P0_P: float = 0                  # Total to Static Pressure Ratio
    P0_ratio: float = 0              # Total-to-total pressure ratio target (design input; may be overwritten in legacy diagnostics)
    P0_ratio_target: float = 0       # Frozen design target for P0_ratio (never overwritten; used for initial guesses)
    Power_Type: PowerType = PowerType.P0_P
    euler_power: float = 0
    Reynolds: float = 0
    eta_poly: float = 0.0            # Optional per-row polytropic efficiency target
    num_blades: int = 0
    
    # Used for loss calculations
    _blade_to_blade_gap: float = 0.025 # Gap between blade in terms of percent chord.
    
    _aspect_ratio: float = 0.9 # 
    _pitch_to_chord: npt.NDArray = field(default_factory=lambda: np.array([0.7])) # Pitch to chord ratio, used to determine number of blades and compute loss 
    
    _axial_chord: float = -1 
    _chord: npt.NDArray = field(default_factory=lambda: np.array([-1.0]))
    _stagger: npt.NDArray = field(default_factory=lambda: np.array([42.0]))
    _te_s: float = 0.08
    _tip_clearance: float = 0 # Clearance as a percentage of span or blade height

    _inlet_to_outlet_pratio: list = field(default_factory=lambda: [0.06,0.95])

    def __post_init__(self):
        if self.shroud_location == 0:
            self.shroud_location = self.hub_location
        # Preserve any user-specified target ratio so later diagnostics can safely overwrite P0_ratio.
        if self.P0_ratio_target == 0 and self.P0_ratio != 0:
            self.P0_ratio_target = self.P0_ratio
    
    @property
    def inlet_to_outlet_pratio(self) -> Tuple[float,float]:
        """This is what is varied by the optimization. 
        The range is between [0 and 1] but you should 

        Returns:
            List[float]: _description_
        """
        return self._inlet_to_outlet_pratio # type: ignore
    
    @inlet_to_outlet_pratio.setter
    def inlet_to_outlet_pratio(self,val:Tuple[float,float]=(0.06,0.7)):
        """Sets the inlet_to_outlet pratio of the blade row 

        Args:
            val (Tuple[float,float], optional): guess value from [0,1], you should not use 0 or 1 though. Defaults to (0.06,0.7).
        """
        self._inlet_to_outlet_pratio = val 
        
    @property
    def blade_to_blade_gap(self) -> float:
        """Returns the blade to blade gap value 

        Returns:
            float: _description_
        """
        return self._blade_to_blade_gap
    
    @blade_to_blade_gap.setter
    def blade_to_blade_gap(self,val:float):
        """Sets the blade to blade gap as a percent. This applies to the next row. So (row1) (row2) if (row1) gap is set to 0.25 then (row2) is offset by 0.25*row1.chord

        Args:
            val (float): percentage of chord to space out the blade 
        """
        self._blade_to_blade_gap = val

    @property
    def aspect_ratio(self):
        return self._aspect_ratio
    
    @aspect_ratio.setter
    def aspect_ratio(self,val:float):
        """Sets the aspect ratio

        Aspect ratio is defined as the height/chord not height/(axial chord)

        Args:
            val (float): new aspect ratio 
        """
        self._aspect_ratio = val

    @property
    def axial_chord(self) -> float:
        """Returns the mean axial chord defined in the x-direction 

        Returns:
            float: Axial Chord
        """
        
        return self._axial_chord
    
    @axial_chord.setter
    def axial_chord(self,val:float):
        self._axial_chord = val
        
        
    @property
    def pitch_to_chord(self) -> npt.NDArray:
        """Gets the pitch to chord ratio 

        Returns:
            float: pitch to chord ratio 
        """
        return self._pitch_to_chord
    
    @pitch_to_chord.setter
    def pitch_to_chord(self,val:float):
        """Set the pitch to chord ratio 

        Args:
            val (float): new pitch to chord ratio. Typically stators are 0.8 to 0.95. Rotors 0.7 to 0.8 
        """
        self._pitch_to_chord = convert_to_ndarray(val)
    
    @property
    def solidity(self) -> npt.NDArray:
        """Inverse of pitch to chord ratio

        Returns:
            float: solidity value
        """
        return 1/self._pitch_to_chord 
    
    @solidity.setter
    def solidity(self,val:float):
        """Inverse of pitch to chord ratio

        Args:
            val (float): sets the inverse of pitch to chord ratio

        Returns:
            float: solidity
        """
        self._pitch_to_chord = 1/convert_to_ndarray(val)

    @property
    def metal_inlet_angle(self) -> npt.NDArray:
        """Blade metal inlet angle (degrees)."""
        return np.degrees(self._beta1_metal)
    
    @property
    def beta1_metal(self) -> npt.NDArray:
        """Backward-compatible alias for metal_inlet_angle."""
        return self.metal_inlet_angle
    
    @property
    def metal_exit_angle(self) -> npt.NDArray:
        """Blade metal exit angle (degrees)."""
        return np.degrees(self._beta2_metal)
    
    @property
    def beta2_metal(self) -> npt.NDArray:
        """Backward-compatible alias for metal_exit_angle."""
        return self.metal_exit_angle
    
    @property
    def stagger(self) -> float:
        """Average stagger angle

        Returns:
            float: stagger angle
        """
        return self._stagger
    
    @stagger.setter
    def stagger(self,val:float):
        """Set the stagger angle in degrees 

        Args:
            val (float): stagger angle. Degrees
        """
        self._stagger = val
    
    @property
    def chord(self) -> npt.NDArray:
        """Chord defined at mean radius

        Returns:
            float: axial chord
        """
        return self.axial_chord / np.cos(np.radians(self.stagger))
    
    @property
    def pitch(self) -> float:
        """Returns the pitch which is the distance from blade to blade

        Returns:
            float: pitch
        """
        return self.pitch_to_chord*self.chord
    
    @property
    def throat(self) -> float:
        """Throat distance

        Returns:
            float: throat 
        """
        if self.row_type == RowType.Stator:
            return self.pitch*np.sin(np.pi/2-self.alpha2.mean())
        else:
            return self.pitch*np.sin(np.pi/2-self.beta2.mean())
    
    _num_blades: float = 0

    @property
    def num_blades(self) -> float:
        """Configured number of blades (set during design/initialization)."""
        return self._num_blades

    @num_blades.setter
    def num_blades(self, val: float) -> None:
        self._num_blades = val
    
    @property
    def camber(self) -> float:
        """Estimates the camber of the blade using a bezier curve. This is not as accurate because thickness is not defined on suction and pressure sides. 

        Returns:
            float: camber length
        """
        if self.row_type == RowType.Stator:
            a2d = Airfoil2D(np.degrees(self.alpha1.mean()),
                        np.degrees(self.alpha2.mean()),
                        self.axial_chord,
                        self.stagger)
        else:
            a2d = Airfoil2D(np.degrees(self.beta1.mean()),
                        np.degrees(self.beta2.mean()),
                        self.axial_chord,
                        self.stagger)
        
        return a2d.camberBezier.get_curve_length()
    
    @property
    def tip_clearance(self):
        """Tip clearance as a percentage of annulus height
        """
        return self._tip_clearance
    
    @tip_clearance.setter
    def tip_clearance(self,val:float):
        """Sets the tip clearance

        Args:
            val (float): tip clearance as a percentage of annulus height.
        """
        self._tip_clearance = val
        
    # Backwards-compatible alias
    @property
    def location(self) -> float:
        return self.hub_location

    @location.setter
    def location(self, val: float) -> None:
        self.hub_location = val
    
    @metal_inlet_angle.setter
    def metal_inlet_angle(self, metal_inlet_angle: List[float], percent: List[float] = []):
        """Sets the leading edge metal angle for the blade (degrees)."""
        arr = np.radians(convert_to_ndarray(metal_inlet_angle))
        if len(percent) != len(metal_inlet_angle):
            percent = np.linspace(0, 1, len(arr)).tolist()  # type: ignore
        self._beta1_metal = arr
        self.beta1_metal_radii = convert_to_ndarray(percent)
        self.beta1_fixed = True
        self.beta1 = self.metal_inlet_angle.copy()
    
    @beta1_metal.setter
    def beta1_metal(self, beta1_metal: List[float], percent: List[float] = []):
        """Backward-compatible alias for metal_inlet_angle setter."""
        self.metal_inlet_angle = beta1_metal
        
    @metal_exit_angle.setter
    def metal_exit_angle(self, metal_exit_angle: List[float], percent: List[float] = []):
        """Sets the trailing edge metal angle for the blade (degrees)."""
        arr = np.radians(convert_to_ndarray(metal_exit_angle))
        if len(percent) != len(metal_exit_angle):
            percent = np.linspace(0, 1, len(arr)).tolist()  # type: ignore
        self._beta2_metal = arr
        self.beta2_metal_radii = convert_to_ndarray(percent)
        self.beta2_fixed = True

        # Apply deviation if defined; deviation_function returns degrees
        deviation_func = getattr(self, "deviation_function", None)
        deviation_rad = 0.0
        if callable(deviation_func):
            try:
                deviation_val = deviation_func(self, None)
                deviation_rad = np.radians(deviation_val)
            except Exception:
                deviation_rad = 0.0

        beta2_effective = self._beta2_metal + deviation_rad
        if self.row_type == RowType.Stator:
            self.alpha2 = beta2_effective.copy()
            self.beta2 = beta2_effective.copy()
        else:
            self.beta2 = beta2_effective.copy()
        self.deviation = np.full_like(self.beta2, deviation_rad)
    
    @beta2_metal.setter
    def beta2_metal(self, beta2_metal: List[float], percent: List[float] = []):
        """Backward-compatible alias for metal_exit_angle setter."""
        self.metal_exit_angle = beta2_metal
        
    @property
    def rpm(self):
        return self._rpm 
    
    @rpm.setter
    def rpm(self,val:float):
        self._rpm = val
        self.omega = self._rpm * np.pi/30 # rev/min * 2pi rads/1 rev * 1 min/60 sec
    
    @property
    def coolant(self):
        return self._coolant
    
    @coolant.setter
    def coolant(self,coolant:Coolant):
        """Add a coolant to the end of the blade row

        Args:
            coolant (Coolant): Coolant
        """
        self._coolant = coolant
    
    @property    
    def loss_model(self):
        return self.loss_function
    
    @loss_model.setter
    def loss_model(self, model:Union[LossBaseClass, Sequence[LossBaseClass]]):
        """Assign one or more loss models that inherit :class:`LossBaseClass`.

        Args:
            model: Either a single loss model or a sequence of models.
        """
        if isinstance(model, LossBaseClass):
            self.loss_function = model
            return

        raise TypeError("Loss models must inherit LossBaseClass.")
    
    @property
    def te_pitch(self):
        """Trailing edge to pitch ratio 

        Returns:
            float: trailing edge to pitch ratio
        """
        return self._te_s
    
    @te_pitch.setter
    def te_pitch(self,val:float):
        """Set the trailing edge to pitch ratio. Typical values from 0.02 to 0.12

        Args:
            val (float): new trailing edge to pitch ratio
        """
        self._te_s = val
    
    def __repr__(self):
        return f"{self.row_type.name} P0:{np.mean(self.P0):0.2f} T0:{np.mean(self.T0):0.2f} P:{np.mean(self.P):0.2f} massflow:{np.mean(self.total_massflow_no_coolant):0.3f}"

    def synchronize_blade_geometry(self) -> None:
        """Couple num_blades, pitch-to-chord/solidity, chord, and stagger.

        Uses mean radius from interpolated streamlines to derive pitch, chord,
        and stagger (axial chord / chord).
        """
        if self.num_blades <= 0 or self.r.size == 0:
            return

        # Pitch from blade count and local radius
        pitch = 2 * np.pi * self.r / self.num_blades

        # Pitch-to-chord (or 1/solidity) may be scalar or spanwise; broadcast it
        ptc = convert_to_ndarray(self.pitch_to_chord)
        if ptc.size == 1:
            ptc = ptc * np.ones_like(self.r, dtype=float)
        else:
            t_src = np.linspace(0, 1, ptc.size)
            ptc = np.interp(self.percent_hub_shroud, t_src, ptc)

        chord = pitch / np.maximum(ptc, 1e-9)
        self._chord = chord
        self._pitch_to_chord = ptc

        axial = self.axial_chord if self.axial_chord > 0 else float(np.mean(chord))
        if self.axial_chord <= 0:
            self.axial_chord = axial

        ratio = np.clip(axial / np.maximum(chord, 1e-9), -1.0, 1.0)
        stagger_rad = np.arccos(ratio)
        # Store stagger distribution in degrees
        self._stagger = np.degrees(stagger_rad)
    
    def to_dict(self):
        
        data = {
            "StageID":self.stage_id,
            "RowType":self.row_type.name,
            "R":self.R,
            "gamma":self.gamma,
            "Cp":self.Cp,
            "Cv":self.Cv,
            "P0_P":self.P0_P,
            "rp":self.rp,
            "total_massflow":self.total_massflow,
            "massflow":self.massflow.tolist(),
            "calculated_massflow":self.calculated_massflow,
            "alpha1":np.degrees(self.alpha1).tolist(),
            "alpha2":np.degrees(self.alpha2).tolist(),
            "beta1":np.degrees(self.beta1).tolist(),
            "beta2":np.degrees(self.beta2).tolist(),
            "beta1_metal":np.degrees(self._beta1_metal).tolist(),
            "beta2_metal":np.degrees(self._beta2_metal).tolist(),
            "Vm":self.Vm.tolist(),
            "Vx":self.Vx.tolist(),
            "Vr":self.Vr.tolist(),
            "Vt":self.Vt.tolist(),
            "U":self.U.tolist(),
            "V":self.V.tolist(),
            "M":self.M.tolist(),
            "M_rel":self.M_rel.tolist(),
            "W":self.W.tolist(),
            "Wt":self.Wt.tolist(),
            "omega":self.omega,
            "P0":self.P0.tolist(),
            "T0":self.T0.tolist(),
            "P0R":self.P0R.tolist(),
            "T0R":self.T0R.tolist(),
            "P":self.P.tolist(),
            "T":self.T.tolist(),
            "rho":self.rho.tolist(),
            "mu":self.mu,
            "Yp":self.Yp,
            "flow_coefficient": self.flow_coefficient,
            "Power":self.power,
            "P0_P": self.P0_P,
            "eta_total":self.eta_total,
            "eta_static":self.eta_static,
            "eta_poly": self.eta_poly,
            "euler_power":self.euler_power,
            "axial_chord":self.axial_chord,
            "aspect_ratio":self.aspect_ratio,
            "num_blades":self.num_blades,
            "total_area": self.total_area,
            "area": self.area.tolist(),
            "radius":self.r.tolist(),
            "x":self.x.tolist(),
            "dx":self.x[-1]-self.x[0],
            "dr":self.r[-1]-self.r[0],
            "mprime":self.mprime[-1],
            "Reynolds":self.Reynolds,
            "axial_chord":self.axial_chord
        }

        return data

#* Some functions related to blade row 
def interpolate_streamline_quantities(row:BladeRow,passage:Passage,num_streamlines:int=3):
    """Interpolate all quantities onto the streamline and allocates variables. 
    Run this after setting some initial conditions 

    Args:
        r_streamline (npt.NDArray): Radii describing the streamline 
        passage (Passage): Passage object describing the geometry of the hub and shroud
        num_streamlines (int): number of streamlines to consider

    Returns:
        (BladeRow): new row object with quantities interpolated
    """
    src_percent = convert_to_ndarray(row.percent_hub_shroud)

    row.cutting_line,_,_ = passage.get_cutting_line(row.location)
    t_span = np.array([0.5]) if num_streamlines <= 1 else np.linspace(0, 1, num_streamlines)
    row.x, row.r = row.cutting_line.get_point(t_span)
    if num_streamlines <= 1:
        streamline_percent_length = np.array([0.5])
        row.total_area = passage.get_area(row.location)
        row.area = np.array([row.total_area])
    else:
        streamline_percent_length = np.sqrt((row.r-row.r[0])**2+(row.x-row.x[0])**2)/row.cutting_line.length
    
    # Flow angles 
    row._beta1_metal = row._beta1_metal.default_factory() if type(row._beta1_metal) == Field else row._beta1_metal

    if row.row_type==RowType.Stator:
        assert type(row._beta2_metal)!=Field,"Stator exit Flow angle must be set"

    row._beta1_metal = interpolate_quantities(row._beta1_metal,row.beta1_metal_radii,streamline_percent_length)
    row._beta2_metal = interpolate_quantities(row._beta2_metal,row.beta2_metal_radii,streamline_percent_length)
    row.beta1_metal_radii = streamline_percent_length
    row.beta2_metal_radii = streamline_percent_length
    row.deviation = streamline_percent_length * 0
    
    row.mprime = interpolate_quantities(row.mprime, src_percent, streamline_percent_length)

    row.alpha1 = safe_interpolate(row.alpha1, src_percent, streamline_percent_length, radians=False)
    row.alpha2 = safe_interpolate(row.alpha2, src_percent, streamline_percent_length, radians=False)
    row.beta1 = safe_interpolate(row.beta1, src_percent, streamline_percent_length, radians=False)
    row.beta2 = safe_interpolate(row.beta2, src_percent, streamline_percent_length, radians=False)
    
    # Velocities 
    row.Vm = interpolate_quantities(row.Vm, src_percent, streamline_percent_length)
    row.Vx = interpolate_quantities(row.Vx, src_percent, streamline_percent_length)
    row.Vt = interpolate_quantities(row.Vt, src_percent, streamline_percent_length)
    row.Vr = interpolate_quantities(row.Vr, src_percent, streamline_percent_length)
    row.V = interpolate_quantities(row.V, src_percent, streamline_percent_length)
    row.V2 = interpolate_quantities(row.V2, src_percent, streamline_percent_length)
    row.M = interpolate_quantities(row.M, src_percent, streamline_percent_length)
    row.M_rel = interpolate_quantities(row.M_rel, src_percent, streamline_percent_length)
    row.U = interpolate_quantities(row.U, src_percent, streamline_percent_length)
    row.W = interpolate_quantities(row.W, src_percent, streamline_percent_length)
    row.Wt = interpolate_quantities(row.Wt, src_percent, streamline_percent_length)

    # Total Quantities
    row.T0 = interpolate_quantities(row.T0, src_percent, streamline_percent_length)
    row.T0_is = interpolate_quantities(row.T0, src_percent, streamline_percent_length) # For Turbines
    row.P0 = interpolate_quantities(row.P0, src_percent, streamline_percent_length)
    row.P0_is = interpolate_quantities(row.P0, src_percent, streamline_percent_length) # For Compressors 
    row.P0_stator_inlet = interpolate_quantities(row.P0_stator_inlet, src_percent, streamline_percent_length)
    
    # Relative Quantities
    row.P0R = interpolate_quantities(row.P0R, src_percent, streamline_percent_length)
    row.P0R_is = interpolate_quantities(row.P0, src_percent, streamline_percent_length)
    row.T0R = interpolate_quantities(row.T0R, src_percent, streamline_percent_length)

    # Static Quantities 
    row.P = interpolate_quantities(row.P, src_percent, streamline_percent_length)
    row.T = interpolate_quantities(row.T, src_percent, streamline_percent_length)
    row.T_is = interpolate_quantities(row.T_is, src_percent, streamline_percent_length)
    row.rho = interpolate_quantities(row.rho, src_percent, streamline_percent_length)
    row.entropy_rise = interpolate_quantities(row.entropy_rise, src_percent, streamline_percent_length)

    row.percent_hub_shroud = streamline_percent_length

    return row


def sutherland(T:Union[float,npt.NDArray]) -> Union[float,npt.NDArray]:    
    """Sutherland viscosity calculation used for reynolds number 

    Args:
        T (float): Temperature in Kelvin

    Returns:
        float: Dynamic Viscosity (mu) in Pa*s
    """
    S = 110.4
    C1 = 1.458E-6
    return C1*T**1.5 / (T+S)

def interpolate_quantities(q:npt.NDArray,r:npt.NDArray,r2:npt.NDArray):
    """Interpolates array q

    Args:
        q (npt.NDArray): quantities defined at radius r 
        r (npt.NDArray): radius where quantities `q` are defined
        r2 (npt.NDArray): new radius to interpolate the quantities to e.g. streamline radius

    Returns:
        npt.NDArray: quantities interpolated onto r2
    """
    if (type(q) == Field):
        q = q.default_factory()
    if (type(r) == Field):
        r = r.default_factory()
    if len(q)==1:
        q2 = np.zeros(shape=r2.shape)
        return q[0]+q2
    else:
        if len(r) != len(q):
            r = np.linspace(0, 1, len(q))
        return interp1d(r,q,kind='linear')(r2)
    
def compute_gas_constants(row:BladeRow,fluid:Optional[Solution]=None) -> None:
    """Updates the Cp, Gamma, and density for a blade row. If fluid is not specified then only density and viscosity is updated. 
    
    Args:
        row (BladeRow): _description_
        fluid (Solution, optional): _description_. Defaults to None.

    Returns:
        (BladeRow): updated row
    """
    if fluid:
        Tm = row.T.mean()
        Pm = row.P.mean()
        # A single ct.Solution is shared across all blade rows (TurbineSpool/CompressorSpool
        # assign br.fluid = self._fluid), so mutating it in place leaks this row's (T,P) into
        # any later property read on the shared object. Snapshot and restore the state around
        # the read so the shared Solution is left untouched.
        _saved_state = fluid.state
        try:
            fluid.TP = Tm,Pm
            row.Cp = fluid.cp
            row.Cv = fluid.cv
        finally:
            fluid.state = _saved_state
        row.R = row.Cp-row.Cv
        row.gamma = row.Cp/row.Cv
    # Use Ideal Gas 
    row.rho = row.P/(row.T*row.R)
    row.mu = sutherland(row.T) # type: ignore
