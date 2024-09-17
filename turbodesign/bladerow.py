from dataclasses import field, Field
from typing import Any, Callable, List, Tuple, Union
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
from .passage import Passage
    

class BladeRow:
    stage_id:int = 0
    row_type: RowType = RowType.Stator
    loss_function:LossBaseClass
    cutting_line:line2D         # Line perpendicular to the streamline
    rp:float = 0.4              # Degree of Reaction
    
    # Fluid
    R: float = 287.15           # Ideal Gas constant J/(Kg K)
    gamma: float = 1.33         # Ratio of Cp/Cv
    Cp: float = 1019            # Cp J/(Kg*K)
    Cv: float = 1019/1.14       # Cv J/(Kg*K)
    _coolant:Coolant = None     # type: ignore # Coolant Fluid
    fluid: composite.Solution = Solution("air.yaml")
    mu:float = 0 
    
    total_massflow:float = 0    # Massflow spool + all upstream cooling flow [kg/s]
    massflow:npt.NDArray = field(default_factory=lambda: np.array([0]))  # Massflow per radii
    total_massflow_no_coolant:float = 0     # Inlet massflow
    # ----------------------------------

    # Streamline Properties 
    axial_location:float = 0 # Where blade row is defined along the hub. 
    percent_hub_shroud: npt.NDArray = field(default_factory=lambda: np.array([0]))    # Percent streamline length from hub to shroud.
    x: npt.NDArray = field(default_factory=lambda: np.array([0]))       # x - coordinates (useful for computing axial chord)
    r: npt.NDArray = field(default_factory=lambda: np.array([0]))       # Radius - coordinates 
    area:float = 0
    # Calculated massflow is the massflow computed after radial eq solver
    calculated_massflow: float = 0
    
    # Row Efficiency (calculated or specified)
    eta_total:float = 0     # Total to Total
    eta_static:float = 0    # Total to static
    stage_loading:float = 0 # stage loading how much work done per stage

    alpha1: npt.NDArray = field(default_factory=lambda: np.array([0]))               # Blade inlet absolute flow angle
    alpha2: npt.NDArray = field(default_factory=lambda: np.array([0]))               # Blade exit absolute flow angle
    
    beta1: npt.NDArray = field(default_factory=lambda: np.array([0]))                 # Blade inlet relative flow angle
    beta2: npt.NDArray = field(default_factory=lambda: np.array([0]))                 # Blade exit relative flow angle
    
    _beta1_metal:npt.NDArray = field(default_factory=lambda: np.array([0]))           # blade inlet metal angle
    beta1_metal_radii:npt.NDArray = field(default_factory=lambda: np.array([0]))      # radii where metal angle is defined
    
    _beta2_metal:npt.NDArray = field(default_factory=lambda: np.array([0]))           # blade exit metal angle
    beta2_metal_radii:npt.NDArray = field(default_factory=lambda: np.array([0]))      # radii where metal angle is defined
    
    beta1_fixed:bool = False    # Geometry already defined. This affects the inlet flow angle
    beta2_fixed:bool = False    # Geometry already defined. This affects the exit flow angle

    # Velocities 
    Vm: npt.NDArray = field(default_factory=lambda: np.array([0]))               # Meridional velocity
    Vx: npt.NDArray = field(default_factory=lambda: np.array([0]))               # Axial Velocity
    Vt: npt.NDArray = field(default_factory=lambda: np.array([0]))               # Tangential Velocity
    Vr:npt.NDArray = field(default_factory=lambda: np.array([0]))                # Radial velocity 
    V: npt.NDArray = field(default_factory=lambda: np.array([0]))                # Absolute Velocity in 3D coordinate system
    V2: npt.NDArray = field(default_factory=lambda: np.array([0]))               # Absolute Velocity in Theta-Axial plane
    M: npt.NDArray = field(default_factory=lambda: np.array([0]))                # Mach Number
    M_rel: npt.NDArray = field(default_factory=lambda: np.array([0]))            # Relative Mach Number
    U: npt.NDArray = field(default_factory=lambda: np.array([0]))                # Peripheral velocity
    W: npt.NDArray = field(default_factory=lambda: np.array([0]))                # Relative Velocity in Theta-Axial plane
    Wt: npt.NDArray = field(default_factory=lambda: np.array([0]))               # Relative Tangential Velocity    

    _rpm: float = 0 
    omega:float = 0 # angular velocity rad/s
    
    P0_stator_inlet = field(default_factory=lambda: np.array([0]))              # Every quantity is an exit quantity, This is used  for efficiency calcs
    T0_stator_inlet = field(default_factory=lambda: np.array([0]))              # Every quantity is an exit quantity, This is used  for efficiency calcs
    P0: npt.NDArray = field(default_factory=lambda: np.array([0]))              # Total Quantities 
    T0: npt.NDArray = field(default_factory=lambda: np.array([0]))          
    T0_is:npt.NDArray = field(default_factory=lambda: np.array([0])) 
    P0R: npt.NDArray = field(default_factory=lambda: np.array([0]))             # Relative Total Pressure (Pa)
    T0R: npt.NDArray = field(default_factory=lambda: np.array([0]))
    
    # Static Quantities
    P: npt.NDArray = field(default_factory=lambda: np.array([0]))
    T: npt.NDArray = field(default_factory=lambda: np.array([0]))
    T_is: npt.NDArray = field(default_factory=lambda: np.array([0]))
    rho: npt.NDArray = field(default_factory=lambda: np.array([0]))

    # Related to streamline curvature
    phi:npt.NDArray = field(default_factory=lambda: np.array([0]))                      # Inclination angle x,r plane. AY td2.f
    rm: npt.NDArray = field(default_factory=lambda: np.array([0]))                      # Curvature
    incli_curve_radii: npt.NDArray = field(default_factory=lambda: np.array([0]))       # radius at which curvature was evaluated
    
    Yp: float = 0                   # Pressure loss
    power:float = 0                 # Watts 
    power_distribution:npt.NDArray  # How power is divided by radius. Example: Equal distribution [0.33 0.33 0.33]. More at Tip [0.2,0.3,0.5]. More at Hub [0.6 0.5 ]
    P0_P:float = 0                  # Total to Static Pressure Ratio 
    Power_Type:PowerType
    euler_power:float = 0

    # Used for loss calculations
    _blade_to_blade_gap:float = 0.025 # Gap between blade in terms of percent chord.
    
    _aspect_ratio:float = 0.9 # 
    _pitch_to_chord:float = 0.7 # Pitch to chord ratio, used to determine number of blades and compute loss 
    
    _axial_chord:float = -1 
    _chord:float = -1 
    _stagger:float = 42
    _te_s:float = 0.08
    _tip_clearance:float = 0 # Clearance as a percentage of span or blade height

    _inlet_to_outlet_pratio = [0.06,0.7]
        
    @property
    def inlet_to_outlet_pratio(self) -> Tuple[float,float]:
        """This is what is varied by the optimization. 
        The range is between [0 and 1] but you should 

        Returns:
            List[float]: _description_
        """
        return self._inlet_to_outlet_pratio
    
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
    def pitch_to_chord(self) -> float:
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
        self._pitch_to_chord = val
    
    @property
    def solidity(self) -> float:
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
        self._pitch_to_chord = 1/val

    @property
    def beta1_metal(self) -> npt.NDArray:
        return np.degrees(self._beta1_metal)
    
    @property
    def beta2_metal(self) -> npt.NDArray:
        return np.degrees(self._beta2_metal)
    
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
    def chord(self) -> float:
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
    
    @property
    def num_blades(self) ->float:
        """returns the number of blades 

        Returns:
            float: number of blades
        """
        return int(2*np.pi*self.r.mean() / self.pitch)
    
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
        
    def __init__(self,axial_location:float,row_type:RowType=RowType.Stator,stage_id:int = 0):
        """Initializes the blade row to be a particular type

        Args:
            axial_location (float): Location of the blade row as a percentage of the total axial length
            row_type (RowType): Specifies the Type. Defaults to RowType.Stator
            power (float, optional): power . Defaults to 0.
            P0_P (float, optional): Total to Static Pressure Ratio
            stage_id (int, optional): ID of the stage so if you have 9 stages, the id could be 9. It's used to separate the stages. Each stage will have it's own unique degree of reaction 
        """
        self.row_type = row_type
        self.axial_location = axial_location 
        self.Yp = 0 # Loss
        self.stage_id = stage_id
    
    @beta1_metal.setter
    def beta1_metal(self,beta1_metal:List[float],percent:List[float]=[]):
        """Sets the leading edge metal angle for the blade

        Args:
            beta1_metal (List[float]): blade leading edge angle
            percent (List[float]): percent location of metal angles from hub to shroud.
        """
        self._beta1_metal = np.radians(convert_to_ndarray(beta1_metal))
        if len(percent) != len(beta1_metal):
            percent = np.linspace(0,1,len(self._beta1_metal)).tolist()
        self.beta1_metal_radii = convert_to_ndarray(percent)
        self.beta1_fixed = True
        self.beta1 = self.beta1_metal.copy()
        
    @beta2_metal.setter
    def beta2_metal(self,beta2_metal:List[float],percent:List[float]=[]):
        """Sets the trailing edge metal angle for the blade

        Args:
            beta2_metal (List[float]): Blade exit metal angle
            percent (List[float]): percent location of metal angles from hub to shroud.

        """
        self._beta2_metal = np.radians(convert_to_ndarray(beta2_metal))
        if len(percent) != len(beta2_metal):
            percent = np.linspace(0,1,len(self._beta2_metal)).tolist()
        self.beta2_metal_radii = convert_to_ndarray(percent)
        self.beta2_fixed = True
        self.beta2 = self._beta2_metal.copy()
        
        if self.row_type == RowType.Stator:
            self.alpha2 = self._beta2_metal.copy()
        
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
    def loss_model(self, model:Callable[[Any], Any]):
        """Add in custom loss model

        Args:
            model (function): custom loss function. Input will be of the format blade row

        Example:
        
        def mylossfunction(row:BladeRow) -> float
            code to do something with machine learning
            return pressure loss
        """
        self.loss_function = model
    
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
            "V":self.V.tolist(),
            "M":self.M.tolist(),
            "M_rel":self.M_rel.tolist(),
            "U":self.U.tolist(),
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
            "Yp":self.Yp,
            "Power":self.power,
            "P0_P": self.P0_P,
            "eta_total":self.eta_total,
            "eta_static":self.eta_static,
            "euler_power":self.euler_power,
            "axial_chord":self.axial_chord,
            "aspect_ratio":self.aspect_ratio,
            "area": self.area
        }

        return data

#* Some functions related to blade row 
def interpolate_streamline_radii(row:BladeRow,passage:Passage,num_streamlines:int=3):
    """Interpolate all quantities onto the streamline and allocates variables. 
    Run this after setting some initial conditions 

    Args:
        r_streamline (npt.NDArray): Radii describing the streamline 
        passage (Passage): Passage object describing the geometry of the hub and shroud
        num_streamlines (int): number of streamlines to consider

    Returns:
        (BladeRow): new row object with quantities interpolated
    """
    row.cutting_line,_,_ = passage.get_cutting_line(row.axial_location)
    row.x,row.r = row.cutting_line.get_point(np.linspace(0,1,num_streamlines))
    streamline_percent_length = np.sqrt((row.r-row.r[0])**2+(row.x-row.x[0])**2)/row.cutting_line.length
    
    # Flow angles 
    row._beta1_metal = row._beta1_metal.default_factory() if type(row._beta1_metal) == Field else row._beta1_metal

    if row.row_type==RowType.Stator:
        assert type(row._beta2_metal)!=Field,"Stator exit Flow angle must be set"

    row._beta1_metal = interpolate_quantities(row._beta1_metal,row.beta1_metal_radii,streamline_percent_length)
    row._beta2_metal = interpolate_quantities(row._beta2_metal,row.beta2_metal_radii,streamline_percent_length)
    row.beta1_metal_radii = streamline_percent_length
    row.beta2_metal_radii = streamline_percent_length
    
    if type(row.percent_hub_shroud) == Field: 
        row.percent_hub_shroud = streamline_percent_length
    else:
        row.percent_hub_shroud = streamline_percent_length # Reset the radii to streamline radii

    row.alpha1 = interpolate_quantities(row.alpha1,row.percent_hub_shroud,streamline_percent_length)
    row.alpha2 = interpolate_quantities(row.alpha2,row.percent_hub_shroud,streamline_percent_length)
    row.beta1 = interpolate_quantities(row.beta1,row.percent_hub_shroud,streamline_percent_length)
    row.beta2 = interpolate_quantities(row.beta2,row.percent_hub_shroud,streamline_percent_length)
    
    # Velocities 
    row.Vm = interpolate_quantities(row.Vm, row.percent_hub_shroud, streamline_percent_length)
    row.Vx = interpolate_quantities(row.Vx, row.percent_hub_shroud, streamline_percent_length)
    row.Vt = interpolate_quantities(row.Vt, row.percent_hub_shroud, streamline_percent_length)
    row.Vr = interpolate_quantities(row.Vr, row.percent_hub_shroud, streamline_percent_length)
    row.V = interpolate_quantities(row.V, row.percent_hub_shroud, streamline_percent_length)
    row.V2 = interpolate_quantities(row.V2, row.percent_hub_shroud, streamline_percent_length)
    row.M = interpolate_quantities(row.M, row.percent_hub_shroud, streamline_percent_length)
    row.M_rel = interpolate_quantities(row.M_rel, row.percent_hub_shroud, streamline_percent_length)
    row.U = interpolate_quantities(row.U, row.percent_hub_shroud, streamline_percent_length)
    row.W = interpolate_quantities(row.W, row.percent_hub_shroud, streamline_percent_length)
    row.Wt = interpolate_quantities(row.Wt, row.percent_hub_shroud, streamline_percent_length)

    # Total Quantities
    row.T0 = interpolate_quantities(row.T0,row.percent_hub_shroud,streamline_percent_length)
    row.T0_is = interpolate_quantities(row.T0,row.percent_hub_shroud,streamline_percent_length)
    row.P0 = interpolate_quantities(row.P0,row.percent_hub_shroud,streamline_percent_length)
    row.P0_stator_inlet = interpolate_quantities(row.P0_stator_inlet,row.percent_hub_shroud,streamline_percent_length)
    
    # Relative Quantities
    row.P0R = interpolate_quantities(row.P0R,row.percent_hub_shroud,streamline_percent_length)
    row.T0R = interpolate_quantities(row.T0R,row.percent_hub_shroud,streamline_percent_length)

    # Static Quantities 
    row.P = interpolate_quantities(row.P,row.percent_hub_shroud,streamline_percent_length)
    row.T = interpolate_quantities(row.T,row.percent_hub_shroud,streamline_percent_length)
    row.T_is = interpolate_quantities(row.T_is,row.percent_hub_shroud,streamline_percent_length)
    row.rho = interpolate_quantities(row.rho,row.percent_hub_shroud,streamline_percent_length)

    if row.row_type == RowType.Inlet:
        row.P0_fun = interp1d(row.percent_hub_shroud,row.P0) 
        row.T0_fun = interp1d(row.percent_hub_shroud,row.T0) 
    elif row.row_type == RowType.Outlet:
        row.P_fun = interp1d(row.percent_hub_shroud,row.P) 

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
        return interp1d(r,q,kind='linear')(r2)
    
def compute_gas_constants(row:BladeRow):
    """Calculates all the gas constants for a row. 
        This should be done if T or P change 
    """
    Tm = row.T.mean()
    Pm = row.P.mean()
    row.fluid.TP = Tm,Pm
    row.Cp = row.fluid.cp
    row.Cv = row.fluid.cv
    row.R = row.Cp-row.Cv
    row.gamma = row.Cp/row.Cv
    row.mu = sutherland(Tm) # type: ignore
    row.rho[:] = row.fluid.density
    # i = 0 
    # for T,P in zip(row.T,row.P):
    #     row.rho[i] = P/(T*row.R)
    #     i+=1
    return row