from dataclasses import dataclass, field
from optparse import Option
from typing import List, Optional, Union

from .enums import RowType
from .bladerow import BladeRow, compute_gas_constants
from .arrayfuncs import convert_to_ndarray, safe_interpolate
import numpy as np 
from cantera import Solution
from .passage import Passage
from scipy.interpolate import interp1d
import numpy.typing as npt 

class Inlet(BladeRow):
    """Station defined at Inlet

    Inherits:
        (BladeRow): Defines the properties of the blade row
    """
    fun: interp1d
    static_defined: bool
    def __init__(self, 
                 hub_location:float=0,
                 shroud_location:Optional[float]=None,
                 alpha:Union[float,List[float]]=[0]):
        """Initializes the inlet station. 
            Uses the beta and exit mach number to predict a value for Vm

        Args:
            location (float): Location as a percentage of hub length
            beta (Union[float,List[float]], optional): Inlet flow angle in relative direction. Defaults to [].

        """
        super().__init__(row_type=RowType.Inlet,hub_location=hub_location,shroud_location=shroud_location,stage_id=-1)
        self.beta1 = convert_to_ndarray([0.0])
        # Default absolute angles to zero to avoid attribute errors during interpolation
        self.alpha1 = convert_to_ndarray([0.0])
        self.alpha2 = convert_to_ndarray(alpha)
        self.beta2 = convert_to_ndarray([0.0])
                   
    
    def init_static(self,P:Union[float,List[float]],T:Union[float,List[float]],M:Union[float,List[float]],percent_radii:Union[float,List[float]]=[0.5]):
        """Initializes the inlet with static quantities at the inlet

        Args:
            P (Union[float,List[float]]): Static Pressure either as a float or array
            T (Union[float,List[float]]): Static Temperature either as a float or array 
            M (Union[float,List[float]]): Mach Number either as a float or array 
            percent_radii (Union[float,List[float]], optional): Percent radii where P,T, and M are defined. Defaults to [0.5].
        """
        self.P = convert_to_ndarray(P)
        self.M = convert_to_ndarray(M)
        self.T = convert_to_ndarray(T)
        self.static_defined = True
        self.percent_hub_shroud = convert_to_ndarray(percent_radii)
        
    def init_total(
        self,
        P0:Union[float,List[float]],
        T0:Union[float,List[float]],
        M:Union[float,List[float]],
        percent_radii:Optional[Union[float,List[float]]]=None
    ):
        """Initializes the inlet with total quantities at the inlet

        Args:
            P0 (Union[float,List[float]]): Total Pressure either as a float or array 
            T0 (Union[float,List[float]]): Total Temperature either as a float or array 
            M (Union[float,List[float]]): Mach Number either as a float or array 
            percent_radii (Optional[Union[float,List[float]]], optional): Percent radii where P0, T0, and M are defined. Defaults to `None`, which uses evenly spaced radii when multiple values exist or `[0.5]` otherwise.
        """
        self.P0 = convert_to_ndarray(P0)
        self.T0 = convert_to_ndarray(T0) 
        self.M = convert_to_ndarray(M)
        if percent_radii is None:
            percent_radii = convert_to_ndarray([0.5]) # type: ignore
        if len(self.M)>1: 
            percent_radii = np.linspace(0,1,len(self.M))  # type: ignore
        self.static_defined = False
        self.percent_hub_shroud = percent_radii
        
    def __interpolate_quantities__(self,num_streamlines:int=5):
        """Initializes the inputs 
        
        Args:
            num_streamlines (int, optional): _description_. Defaults to 5.
            IsCompressor (bool, optional): This is if static pressure is defined at the inlet and total pressure at the outlet. Defaults to False.
        """
        dst = np.array([0.5]) if num_streamlines <= 1 else np.linspace(0,1,num_streamlines)
        self.M = safe_interpolate(self.M, self.percent_hub_shroud, dst)
        if self.static_defined: # This comes from the initialization
            self.P = safe_interpolate(self.P, self.percent_hub_shroud, dst)
        else:
            self.P0 = safe_interpolate(self.P0, self.percent_hub_shroud, dst)
        self.T0 = safe_interpolate(self.T0, self.percent_hub_shroud, dst)
        # Angles: default to 0 if unspecified
        self.beta1 = safe_interpolate(self.beta1, self.percent_hub_shroud, dst, radians=True)
        self.beta2 = safe_interpolate(self.beta2, self.percent_hub_shroud, dst, radians=True)
        self.alpha1 = safe_interpolate(self.alpha1, self.percent_hub_shroud, dst, radians=True)
        self.alpha2 = safe_interpolate(self.alpha2, self.percent_hub_shroud, dst, radians=True)
        
    def __initialize_fluid__(self,fluid:Optional[Solution]=None,R:float=287.15,gamma:float=1.4,Cp:float=1024):
        """Initialize the inlet using the fluid. This function should be called by a class that inherits from spool

        Args:
            fluid (Solution, optional): Cantera fluid object. Defaults to None.
            R (float, optional): Ideal Gas Constant. Defaults to 287.15 J/(Kg K) for air
            gamma (float, optional): _description_. Defaults to 1.4.
            Cp (float, optional): _description_. Defaults to 1024 J/(Kg K).
        """
        self.loss_function = None
        
        if fluid:
            fluid.TP = self.T0.mean(),self.P0.mean()
            self.gamma = fluid.cp/fluid.cv
            if self.static_defined:
                self.P0 = self.P * (1+(self.gamma-1)/2 * self.M**2) ** (self.gamma/(self.gamma-1))
                self.T0 = self.T * (1+(self.gamma-1)/2 * self.M**2)
            else:
                self.P = self.P0 * 1/(1 + (self.gamma-1) * self.M**2)**(self.gamma/(self.gamma-1))
            self.T = self.T0 * 1/(1 + (self.gamma-1) * self.M**2)
            fluid.TP = self.T.mean(),self.P.mean()
            self.rho = convert_to_ndarray([fluid.density])
        else:
            self.Cp = Cp
            self.gamma = gamma
            self.R = R
            self.T = self.T0 * 1/(1 + (self.gamma-1) * self.M**2)
            if self.static_defined:
                self.P0 = self.P * (1+(self.gamma-1)/2 * self.M**2) ** (self.gamma/(self.gamma-1)) 
                self.T0 = self.T * (1+(self.gamma-1)/2 * self.M**2)
            else:
                self.P = self.P0 * 1/(1 + (self.gamma-1) * self.M**2)**(self.gamma/(self.gamma-1))
            self.rho = self.P/(self.R*self.T)

        self.beta1_metal = [0] 
        self.beta2_metal = [0]
        if len(self.percent_hub_shroud) == 1:
            self.percent_hub_shroud = np.linspace(0,1,2)
            self.P0 = self.percent_hub_shroud*0+self.P0[0]
            self.T0 = self.percent_hub_shroud*0+self.T0[0]
        self.P0_fun = interp1d(self.percent_hub_shroud,self.P0)
        self.T0_fun = interp1d(self.percent_hub_shroud,self.T0)
        self.mprime = [0] # type: ignore
        
    def __initialize_velocity__(self,passage:Passage,num_streamlines:int):
        """Initialize velocity calculations. Assumes streamlines and inclination angles have been calculated 
            Call this before performing calculations
            
        Args:
            passage (Passage): Passage object
            num_streamlines (int): number of streamlines
        
        """
        # Perform Calculations on Velocity 
        Vm_prev = 0; Vm_err = 0 

        cutline,_,_ = passage.get_cutting_line(t_hub=self.location,t_shroud=self.shroud_location)
        t_span = np.array([0.5]) if num_streamlines <= 1 else np.linspace(0,1,num_streamlines)
        self.x,self.r = cutline.get_point(t_span)
        for _ in range(2):
            T0_T = (1+(self.gamma-1)/2 * self.M**2)
            
            self.Vm = self.M**2 * self.gamma*self.R*self.T0/T0_T \
                        / (1+np.cos(self.phi)**2 * np.tan(self.alpha2)**2)

            self.Vm = np.sqrt(self.Vm)
            self.T = self.T0/T0_T
            self.P = self.P0/(T0_T)**(self.gamma/(self.gamma-1))
            self.rho = self.P/(self.R*self.T)
            
            self.Vx = self.Vm * np.cos(self.phi)
            self.Vt = self.Vm * np.cos(self.phi) * np.tan(self.alpha2)
            self.V = np.sqrt(self.Vm**2 + self.Vt**2)        
            self.Vr = self.Vm * np.sin(self.phi) 
            
            compute_gas_constants(self)
            rho_mean = self.rho.mean()
            Vm_tube = np.zeros(max(len(self.massflow)-1, 1))
            if len(self.massflow) <= 1:
                Vm_tube[0] = float(self.Vm.mean())
            # Compute tube-averaged Vm from massflow differences
            for i in range(1, len(self.massflow)):
                tube_massflow = self.massflow[i]-self.massflow[i-1]
                rho_bar = rho_mean if len(self.rho) == 1 else 0.5 * (self.rho[i] + self.rho[i-1])
                if np.abs((self.x[-1]-self.x[0]))<1E-5: # Axial Machines
                    area = np.pi*(self.r[i]**2-self.r[i-1]**2)
                else:   # Radial Machines
                    dx = self.x[i]-self.x[i-1]
                    S = (self.r[i]-self.r[i-1])
                    C = np.sqrt(1+((self.r[i]-self.r[i-1])/dx)**2)
                    area = 2*np.pi*C*(S/2*dx**2+self.r[i-1]*dx)
                Vm_tube[i-1] = tube_massflow/(rho_bar*area + 1e-12)

            # Recover per-streamline Vm; handle single-streamline as meanline
            if len(self.Vm) <= 1:
                self.Vm = np.array([Vm_tube[0] if len(Vm_tube) else rho_mean*0])
            else:
                self.Vm[0] = Vm_tube[0]
                for i in range(1, len(self.Vm)):
                    self.Vm[i] = 2 * Vm_tube[i-1] - self.Vm[i-1]
            
            self.M = self.V /np.sqrt(self.gamma*self.R*self.T)
            Vm_err = np.max(abs(self.Vm-Vm_prev)/self.Vm)
            Vm_prev = self.Vm
            if Vm_err < 1E-4:
                break
        
        if num_streamlines <= 1:
            Area = passage.get_area(self.location)
        else:
            Area = 0
            for j in range(1,num_streamlines):
                if np.abs((self.x[j]-self.x[j-1]))<1E-12: # Axial Machines  
                    Area += np.pi*(self.r[j]**2-self.r[j-1]**2)
                else:   # Radial Machines
                    dx = self.x[j]-self.x[j-1]
                    S = (self.r[j]-self.r[j-1])
                    C = np.sqrt(1+((self.r[j]-self.r[j-1])/dx)**2)
                    Area += 2*np.pi*C*(S/2*dx**2+self.r[j-1]*dx)
        self.calculated_massflow = self.rho.mean()*self.Vm.mean() * Area


    def get_total_pressure(self,percent_hub_shroud:Union[float,npt.NDArray]):
        """Returns the static pressure at a certain percent hub_shroud

        Args:
            percent_hub_shroud (Union[float,npt.NDArray]): _description_

        Returns:
            _type_: _description_
        """
        if type(percent_hub_shroud) == float:
            return float(self.P0_fun(percent_hub_shroud))
        else:
            return self.P0_fun(percent_hub_shroud)
    
    def get_total_temperature(self,percent_hub_shroud:Union[float,npt.NDArray]):
        """Returns the static pressure at a certain percent hub_shroud

        Args:
            percent_hub_shroud (Union[float,npt.NDArray]): _description_

        Returns:
            _type_: _description_
        """
        if type(percent_hub_shroud) == float:
            return float(self.T0_fun(percent_hub_shroud))
        else:
            return self.T0_fun(percent_hub_shroud)
