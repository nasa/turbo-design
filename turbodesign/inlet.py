from dataclasses import dataclass, field
from typing import List, Union
from .enums import RowType
from .bladerow import BladeRow, compute_gas_constants, interpolate_quantities
from .arrayfuncs import convert_to_ndarray
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
    
    def __init__(self,M:float,T0:Union[float,List[float]],
                 P0:Union[float,List[float]],
                 location:float=0,
                 beta:Union[float,List[float]]=[0],
                 percent_radii:Union[float,List[float]]=[0.5]):
        """Initializes the inlet station. 
            Uses the beta and exit mach number to predict a value for Vm

        Args:
            M (float): Mach number at the inlet plane
            T0 (Union[float,List[float]]): Total Temperature Array
            P0 (Union[float,List[float]]): Total Pressure Array
            percent_radii (Union[float,List[float]]): Radius where total pressure and temperature are defined
            location (float): Location as a percentage of hub length
            beta (Union[float,List[float]], optional): Inlet flow angle in relative direction. Defaults to [].

        """
        super().__init__(row_type=RowType.Inlet,location=location,stage_id=-1)
        self.beta1 = convert_to_ndarray(beta)
        self.M = convert_to_ndarray(M)
        self.T0 = convert_to_ndarray(T0)
        self.P0 = convert_to_ndarray(P0)
        self.percent_hub_shroud = convert_to_ndarray(percent_radii)
   
    def initialize_inputs(self,num_streamlines:int=5):
        self.M = interpolate_quantities(self.M, self.percent_hub_shroud, np.linspace(0,1,num_streamlines))
        self.P0 = interpolate_quantities(self.P0,self.percent_hub_shroud, np.linspace(0,1,num_streamlines))
        self.T0 = interpolate_quantities(self.T0,self.percent_hub_shroud, np.linspace(0,1,num_streamlines)) 
        # if it's inlet alpha and beta are the same, relative flow angle = absolute. 
        self.beta1 = interpolate_quantities(self.beta1,self.percent_hub_shroud, np.linspace(0,1,num_streamlines)) 
        self.beta2 = np.radians(convert_to_ndarray(self.beta1))
        self.alpha1 = np.radians(convert_to_ndarray(self.beta1))         
        
    def initialize_fluid(self,fluid:Solution=None,R:float=287.15,gamma:float=1.4,Cp:float=1024):
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
            self.T = self.T0 * 1/(1 + (self.gamma-1) * self.M**2)
            self.P = self.P0 * 1/(1 + (self.gamma-1) * self.M**2)**(self.gamma/(self.gamma-1))
            fluid.TP = self.T.mean(),self.P.mean()
            self.rho = convert_to_ndarray([fluid.density])
        else:
            self.Cp = Cp
            self.gamma = gamma
            self.R = R
            self.T = self.T0 * 1/(1 + (self.gamma-1) * self.M**2)
            self.P = self.P0 * 1/(1 + (self.gamma-1) * self.M**2)**(self.gamma/(self.gamma-1))
            self.rho = self.P/(self.R*self.T)

        self.rpm = 0
        self.beta1_metal = [0] 
        self.beta2_metal = [0]
        if len(self.percent_hub_shroud) == 1:
            self.percent_hub_shroud = np.linspace(0,1,2)
            self.P0 = self.percent_hub_shroud*0+self.P0[0]
            self.T0 = self.percent_hub_shroud*0+self.T0[0]
        self.P0_fun = interp1d(self.percent_hub_shroud,self.P0)
        self.T0_fun = interp1d(self.percent_hub_shroud,self.T0)
        self.mprime = [0]
        
    def initialize_velocity(self,passage:Passage,num_streamlines:int):
        """Initialize velocity calculations. Assumes streamlines and inclination angles have been calculated 
            Call this before performing calculations
            
        Args:
            passage (Passage): Passage object
            num_streamlines (int): number of streamlines
        
        """
        # Perform Calculations on Velocity 
        Vm_prev = 0; Vm_err = 0 

        cutline,_,_ = passage.get_cutting_line(self.location)
        self.x,self.r = cutline.get_point(np.linspace(0,1,num_streamlines))
        for _ in range(10):
            T0_T = (1+(self.gamma-1)/2 * self.M**2)

            self.Vm = self.M**2 * self.gamma*self.R*self.T0/T0_T \
                        / (1+np.cos(self.phi)**2 * np.tan(self.alpha1)**2)

            self.Vm = np.sqrt(self.Vm)
            self.T = self.T0/T0_T
            self.P = self.P0/(T0_T)**(self.gamma/(self.gamma-1))
            self.rho = self.P/(self.R*self.T)
            
            self.Vx = self.Vm * np.cos(self.phi)
            self.Vt = self.Vm * np.cos(self.phi) * np.tan(self.beta1)
            self.V = np.sqrt(self.Vm**2 + self.Vt**2)        
            self.Vr = self.Vm * np.sin(self.phi) 
            
            compute_gas_constants(self)
            rho_mean = self.rho.mean()
            for i in range(len(self.massflow)-1):    
                tube_massflow = self.massflow[i+1]-self.massflow[i]
                if np.abs((self.x[-1]-self.x[0]))<1E-5: # Axial Machines
                    self.Vm[i+1] = tube_massflow/(rho_mean*np.pi*(self.r[i+1]**2-self.r[i]**2))
                else:   # Radial Machines
                    dx = self.x[i]-self.x[i-1]
                    S = (self.r[i]-self.r[i-1])
                    C = np.sqrt(1+((self.r[i]-self.r[i-1])/dx)**2)
                    area = 2*np.pi*C*(S/2*dx**2+self.r[i-1]*dx)
                    self.Vm[i+1] = tube_massflow/(rho_mean*area)
            self.Vm[0] = 1/(len(self.Vm)-1)*self.Vm[1:].sum()
            
            self.M = self.Vm /np.sqrt(self.gamma*self.R*self.T)
            Vm_err = np.max(abs(self.Vm-Vm_prev)/self.Vm)
            Vm_prev = self.Vm
            if Vm_err < 1E-4:
                break
        
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