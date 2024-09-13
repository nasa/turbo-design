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
    def __init__(self,M:float,T0:Union[float,List[float]],P0:Union[float,List[float]],percent_radii:Union[float,List[float]],fluid:Solution,axial_location:float=0,beta:Union[float,List[float]]=[0]):
        """Initializes the inlet station. 
            Uses the beta and exit mach number to predict a value for Vm

        Args:
            M (float): Mach number at the inlet plane
            beta (Union[float,List[float]]): exit relative flow angle
            T0 (Union[float,List[float]]): Total Temperature Array
            P0 (Union[float,List[float]]): Total Pressure Array
            percent_radii (Union[float,List[float]]): Radius where total pressure and temperature are defined
            fluid (ct.Solution): Cantera mixture
            axial_location (float): Axial Location as a percentage of hub length
            beta (Union[float,List[float]], optional): Inlet flow angle in relative direction. Defaults to [].
        """
        super().__init__(row_type=RowType.Inlet,axial_location=axial_location,stage_id=-1)
        self.loss_function = None
        self.beta1 = convert_to_ndarray(beta)
        self.M = convert_to_ndarray(M)
        self.T0 = convert_to_ndarray(T0)
        self.P0 = convert_to_ndarray(P0)
        self.percent_hub_shroud = convert_to_ndarray(percent_radii)
        # if it's inlet alpha and beta are the same, relative flow angle = absolute. 
        self.beta2 = np.radians(convert_to_ndarray(beta))
        self.alpha1 = np.radians(convert_to_ndarray(beta))
        fluid.TP = self.T0.mean(),self.P0.mean()
        self.gamma = fluid.cp/fluid.cv
        
        self.T = self.T0 * 1/(1 + (self.gamma-1) * self.M**2)
        self.P = self.P0 * 1/(1 + (self.gamma-1) * self.M**2)**(self.gamma/(self.gamma-1))
        fluid.TP = self.T.mean(),self.P.mean()
        self.rho = convert_to_ndarray([fluid.density])
        self.fluid = fluid
        self.rpm = 0
        
        self.beta1_metal = [0] 
        self.beta2_metal = [0]
        self.P0_fun = interp1d(self.percent_hub_shroud,P0)
        self.T0_fun = interp1d(self.percent_hub_shroud,T0)
        
        
    def initialize_velocity(self,passage:Passage,num_streamlines:int):
        """Initialize velocity calculations. Assumes streamlines and inclination angles have been calculated 
     
        """
        # Perform Calculations on Velocity 
        Vm_prev = 0; Vm_err = 0 
        t,x,radius = passage.get_streamline(self.percent_hub_shroud)
        radius = radius[0]

        cutline,_,_ = passage.get_cutting_line(self.axial_location)
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
            
            self = compute_gas_constants(self)
            rho_mean = self.rho.mean()
            for i in range(len(self.massflow)-1):    
                tube_massflow = self.massflow[i+1]-self.massflow[i]
                if np.abs((self.x[i]-self.x[i-1]))<1E-12: # Axial Machines
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