from dataclasses import dataclass, field
from enum import Enum
from typing import List, Optional, Union

from .enums import RowType
from .bladerow import BladeRow
from .arrayfuncs import convert_to_ndarray
import numpy as np
import numpy.typing as npt 
from scipy.interpolate import interp1d

class OutletType(Enum):
    static_pressure = 1
    total_pressure = 2
    massflow_static_pressure = 3 

class Outlet(BladeRow):
    P_fun:interp1d
    P0_fun:interp1d
    num_streamlines:int
    outlet_type: OutletType = OutletType.static_pressure
    percent_hub_shroud:npt.NDArray

    def __init__(self,num_streamlines:int=3,location:float=1):
        """Initialize the outlet with streamlines and a location as a percentage along the hub

        Args:
            num_streamlines (int, optional): _description_. Defaults to 3.
            location (float, optional): Location as percentage along hub curve. Defaults to 1.
        """
        super().__init__(hub_location=location, row_type=RowType.Outlet, stage_id=-1)
        self.loss_function = None
        self.location = location
        self.num_streamlines = num_streamlines
        # Evenly spaced from hub (0) to shroud (1)
        self.percent_hub_shroud = np.linspace(0, 1, self.num_streamlines)
        
    def init_static(self,P:Union[List[float],float],percent_radii:Union[List[float],float],massflow:float | None = None):
        """Initialize turbine inputs

        Args:
            P (float): Exit static pressure [Pa]
            percent_radii (Union[List[float],float]): percent radii where the exit static pressure is defined.
            massflow (float | None): initialize with a massflow. Note this is only valid if you set the solver to solve to MassflowConstraint.AngleMatch
        """
        self.percent_hub_shroud = convert_to_ndarray(percent_radii)
        if len(self.percent_hub_shroud)==1:
            self.percent_hub_shroud = np.linspace(0,1,self.num_streamlines)
        self.P = convert_to_ndarray(P)
        self.P = self.P[0]+0*self.percent_hub_shroud*0
        self.P_fun = interp1d(self.percent_hub_shroud,self.P)
        
        if massflow is not None:
            self.total_massflow = massflow
            self.outlet_type = OutletType.massflow_static_pressure
        else:
            self.outlet_type = OutletType.static_pressure    
        
    def init_total(self,P0:Union[List[float],float],percent_radii:Union[List[float],float]):
        """Initialize compressor inputs

        Args:
            P0 (Union[List[float],float]): Exit Total Pressure (this will be matched)
            percent_radii (Union[List[float],float]): percent radii where exit total pressure is defined
        """

        self.percent_hub_shroud = convert_to_ndarray(percent_radii)
        self.IsCompressor = True
        self.P0 = convert_to_ndarray(P0)
        self.P0 = self.P0[0]+0*self.percent_hub_shroud*0
        self.P0_fun = interp1d(self.percent_hub_shroud,self.P0)
        self.IsCompressor = True
        self.outlet_type = OutletType.total_pressure
        self.P = self.P0 # Do this first but we will adjust
        
    def transfer_quantities(self,upstream:BladeRow):
        """Transfer quantities from upstream row to outlet while maintaining the outlet static pressure

        Args:
            upstream (BladeRow): Upstream row, for turbines this is a rotor, for compressors this is a stator.
        """
        self.__dict__ = upstream.__dict__.copy() # Copies P and hub shroud percentage
        if self.outlet_type == OutletType.static_pressure:
            self.P_fun = interp1d(self.percent_hub_shroud,self.P)
        else: # Compressor
            self.P0_fun = interp1d(self.percent_hub_shroud,self.P0)
        self.row_type = RowType.Outlet
        
    
    def get_static_pressure(self,percent_hub_shroud:Union[float,npt.NDArray]):
        """Returns the static pressure at a certain percent hub_shroud

        Args:
            percent_hub_shroud (Union[float,npt.NDArray]): value or array from 0 to 1 where you want the exit static pressure

        Returns:
            npt.NDArray: Returns an array of static pressure
        """
        if type(percent_hub_shroud) == float:
            return float(self.P_fun(percent_hub_shroud))
        else:
            return self.P_fun(percent_hub_shroud)
    
    def get_total_pressure(self, percent_hub_shroud:Union[float,npt.NDArray]):
        """Returns the total pressure at a certain percent hub shroud 

        Args:
            percent_hub_shroud (Union[float,npt.NDArray]): value or array from 0 to 1 where you want the exit total pressure

        Returns:
            npt.NDArray: Returns an array of total pressure
        """
        if type(percent_hub_shroud) == float:
            return float(self.P0_fun(percent_hub_shroud))
        else:
            return self.P0_fun(percent_hub_shroud)
        
