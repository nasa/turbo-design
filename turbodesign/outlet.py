from dataclasses import dataclass, field
from typing import List, Union
from .enums import RowType
from .bladerow import BladeRow, compute_gas_constants, interpolate_quantities
from .arrayfuncs import convert_to_ndarray
import numpy as np
import numpy.typing as npt 
import copy 
from scipy.interpolate import interp1d


class Outlet(BladeRow):
    P_fun:interp1d
    
    def __init__(self,P:Union[float,List[float]],percent_radii:List[float],num_streamlines:int=3):
        """Initialize the outlet

        Args:
            P (Union[float,List[float]]): List of static pressure profile at outlet
            percent_radii (List[float]): Percent Radius from 0 to 1 or just put 0.5 and it uses all of it 
        """
        self.P = convert_to_ndarray(P)
        self.percent_hub_shroud = convert_to_ndarray(percent_radii)
        if len(self.percent_hub_shroud)==1:
            self.percent_hub_shroud = np.arange(0,1,num_streamlines)
        self.P_fun = interp1d(self.percent_hub_shroud,self.P)
        self.row_type = RowType.Outlet
        self.loss_function = None
    
    
    def transfer_quantities(self,upstream:BladeRow):
        """Transfer quantities from upstream row to outlet while maintaining the outlet static pressure

        Args:
            upstream (BladeRow): Upstream row, for turbines this is a rotor, for compressors this is a stator. 
        """
        self.__dict__ = upstream.__dict__.copy() # Copies P and hub shroud percentage  
        self.P_fun = interp1d(self.percent_hub_shroud,self.P)
        self.row_type = RowType.Outlet
        
    
    def get_static_pressure(self,percent_hub_shroud:Union[float,npt.NDArray]):
        """Returns the static pressure at a certain percent hub_shroud

        Args:
            percent_hub_shroud (Union[float,npt.NDArray]): _description_

        Returns:
            _type_: _description_
        """
        if type(percent_hub_shroud) == float:
            return float(self.P_fun(percent_hub_shroud))
        else:
            return self.P_fun(percent_hub_shroud)
        
        