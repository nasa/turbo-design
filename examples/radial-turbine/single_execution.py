from typing import Tuple
from turbodesign import PassageType
from turbodesign import TurbineSpool, Inlet, RowType, BladeRow, Passage, Outlet
from turbodesign.enums import MassflowConstraint
from turbodesign.coolant import Coolant
from turbodesign.loss.turbine import FixedPressureLoss
from cantera import Solution
from scipy.optimize import minimize_scalar
from scipy.interpolate import pchip
import matplotlib.pyplot as plt
import cantera as ct
from endwall import build_endwalls
import numpy.typing as npt 

def radial_turbine(P0:float,P0_P:float,
                   T0:float,Design_RPM:float,
                   hub:npt.NDArray,shroud:npt.NDArray,
                   blade_position:Tuple[float,float],
                   alpha2:float=50,beta3:float=-50,
                   massflow:float=1):
    """Performs a 1D Analysis of a radial turbine

    Args:
        P0 (float): Inlet Total Pressure
        P0_P (float): Total Pressure to static pressure ratio
        T0 (float): Inlet Total Temperature
        Design_RPM (float): Rotation Rate Rev/min
        inlet_hub_shroud_ratio (float): Inlet hub to shroud ratio
        outlet_hub_shroud_ratio (float): exit hub to shroud ratio
        radius (float): hub radius in meters
        rhub_out (float): outlet hub radius
        alpha2 (float): outlet metal angle stator
        beta3 (float): outlet metal angle rotor
        massflow (float): massflow rate in kg/s

    """
    
    # import matplotlib.pyplot as plt
    # plt.plot(hub[:,0],hub[:,1],shroud[:,0],shroud[:,1],color='k',label='Passage')
    # plt.axis('equal')
    # plt.ylim([-max(hub[:,1])*1.1,1.1*max(hub[:,1])])
    # plt.legend()
    # plt.savefig('passage.jpg',dpi=150)
    
    
    # spool.plot()
    

