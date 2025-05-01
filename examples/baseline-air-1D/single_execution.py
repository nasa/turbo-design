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
    
    passage = Passage(hub[:,0],hub[:,1],
                    shroud[:,0],shroud[:,1],
                    passageType=PassageType.Centrifugal)
    # passage.get_streamline(0)
    # passage.plot_cuts()

    #%% Design Conditions 
    P = P0/P0_P         # Outlet Static Pressure [Pascal]


    #%% Defining the Inlet
    inlet = Inlet(M=0.1,
                    P0=[P0],
                    T0=[T0],
                    beta=[0],
                    percent_radii=0.5,
                    location=0)
    
    outlet = Outlet(P=P,percent_radii=0.5,num_streamlines=5)

    stator = BladeRow(row_type=RowType.Stator, location=blade_position[0])
    stator.R = 287.15
    stator.gamma = 1.35
    stator.Cp = stator.gamma*stator.R/(stator.gamma-1)
    
    rotor = BladeRow(row_type=RowType.Rotor, location=blade_position[1])
    rotor.R = 287.15
    rotor.gamma = 1.35
    rotor.Cp = stator.gamma*stator.R/(stator.gamma-1)
    
    # if coolant has 0 massflow then it isn't used
    stator.coolant = Coolant(T0=T0*0.555556,P0=5E5,Cp=900,massflow_percentage=0) 
    rotor.coolant = Coolant(T0=T0*0.555556,P0=5E5,Cp=900,massflow_percentage=0)

    # Add in turning angles
    stator.beta2_metal = [alpha2,alpha2,alpha2] # Angle, hub,mean,tip
    stator.loss_model = FixedPressureLoss(0.0)

    rotor.beta2_metal = [beta3,beta3,beta3] # Angle, hub,mean,tip
    rotor.loss_model = FixedPressureLoss(0.15669278543371953)

    #%% Initialize the Spool
    spool = TurbineSpool(passage=passage,
                rpm=Design_RPM, 
                num_streamlines=3,
                massflow=massflow,
                rows=[inlet,stator,rotor,outlet],
                fluid=None)

    spool.adjust_streamlines = False 

    spool.massflow_constraint = MassflowConstraint.BalanceMassFlow
    
    spool.solve() # This also initializes streamlines
    spool.plot_velocity_triangles()
    spool.export_properties("output.json")
    # spool.plot()
    

