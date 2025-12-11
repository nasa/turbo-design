'''
    1D meanline example from 
    Turbine Passage Design Methodology to Minimize Entropy Production—A Two-Step Optimization Strategy. 
    https://doi.org/10.3390/e21060604 
    
    In this example the blade exit angles are fixed and only degree of reaction changes between the stage to match the massflow
'''
#%% Import Library
from pathlib import Path
from turbodesign import PassageType
from turbodesign.row_factory import make_rotor_row, make_stator_row
from turbodesign import TurbineSpool, Inlet, RowType, BladeRow, Passage,Outlet
from turbodesign.enums import MassflowConstraint
from turbodesign import Coolant
from turbodesign.loss.turbine import FixedPressureLoss, KackerOkapuu
import numpy as np 
from cantera import Solution

#%% Define the Passage 
# Geometry from OptTurb
rmean = 0.389
H1 = 0.04  
H2 = 1.159*H1
H3 = 1.317*H2
cax = (H1+H2+H3)/3 
                        # Inlet, Stator Inlet, Stator Exit, Rotor Exit
rhub = [rmean-H1/2,rmean-H1/2,rmean-H2/2,rmean-H3/2]
rshroud = [rmean+H1/2,rmean+H1/2,rmean+H2/2,rmean+H3/2]
xhub = np.array([0, cax, 2*cax, 3*cax])
xshroud = np.array([0, cax, 2*cax, 3*cax])
axial_len = xhub[-1]-xhub[0]

passage = Passage(xhub,rhub,
                 xshroud,rshroud,
                 passageType=PassageType.Axial)
#%% Design Conditions 
Design_RPM = 7500
massflow = 35.9 # kg/s
P0 = 500000     # Pascal 
T0 = 676.3      # Kelvin

# Fluid
fluid = Solution('air.yaml')
fluid.TP = T0, P0 # Use pascal for cantera
print(f"Coefficient of Pressure [J/Kg] {fluid.cp:0.4f}")


#%% Defining the Inlet
inlet = Inlet(beta=[0], hub_location=0)
inlet.init_total(
    P0=[P0],
    T0=[T0],
    M=[0.2],
    percent_radii=[0.5],
)
outlet = Outlet(P=P0/3.96,percent_radii=0.5,num_streamlines=3)

#%% Define Blade Rows 
# Axial location is a percentage along the hub where row exit is defined
stator1 = make_stator_row(row_type=RowType.Stator, hub_location=2*cax/axial_len)
rotor1 = make_rotor_row(row_type=RowType.Rotor, hub_location=3*cax/axial_len)

stator1.axial_chord = cax # Set an axial chord
rotor1.axial_chord = cax

# Coolant Definition: Use Kelvin and Pascal
stator1.coolant = Coolant(T0=616*0.555556, P0= 50.6 * 6894.76, massflow_percentage=0,Cp=fluid.cp) 
rotor1.coolant = Coolant(T0=622*0.555556, P0=50.3 * 6894.76,massflow_percentage=0,Cp=fluid.cp)

# Add in turning angles
stator1.beta2_metal = [73,73,73] # Angle, hub,mean,tip
stator1.loss_model = KackerOkapuu()
rotor1.loss_model = KackerOkapuu()
rotor1.beta2_metal = [-67.6,-67.6,-67.6] # Angle, hub,mean,tip

#%% Initialize the TurbineSpool
spool = TurbineSpool(passage=passage,
            massflow=massflow,
            inlet=inlet,
            outlet=outlet,
            rows=[stator1,rotor1],
            rpm=Design_RPM,
            num_streamlines=3)
spool.fluid = fluid
spool.massflow_constraint = MassflowConstraint.PressureBalance # Fixes the exit angle and changes degree of reaction
# spool.plot_geometry()
spool.solve() # This also initializes streamlines
export_path = Path(__file__).resolve().parent / "optturb.json"
spool.export_properties(str(export_path))
spool.plot()
spool.plot_velocity_triangles()
print('check')
