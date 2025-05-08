'''
    1D meanline example from 
    Turbine Passage Design Methodology to Minimize Entropy Production—A Two-Step Optimization Strategy. 
    https://doi.org/10.3390/e21060604 
    
    In this example the blade exit angles are fixed and only degree of reaction changes between the stage to match the massflow
'''
#%% Import Library
from turbodesign import PassageType
from turbodesign import TurbineSpool, Inlet, RowType, BladeRow, Passage, Outlet
from turbodesign.enums import MassflowConstraint
from turbodesign.coolant import Coolant
from turbodesign.loss.turbine import FixedPressureLoss
from pyturbo.helper import line2D
import numpy as np 
from cantera import Solution

#! Check the Cp!! That could be why massflow is not converging
#%% Define the Passage 

hub_shroud = np.array([[-3, 6, -3, 7.5],
                       [ 7, 6,  7, 7.5]]) * 0.0254
hub = hub_shroud[:,:2]
shroud = hub_shroud[:,2:]
row1X = np.array([[-2.5, 6, -2.5, 7.5],[0.358, 6, 0.358, 7.5]]) * 0.0254 # Inch to meter
row2X = np.array([[0.358, 6, 0.358, 8],[2.68, 6, 2.68, 7.5]]) * 0.0254 
row3X = np.array([[2.68, 6, 2.68, 7.5], [6, 4, 6, 7.5]]) * 0.0254 
axial_len = row3X[-1,0] - row1X[0,0]
cax1 = row1X[-1,0]-row1X[0,0]
cax2 = row2X[-1,0]-row2X[0,0]
cax3 = row3X[-1,0]-row3X[0,0]

hub = line2D(hub[0,:],hub[1,:])
xhub,rhub = hub.get_point(np.linspace(0,1,100))

shroud = line2D(shroud[0,:],shroud[1,:])
xshroud,rshroud = shroud.get_point(np.linspace(0,1,100))


passage = Passage(xhub.flatten(),rhub.flatten(),
                 xshroud.flatten(),rshroud.flatten(),
                 passageType=PassageType.Axial)
#%% Design Conditions 
Design_RPM = -9000
massflow =9    # kg/s, Initial guess
P0 = 413.617*1000    # Pascal 
T0 = 555.5        # Kelvin

# Fluid
air = Solution('air.yaml')
air.TP = T0, P0 # Use pascal for cantera
print(f"Coefficient of Pressure [J/Kg] {air.cp:0.4f}")

#%% Defining the Inlet
inlet = Inlet(M=0.02, 
                 P0=[P0],
                 T0=[T0], 
                 beta=[0],
                 percent_radii=0.5,
                 location=0)

outlet = Outlet(P=206.799*1000,percent_radii=0.5,num_streamlines=5)

#%% Define Blade Rows 
# Axial location is a percentage along the hub where row exit is defined
stator1 = BladeRow(row_type=RowType.Stator,location=cax1/axial_len)
rotor1 = BladeRow(row_type=RowType.Rotor, location=(cax1+cax2)/axial_len)
stator2 = BladeRow(row_type=RowType.Stator,location=1)

# stator1.gamma = 1.38
# stator1.Cp = 1042.8
# stator1.R = 287.15

# rotor1.gamma = 1.38
# rotor1.Cp = 1042.8
# rotor1.R = 287.15

# stator2.gamma = 1.38
# stator2.Cp = 1042.8
# stator2.R = 287.15

stator1.axial_chord = cax1 # Set an axial chord
rotor1.axial_chord = cax2
stator2.axial_chord = cax3 

stator1.stage_id = 0; rotor1.stage_id = 0
stator2.stage_id = 0; 

# Coolant Definition: Use Kelvin and Pascal. Coolant only needs P0, T0, massflow, and Cp
stator1.coolant = Coolant(T0=T0*0.5, P0 = P0 * 6894.76, massflow_percentage=0, Cp=air.cp) 
rotor1.coolant = Coolant(T0*0.5, P0 = P0 * 6894.76, massflow_percentage=0, Cp=air.cp)
stator2.coolant = Coolant(T0=T0*0.5, P0 = P0 * 6894.76, massflow_percentage=0, Cp=air.cp) 

# Add in turning angles
stator1.beta2_metal = [-67.1,-67.1,-67.1,-67.1,-67.1]                  # Alpha2
# Including beta1_metal angles could be useful for loss calculations if the custom loss function uses inlet metal angle
rotor1.beta1_metal = [-30,-30,-30,-30,-30]              # Beta2
rotor1.beta2_metal = [62.7,62.7,62.7,62.7,62.7]         # Beta25

stator2.beta2_metal = [-65.1,-65.1,-65.1,-65.1,-65.1]   # Alpha3

# These are all guessed values 
stator1.loss_model = FixedPressureLoss(0.048)
rotor1.loss_model = FixedPressureLoss(0.108)
stator2.loss_model = FixedPressureLoss(0.082)
stator1.inlet_to_outlet_pratio = (0.01,0.9)
rotor1.inlet_to_outlet_pratio = (0.01,0.9)
stator2.inlet_to_outlet_pratio = (0.01,0.9)

#%% Initialize the Spool
spool = TurbineSpool(passage=passage,
            rpm=Design_RPM, 
            num_streamlines=5, 
            massflow=massflow, 
            fluid=air,
            rows=[inlet,stator1,rotor1,stator2,outlet])
spool.massflow_constraint = MassflowConstraint.BalanceMassFlow # Fixes the exit angle and changes degree of reaction
# spool.plot_geometry()
spool.solve() # This also initializes streamlines
spool.export_properties("3RowSteady.json")
# spool.plot()
spool.plot_velocity_triangles()
print('check')