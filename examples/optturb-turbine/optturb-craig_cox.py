'''
    GEE3HP Turbine
    2 stage cooled turbine

'''
#%% Import Library
import sys
sys.path.insert(0,'../../')
from td3 import PassageType
from td3 import TurbineSpool, Inlet, RowType, BladeRow, Passage, Outlet
from td3.enums import MassflowConstraint
from td3.coolant import Coolant
from td3.loss.turbine import CraigCox
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
xhub = np.array([-cax, 0.0, cax, 2*cax])
xshroud = np.array([-cax, 0.0, cax, 2*cax])
axial_len = xhub[-1]-xhub[0]

passage = Passage(xhub,rhub,
                 xshroud,rshroud,
                 passageType=PassageType.Axial)

#%% Design Conditions 
Design_RPM = 7500
P0 = 500000     # Pascal 
T0 = 676.3      # Kelvin

# Fluid 
fluid = Solution('air.yaml')
fluid.TP = T0, P0 # Use pascal for cantera
print(f"Coefficient of Pressure [J/Kg] {fluid.cp:0.4f}")

#%% Defining the Inlet
inlet = Inlet(M=0.2, 
                 P0=[P0], 
                 T0=[T0], 
                 beta=[0], 
                 fluid=fluid, 
                 percent_radii=0.5, 
                 axial_location=0)
outlet = Outlet(P=P0/3.96,percent_radii=0.5,num_streamlines=3)

stator1 = BladeRow(row_type=RowType.Stator,axial_location=2*cax/axial_len)
rotor1 = BladeRow(row_type=RowType.Rotor,axial_location=3*cax/axial_len)

stator1.axial_chord = cax # Set an axial chord. Used for plotting
rotor1.axial_chord = cax

# Coolant Definition: Use Kelvin and Pascal
stator1.coolant = Coolant(fluid, T0=616*0.555556, P0= 50.6 * 6894.76, massflow_percentage=0) 
rotor1.coolant = Coolant(fluid, 622*0.555556, 50.3 * 6894.76,massflow_percentage=0)

# Add in turning angles
stator1.beta2_metal = [73,73,73] # Angle, hub,mean,tip
stator1.loss_model = CraigCox()
rotor1.loss_model = CraigCox()
rotor1.beta2_metal = [-67.6,-67.6,-67.6] # Angle, hub,mean,tip

# Coolant: Use Kelvin and Pascal
massflow = 35.9 # kg/s, guessed value

spool = TurbineSpool(passage=passage,
            rpm=Design_RPM, 
            num_streamlines=3, 
            massflow=massflow, 
            rows=[inlet,stator1,rotor1,outlet])
spool.fluid = fluid
spool.massflow_constraint = MassflowConstraint.BalanceMassFlow # Fixes the exit angle and changes degree of reaction

# spool.plot_geometry()
spool.solve() # This also initializes streamlines
spool.export_properties("optturb.json")
spool.plot()
spool.plot_velocity_triangles()
print('check')