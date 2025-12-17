'''
    1D meanline example from 
    Turbine Passage Design Methodology to Minimize Entropy Production—A Two-Step Optimization Strategy. 
    https://doi.org/10.3390/e21060604

    In this example the blade exit angles are allowed to change 
'''
#%% Import Library
from turbodesign import PassageType
from turbodesign.row_factory import make_rotor_row, make_stator_row
from turbodesign import TurbineSpool, Inlet, RowType, BladeRow, Passage, Outlet
from turbodesign.enums import MassflowConstraint
from turbodesign.coolant import Coolant
from turbodesign.loss import FixedPressureLoss
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
# Design Conditions 
Design_RPM = 7500
power = 3.64E6  # Watts
massflow = 35.9 # kg/s; Guessed value to start the calculation
P0 = 500000     # Pascal 
T0 = 676.3      # Kelvin

# Fluid
fluid = Solution('air.yaml')
fluid.TP = T0, P0 # Use pascal for cantera
print(f"Coefficient of Pressure [J/Kg] {fluid.cp:0.4f}")

#%% Defining the Inlet/Outlet
inlet = Inlet(hub_location=0, alpha=[0])
inlet.init_total(
    P0=[P0],
    T0=[T0],
    M=[0.2],
    percent_radii=[0.5],
)
outlet = Outlet(num_streamlines=3)
outlet.init_static(P=P0 / 3.96, percent_radii=[0.5])

#%% Define Blade Rows 
# Axial location is a percentage along the hub where row exit is defined
stator1 = make_stator_row(hub_location=2 * cax / axial_len)
rotor1 = make_rotor_row(hub_location=3 * cax / axial_len)
stator1.axial_chord = cax
rotor1.axial_chord = cax
rotor1.rp = 0.3924 # Degree of Reaction guessed value 
# Coolant: Use Kelvin and Pascal
stator1.coolant = Coolant(T0=616*0.555556, P0=50.6*6894.76, massflow_percentage=0)
rotor1.coolant = Coolant(T0=622*0.555556, P0=50.3*6894.76, massflow_percentage=0)

# Add in turning angles
stator1.beta2_metal = [72,72,73] # Angle, hub,mean,tip
rotor1.beta2_metal = [-67.6,-67.6,-67.6] # Angle, hub,mean,tip
# Loss Values for Stator and Rotor 
stator1.loss_model = FixedPressureLoss(0.118)
rotor1.loss_model = FixedPressureLoss(0.214)
#%% Initialize the TurbineSpool
spool = TurbineSpool(
    passage=passage,
    massflow=massflow,
    inlet=inlet,
    outlet=outlet,
    rows=[stator1, rotor1],
    rpm=Design_RPM,
    num_streamlines=3,
)
spool.fluid = fluid
spool.massflow_constraint = MassflowConstraint.AngleMatch # changes the exit angle
# spool.plot_geometry()
spool.solve() # This also initializes streamlines
spool.export_properties("optturb.json")
spool.plot()
spool.plot_velocity_triangles()
print('check')
