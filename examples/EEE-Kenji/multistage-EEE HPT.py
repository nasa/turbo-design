'''
    EEE HPT example from
    
    In this example the blade exit angles are fixed and only degree of reaction changes between the stage to match the massflow
'''
#%% Import Library
import sys
from typing import List
sys.path.insert(0,'../../')
from td3 import PassageType
from td3 import TurbineSpool, Inlet, RowType, BladeRow, Passage, Outlet
from td3.enums import MassflowConstraint
from td3.coolant import Coolant
from td3.loss.turbine import FixedPressureLoss
import numpy as np 
from cantera import Solution
import pickle, math
import numpy.typing as npt
import matplotlib.pyplot as plt
from pyturbo.helper import line2D
from scipy.interpolate import interp1d

# Convert to meters 
def convert_to_meters(array:List[npt.NDArray]):
    new_array = []
    for a in array:
        new_array.append(a*0.0254)
    return new_array

def axial_gap(blade1:List[npt.NDArray],blade2:List[npt.NDArray]):
    return min(blade2[1][:,0])-max(blade1[0][:,0])

def axial_chord(blade:List[npt.NDArray]):
    return max(blade[0][:,0])-min(blade[0][:,0])

def exit_flow_angle(blade:List[npt.NDArray]):
    beta = []
    for section in blade:
        n = math.floor(section.shape[0]/2)
        p1 = section[n-9,:]
        p3 = section[n+9,:]
        l = line2D(p1,p3)
        beta.append(np.degrees(np.arctan2(l.dx,-l.dy)))
    return beta

def plot(hub:npt.NDArray,shroud:npt.NDArray,blades:List[npt.NDArray]):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.plot3D(hub[:,0],hub[:,1])
    ax.plot3D(shroud[:,0],shroud[:,1])
    for blade in blades:
        for section in blade:
            ax.plot3D(section[:,0],section[:,2],section[:,1])
    plt.axis('scaled')
    plt.show()
    
def plot_profile(blade):
    n = blade[0].shape[0]
    plt.figure(num=0,clear=True)
    plt.plot(blade[0][:,0],blade[0][:,1],'.')
    plt.plot(blade[0][math.floor(n/2),0],blade[0][math.floor(n/2),1],'or')
    plt.axis('scaled')
    plt.show()
    
hub_shroud = pickle.load(open('hub_shroud.pkl','rb'))       # Units are in in inches
stator_rotor = pickle.load(open('stator_rotor.pkl','rb'))   

hub = hub_shroud['Hub']*0.0254; shroud = hub_shroud['Shroud']*0.0254
hub2 = shroud*0
hub2[:,0] = shroud[:,0]
hub2[:,1] = interp1d(hub[:,0],hub[:,1])(shroud[:,0])

stator1 = convert_to_meters(stator_rotor['Stator1'])
stator2 = convert_to_meters(stator_rotor['Stator2'])
rotor1 = convert_to_meters(stator_rotor['Rotor1'])
rotor2 = convert_to_meters(stator_rotor['Rotor2'])
gap1 = axial_gap(stator1,rotor1)
gap2 = axial_gap(rotor1,stator2)
gap3 = axial_gap(stator2,rotor2)

stator1_cax = axial_chord(stator1)
rotor1_cax = axial_chord(rotor1)
stator2_cax = axial_chord(stator2)
rotor2_cax = axial_chord(rotor2)

# plot_profile(stator1)

# plot(hub,shroud,[stator1,rotor1,stator2,rotor2])

#%% Define the Passage 
# Geometry from OptTurb
passage = Passage(hub2[:,0],hub2[:,1],
                  shroud[:,0],shroud[:,1],
                  passageType=PassageType.Axial)

# #%% Design Conditions 
Design_RPM = 13119
massflow = 69.1 # kg/s
P0 = 3163392     # Pascal 
T0 = 1810      # Kelvin

# Fluid
fluid = Solution('air.yaml')
fluid.TP = T0, P0 # Use pascal for cantera
print(f"Coefficient of Pressure [J/Kg] {fluid.cp:0.4f}")

#%% Defining the Inlet
inlet = Inlet(M=0.1, 
                 P0=[P0],
                 T0=[T0], 
                 beta=[0], 
                 fluid=fluid, 
                 percent_radii=0.5,
                 axial_location=0)
outlet = Outlet(P=829869.3,percent_radii=0.5,num_streamlines=3)

hub_len = max(hub[:,0])-min(hub[:,0])
#%% Define Blade Rows 
# Axial location is a percentage along the hub where row exit is defined
stator1 = BladeRow(row_type=RowType.Stator,axial_location=(min(stator1[0][:,0]) - min(hub[:,0]))/hub_len, stage_id=0)
rotor1 = BladeRow(row_type=RowType.Rotor, axial_location=(min(rotor1[0][:,0]) - min(hub[:,0]))/hub_len,stage_id=0)
stator2 = BladeRow(row_type=RowType.Stator,axial_location=(min(stator2[0][:,0]) - min(hub[:,0]))/hub_len,stage_id=1)
rotor2 = BladeRow(row_type=RowType.Rotor, axial_location=(min(rotor2[0][:,0]) - min(hub[:,0]))/hub_len,stage_id=1)

stator1.axial_chord = stator1_cax # Set an axial chord
rotor1.axial_chord = rotor1_cax
stator2.axial_chord = stator2_cax 
rotor2.axial_chord = rotor2_cax

# Coolant Definition: Use Kelvin and Pascal
stator1.coolant = Coolant(fluid, T0=T0*0.555556, P0= P0,massflow_percentage=0) 
rotor1.coolant = Coolant(fluid, T0*0.555556, P0=P0,massflow_percentage=0)
stator2.coolant = Coolant(fluid, T0=T0*0.555556, P0=P0,massflow_percentage=0) 
rotor2.coolant = Coolant(fluid, T0*0.555556, P0=P0,massflow_percentage=0)

# Add in turning angles
stator1.beta2_metal = [72.2,72.2,72.2]        # Angle, hub,mean,tip
rotor1.beta2_metal = [-63.4,-63.4,-63.4] 
stator2.beta2_metal = [68.5,68.5,68.5]
rotor2.beta2_metal = [-63.5,-63.5,-63.5]

# These are all guessed values 
stator1.loss_model = FixedPressureLoss(0.06)
rotor1.loss_model = FixedPressureLoss(0.06)
stator2.loss_model = FixedPressureLoss(0.08)
rotor2.loss_model = FixedPressureLoss(0.09)

# inlet_to_outlet_pratio tells the optimize to guess within this range
stator1.inlet_to_outlet_pratio = [0.07, 0.15]
rotor1.inlet_to_outlet_pratio = [0.15, 0.7]
stator2.inlet_to_outlet_pratio = [0.10, 0.9]
rotor2.inlet_to_outlet_pratio = [0.15, 0.5]

#%% Initialize the Spool
spool = TurbineSpool(passage=passage,
            rpm=Design_RPM, 
            num_streamlines=3, 
            massflow=massflow, 
            rows=[inlet,stator1,rotor1,stator2,rotor2,outlet])
spool.fluid = fluid
spool.massflow_constraint = MassflowConstraint.BalanceMassFlow # Fixes the exit angle and changes degree of reaction
# spool.plot_geometry()
spool.solve() # This also initializes streamlines
spool.export_properties("EEE-HPT.json")
spool.plot()
spool.plot_velocity_triangles()
print('check')