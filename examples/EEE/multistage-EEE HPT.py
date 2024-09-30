'''
    EEE HPT example from
    
    In this example the blade exit angles are fixed and only degree of reaction changes between the stage to match the massflow
'''
#%% Import Library
from typing import List
from turbodesign import PassageType
from turbodesign import TurbineSpool, Inlet, RowType, BladeRow, Passage, Outlet
from turbodesign.enums import MassflowConstraint
from turbodesign.coolant import Coolant
from turbodesign.loss.turbine import FixedPressureLoss
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
    
hub_shroud = pickle.load(open('examples/EEE/hub_shroud.pkl','rb'))       # Units are in in inches
stator_rotor = pickle.load(open('examples/EEE/stator_rotor.pkl','rb'))   

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

# #%% Design Conditions from Table I "EEE_HPT_report.pdf"
Design_RPM = 12680
massflow = 69.1 # kg/s
P0 = 60E5    # Pascal 
T0 = 1588      # Kelvin 
Me = 0.421
PT4_PT42 = 4.933
PT42 = P0/PT4_PT42
Pexit = PT42/(1+(1.4-1)/2*Me**2)**(1.4/0.4)
# Fluid
fluid = Solution('air.yaml')
fluid.TP = T0, P0 # Use pascal for cantera
print(f"Coefficient of Pressure [J/Kg] {fluid.cp:0.4f}")

#%% Defining the Inlet
hub_len = max(hub[:,0])-min(hub[:,0])
inlet = Inlet(M=0.1, 
                 P0=[P0],
                 T0=[T0], 
                 beta=[0], 
                 fluid=fluid, 
                 percent_radii=0.5,
                 axial_location=(min(stator1[0][:,0]) - min(hub[:,0]))/hub_len)
outlet = Outlet(P=Pexit,percent_radii=0.5,num_streamlines=3)

#%% Define Blade Rows 
# Axial location is a percentage along the hub where row exit is defined
stator1 = BladeRow(row_type=RowType.Stator,axial_location=(max(stator1[0][:,0]) - min(hub[:,0]))/hub_len, stage_id=0)
rotor1 = BladeRow(row_type=RowType.Rotor, axial_location=(max(rotor1[0][:,0]) - min(hub[:,0]))/hub_len,stage_id=0)
stator2 = BladeRow(row_type=RowType.Stator,axial_location=(max(stator2[0][:,0]) - min(hub[:,0]))/hub_len,stage_id=1)
rotor2 = BladeRow(row_type=RowType.Rotor, axial_location=(max(rotor2[0][:,0]) - min(hub[:,0]))/hub_len,stage_id=1)

stator1.axial_chord = stator1_cax # Set an axial chord
rotor1.axial_chord = rotor1_cax
stator2.axial_chord = stator2_cax 
rotor2.axial_chord = rotor2_cax

# Coolant Definition: Use Kelvin and Pascal
stator1.coolant = Coolant(fluid, T0=T0*0.555556, P0= P0,massflow_percentage=0) 
rotor1.coolant = Coolant(fluid, T0*0.555556, P0=P0,massflow_percentage=0)
stator2.coolant = Coolant(fluid, T0=T0*0.555556, P0=P0,massflow_percentage=0) 
rotor2.coolant = Coolant(fluid, T0*0.555556, P0=P0,massflow_percentage=0)

# Add in turning angles comes from "HPT_EEE_report.pdf" table III
stator1.beta2_metal = [75.4, 74.2, 73.1]        # Angle, hub,mean,tip
rotor1.beta2_metal = [-64.4,-66.9,-65.6] 
stator2.beta2_metal = [69,69,69]
rotor2.beta2_metal = [-59.8,-59.8,-59.9]

# These are all guessed values 
stator1.loss_model = FixedPressureLoss(0.03)
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
spool.export_properties("examples/EEE/EEE-HPT.json")
spool.plot()
spool.plot_velocity_triangles()
print('check')