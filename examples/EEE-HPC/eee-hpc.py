from turbodesign import Spool, Inlet, RowType, BladeRow, Passage, Outlet, PassageType
from turbodesign.enums import MassflowConstraint
from turbodesign import Coolant
from turbodesign.loss import FixedPressureLoss
import numpy as np
from cantera import Solution
import pickle
# Geometry Import 

e3_hpc = pickle.load(open('e3_hpc_processed.pkl','rb'))
hub = e3_hpc['hub']
shroud = e3_hpc['shroud']

# TurboDesign Setup - Climb
P0_Ratio = 23
P0 = 59641.8        # Total Pressure [Pa]
T0 = 1587           # Total Temperature [K]
P = P0              # Static Pressure [Pa] @ mach = 0
n_streamlines = 12

# Fluid
fluid = Solution('air.yaml')
fluid.TP = T0, P0 # Use pascal for cantera
print(f"Coefficient of Pressure [J/Kg] {fluid.cp:0.4f}")
#%% Defining the Inlet
inlet = Inlet(hub_location=0, shroud_location=0,beta=[0])
inlet.init_static(M=0, P=[P0], T=[T0])
outlet = Outlet(num_streamlines=n_streamlines)
outlet.init_total(P0=P0*P0_Ratio,percent_radii=[0.5])
#%% Define Blade Rows, processed data is already in mm, hub is already in mm 
cax_arr = [] 
for blade in e3_hpc['blades']: # There should be 21 blades starting with igv and moving into rotor stator pairs
    cax = (blade[0,:,0].max()-blade[0][0,:,0].min())/1000
    cax_arr.append(cax)
    
hub_exit_locations = []; shroud_exit_locations = []
# Get the exit locations
for i in range(1,len(e3_hpc['rotor_stator'])):
    prev_blade = e3_hpc['rotor_stator'][i-1]
    blade = e3_hpc['rotor_stator'][i]
    
    hub_exit = ((prev_blade[0,:,0].max() + blade[0,:,0].min())/2  - hub[:,0].min()) / (hub[:,0].max() - hub[:,0].min())
    hub_exit_locations.append(hub_exit)
    
    shroud_exit = ((prev_blade[-1,:,0].max() + blade[-1,:,0].min())/2  - shroud[:,0].min()) / (shroud[:,0].max() - shroud[:,0].min())
    shroud_exit_locations.append(shroud_exit)

lastblade = e3_hpc['rotor_stator'][-1]
hub_exit_locations.append((lastblade[0,:,0].max()  - hub[:,0].min()) / (hub[:,0].max() - hub[:,0].min()))
shroud_exit_locations.append((lastblade[-1,:,0].max()  - shroud[:,0].min()) / (shroud[:,0].max() - shroud[:,0].min()))
    
beta_exit_flow = [73.6,-67.2,69.5,-63.9]
P0_Loss = [0.057,0.088,0.069,0.014]         # (P01-P02)/(P01-P2)

# Axial location is a percentage along the hub where row exit is defined
IGV1 = BladeRow(row_type=RowType.Stator, hub_location=hub_exit_locations[0],shroud_location=shroud_exit_locations[0],stage_id=1)
rotor1 = BladeRow(row_type=RowType.Rotor, hub_location=hub_exit_locations[1],shroud_location=shroud_exit_locations[1],stage_id=1)
stator1 = BladeRow(row_type=RowType.Stator, hub_location=hub_exit_locations[2],shroud_location=shroud_exit_locations[2],stage_id=1)

rotor2 = BladeRow(row_type=RowType.Rotor, hub_location=hub_exit_locations[3],shroud_location=shroud_exit_locations[3],stage_id=2)
stator2 = BladeRow(row_type=RowType.Stator, hub_location=hub_exit_locations[4],shroud_location=shroud_exit_locations[4],stage_id=2)

rotor3 = BladeRow(row_type=RowType.Rotor, hub_location=hub_exit_locations[5],shroud_location=shroud_exit_locations[5],stage_id=3)
stator3 = BladeRow(row_type=RowType.Stator, hub_location=hub_exit_locations[6],shroud_location=shroud_exit_locations[6],stage_id=3)

rotor4 = BladeRow(row_type=RowType.Rotor, hub_location=hub_exit_locations[7],shroud_location=shroud_exit_locations[7],stage_id=4)
stator4 = BladeRow(row_type=RowType.Stator, hub_location=hub_exit_locations[8],shroud_location=shroud_exit_locations[8],stage_id=4)

rotor5 = BladeRow(row_type=RowType.Rotor, hub_location=hub_exit_locations[9],shroud_location=shroud_exit_locations[9],stage_id=5)
stator5 = BladeRow(row_type=RowType.Stator, hub_location=hub_exit_locations[10],shroud_location=shroud_exit_locations[10],stage_id=5)

rotor6 = BladeRow(row_type=RowType.Rotor, hub_location=hub_exit_locations[11],shroud_location=shroud_exit_locations[11],stage_id=6)
stator6 = BladeRow(row_type=RowType.Stator, hub_location=hub_exit_locations[12],shroud_location=shroud_exit_locations[12],stage_id=6)

rotor7 = BladeRow(row_type=RowType.Rotor, hub_location=hub_exit_locations[13],shroud_location=shroud_exit_locations[13],stage_id=7)
stator7 = BladeRow(row_type=RowType.Stator, hub_location=hub_exit_locations[14],shroud_location=shroud_exit_locations[14],stage_id=7)

rotor8 = BladeRow(row_type=RowType.Rotor, hub_location=hub_exit_locations[15],shroud_location=shroud_exit_locations[15],stage_id=8)
stator8 = BladeRow(row_type=RowType.Stator, hub_location=hub_exit_locations[16],shroud_location=shroud_exit_locations[16],stage_id=8)

rotor9 = BladeRow(row_type=RowType.Rotor, hub_location=hub_exit_locations[17],shroud_location=shroud_exit_locations[17],stage_id=9)
stator9 = BladeRow(row_type=RowType.Stator, hub_location=hub_exit_locations[18],shroud_location=shroud_exit_locations[18],stage_id=9)

rotor10 = BladeRow(row_type=RowType.Rotor, hub_location=hub_exit_locations[19],shroud_location=shroud_exit_locations[19],stage_id=10)
stator10 = BladeRow(row_type=RowType.Stator, hub_location=hub_exit_locations[20],shroud_location=shroud_exit_locations[20],stage_id=10)
 
 # Set an axial chord and number of blades (solidity)
IGV1.axial_chord = cax_arr[0]; IGV1.num_blades = 32
rotor1.axial_chord = cax_arr[1]; rotor1.num_blades = 28 
stator1.axial_chord = cax_arr[2]; stator1.num_blades = 50

rotor2.axial_chord = cax_arr[3]; rotor2.num_blades = 38
stator2.axial_chord = cax_arr[4]; stator2.num_blades = 68

rotor3.axial_chord = cax_arr[5]; rotor3.num_blades = 50
stator3.axial_chord = cax_arr[6]; stator3.num_blades = 83

rotor4.axial_chord = cax_arr[7]; rotor4.num_blades = 60
stator4.axial_chord = cax_arr[8]; stator4.num_blades = 92

rotor5.axial_chord = cax_arr[9]; rotor5.num_blades = 70 
stator5.axial_chord = cax_arr[10]; stator5.num_blades = 110

rotor6.axial_chord = cax_arr[11]; rotor6.num_blades = 80
stator6.axial_chord = cax_arr[12]; stator6.num_blades = 120

rotor7.axial_chord = cax_arr[13]; rotor7.num_blades = 82 
stator7.axial_chord = cax_arr[14]; stator7.num_blades = 112

rotor8.axial_chord = cax_arr[15]; rotor8.num_blades = 84
stator8.axial_chord = cax_arr[16]; stator8.num_blades = 104

rotor9.axial_chord = cax_arr[17]; rotor9.num_blades = 88
stator9.axial_chord = cax_arr[18]; stator9.num_blades = 118

rotor10.axial_chord = cax_arr[19]; rotor10.num_blades = 95
stator10.axial_chord = cax_arr[20]; stator10.num_blades = 140 

gamma = 1.4
Cp = gamma/(gamma-1) * 287.15

rows = [IGV1,rotor1,stator1,rotor2,stator2,rotor3,stator3,rotor4,stator4,
        rotor5,stator5,rotor6,stator6,rotor7,stator7,rotor8,stator8,
        rotor9,stator9,rotor10,stator10]
rows.insert(0,inlet)
rows.append(outlet)

for row in rows:
    row.gamma = gamma
    row.Cp = Cp
    row.coolant = Coolant(T0=293,P0=101325,massflow_percentage=0)
    
# Add in turning angles
IGV1.beta2_metal = [0 for _ in range(n_streamlines)]        # Angle, hub,mean,tip
rotor1.beta2_metal = [65.76, 64.2, 62.77, 61.16, 59.76, 58.45, 57.28, 56.3, 55.64, 55.34, 56.17, 57.08]
stator1.beta2_metal = [57.08, 50.44, 47.32, 45.75, 45.12, 44.83, 44.66, 44.69, 44.86, 45.56, 47.13, 48.32]

rotor2.beta2_metal =  [65.17, 64.43, 63.47, 62.04, 60.63, 59.25, 57.86, 56.47, 55.13, 53.86, 52.48, 51.52]
stator2.beta2_metal =  [58.25, 51.41, 47.93, 46.06, 45.56, 45.23, 44.94, 44.82, 44.91, 46.02, 48.84, 51.02]

rotor3.beta2_metal =  [64.32, 63.66, 62.87, 61.83, 60.70, 59.50, 58.26, 57.03, 55.76, 54.46, 52.81, 51.72]
stator3.beta2_metal = [59.81, 51.18, 47.15, 45.23, 44.84, 44.81, 44.88, 45.02, 45.46, 47.24, 51.8, 55.39]

rotor4.beta2_metal = [63.30, 62.67, 61.96, 60.87, 59.74, 58.65, 57.60, 56.55, 55.51, 54.39, 52.78, 51.93]
stator4.beta2_metal = [60.88, 52.16, 47.94, 45.82, 45.51, 45.55, 45.57, 45.64, 45.88, 47.65, 53.19, 57.95]

rotor5.beta2_metal = [62.72, 62.08, 61.48, 60.56, 59.6, 58.67, 57.76, 56.88, 56.02, 55.05, 53.45, 52.33]
stator5.beta2_metal = [62.45, 53.7, 49.42, 47.09, 46.63, 46.65, 46.69, 46.76, 47.14, 48.9, 54.39, 59.67] 

rotor6.beta2_metal = [61.92, 61.62, 61.25, 60.59, 59.84, 59.08, 58.33, 57.6, 56.88, 56.03, 54.72, 53.84]
stator6.beta2_metal = [58.18, 52.16, 47.93, 45.19, 44.41, 44.33, 44.31, 44.31, 44.75, 46.65, 51.69, 56.37]

rotor7.beta2_metal = [61.99, 61.66, 61.31, 60.75, 60.11, 59.42, 58.74, 58.08, 57.41, 56.59, 55.34, 54.51]
stator7.beta2_metal = [59.45, 53.27, 48.64, 45.23, 43.71, 43.33, 43.31, 43.47, 44.28, 46.60, 51.64, 55.93]

rotor8.beta2_metal = [61.97, 61.81, 61.47, 60.94, 60.34, 59.72, 59.12, 58.59, 58.03, 57.3, 56.27, 55.56]
stator8.beta2_metal = [59.99, 56.43, 51.45, 47.44, 45.24, 44.3, 44.17, 44.55, 45.8, 48.26, 53.05, 57.06]
rotor9.beta2_metal = [63.49, 63.23, 62.91, 62.4, 61.86, 61.31, 60.78, 60.28, 59.77, 59.18, 58.29, 57.72]
stator9.beta2_metal = [62.42, 58.47, 53.20, 48.92, 46.36, 45.18, 45.01, 45.65, 47.20, 49.95, 54.66, 58.24]
rotor10.beta2_metal = [64.50, 64.25, 63.96, 63.48, 62.99, 62.50, 62.03, 61.56, 61.12, 60.57, 59.89, 59.43]
stator10.beta2_metal = [61.99, 58.41, 53.45, 49.33, 46.9, 45.76, 45.66, 46.61, 48.31, 51.07, 55.74, 59.26]

# These are all guessed values
IGV1.loss_model = FixedPressureLoss(0.0509)
rotor1.loss_model = FixedPressureLoss(0.087)
stator2.loss_model = FixedPressureLoss(0.0735)
rotor2.loss_model = FixedPressureLoss(0.07225)
stator3.loss_model = FixedPressureLoss(0.06675)
rotor3.loss_model = FixedPressureLoss(0.06675)
stator3.loss_model = FixedPressureLoss(0.0636)
rotor3.loss_model = FixedPressureLoss(0.0619)

rotor4.loss_model = FixedPressureLoss(0.0559)
stator4.loss_model = FixedPressureLoss(0.8024)

rotor5.loss_model = FixedPressureLoss(0.089475) 
stator5.loss_model = FixedPressureLoss(0.056125)

rotor6.loss_model = FixedPressureLoss(0.0564)
stator6.loss_model = FixedPressureLoss(0.0565)

rotor7.loss_model = FixedPressureLoss(0.0635)
stator7.loss_model = FixedPressureLoss(0.0616)

rotor8.loss_model = FixedPressureLoss(0.06696)
stator8.loss_model = FixedPressureLoss(0.0673)

rotor9.loss_model = FixedPressureLoss(0.0698)
stator9.loss_model = FixedPressureLoss(0.0703)

rotor10.loss_model = FixedPressureLoss(0.0728)
stator10.loss_model = FixedPressureLoss(0.0935)

hub_m = hub/1000
shroud_m = shroud/1000

passage = Passage(hub_m[:,0],hub_m[:,1],
                 shroud_m[:,0],shroud_m[:,1],
                 passageType=PassageType.Axial) # type: ignore

spool = Spool(passage=passage,
            rpm=12400, 
            num_streamlines=n_streamlines, 
            massflow=20, 
            fluid=None,
            rows=rows)
spool.massflow_constraint = MassflowConstraint.BalanceMassFlow # Fixes the exit angle and changes degree of reaction
# spool.plot_geometry()
spool.adjust_streamlines = False
spool.solve() # This also initializes streamlines
spool.export_properties("eee_results.json")
spool.plot()
spool.plot_velocity_triangles()