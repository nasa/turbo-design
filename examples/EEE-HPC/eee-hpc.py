from turbodesign import TurbineSpool, Inlet, RowType, BladeRow, Passage, Outlet, PassageType
from turbodesign.enums import MassflowConstraint
from turbodesign import Coolant
from turbodesign.loss.fixedpressureloss import FixedPressureLoss
import numpy as np
from cantera import Solution
from pathlib import Path
import pandas as pd
import pickle
from read_overall_data import load_overall_data, load_blade_counts
# Geometry Import 

e3_hpc = pickle.load(open(Path(__file__).resolve().parent / 'e3_hpc_processed.pkl','rb'))
hub = e3_hpc['hub']
shroud = e3_hpc['shroud']

excel_data = load_overall_data(Path(__file__).resolve().parent / 'E3_HPC_Overall_Data.xlsx')
blade_counts = load_blade_counts(Path(__file__).resolve().parent / 'E3_HPC_Overall_Data.xlsx')

P0 = excel_data['inlet']["Inlet Pt"].mean()
T0 = excel_data['inlet']["TT Exit"].mean()
n_streamlines = 12
P0_Ratio = 1.0  # placeholder until defined from data

# Fluid
fluid = Solution('air.yaml')
fluid.TP = T0, P0 # Use pascal for cantera

print(f"Coefficient of Pressure [J/Kg] {fluid.cp:0.4f}")
#%% Defining the Inlet
inlet = Inlet(hub_location=0, shroud_location=0,beta=[0])
inlet.init_total(P0=excel_data['inlet']["Inlet Pt"].to_numpy().tolist(), T0=excel_data['inlet']["TT Exit"].to_numpy().tolist(),M=excel_data['inlet']["Ma Inlet"].to_numpy().tolist())
outlet = Outlet(num_streamlines=n_streamlines)
outlet.init_total(P0=excel_data['stator10']['PT Exit'].mean(),percent_radii=[0.5])
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

# Axial location is a percentage along the hub where row exit is defined
IGV1 = BladeRow(row_type=RowType.Stator, hub_location=hub_exit_locations[0],shroud_location=shroud_exit_locations[0],stage_id=1)
IGV1.num_blades = blade_counts.get("igv", IGV1.num_blades)

rotor1 = BladeRow(row_type=RowType.Rotor, hub_location=hub_exit_locations[1],shroud_location=shroud_exit_locations[1],stage_id=1)
rotor1.num_blades = blade_counts.get("rotor1", rotor1.num_blades)

stator1 = BladeRow(row_type=RowType.Stator, hub_location=hub_exit_locations[2],shroud_location=shroud_exit_locations[2],stage_id=1)
stator1.num_blades = blade_counts.get("stator1", stator1.num_blades)

rotor2 = BladeRow(row_type=RowType.Rotor, hub_location=hub_exit_locations[3],shroud_location=shroud_exit_locations[3],stage_id=2)
rotor2.num_blades = blade_counts.get("rotor2", rotor2.num_blades)

stator2 = BladeRow(row_type=RowType.Stator, hub_location=hub_exit_locations[4],shroud_location=shroud_exit_locations[4],stage_id=2)
stator2.num_blades = blade_counts.get("stator2", stator2.num_blades)

rotor3 = BladeRow(row_type=RowType.Rotor, hub_location=hub_exit_locations[5],shroud_location=shroud_exit_locations[5],stage_id=3)
rotor3.num_blades = blade_counts.get("rotor3", rotor3.num_blades)

stator3 = BladeRow(row_type=RowType.Stator, hub_location=hub_exit_locations[6],shroud_location=shroud_exit_locations[6],stage_id=3)
stator3.num_blades = blade_counts.get("stator3", stator3.num_blades)

rotor4 = BladeRow(row_type=RowType.Rotor, hub_location=hub_exit_locations[7],shroud_location=shroud_exit_locations[7],stage_id=4)
rotor4.num_blades = blade_counts.get("rotor4", rotor4.num_blades)
stator4 = BladeRow(row_type=RowType.Stator, hub_location=hub_exit_locations[8],shroud_location=shroud_exit_locations[8],stage_id=4)
stator4.num_blades = blade_counts.get("stator4", stator4.num_blades)

rotor5 = BladeRow(row_type=RowType.Rotor, hub_location=hub_exit_locations[9],shroud_location=shroud_exit_locations[9],stage_id=5)
rotor5.num_blades = blade_counts.get("rotor5", rotor5.num_blades)
stator5 = BladeRow(row_type=RowType.Stator, hub_location=hub_exit_locations[10],shroud_location=shroud_exit_locations[10],stage_id=5)
stator5.num_blades = blade_counts.get("stator5", stator5.num_blades)

rotor6 = BladeRow(row_type=RowType.Rotor, hub_location=hub_exit_locations[11],shroud_location=shroud_exit_locations[11],stage_id=6)
rotor6.num_blades = blade_counts.get("rotor6", rotor6.num_blades)
stator6 = BladeRow(row_type=RowType.Stator, hub_location=hub_exit_locations[12],shroud_location=shroud_exit_locations[12],stage_id=6)
stator6.num_blades = blade_counts.get("stator6", stator6.num_blades)

rotor7 = BladeRow(row_type=RowType.Rotor, hub_location=hub_exit_locations[13],shroud_location=shroud_exit_locations[13],stage_id=7)
rotor7.num_blades = blade_counts.get("rotor7", rotor7.num_blades)
stator7 = BladeRow(row_type=RowType.Stator, hub_location=hub_exit_locations[14],shroud_location=shroud_exit_locations[14],stage_id=7)
stator7.num_blades = blade_counts.get("stator7", stator7.num_blades)

rotor8 = BladeRow(row_type=RowType.Rotor, hub_location=hub_exit_locations[15],shroud_location=shroud_exit_locations[15],stage_id=8)
rotor8.num_blades = blade_counts.get("rotor8", rotor8.num_blades)
stator8 = BladeRow(row_type=RowType.Stator, hub_location=hub_exit_locations[16],shroud_location=shroud_exit_locations[16],stage_id=8)
stator8.num_blades = blade_counts.get("stator8", stator8.num_blades)

rotor9 = BladeRow(row_type=RowType.Rotor, hub_location=hub_exit_locations[17],shroud_location=shroud_exit_locations[17],stage_id=9)
rotor9.num_blades = blade_counts.get("rotor9", rotor9.num_blades)
stator9 = BladeRow(row_type=RowType.Stator, hub_location=hub_exit_locations[18],shroud_location=shroud_exit_locations[18],stage_id=9)
stator9.num_blades = blade_counts.get("stator9", stator9.num_blades)

rotor10 = BladeRow(row_type=RowType.Rotor, hub_location=hub_exit_locations[19],shroud_location=shroud_exit_locations[19],stage_id=10)
rotor10.num_blades = blade_counts.get("rotor10", rotor10.num_blades)
stator10 = BladeRow(row_type=RowType.Stator, hub_location=hub_exit_locations[20],shroud_location=shroud_exit_locations[20],stage_id=10)
stator10.num_blades = blade_counts.get("stator10", stator10.num_blades)

# Set axial chord from geometry
IGV1.axial_chord = cax_arr[0]
rotor1.axial_chord = cax_arr[1]
stator1.axial_chord = cax_arr[2]

rotor2.axial_chord = cax_arr[3]
stator2.axial_chord = cax_arr[4]

rotor3.axial_chord = cax_arr[5]
stator3.axial_chord = cax_arr[6]

rotor4.axial_chord = cax_arr[7]
stator4.axial_chord = cax_arr[8]

rotor5.axial_chord = cax_arr[9]
stator5.axial_chord = cax_arr[10]

rotor6.axial_chord = cax_arr[11]
stator6.axial_chord = cax_arr[12]

rotor7.axial_chord = cax_arr[13]
stator7.axial_chord = cax_arr[14]

rotor8.axial_chord = cax_arr[15]
stator8.axial_chord = cax_arr[16]

rotor9.axial_chord = cax_arr[17]
stator9.axial_chord = cax_arr[18]

rotor10.axial_chord = cax_arr[19]
stator10.axial_chord = cax_arr[20]

# Metal exit angles pulled from Excel Beta column
def _set_beta2_metal(row: BladeRow, key: str):
    df = excel_data.get(key)
    if df is None or "Beta" not in df.columns:
        return
    beta_vals = pd.to_numeric(df["Beta"], errors="coerce").tolist()
    if len(beta_vals) >= n_streamlines:
        row.beta2_metal = beta_vals[:n_streamlines]

_set_beta2_metal(IGV1, "inlet")
_set_beta2_metal(rotor1, "rotor1")
_set_beta2_metal(stator1, "stator1")
_set_beta2_metal(rotor2, "rotor2")
_set_beta2_metal(stator2, "stator2")
_set_beta2_metal(rotor3, "rotor3")
_set_beta2_metal(stator3, "stator3")
_set_beta2_metal(rotor4, "rotor4")
_set_beta2_metal(stator4, "stator4")
_set_beta2_metal(rotor5, "rotor5")
_set_beta2_metal(stator5, "stator5")
_set_beta2_metal(rotor6, "rotor6")
_set_beta2_metal(stator6, "stator6")
_set_beta2_metal(rotor7, "rotor7")
_set_beta2_metal(stator7, "stator7")
_set_beta2_metal(rotor8, "rotor8")
_set_beta2_metal(stator8, "stator8")
_set_beta2_metal(rotor9, "rotor9")
_set_beta2_metal(stator9, "stator9")
_set_beta2_metal(rotor10, "rotor10")
_set_beta2_metal(stator10, "stator10")

gamma = 1.4
Cp = gamma/(gamma-1) * 287.15

rows = [IGV1,rotor1,stator1,rotor2,stator2,rotor3,stator3,rotor4,stator4,
        rotor5,stator5,rotor6,stator6,rotor7,stator7,rotor8,stator8,
        rotor9,stator9,rotor10,stator10]

for row in [inlet, *rows, outlet]:
    row.gamma = gamma
    row.Cp = Cp
    row.coolant = Coolant(T0=293,P0=101325,massflow_percentage=0)

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

spool = TurbineSpool(
            passage=passage,
            massflow=20,
            inlet=inlet,
            outlet=outlet,
            rows=rows,
            rpm=12400,
            num_streamlines=n_streamlines,
            fluid=None)
spool.massflow_constraint = MassflowConstraint.PressureBalance # Fixes the exit angle and changes degree of reaction
# spool.plot_geometry()
spool.adjust_streamlines = False
spool.solve() # This also initializes streamlines
spool.export_properties("eee_results.json")
spool.plot()
spool.plot_velocity_triangles()
