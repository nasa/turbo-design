"""Background

EEE report (the xlsx) was extracted from a publication in 1982 and it was scanned who knows when. The document looks old and hard to read. A lot of this work looks to be good but without the code used to generate the E3 spreadsheet or CFD results from the 1980s if they even have that, all this falls under the category of 'trust me bro'

Todo: 
1. Run a simulation of EEE and do the calculations then modify this file or create a new one with the latest data. 
2. Match the 1D with CFD to validate the math. 

"""

from turbodesign import TurbineSpool, Inlet, RowType, BladeRow, Passage, Outlet, PassageType
from turbodesign.compressor_spool import CompressorSpool
from turbodesign.enums import MassflowConstraint
from turbodesign import Coolant
from turbodesign.loss.fixedpressureloss import FixedPressureLoss
from turbodesign.deviation.fixed_deviation import FixedDeviation
import numpy as np
from cantera import Solution
from pathlib import Path
import pandas as pd
import pickle
from read_overall_data import load_overall_data, load_blade_counts, compute_entropy_rise
# Geometry Import 

e3_hpc = pickle.load(open(Path(__file__).resolve().parent / 'e3_hpc_processed.pkl','rb'))
hub = e3_hpc['hub']
shroud = e3_hpc['shroud']

excel_data = load_overall_data(
    Path(__file__).resolve().parent / 'E3_HPC_Overall_Data.xlsx',
    sheet_name=None,
    loss_sheet_name="Detailed Report Data",
    convert_units=True,
)
blade_counts = load_blade_counts(Path(__file__).resolve().parent / 'E3_HPC_Overall_Data.xlsx')


P0 = excel_data['inlet']["Inlet Pt"].mean()
T0 = excel_data['inlet']["TT Exit"].mean()
n_streamlines = 12
P0_Ratio = 1.0  # placeholder until defined from data

# Fluid
fluid = Solution('air.yaml')
fluid.TP = T0, P0 # Use pascal for cantera
entropy_calcs = compute_entropy_rise(Path(__file__).resolve().parent / 'E3_HPC_Overall_Data.xlsx',fluid=fluid)

print(f"Coefficient of Pressure [J/Kg] {fluid.cp:0.4f}")
#%% Defining the Inlet
inlet = Inlet(hub_location=0, shroud_location=0,beta=[0])
inlet.init_total(P0=excel_data['inlet']["Inlet Pt"].to_numpy().tolist(), T0=excel_data['inlet']["TT Exit"].to_numpy().tolist(),M=excel_data['inlet']["Ma Inlet"].to_numpy().tolist())
outlet = Outlet(num_streamlines=n_streamlines)
outlet.init_total(P0=excel_data['stator10']['PT Exit'].mean(),percent_radii=[0.5])
#%% Define Blade Rows, processed data is already in mm, hub is already in mm 
cax_arr = [] 
# There should be 21 blades starting with igv and moving into rotor stator pairs
for blade in e3_hpc['blades']:
    ss = blade[0]; ps = blade[1]    
    cax1 = (ss[0,:,0].max()-ss[0,:,0].min())/1000
    cax2 = (ps[0,:,0].max()-ps[0,:,0].min())/1000
    cax_arr.append(max([cax1,cax2]))
    
hub_exit_locations = []; shroud_exit_locations = []

# Get the exit locations as a percentage along the hub and shroud 
for i in range(1,len(e3_hpc['blades'])):
    prev_blade = e3_hpc['blades'][i-1]
    blade = e3_hpc['blades'][i]

    hub_exit = ((prev_blade[0][0,:,0].max() + blade[0][0,:,0].min())/2  - hub[:,0].min()) / (hub[:,0].max() - hub[:,0].min())
    hub_exit_locations.append(hub_exit)
    
    shroud_exit = ((prev_blade[0][-1,:,0].max() + blade[0][-1,:,0].min())/2  - shroud[:,0].min()) / (shroud[:,0].max() - shroud[:,0].min())
    shroud_exit_locations.append(shroud_exit)

lastblade = e3_hpc['blades'][-1]
hub_exit_locations.append((lastblade[0][0,:,0].max()  - hub[:,0].min()) / (hub[:,0].max() - hub[:,0].min()))
shroud_exit_locations.append((lastblade[0][-1,:,0].max()  - shroud[:,0].min()) / (shroud[:,0].max() - shroud[:,0].min()))

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

# Metal inlet/exit angles pulled from Excel Beta columns
def set_beta_metal(row: BladeRow, key: str):
    df = excel_data.get(key)
    if df is None:
        return
    beta_exit_col = None
    beta_inlet_col = None
    if "Beta.1" in df.columns:
        beta_exit_col = "Beta.1"
        if "Beta" in df.columns:
            beta_inlet_col = "Beta"
    elif "Beta" in df.columns:
        beta_exit_col = "Beta"

    def _prepare_beta(values: pd.Series) -> list[float]:
        beta_vals = pd.to_numeric(values, errors="coerce").dropna().tolist()
        trimmed = beta_vals[:n_streamlines]
        if row.row_type == RowType.Rotor:
            return [-abs(val) for val in trimmed]
        return trimmed

    if beta_exit_col:
        exit_vals = _prepare_beta(df[beta_exit_col])
        if exit_vals:
            row.beta2_metal = exit_vals
    if beta_inlet_col:
        inlet_vals = _prepare_beta(df[beta_inlet_col])
        if inlet_vals:
            row.beta1_metal = inlet_vals

set_beta_metal(IGV1, "inlet")
set_beta_metal(rotor1, "rotor1")
set_beta_metal(stator1, "stator1")
set_beta_metal(rotor2, "rotor2")
set_beta_metal(stator2, "stator2")
set_beta_metal(rotor3, "rotor3")
set_beta_metal(stator3, "stator3")
set_beta_metal(rotor4, "rotor4")
set_beta_metal(stator4, "stator4")
set_beta_metal(rotor5, "rotor5")
set_beta_metal(stator5, "stator5")
set_beta_metal(rotor6, "rotor6")
set_beta_metal(stator6, "stator6")
set_beta_metal(rotor7, "rotor7")
set_beta_metal(stator7, "stator7")
set_beta_metal(rotor8, "rotor8")
set_beta_metal(stator8, "stator8")
set_beta_metal(rotor9, "rotor9")
set_beta_metal(stator9, "stator9")
set_beta_metal(rotor10, "rotor10")
set_beta_metal(stator10, "stator10")


# Assign loss models from Excel Loss column where available
def set_loss_model(row: BladeRow, key: str):
    df = excel_data.get(key)
    if df is None or "Loss" not in df.columns:
        return
    loss_vals = pd.to_numeric(df["Loss"], errors="coerce").to_numpy()
    if loss_vals.size == 0:
        return
    row.loss_model = FixedPressureLoss(loss_vals)

set_loss_model(IGV1, "inlet")
set_loss_model(rotor1, "rotor1")
set_loss_model(stator1, "stator1")
set_loss_model(rotor2, "rotor2")
set_loss_model(stator2, "stator2")
set_loss_model(rotor3, "rotor3")
set_loss_model(stator3, "stator3")
set_loss_model(rotor4, "rotor4")
set_loss_model(stator4, "stator4")
set_loss_model(rotor5, "rotor5")
set_loss_model(stator5, "stator5")
set_loss_model(rotor6, "rotor6")
set_loss_model(stator6, "stator6")
set_loss_model(rotor7, "rotor7")
set_loss_model(stator7, "stator7")
set_loss_model(rotor8, "rotor8")
set_loss_model(stator8, "stator8")
set_loss_model(rotor9, "rotor9")
set_loss_model(stator9, "stator9")
set_loss_model(rotor10, "rotor10")
set_loss_model(stator10, "stator10")

# Assign deviation models from Excel Deviation column where available
def set_deviation_model(row: BladeRow, key: str):
    df = excel_data.get(key)
    if df is None or "Deviation" not in df.columns:
        return
    dev_vals = pd.to_numeric(df["Deviation"], errors="coerce").to_numpy()
    if dev_vals.size == 0:
        return
    row.deviation_function = FixedDeviation(dev_vals)

set_deviation_model(IGV1, "inlet")      # This sets the deviation value according to the spreadsheet
set_deviation_model(rotor1, "rotor1")
set_deviation_model(stator1, "stator1")
set_deviation_model(rotor2, "rotor2")
set_deviation_model(stator2, "stator2")
set_deviation_model(rotor3, "rotor3")
set_deviation_model(stator3, "stator3")
set_deviation_model(rotor4, "rotor4")
set_deviation_model(stator4, "stator4")
set_deviation_model(rotor5, "rotor5")
set_deviation_model(stator5, "stator5")
set_deviation_model(rotor6, "rotor6")
set_deviation_model(stator6, "stator6")
set_deviation_model(rotor7, "rotor7")
set_deviation_model(stator7, "stator7")
set_deviation_model(rotor8, "rotor8")
set_deviation_model(stator8, "stator8")
set_deviation_model(rotor9, "rotor9")
set_deviation_model(stator9, "stator9")
set_deviation_model(rotor10, "rotor10")
set_deviation_model(stator10, "stator10")

gamma = 1.4
Cp = gamma/(gamma-1) * 287.15

rows = [IGV1,rotor1,stator1,rotor2,stator2,rotor3,stator3,rotor4,stator4,
        rotor5,stator5,rotor6,stator6,rotor7,stator7,rotor8,stator8,
        rotor9,stator9,rotor10,stator10]

hub_m = hub/1000
shroud_m = shroud/1000

passage = Passage(hub_m[:,0],hub_m[:,1],
                 shroud_m[:,0],shroud_m[:,1],
                 passageType=PassageType.Axial) # type: ignore

spool = CompressorSpool(
            passage=passage,
            massflow=54,
            inlet=inlet,
            outlet=outlet,
            rows=rows,
            rpm=12400,
            num_streamlines=n_streamlines,
            fluid=fluid)


spool.massflow_constraint = MassflowConstraint.PressureBalance # Fixes the exit angle and changes degree of reaction
# spool.plot_geometry()
spool.adjust_streamlines = False
spool.solve() # This also initializes streamlines
spool.export_properties("E3-HPC-Results.json")
spool.plot()
spool.plot_velocity_triangles()
