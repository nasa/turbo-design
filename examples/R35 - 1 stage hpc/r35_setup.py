"""Background

Rotor 35 1 stage HPC 

"""

from turbodesign import TurbineSpool, Inlet, RowType, BladeRow, Passage, Outlet, PassageType
from turbodesign.row_factory import make_blade_row
from turbodesign.enums import MassflowConstraint
from turbodesign import read_agf
from turbodesign.loss.fixedpressureloss import FixedPressureLoss
from turbodesign.deviation.fixed_deviation import FixedDeviation
import numpy as np
from cantera import Solution
from pathlib import Path
import pandas as pd
import pickle
# Geometry Import 

blade_counts = [36, 46]

rotor = read_agf(str(Path(__file__).resolve().parent / 'R35.agf'))
stator = read_agf(str(Path(__file__).resolve().parent / 'S35.agf'))
P0 = rotor['inlet'].ptin * 6894.76 # Pa
T0 = rotor['inlet'].ttin/1.8 # K
M = rotor['inlet'].machin
P03 = 2.813195e+01 * 6894.76 # Pa 
gamma = 1.4
Cp = gamma/(gamma-1) * 287.15

n_streamlines = 5

# Fluid
fluid = Solution('air.yaml')
fluid.TP = T0, P0 # Use pascal for cantera

print(f"Coefficient of Pressure [J/Kg] {fluid.cp:0.4f}")
#%% Defining the Inlet
inlet = Inlet(hub_location=0, shroud_location=0,beta=[0])
inlet.init_total(P0=P0, T0=T0,M=M)
outlet = Outlet(num_streamlines=n_streamlines)
outlet.init_total(P0=P03,percent_radii=[0.5])
#%% Define Blade Rows, processed data is already in mm, hub is already in mm 
# Note: You dont really need this to run the code but it could be helpful to some calculations. You do need to know where the blades are placed in the passage 
cax_rotor = np.max(rotor['sections'][0][:,0]) - np.min(stator['sections'][0][:,0]) 
cax_stator = np.max(stator['sections'][0][:,0]) - np.min(stator['sections'][0][:,0]) 
# There should be 21 blades starting with igv and moving into rotor stator pairs
hub = stator['hub']; shroud = stator['shroud'] # rotor has same exact data 
hub_exit_locations = [np.max(rotor['sections'][0][:,0]), np.max(stator['sections'][0][:,0])]; 
shroud_exit_locations = [np.max(rotor['sections'][-1][:,0]), np.max(stator['sections'][-1][:,0])]; 


# Axial location is a percentage along the hub where row exit is defined
rotor1 = make_blade_row(row_type=RowType.Rotor, hub_location=hub_exit_locations[1],shroud_location=shroud_exit_locations[1],stage_id=1)
rotor1.num_blades = blade_counts[0]

stator1 = make_blade_row(row_type=RowType.Stator, hub_location=hub_exit_locations[2],shroud_location=shroud_exit_locations[2],stage_id=1)
stator1.num_blades = blade_counts[1]

# Set axial chord from geometry
rotor1.axial_chord = cax_rotor
stator1.axial_chord = cax_stator
rotor1.beta2_metal = 
stator1.beta2_metal = 
RPM = rotor['outlet'].rpm

rows = [rotor1,stator1]

hub_m = hub*0.0254
shroud_m = shroud*0.0254

passage = Passage(hub_m[:,0],hub_m[:,1],
                 shroud_m[:,0],shroud_m[:,1],
                 passageType=PassageType.Axial) # type: ignore

spool = TurbineSpool(
            passage=passage,
            massflow=20,
            inlet=inlet,
            outlet=outlet,
            rows=rows,
            rpm=RPM,
            num_streamlines=n_streamlines,
            fluid=fluid)


spool.massflow_constraint = MassflowConstraint.PressureBalance # Fixes the exit angle and changes degree of reaction
# spool.plot_geometry()
spool.adjust_streamlines = False
spool.solve() # This also initializes streamlines
spool.export_properties("R35-Results.json")
spool.plot()
spool.plot_velocity_triangles()
