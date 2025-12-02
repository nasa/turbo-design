from turbodesign import TurbineSpool, Inlet, RowType, BladeRow, Passage, Outlet, PassageType, Coolant
from turbodesign.enums import MassflowConstraint
from turbodesign.loss.fixedpressureloss import FixedPressureLoss
import numpy as np
from cantera import Solution
import pickle
from pathlib import Path
from scipy.interpolate import PchipInterpolator
# Geometry Import 
from get_ss_ps import split_airfoil_by_angle_distance, resample_curve, plot_blade

blades = pickle.load(open(Path(__file__).resolve().parent / 'stator_rotor.pkl','rb'))

processed_data = []
npts = 400
for i,blade in enumerate(blades):
    nsections = len(blade)
    n,m = blade[0].shape
    ss = np.zeros(shape=(nsections,npts,m))
    ps = np.zeros(shape=(nsections,npts,m))
    for section_index in range(nsections):
        ss_temp, ps_temp = split_airfoil_by_angle_distance(blade[section_index],npts=npts)
        # plot_blade(ss,ps,'test')
        ss[section_index,:,:] = ss_temp
        ps[section_index,:,:] = ps_temp
    ss = ss[np.argsort(ss[:,0,2]),:,:]
    ps = ps[np.argsort(ps[:,0,2]),:,:]
    ss *= 25.4 # convert to mm
    ps *= 25.4 # convert to mm
    print(f'Blade-{i} has {nsections} sections')
    processed_data.append((ss, ps))


data_dir = Path(__file__).resolve().parent

hub_shroud = pickle.load(open(data_dir / 'hub_shroud.pkl','rb'))
x = np.linspace(hub_shroud['Hub'][:,0].min(),hub_shroud['Hub'][:,0].max(),200)

hub = np.vstack([x,PchipInterpolator(hub_shroud['Hub'][:,0],hub_shroud['Hub'][:,1])(x)]).transpose()
hub *= 25.4          # Convert in to mm

shroud = np.vstack([x,PchipInterpolator(hub_shroud['Shroud'][:,0],hub_shroud['Shroud'][:,1])(x)]).transpose()
shroud *= 25.4 

# TurboDesign Setup 
nblades = [46,76,48,70] # Vanes, Rotors, Vanes, Rotors
P0 = 1257450        # Pa
T0 = 1587           # K
P = 230.295*1000    # Pa
n_streamlines = 5

# Fluid
fluid = Solution('air.yaml')
fluid.TP = T0, P0 # Use pascal for cantera
print(f"Coefficient of Pressure [J/Kg] {fluid.cp:0.4f}")
#%% Defining the Inlet
inlet = Inlet(beta=[0,0], hub_location=0)
inlet.init_total(P0=[P0,P0],T0=[T0,T0],M=0.1)
outlet = Outlet(num_streamlines=n_streamlines)
outlet.init_static(P=P,percent_radii=[0.5])
#%% Define Blade Rows, processed data is already in mm, hub is already in mm 
cax1 = ( processed_data[0][0][0,:,0].max()-processed_data[0][0][0,:,0].min() )/1000
cax2 = ( processed_data[1][0][0,:,0].max()-processed_data[1][0][0,:,0].min() )/1000
cax3 = ( processed_data[2][0][0,:,0].max()-processed_data[2][0][0,:,0].min() )/1000
cax4 = ( processed_data[3][0][0,:,0].max()-processed_data[3][0][0,:,0].min() )/1000
location1 = (processed_data[0][0][0,:,0].max() - hub[:,0].min()) / (hub[:,0].max()-hub[:,0].min())  # Exit Locations
location2 = (processed_data[1][0][0,:,0].max() - hub[:,0].min()) / (hub[:,0].max()-hub[:,0].min())
location3 = (processed_data[2][0][0,:,0].max() - hub[:,0].min()) / (hub[:,0].max()-hub[:,0].min())
location4 = (processed_data[3][0][0,:,0].max() - hub[:,0].min()) / (hub[:,0].max()-hub[:,0].min())
beta_exit_flow = [73.6,-67.2,69.5,-63.9]
P0_Loss = [0.057,0.088,0.069,0.014]         # (P01-P02)/(P01-P2)

# Axial location is a percentage along the hub where row exit is defined
stator1 = BladeRow(row_type=RowType.Stator, hub_location=location1)
rotor1 = BladeRow(row_type=RowType.Rotor, hub_location=location2)
stator2 = BladeRow(row_type=RowType.Stator, hub_location=location3)
rotor2 = BladeRow(row_type=RowType.Rotor, hub_location=location4)

stator1.axial_chord = cax1 # Set an axial chord
rotor1.axial_chord = cax2
stator2.axial_chord = cax3
rotor2.axial_chord = cax4

inlet.gamma = 1.30
stator1.gamma = 1.30
rotor1.gamma = 1.32
stator2.gamma = 1.33
rotor2.gamma = 1.39

inlet.Cp = inlet.gamma/(inlet.gamma-1) * 287.15
stator1.Cp = stator1.gamma/(stator1.gamma-1) * 287.15
rotor1.Cp = rotor1.gamma/(rotor1.gamma-1) * 287.15
stator2.Cp = stator2.gamma/(stator2.gamma-1) * 287.15
rotor2.Cp = rotor2.gamma/(rotor2.gamma-1) * 287.15

stator1.stage_id = 0; rotor1.stage_id = 0
stator2.stage_id = 1; rotor2.stage_id = 1

# Coolant Definition: Use Kelvin and Pascal. Coolant only needs P0, T0, massflow, and Cp
stator1.coolant = Coolant(T0=T0*0.6, P0= P0, massflow_percentage=0,Cp=fluid.cp) 
rotor1.coolant = Coolant(T0=T0*0.6, P0=P0,massflow_percentage=0,Cp=fluid.cp)
stator2.coolant = Coolant(T0=T0*0.6, P0=P0, massflow_percentage=0,Cp=fluid.cp) 
rotor2.coolant = Coolant(T0=T0*0.6, P0=P0,massflow_percentage=0,Cp=fluid.cp)

# Add in turning angles
stator1.beta2_metal = [beta_exit_flow[0] for _ in range(n_streamlines)]        # Angle, hub,mean,tip
rotor1.beta2_metal = [beta_exit_flow[1] for _ in range(n_streamlines)]
stator2.beta2_metal = [beta_exit_flow[2] for _ in range(n_streamlines)]
rotor2.beta2_metal = [beta_exit_flow[3] for _ in range(n_streamlines)]


# These are all guessed values
stator1.loss_model = FixedPressureLoss(P0_Loss[0])  # type: ignore
rotor1.loss_model = FixedPressureLoss(P0_Loss[1])   # type: ignore
stator2.loss_model = FixedPressureLoss(P0_Loss[2])  # type: ignore
rotor2.loss_model = FixedPressureLoss(P0_Loss[3])   # type: ignore
hub_m = hub/1000
shroud_m = shroud/1000
passage = Passage(hub_m[:,0],hub_m[:,1],
                 shroud_m[:,0],shroud_m[:,1],
                 passageType=PassageType.Axial) # type: ignore

spool = TurbineSpool(passage=passage,
            massflow=20,
            inlet=inlet,
            outlet=outlet,
            rows=[stator1,rotor1,stator2,rotor2],
            rpm=12400,
            num_streamlines=n_streamlines,
            fluid=None)
spool.massflow_constraint = MassflowConstraint.PressureBalance # Fixes the exit angle and changes degree of reaction
# spool.plot_geometry()
spool.adjust_streamlines = False
spool.solve() # This also initializes streamlines
spool.export_properties(str(data_dir / "eee_results_P0.json"))
spool.plot()
spool.plot_velocity_triangles()
