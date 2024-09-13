'''
    GEE3HP Turbine
    2 stage cooled turbine

'''

import sys
sys.path.insert(0,'../')
from td3 import CoolingType, Units, Spool, Inlet, RowType, BladeRow
from td3.coolant import Coolant
from td3.loss
import numpy as np 
from cantera import Solution

coolingType = CoolingType.Cooled_EFFICIENCY_MF_T0_P0 
#%% Initialize the Spool
# Geometry - From TD2
annulus_radii = np.array([[12.6, 12.6, 12.8, 12.73, 12.29, 12.25, 12.25],
                          [14.75, 14.75, 14.44, 14.44, 15.04, 15.04, 15.04]])*0.0254 # Convert to meters 
station_axial_loc = np.array([-1, 0.0, 1.4, 3, 6, 7.8, 9])*0.0254

# Design Conditions 
Design_RPM = 8284
power=3663 * 745.7 # Watts
massflow = 23.3 * 0.453592 # kg/s
P0 = 50*6894.76     # Pascal 
T0 = 1291*0.55      # Kelvin

# Fluid 
fluid = Solution('air.yaml')
fluid.TP = T0, P0 # Use pascal for cantera
print(f"Coefficient of Pressure [J/Kg] {fluid.cp:0.4f}")

# Coolant: Use Kelvin and Pascal
stator_coolant1 = Coolant(fluid, T0=616*0.555556, P0= 50.6 * 6894.76, massflow_percentage=0.0941) 
rotor_coolant1 = Coolant(fluid, 622*0.555556, 50.3 * 6894.76,massflow_percentage=0.07)
stator_coolant2 = Coolant(fluid, 640*0.555556, 23.7 * 6894.76,massflow_percentage=0.0245)
rotor_coolant2 = Coolant(fluid, 640*0.555556, 23.7 * 6894.76,massflow_percentage=0.0)

cooling_type = CoolingType.Cooled_EFFICIENCY_MF_T0_P0
mean_radius = (annulus_radii[0,0]+annulus_radii[1,0])/2 # Mean radius

station1 = Inlet(M=0.4,P0=[P0], T0=[T0], beta=[0], fluid=fluid, radii=mean_radius)

station2 = BladeRow(RowType.Stator, power=0)
station3 = BladeRow(RowType.Rotor, power=power*0.57)
station4 = BladeRow(RowType.Stator,power=0)
station5 = BladeRow(RowType.Rotor,power=power*0.43)

station2.add_coolant(stator_coolant1)
station3.add_coolant(rotor_coolant1)
station4.add_coolant(stator_coolant2)
station5.add_coolant(rotor_coolant2)

# Add in turning angles
station2.set_blade_trailing_edge_geometry([73.1, 74.2, 75.4]) # Angle, hub,mean,tip
station4.set_blade_trailing_edge_geometry([69, 69, 69])

spool = Spool(annulus_radii=annulus_radii,
            x_station=station_axial_loc,
            rpm=Design_RPM, 
            num_streamlines=3, 
            massflow=massflow, 
            rows=[station1,station2,station3,station4,station5])
spool.fluid = fluid

# spool.plot_geometry()
spool.solve() # This also initializes streamlines
print('check')