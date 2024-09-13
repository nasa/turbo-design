'''
    GEE3HP Turbine
    2 stage cooled turbine

'''

import sys
sys.path.insert(0,'../../')
from td3 import CoolingType, LossType, Units, TurbineSpool, Inlet, RowType, BladeRow
from td3.coolant import Coolant
from td3.loss.turbine import Traupel
import numpy as np 
from cantera import Solution

lossType = LossType.TD2_Pressure_Loss
coolingType = CoolingType.Cooled_EFFICIENCY_MF_T0_P0 
#%% Initialize the Spool
# Geometry - From TD2
rmean = 0.389
H1 = 0.063
H2 = 1.159*H1
H3 = 1.317*H2
cax = (H1+H2+H3)/3 
                        # Inlet, Stator Inlet, Stator Exit, Rotor Exit
rhub = [rmean-H1/2,rmean-H1/2,rmean-H2/2,rmean-H3/2]
rshroud = [rmean+H1/2,rmean+H1/2,rmean+H2/2,rmean+H3/2]
xhub = np.array([-cax, 0.0, cax, 2*cax])
xshroud = np.array([-cax, 0.0, cax, 2*cax])

# Design Conditions 
Design_RPM = 7500
power = 5.74E6 # Watts
massflow = 35.9 # kg/s
P0 = 500000     # Pascal 
T0 = 676.3      # Kelvin

# Fluid 
fluid = Solution('air.yaml')
fluid.TP = T0, P0 # Use pascal for cantera
print(f"Coefficient of Pressure [J/Kg] {fluid.cp:0.4f}")

# Coolant: Use Kelvin and Pascal


cooling_type = CoolingType.Cooled_EFFICIENCY_MF_T0_P0

station1 = Inlet(M=0.4,P0=[P0], T0=[T0], beta=[0], fluid=fluid, percent_radii=0.5)
station2 = BladeRow(RowType.Stator, power=0)
station3 = BladeRow(RowType.Rotor, power=power)

station2.coolant = Coolant(fluid, T0=616*0.555556, P0= 50.6 * 6894.76, massflow_percentage=0) 
station3.coolant = Coolant(fluid, 622*0.555556, 50.3 * 6894.76,massflow_percentage=0)

# Add in turning angles
station2.beta2_metal = [73,73,73] # Angle, hub,mean,tip
station2.loss_model = Traupel()
station3.loss_model = Traupel()

spool = TurbineSpool(rhub=rhub,xhub=xhub,
            rshroud=rshroud,xshroud=xshroud,
            rpm=Design_RPM, 
            num_streamlines=3, 
            massflow=massflow, 
            rows=[station1,station2,station3])
spool.fluid = fluid

# spool.plot_geometry()
spool.solve() # This also initializes streamlines
spool.export_properties("optturb.json")
spool.plot()
spool.plot_velocity_triangles()



print('check')