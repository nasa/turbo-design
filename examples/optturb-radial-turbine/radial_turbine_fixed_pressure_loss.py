'''
    Radial Turbine Example
    
    In this example the blade exit angles are fixed and only degree of reaction changes between the rows to match the massflow
'''
#%% Import Library
import sys
sys.path.insert(0,'../../')
from td3 import PassageType
from td3 import TurbineSpool, Inlet, RowType, BladeRow, Passage, Outlet
from td3.enums import MassflowConstraint
from td3.coolant import Coolant
from td3.loss.turbine import FixedPressureLoss
from pyturbo.helper import bezier
import numpy as np 
from cantera import Solution
from scipy.optimize import minimize_scalar
from scipy.interpolate import pchip
import matplotlib.pyplot as plt

#%% Define the Passage 
# Shroud is defined using a thickness offset from the hub to construct a spline
rhub_ctrl_pts = [0.12,0.10,0.085,
                 0.06,0.04,
                 0.0235, 0.0235,0.0235]

xhub_ctrl_pts = [0.0, 0.0, 0.0,
                 0.02,0.05,
                 0.08,0.12,0.13]

dr = [0.008, 0.008, 0.008, 
      0.015, 0.02,
      0.025,0.025,0.025]
t = [0, 0.1, 0.2,
     0.4, 0.6,
     0.92, 0.98, 1.0]

hub = bezier(xhub_ctrl_pts,rhub_ctrl_pts)
shroud_dh = bezier(t,dr)

def r2(x:float,x1:float,r1:float,slope:float):
    return slope*(x-x1)+r1

def dh_error(x2:float,x1:float,r1:float,dx:float,dr:float,h:float):
    slope = -dx/dr
    r2_guess = r2(x2,x1,r1,slope)
    return np.abs(h-np.sqrt((x1-x2)**2+(r1-r2_guess)**2))

# Build Shroud
npts = 30
xhub,rhub = hub.get_point(np.linspace(0,1,npts))
dx_pts = np.gradient(xhub, np.linspace(0,1,npts))
dr_pts = np.gradient(rhub, np.linspace(0,1,npts))
_, h_pts = shroud_dh.get_point(np.linspace(0,1,npts))
xshroud = xhub*0
rshroud = xhub*0; i = 0
for dx,dr,x1,r1,h in zip(dx_pts,dr_pts,xhub,rhub,h_pts): 
    if abs(dx/dr) >20:
        xshroud[i] = x1
        rshroud[i] = r1+h
    else:
        res = minimize_scalar(dh_error,bounds=[x1,x1+1.5*h],args=(x1,r1,dx,dr,h))
        if r2(res.x,x1,r1,-dx/dr)<r1:
            res = minimize_scalar(dh_error,bounds=[x1-1.5*h,x1],args=(x1,r1,dx,dr,h))
        
        xshroud[i] = res.x
        rshroud[i] = r2(xshroud[i],x1,r1,-dx/dr)
        h_check = np.sqrt((x1-xshroud[i])**2+(r1-rshroud[i])**2)
        # print(f"h = {h} h_check = {h_check}")
        
    i+=1
# plt.figure(num=1,clear=True)
# plt.plot(xhub,rhub)
# plt.plot(xshroud,rshroud,'.')
# plt.axis('scaled')
# plt.show()

passage = Passage(xhub,rhub,
                 xshroud,rshroud,
                 passageType=PassageType.Centrifugal)

# passage.plot_cuts([0,1])
# passage.plot_cuts([0,0.3,0.5,1])

#%% Design Conditions 
Design_RPM = 1000*60
massflow = 0.5 #    # Guessed massflow [kg/s], doesn't matter because code will adjust to match the massflow defined by P0_P 
P0_P = 2            # Total to static Pressure ratio for the entire row
P = 101325          # Outlet Static Pressure [Pascal]
P0 = P * P0_P
T0 = 1000           # Kelvin

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
outlet = Outlet(P=P,percent_radii=0.5,num_streamlines=3)

#%% Define Blade Rows 
# Axial location is a percentage along the hub where row exit is defined
stator0 = BladeRow(row_type=RowType.Stator, axial_location=0.0)
rotor80 = BladeRow(row_type=RowType.Rotor, axial_location=0.80)
rotor100 = BladeRow(row_type=RowType.Rotor, axial_location=1.00)

# Coolant Definition: Use Kelvin and Pascal
stator0.coolant = Coolant(fluid, T0=T0*0.555556,P0=5E5,massflow_percentage=0) 
rotor80.coolant = Coolant(fluid, T0=T0*0.555556,P0=5E5,massflow_percentage=0)
rotor100.coolant = Coolant(fluid, T0=T0*0.555556,P0=5E5,massflow_percentage=0)


# Add in turning angles
stator0.beta2_metal = [60,60,60] # Angle, hub,mean,tip
stator0.loss_model = FixedPressureLoss(0.0)

rotor80.beta2_metal = [-30,-30,-30] # Angle, hub,mean,tip
rotor80.loss_model = FixedPressureLoss(0.15)

rotor100.beta2_metal = [-55,-55,-55] # Angle, hub,mean,tip
rotor100.loss_model = FixedPressureLoss(0.10)
#%% Initialize the Spool
spool = TurbineSpool(passage=passage,
            rpm=Design_RPM, 
            num_streamlines=3, 
            massflow=massflow,
            rows=[inlet,stator0,rotor80,rotor100,outlet])

spool.fluid = fluid
spool.massflow_constraint = MassflowConstraint.BalanceMassFlow
# MassflowConstraint.MatchMassFlow 
# # Fixes the exit angle and changes degree of reaction
# spool.plot_geometry()
spool.solve() # This also initializes streamlines

spool.export_properties("optturb.json")
spool.plot()
spool.plot_velocity_triangles()
print('check')