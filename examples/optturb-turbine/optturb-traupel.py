"""
    GEE3HP Turbine
    2 stage cooled turbine

"""
from turbodesign import PassageType, TurbineSpool, Inlet, RowType, BladeRow, Passage, Outlet
from turbodesign.row_factory import make_rotor_row, make_stator_row
from turbodesign.coolant import Coolant
from turbodesign.loss.turbine import Traupel
import numpy as np 
from cantera import Solution

#%% Initialize the TurbineSpool
# Geometry - From TD2
rmean = 0.389
H1 = 0.063
H2 = 1.159 * H1
H3 = 1.317 * H2
cax = (H1 + H2 + H3) / 3
                        # Inlet, Stator Inlet, Stator Exit, Rotor Exit
rhub = [rmean - H1 / 2, rmean - H1 / 2, rmean - H2 / 2, rmean - H3 / 2]
rshroud = [rmean + H1 / 2, rmean + H1 / 2, rmean + H2 / 2, rmean + H3 / 2]
xhub = np.array([-cax, 0.0, cax, 2 * cax])
xshroud = np.array([-cax, 0.0, cax, 2 * cax])
axial_len = xhub[-1] - xhub[0]

passage = Passage(
    xhub,
    rhub,
    xshroud,
    rshroud,
    passageType=PassageType.Axial,
)

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


station1 = Inlet(hub_location=0, alpha=[0])
station1.init_total(P0=[P0], T0=[T0], M=[0.4], percent_radii=[0.5])

station2 = make_stator_row(hub_location=2 * cax / axial_len)
station3 = make_rotor_row(hub_location=3 * cax / axial_len)
station3.power = power

outlet = Outlet(num_streamlines=3)
outlet.init_static(P=P0 / 3.96, percent_radii=[0.5])


station2.coolant = Coolant(T0=616*0.55, P0=50.6*6894.76, massflow_percentage=0, Cp=1012)
station3.coolant = Coolant(T0=616*0.55, P0=50.6*6894.76, massflow_percentage=0, Cp=1012)

# Add in turning angles
station2.beta2_metal = [73,73,73] # Angle, hub,mean,tip
station2.loss_model = Traupel()
station3.loss_model = Traupel()

spool = TurbineSpool(
    passage=passage,
    massflow=massflow,
    inlet=station1,
    outlet=outlet,
    rows=[station2, station3],
    rpm=Design_RPM,
    num_streamlines=3,
)
spool.fluid = fluid

# spool.plot_geometry()
spool.solve() # This also initializes streamlines
spool.export_properties("optturb.json")
spool.plot()
spool.plot_velocity_triangles()



print('check')
