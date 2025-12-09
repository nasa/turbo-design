
from ast import Pass
from math import radians
from multiprocessing.heap import Arena
from turbodesign import Inlet, RowType, BladeRow, Passage, Outlet, PassageType
from turbodesign.compressor_spool import CompressorSpool
from turbodesign.enums import MassflowConstraint
from turbodesign import Coolant
from turbodesign.loss import fixedpressureloss
from turbodesign.loss.fixedpressureloss import FixedPressureLoss
from turbodesign.deviation.fixed_deviation import FixedDeviation
from cantera import Solution
from pathlib import Path
import numpy as np 


def MFP(M:float,gamma:float=1.4, R:float=287.15):
    expo = -(gamma+1)/(2*(gamma-1))
    return np.sqrt(gamma/R) * M * (1+(gamma-1)/2 * M*M)**expo


# Knowns
T01 = 518.7 # R
P01 = 14.7 # psia
omega = 1000 # rad/s
r = 12 # in
alpha1 = 40 # deg
massflow = 50 # lbm/s
M1 = 0.7 
M3 = M1 
u2_u1 = 1.1 

T01_K = 518.7 / 1.8
P01_Pa = P01 * 6894.76
massflow_kg_s = massflow * 0.453592

Area = massflow_kg_s * np.sqrt(T01_K)/(P01_Pa*np.cos(np.radians(alpha1)) * MFP(M1))

# print(MFP(M1))
# print(MFP(M1,1.4,1716/144))
# print(Area_throat*39.3701**2)

rmean = 12 * 0.0254 # convert inch to meter
h = Area / (np.pi * 4*rmean)
cax = 1 * 0.0254                # Assumed axial chord of 1 inch 
xhub_arr = [0,cax,2*cax]
xshroud_arr = [0,cax,2*cax]
rhub2 = rmean*u2_u1 - h   #
rhub_arr = [rmean-h, rhub2, rhub2]  # Inlet Exit, Rotor Exit, Stator Exit
rshroud_arr = [rmean+h, rmean+h, rmean+h]

passage = Passage(xhub_arr,rhub_arr,xshroud_arr,rshroud_arr,passageType=PassageType.Axial)
inlet = Inlet(hub_location=0)
inlet.alpha2 = [40]
inlet.init_total(P01_Pa,T01_K,M=M1)

rotor = BladeRow(hub_location=0,row_type=RowType.Rotor)
rotor.beta2_metal = [23.87]
rotor.loss_function = FixedPressureLoss(0)
stator = BladeRow(hub_location=0,row_type=RowType.Stator)
stator.beta2_metal = [40] 
stator.loss_function = FixedPressureLoss(0)
outlet = Outlet()
outlet.init_total(1.3*P01_Pa,0.5)

spool = CompressorSpool(passage,massflow_kg_s,inlet,outlet,[rotor,stator],rpm=omega*30/np.pi)
spool.solve()