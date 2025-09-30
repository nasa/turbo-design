from typing import List, Optional, Tuple
import numpy as np
import numpy.typing as npt
from .bladerow import BladeRow, compute_gas_constants
from .enums import RowType, LossType
from scipy.integrate import trapezoid
from .passage import Passage
from .isentropic import IsenP

# For compressors rotor calc is called first
def compressor_rotor_calc(row:BladeRow,upstream:BladeRow,calculate_vm:bool=True):
    upstream.T0R = upstream.T+upstream.W**2/(2*upstream.Cp)
    upstream.P0R = upstream.P * (upstream.T0R/upstream.T)**((upstream.gamma)/(upstream.gamma-1))      
    upstream.M_rel = upstream.W/np.sqrt(upstream.gamma*upstream.R*upstream.T)
    upstream.U = upstream.rpm*np.pi/30 * upstream.r
    upstream.Wt = upstream.Vt - upstream.U
    upstream.W = np.sqrt(upstream.Vx**2 + upstream.Wt**2 + upstream.Vr**2)
    
    row.U = row.omega*row.r
    pass 