from typing import List, Tuple
from .radeq import radeq, radeq_normalized
from .enums import LossType, RowType, PowerType, MassflowConstraint
from .bladerow import BladeRow
from .td_math import compute_gas_constants
from .td_math import compute_quantities, compute_power, compute_massflow
import numpy.typing as npt
import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar
from .passage import Passage     
    
def adjust_streamlines(blade_rows:List[BladeRow],passage:Passage):
    """Adjust the streamlines to evenly divide the massflow

    Args:
        blade_rows (List[BladeRow]): List of blade rows
        passage (Passage): passage object describing the hub and shroud 

    """
    for _,row in enumerate(blade_rows):
        massflow_fraction =  np.linspace(0,1,len(row.percent_hub_shroud))
        row.total_massflow = row.massflow[-1]
        ideal_massflow_fraction = row.massflow[-1] * massflow_fraction
        
        new_percent_streamline = interp1d(row.massflow,row.percent_hub_shroud)(ideal_massflow_fraction[1:-1])
        row.percent_hub_shroud[1:-1] = new_percent_streamline

        cut_line, thub,_ = passage.get_cutting_line(row.axial_location)
        row.x,row.r = cut_line.get_point(row.percent_hub_shroud)
        # Radii may have shifted, recompute Ay and rm
        for i,tr in enumerate(row.percent_hub_shroud):
            t_streamline, x_streamline, r_streamline = passage.get_streamline(tr)                
            phi, rm, r = passage.streamline_curvature(x_streamline,r_streamline)
            row.phi[i] = float(interp1d(t_streamline,phi)(row.axial_location))
            row.rm[i] = float(interp1d(t_streamline,rm)(row.axial_location))
            row.r[i] = float(interp1d(t_streamline,r)(row.axial_location))