from dataclasses import dataclass, field
from .enums import RowType
from bladerow import BladeRow
import numpy.typing as npt 

@dataclass
class Coolant:
    T0:float = field(default=900)                               # Kelvin
    P0:float = field(default=50*101325)                         # Pascal
    massflow_percentage:float = field(default=0.03)     # Fraction of total massflow going through compressor
    Cp:float = field(default=1000)                              # J/K
    

def T0_coolant_weighted_average(row:BladeRow) -> npt.NDArray:
    """Calculate the new weighted Total Temperature array considering coolant

    Args:
        coolant (Coolant): Coolant
        massflow (np.ndarray): massflow mainstream

    Returns:
        float: Total Temperature drop
    """
    
    massflow = row.massflow
    total_massflow_no_coolant = row.total_massflow_no_coolant
    Cp = row.Cp
    
    Cpc = row.coolant.Cp
    T0c = row.coolant.T0
    massflow_coolant = row.coolant.massflow_percentage*total_massflow_no_coolant*row.massflow[1:]/row.massflow[-1] 
    if massflow_coolant.mean()>0:
        if row.row_type == RowType.Stator:
            T0= row.T0
            dT0 = T0.copy() * 0 
            T0_new = (massflow[1:]*Cp*T0[1:] + massflow_coolant*Cpc*T0c) \
                        /(massflow[1:]*Cp + massflow_coolant*Cpc)
            dT0[1:] = T0_new - row.T0[1:]
            dT0[0] = dT0[1]
        else:
            T0R = row.T0R
            T0R_new = T0R.copy()
            Cp = row.Cp
            T0R_new[1:] = (massflow[1:]*Cp*T0R[1:] + massflow_coolant*Cpc*T0c) \
                        /(massflow[1:]*Cp + massflow_coolant*Cpc)
            T0R_new[0] = T0R_new[1]
            
            T = T0R_new - row.W**2/(2*Cp)   # Dont change the velocity triangle but adjust the static temperature 
            T0_new = T+row.V**2/(2*Cp)      # Use new static temperature to calculate the total temperature 
            dT0 = T0_new - row.T0
        return dT0
    else:
        return row.T0*0