from ...bladerow import BladeRow, sutherland
from ...enums import RowType, LossType
from ..losstype import LossBaseClass
import numpy.typing as npt 

class FixedEfficiency(LossBaseClass):
    efficiency:float
    
    def __init__(self,efficiency:float):
        """Fixed Efficiency Loss 
        """
        super().__init__(LossType.Enthalpy)
        self.efficiency = efficiency
    
    
    def __call__(self,row:BladeRow, upstream:BladeRow) -> npt.NDArray:
        """Fixed efficiency loss 
        
        Args:
            row (BladeRow): Blade row being evaluated.
            upstream (BladeRow): Upstream blade row (unused, kept for API parity).

        Returns:
            numpy.ndarray: Spanwise efficiency array; zeros for stators, fixed value for rotors.
        """
        if row.row_type == RowType.Stator:
            return row.r*0 
        else:
            return self.efficiency + row.r*0
        
