from ...bladerow import BladeRow, sutherland
from ...enums import RowType, LossType
from ..losstype import LossBaseClass

class FixedEfficiency(LossBaseClass):
    efficiency:float
    
    def __init__(self,efficiency:float):
        """Fixed Efficiency Loss 
        """
        super().__init__(LossType.Enthalpy)
        self.efficiency = efficiency
    
    
    def __call__(self,row:BladeRow, upstream:BladeRow) -> float:
        """Fixed efficiency loss 
        
        Args:
            upstream (BladeRow): Upstream blade row
            row (BladeRow): downstream blade row

        Returns:
            float: Stage Efficiency
        """
        if row.row_type == RowType.Stator:
            return 0
        else:
            return self.efficiency
        