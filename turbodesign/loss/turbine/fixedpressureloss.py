from ...bladerow import BladeRow
from ..losstype import LossBaseClass
from ...enums import LossType

class FixedPressureLoss(LossBaseClass):
    pressure_loss:float
    
    def __init__(self,pressure_loss:float):
        """Fixed Pressure Loss
        """
        super().__init__(LossType.Pressure)
        self.pressure_loss = pressure_loss
    
    
    def __call__(self,row:BladeRow, upstream:BladeRow) -> float:
        """Outputs the fixed Pressure Loss
        
        Args:
            upstream (BladeRow): Upstream blade row
            row (BladeRow): downstream blade row

        Returns:
            float: Pressure Loss 
        """
        return self.pressure_loss