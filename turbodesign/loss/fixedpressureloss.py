from .losstype import LossBaseClass
from ..enums import LossType
import numpy.typing as npt 
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from ..bladerow import BladeRow  # for type hints only

class FixedPressureLoss(LossBaseClass):
    pressure_loss:float
    
    def __init__(self,pressure_loss:float):
        """Fixed Pressure Loss
        """
        super().__init__(LossType.Pressure)
        self.pressure_loss = pressure_loss
    
    
    def __call__(self, row: "BladeRow", upstream: "BladeRow") -> npt.NDArray:
        """Outputs the fixed Pressure Loss
        
        Args:
            upstream (BladeRow): Upstream blade row
            row (BladeRow): downstream blade row

        Returns:
            float: Pressure Loss
        """
        Yp = row.r*0+self.pressure_loss
        return Yp