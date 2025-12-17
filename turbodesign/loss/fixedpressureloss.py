from turbodesign.arrayfuncs import convert_to_ndarray
from .losstype import LossBaseClass
from ..enums import LossType
import numpy.typing as npt 
import numpy as np
from scipy.interpolate import interp1d
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from ..bladerow import BladeRow  # for type hints only

class FixedPressureLoss(LossBaseClass):
    """Fixed pressure loss coefficient (scalar or spanwise array)."""

    pressure_loss: npt.NDArray
    
    def __init__(self, pressure_loss: float | npt.ArrayLike):
        super().__init__(LossType.Pressure)
        self.pressure_loss = convert_to_ndarray(pressure_loss)
    
    def __call__(self, row: "BladeRow", upstream: "BladeRow") -> npt.NDArray:
        """Outputs the fixed pressure loss."""
        loss = self.pressure_loss
        if loss.size == 1:
            loss = loss * np.ones_like(row.r) # type: ignore
        elif loss.shape != row.r.shape:
            if len(row.r) == 1:
                return loss.mean()
            else:
                return interp1d(np.linspace(0,1,len(loss)),loss)(row.percent_hub_shroud)
        return loss
