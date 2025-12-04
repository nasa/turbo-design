from typing import Optional, TYPE_CHECKING

import numpy as np
import numpy.typing as npt
from scipy.interpolate import interp1d

from ..arrayfuncs import convert_to_ndarray
from .deviation_base import DeviationBaseClass

if TYPE_CHECKING:
    from ..bladerow import BladeRow  # for type hints only


class FixedDeviation(DeviationBaseClass):
    """Fixed deviation angle (scalar or spanwise distribution)."""

    t: Optional[npt.ArrayLike]
    deviation_angles: npt.NDArray
    
    def __init__(self, deviation_angles: float | npt.ArrayLike, t: Optional[npt.ArrayLike] = None):
        """
        Args:
            deviation_angles: deviation angle in degrees; either scalar or an array
                along the span (0 at hub to 1 at shroud).
            t: optional spanwise coordinates corresponding to ``deviation_angles``.
                If not provided and ``deviation_angles`` is array-like, a uniform
                spanwise distribution is assumed.
        """
        super().__init__()
        self.deviation_angles = convert_to_ndarray(deviation_angles)
        self.t = t
    
    def __call__(self, row: "BladeRow", upstream: "BladeRow") -> npt.NDArray:
        """Outputs the fixed deviation distribution."""
        deviation_angles = self.deviation_angles
        if deviation_angles.size == 1:
            deviation = deviation_angles * np.ones_like(row.r)
        else:
            if self.t is None:
                self.t = np.linspace(0, 1, len(self.deviation_angles))
            deviation = interp1d(self.t, self.deviation_angles, bounds_error=True)(row.percent_hub_shroud)
        return deviation
