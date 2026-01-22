from typing import Optional, TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from ..losstype import LossBaseClass
from ...enums import LossType, RowType

if TYPE_CHECKING:
    from ...bladerow import BladeRow  # for type hints only


class DiffusionLoss(LossBaseClass):
    """Pressure-loss model based on diffusion factor.

    Computes a spanwise diffusion factor and maps it to a pressure-loss
    coefficient Yp via a simple linear ramp above a threshold.
    """

    def __init__(self, df_limit: float = 0.45, df_knee: float = 0.6, df_max: float = 0.9, yp_at_max: float = 0.08):
        """
        Args:
            df_limit: Diffusion factor below which loss is negligible.
            df_knee: Diffusion factor where loss starts ramping up noticeably.
            df_max: Diffusion factor at which Yp reaches ``yp_at_max``.
            yp_at_max: Pressure-loss coefficient when ``df_max`` is reached.
        """
        super().__init__(LossType.Pressure)
        self.df_limit = df_limit
        self.df_knee = df_knee
        self.df_max = df_max
        self.yp_at_max = yp_at_max

    def __call__(self, row: "BladeRow", upstream: "BladeRow") -> npt.NDArray:
        """Return pressure-loss coefficient Yp based on diffusion factor."""
        # Choose absolute vs relative velocities depending on row type
        if row.row_type == RowType.Rotor:
            V1 = upstream.W
            V2 = row.W
            Vt1 = upstream.Wt
            Vt2 = row.Wt
            U = upstream.U
        else:
            V1 = upstream.V
            V2 = row.V
            Vt1 = upstream.Vt
            Vt2 = row.Vt
            U = 0.0

        V1_mag = np.maximum(np.abs(V1), 1e-6)
        V2_mag = np.abs(V2)
        dVt = Vt2 - Vt1

        df = 1.0 - V2_mag / V1_mag + dVt / np.maximum(np.abs(U), 1e-6)

        # Map diffusion factor to Yp
        Yp = np.zeros_like(row.percent_hub_shroud, dtype=float)
        ramp = (df - self.df_limit) / max(self.df_knee - self.df_limit, 1e-6)
        ramp = np.clip(ramp, 0.0, 1.0)
        Yp = ramp * (self.yp_at_max * np.clip((df - self.df_limit) / max(self.df_max - self.df_limit, 1e-6), 0.0, 1.0))
        return Yp