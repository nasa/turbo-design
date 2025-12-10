from turbodesign.arrayfuncs import convert_to_ndarray
from .losstype import LossBaseClass
from ..enums import LossType
import numpy.typing as npt
import numpy as np
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..bladerow import BladeRow  # for type hints only


class FixedPolytropicEfficiency(LossBaseClass):
    """Return a fixed polytropic efficiency target for a row (η_poly)."""

    eta_poly: npt.NDArray

    def __init__(self, eta_poly: float | npt.ArrayLike):
        super().__init__(LossType.Polytropic)
        self.eta_poly = convert_to_ndarray(eta_poly)

    def __call__(self, row: "BladeRow", upstream: "BladeRow") -> npt.NDArray:  # noqa: ARG002
        eta = self.eta_poly
        if eta.size == 1:
            eta = eta * np.ones_like(row.r)  # type: ignore[arg-type]
        elif eta.shape != row.r.shape:
            eta = np.asarray(eta).reshape(row.r.shape)  # type: ignore[attr-defined]
        return eta
