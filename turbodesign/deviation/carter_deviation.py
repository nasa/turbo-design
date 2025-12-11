from typing import Optional, TYPE_CHECKING

import numpy as np
import numpy.typing as npt

from ..arrayfuncs import convert_to_ndarray
from ..enums import RowType
from .deviation_base import DeviationBaseClass

if TYPE_CHECKING:
    from ..bladerow import BladeRow  # for type hints only


class CarterDeviation(DeviationBaseClass):
    """Carter-style exit deviation for axial compressor blades.

    Uses Mattingly's form:
        gamma_e = (4 * alpha_exit * sqrt(sigma) - gamma_i) / (4 * sqrt(sigma) - 1)
        deviation = alpha_exit - gamma_e

    where gamma_i/gamma_e are inlet/exit metal (or flow) angles and sigma is solidity.
    """

    def __init__(self, gamma_inlet: Optional[float | npt.ArrayLike] = None, gamma_exit: Optional[float | npt.ArrayLike] = None, alpha_exit: Optional[float | npt.ArrayLike] = None):
        """
        Args:
            gamma_inlet: Optional inlet metal/camber angle array (radians). If not
                provided, defaults to ``row.alpha1`` (stator) or ``row.beta1`` (rotor).
            gamma_exit: Optional exit metal/camber angle array (radians). If not
                provided, defaults to ``row.alpha2`` (stator) or ``row.beta2`` (rotor).
            alpha_exit: Optional flow exit angle array (radians). If not provided,
                defaults to ``row.alpha2`` for stators and ``row.beta2`` for rotors.
        """
        super().__init__()
        self.gamma_inlet = gamma_inlet
        self.gamma_exit = gamma_exit
        self.alpha_exit = alpha_exit

    def _spanwise(self, value: float | npt.ArrayLike, row: "BladeRow") -> npt.NDArray:
        """Convert scalar/array input to spanwise distribution on row grid."""
        arr = convert_to_ndarray(value)
        if arr.size == 1:
            return arr * np.ones_like(row.percent_hub_shroud, dtype=float)
        t_src = np.linspace(0, 1, arr.size)
        return np.interp(row.percent_hub_shroud, t_src, arr)

    def __call__(self, row: "BladeRow", upstream: "BladeRow") -> npt.NDArray:  # noqa: ARG002
        """Compute deviation (radians) along the span for the supplied row."""
        gamma_i = self.gamma_inlet
        gamma_e = self.gamma_exit
        alpha_exit = self.alpha_exit

        if gamma_i is None:
            gamma_i = getattr(row, "gamma_inlet", None)
        if gamma_e is None:
            gamma_e = getattr(row, "gamma_exit", None)
        if alpha_exit is None:
            alpha_exit = row.alpha2 if row.row_type == RowType.Stator else row.beta2

        # Default to blade flow angles if no explicit camber/metal angles are provided.
        if gamma_i is None:
            gamma_i = row.alpha1 if row.row_type == RowType.Stator else row.beta1
        if gamma_e is None:
            gamma_e = row.alpha2 if row.row_type == RowType.Stator else row.beta2

        gamma_i_span = self._spanwise(gamma_i, row)
        alpha_exit_span = self._spanwise(alpha_exit, row)

        sigma = getattr(row, "solidity", 0.0) or 0.0
        sigma = max(float(sigma), 1e-6)  # avoid divide-by-zero
        k = np.sqrt(sigma)

        if gamma_e is not None:
            gamma_e_span = self._spanwise(gamma_e, row)
        else:
            gamma_e_span = (4.0 * alpha_exit_span * k - gamma_i_span) / np.maximum(4.0 * k - 1.0, 1e-6)

        deviation = alpha_exit_span - gamma_e_span
        return deviation
