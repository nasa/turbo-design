from typing import Tuple

import numpy as np
import numpy.typing as npt

from .bladerow import BladeRow
from .enums import RowType

def compute_streamline_areas(row: BladeRow) -> Tuple[float, npt.NDArray]:
    """Compute total annulus area and individual streamline cross sections."""
    total_area = 0.0
    streamline_area = np.zeros(len(row.percent_hub_shroud))
    for j in range(1, len(row.percent_hub_shroud)):
        if np.abs((row.x[j] - row.x[j - 1])) < 1e-5:  # Axial machines
            delta = np.pi * (row.r[j] ** 2 - row.r[j - 1] ** 2)
            streamline_area[j] = delta
            total_area += delta
        else:  # Radial machines
            dx = row.x[j] - row.x[j - 1]
            S = row.r[j] - row.r[j - 1]
            C = np.sqrt(1 + ((row.r[j] - row.r[j - 1]) / dx) ** 2)
            streamline_area[j] = 2 * np.pi * C * (S / 2 * dx ** 2 + row.r[j - 1] * dx)
            total_area += streamline_area[j]
    return total_area, streamline_area


def compute_massflow(row: BladeRow) -> None:
    """Populate row.massflow, total_massflow, and area fields."""
    massflow_fraction = np.linspace(0, 1, len(row.percent_hub_shroud))
    massflow = row.percent_hub_shroud * 0
    total_area, streamline_area = compute_streamline_areas(row)
    for j in range(1, len(row.percent_hub_shroud)):
        Vm = (row.Vm[j] + row.Vm[j - 1]) / 2
        rho = (row.rho[j] + row.rho[j - 1]) / 2
        massflow[j] = Vm * rho * streamline_area[j] * (1 - row.blockage) + massflow[j - 1]

    row.total_massflow_no_coolant = massflow[-1]
    if row.coolant is not None:
        # account for coolant as a fraction of inlet flow
        massflow += massflow_fraction * row.coolant.massflow_percentage * row.total_massflow_no_coolant
    row.massflow = massflow
    row.calculated_massflow = massflow[-1]
    row.total_massflow = massflow[-1]
    row.total_area = total_area
    row.area = streamline_area


def compute_power(row: BladeRow, upstream: BladeRow) -> None:
    """Calculate power and efficiencies for a blade row using upstream reference."""
    if row.row_type == RowType.Stator:
        row.power = 0
        row.eta_static = 0
        row.eta_total = 0
        row.stage_loading = 0
        row.euler_power = 0
        row.T_is = 0 * row.T0
        row.T0_is = 0 * row.T0  # Make it an array
    else:
        P0_P = (upstream.P0 / row.P).mean()
        row.P0_ratio = (row.P0 / upstream.P0).mean()
        row.T_is = upstream.T0 * (1 / P0_P) ** ((row.gamma - 1) / row.gamma)
        a = np.sqrt(row.gamma * row.R * row.T_is)
        row.T0_is = row.T_is * (1 + (row.gamma - 1) / 2 * (row.V / a) ** 2)

        row.power = row.massflow[-1] * (row.Cp * (upstream.T0 - row.T0)).mean()
        row.eta_static = row.power / (row.massflow[-1] * row.Cp * (upstream.T0.mean() - row.T_is.mean()))
        row.eta_total = (upstream.T0.mean() - row.T0.mean()) / (upstream.T0.mean() - row.T0_is.mean())
        row.stage_loading = row.Cp * (upstream.T0.mean() - row.T0.mean()) / row.U.mean() ** 2
        row.euler_power = row.massflow[-1] * (upstream.U * upstream.Vt - row.U * row.Vt).mean()
