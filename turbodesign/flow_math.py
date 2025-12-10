from typing import Tuple

import numpy as np
import numpy.typing as npt

from .bladerow import BladeRow


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
