from typing import Tuple

import numpy as np
import numpy.typing as npt

from .bladerow import BladeRow
from .enums import RowType

def compute_streamline_areas(row: BladeRow) -> Tuple[float, npt.NDArray]:
    """Compute total annulus area and individual streamline cross-sectional areas.

    Calculates the total annulus area and the cross-sectional area for each streamtube
    based on the radial (r) and axial (x) coordinates of the blade row. Handles both
    axial machines (constant x) and radial machines (varying x).

    Args:
        row: BladeRow object containing percent_hub_shroud, x, r coordinates

    Returns:
        tuple: (total_area, streamline_area) where
            - total_area (float): Total annulus cross-sectional area [m²]
            - streamline_area (ndarray): Array of streamtube areas [m²] matching row.r shape
    """
    total_area = 0.0
    streamline_area = np.zeros(len(row.percent_hub_shroud))
    if len(row.percent_hub_shroud) <= 1:
        if hasattr(row, "total_area") and row.total_area:
            total_area = float(row.total_area)
            streamline_area = np.array([total_area])
        return total_area, streamline_area
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
    """Calculate massflow distribution across streamlines and populate row attributes.

    Computes the cumulative massflow through each streamtube based on density, meridional
    velocity, and streamtube cross-sectional areas. Accounts for blockage and optional
    coolant injection. Updates row attributes in-place.

    Args:
        row: BladeRow object with Vm, rho, percent_hub_shroud, blockage defined

    Returns:
        None. Updates the following row attributes in-place:
            - row.massflow: Cumulative massflow array [kg/s]
            - row.total_massflow: Total massflow including coolant [kg/s]
            - row.total_massflow_no_coolant: Massflow without coolant [kg/s]
            - row.calculated_massflow: Final massflow value [kg/s]
            - row.total_area: Total annulus area [m²]
            - row.area: Streamtube areas array [m²]
    """
    n = len(row.percent_hub_shroud)
    massflow_fraction = np.linspace(0, 1, n)
    total_area, streamline_area = compute_streamline_areas(row)

    if n <= 1:
        Vm = float(row.Vm[0]) if len(row.Vm) else 0.0
        rho = float(row.rho[0]) if len(row.rho) else 0.0
        mass = Vm * rho * (total_area if total_area else 0.0) * (1 - row.blockage)
        massflow = np.array([mass])
        row.total_massflow_no_coolant = mass
        if row.coolant is not None:
            massflow += row.coolant.massflow_percentage * mass
        row.massflow = massflow
        row.calculated_massflow = massflow[-1]
        row.total_massflow = massflow[-1]
        row.total_area = total_area
        row.area = streamline_area
        return

    massflow = np.zeros_like(row.percent_hub_shroud, dtype=float)
    for j in range(1, n):
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


def compute_power(row: BladeRow, upstream: BladeRow | None = None, downstream: BladeRow | None = None, is_compressor: bool | None = None) -> None:
    """Calculate power and efficiencies for a blade row (compressor or turbine).

    Args:
        row: The blade row being evaluated.
        upstream: Upstream reference row (default for turbine-style calculations).
        downstream: Downstream reference row (optional; useful for compressor-style staging).
        is_compressor: Force compressor sign convention when True (power added to flow);
            when False, assumes turbine (power extracted). If None, infers from P0 gain.
    """
    ref = upstream if upstream is not None else downstream
    if ref is None:
        return

    mdot = row.massflow[-1] if getattr(row, "massflow", np.array([])).size else getattr(row, "total_massflow", 0.0)

    if row.row_type == RowType.Stator:
        row.power = 0.0
        row.eta_static = 0.0
        row.eta_total = 0.0
        row.stage_loading = 0.0
        row.euler_power = 0.0
        row.T_is = 0 * row.T0
        row.T0_is = 0 * row.T0  # Make it an array
    else:
        # Preserve any user-configured target ratio. compute_power historically overwrote row.P0_ratio,
        # which makes it hard to treat P0_ratio as a design input elsewhere.
        if getattr(row, "P0_ratio_target", 0.0) == 0 and getattr(row, "P0_ratio", 0.0) != 0:
            row.P0_ratio_target = row.P0_ratio

        P0_P = (ref.P0 / row.P).mean()
        P0_ratio_actual = (row.P0 / ref.P0).mean()
        row.P0_ratio = P0_ratio_actual
        setattr(row, "P0_ratio_actual", float(P0_ratio_actual))
        row.T_is = ref.T0 * (1 / P0_P) ** ((row.gamma - 1) / row.gamma)
        row.T0_is = ref.T0 * (row.P0 / ref.P0) ** ((row.gamma - 1) / row.gamma)

        comp_mode = is_compressor
        if comp_mode is None:
            comp_mode = bool(np.mean(row.P0) > np.mean(ref.P0))

        if comp_mode:
            deltaT = row.T0.mean() - ref.T0.mean()
            row.power = mdot * row.Cp * deltaT
            denom_static = max(row.T.mean() - ref.T0.mean(), 1e-9)
            denom_total = max(row.T0.mean() - ref.T0.mean(), 1e-9)
            row.eta_static = (row.T_is.mean() - ref.T0.mean()) / denom_static
            row.eta_total = (row.T0_is.mean() - ref.T0.mean()) / denom_total
        else:
            deltaT = ref.T0.mean() - row.T0.mean()
            row.power = mdot * row.Cp * deltaT
            row.eta_static = row.power / (mdot * row.Cp * (ref.T0.mean() - row.T_is.mean()))
            # Entropy-based total-total efficiency:  η = w / (w + T_exit·Δs)
            # The standard isentropic formula η = ΔT0/(T01−T0_is) uses the
            # absolute P0 ratio which, for radial machines with large radius
            # change, is dominated by the frame change and barely reflects the
            # relative-frame loss — giving η ≈ 1 even with significant Yp.
            # The entropy-based definition always isolates the irreversibility.
            if np.mean(ref.P0R) > 0 and np.mean(row.P0R) > 0 and deltaT > 0:
                ds = row.R * np.log(np.mean(ref.P0R) / np.mean(row.P0R))
                w_per_mass = row.Cp * deltaT
                row.eta_total = w_per_mass / (w_per_mass + row.T.mean() * max(ds, 0.0))
            else:
                row.eta_total = (ref.T0.mean() - row.T0.mean()) / max(ref.T0.mean() - row.T0_is.mean(), 1e-9)
        
        row.stage_loading = row.Cp * (ref.T0.mean() - row.T0.mean()) / max(row.U.mean() ** 2, 1e-9)
        if is_compressor:
            row.stage_loading *= -1 # Stage_loading will be negative 
        row.euler_power = mdot * (ref.U * ref.Vt - row.U * row.Vt).mean()
        row.flow_coefficient = abs(float(np.mean(row.Vm / row.U)))
