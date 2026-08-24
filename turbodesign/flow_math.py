from typing import Dict, Optional, Sequence, Tuple

import numpy as np
import numpy.typing as npt

from . import isentropic
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
            dr = row.r[j] - row.r[j - 1]
            dl = np.sqrt(dx**2 + dr**2)
            # Signed area: sign follows dx to maintain massflow sign convention
            sign = -1.0 if dx < 0 else 1.0
            streamline_area[j] = sign * np.pi * (row.r[j] + row.r[j - 1]) * dl
            total_area += streamline_area[j]
    return total_area, streamline_area


def radii_for_area(area: float, mean_radius: float) -> Tuple[float, float]:
    """Hub/shroud radii for an annulus of the given area at a fixed mean radius.

    The thin-annulus complement to `isentropic.area_for_massflow`: MFP-based
    sizing gives you a scalar area, not a passage - this is the last step
    that turns "how much area" into "where the walls are," e.g. for sizing a
    component's inlet before any blade row has been solved.

    Args:
        area: Annulus cross-sectional area [m^2].
        mean_radius: Radius to hold fixed while the annulus grows/shrinks [m].

    Returns:
        Tuple of (r_hub, r_shroud) [m].
    """
    half_height = area / (4.0 * np.pi * mean_radius)
    return mean_radius - half_height, mean_radius + half_height


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


def reset_rotor_power(rows: Sequence[BladeRow]) -> None:
    """Zero `power`/`power_mean` on every rotor row.

    `initialize()` seeds its next T0 guess from the *previous* solve's
    `row.power`. Call this before re-solving at a new operating point, or a
    stale power from a very different massflow/rpm can slow convergence.
    """
    for row in rows:
        if row.row_type == RowType.Rotor:
            row.power = 0.0
            row.power_mean = 0.0


def update_choke_diagnostics(row: BladeRow) -> None:
    """Populate row.mass_flow_function / row.choke_margin / row.choke_margin_min.

    Uses M_rel for rotors, M for stators/IGVs. No-op if Mach isn't solved yet.
    """
    M = row.M_rel if row.row_type == RowType.Rotor else row.M
    M_arr = np.atleast_1d(np.asarray(M, dtype=float))
    if not np.any(M_arr):
        return
    gamma = float(np.mean(row.gamma)) if np.size(row.gamma) else 1.4
    m_tilde = np.atleast_1d(np.asarray(isentropic.mass_flow_function(M_arr, gamma), dtype=float))
    m_tilde_max = isentropic.mass_flow_function_max(gamma)
    margin = 1.0 - m_tilde / m_tilde_max
    row.mass_flow_function = m_tilde
    row.choke_margin = margin
    row.choke_margin_min = float(np.min(margin))


def row_flow_capacity(row: BladeRow, massflow: Optional[float] = None) -> Dict[str, float]:
    """Compute the non-dimensional flow-capacity check for a single row.

    Args:
        row: BladeRow with total_area, blockage, gamma, R, and P0/T0 (or
            P0R/T0R for rotors) already populated.
        massflow: Target massflow [kg/s]. Defaults to row.total_massflow.

    Returns:
        Dict with m_tilde_required, m_tilde_max, min_area, area, massflow,
        P0, T0, gamma, R, blockage, frame, and feasible (bool).
    """
    mdot = float(massflow) if massflow is not None else float(np.mean(row.total_massflow))
    gamma = float(np.mean(row.gamma)) if np.size(row.gamma) else 1.4
    R = float(np.mean(row.R)) if np.size(row.R) else 287.0
    blockage = float(row.blockage)
    area = float(row.total_area)

    frame = "absolute"
    P0 = float(np.mean(row.P0))
    T0 = float(np.mean(row.T0))
    if row.row_type == RowType.Rotor:
        P0R = float(np.mean(row.P0R))
        T0R = float(np.mean(row.T0R))
        if P0R > 0 and T0R > 0:
            P0, T0, frame = P0R, T0R, "relative"
        else:
            frame = "absolute (relative unavailable)"

    m_tilde_required = float(isentropic.mass_flow_function_required(mdot, P0, T0, area, gamma, R, blockage))
    m_tilde_max = isentropic.mass_flow_function_max(gamma)
    min_area = float(isentropic.min_area_for_massflow(mdot, P0, T0, gamma, R, blockage))

    return {
        "m_tilde_required": m_tilde_required,
        "m_tilde_max": m_tilde_max,
        "min_area": min_area,
        "area": area,
        "massflow": mdot,
        "P0": P0,
        "T0": T0,
        "gamma": gamma,
        "R": R,
        "blockage": blockage,
        "frame": frame,
        "feasible": m_tilde_required <= m_tilde_max,
    }


def explain_infeasible_massflow(row: BladeRow, target_massflow: float, P0: npt.NDArray, T0: npt.NDArray, area: Optional[float] = None) -> Optional[str]:
    """Return an actionable message if `row` cannot pass `target_massflow` at
    (P0, T0) and its (given or current) area/blockage, else None.

    Translates a bounded Mach-solve `RuntimeError` into a feasibility
    diagnosis at the point of failure, before `assert_flow_capacity`'s own
    upfront check gets a chance to run.

    Args:
        row: The row whose massflow capacity is in question.
        target_massflow: Massflow the row is being asked to pass [kg/s].
        P0: Total pressure to evaluate capacity at [Pa].
        T0: Total temperature to evaluate capacity at [K].
        area: Annulus area override [m^2]. Defaults to `row.total_area`
            (which may still be stale/zero at this point in the solve).

    Returns:
        A multi-line diagnostic string if infeasible, else None.
    """
    if area is None:
        area = float(row.total_area) if np.size(row.total_area) else 0.0
    if area <= 0 or target_massflow <= 0:
        return None
    gamma = float(np.mean(row.gamma)) if np.size(row.gamma) else 1.4
    R = float(np.mean(row.R)) if np.size(row.R) else 287.0
    blockage = float(row.blockage)
    P0_mean = float(np.mean(P0))
    T0_mean = float(np.mean(T0))
    if P0_mean <= 0 or T0_mean <= 0:
        return None

    m_tilde_required = float(isentropic.mass_flow_function_required(target_massflow, P0_mean, T0_mean, area, gamma, R, blockage))
    m_tilde_max = isentropic.mass_flow_function_max(gamma)
    if m_tilde_required <= m_tilde_max:
        return None

    min_area = float(isentropic.min_area_for_massflow(target_massflow, P0_mean, T0_mean, gamma, R, blockage))
    area_eff = area * (1.0 - blockage)
    ratio = min_area / area_eff if area_eff > 0 else float("inf")
    return (
        f"blade row {row.id} ({row.row_type.name}, stage {row.stage_id}) cannot pass the "
        f"requested massflow. Required flow function m~ = {m_tilde_required:.4f} exceeds the "
        f"choked limit m~_max = {m_tilde_max:.4f} for gamma = {gamma:.3f} - no Mach number can "
        f"satisfy this.\n"
        f"  target massflow : {target_massflow:.4f} kg/s\n"
        f"  P0 = {P0_mean:.1f} Pa, T0 = {T0_mean:.1f} K, gamma = {gamma:.3f}, R = {R:.1f} J/(kg K)\n"
        f"  annulus area    : {area:.6f} m^2 (blockage {blockage:.3f}) -> effective {area_eff:.6f} m^2\n"
        f"  minimum area    : {min_area:.6f} m^2  ({ratio:.2f}x larger needed)\n"
        "Fix by: enlarging the annulus at this station, lowering the target massflow, "
        "raising inlet P0, or lowering inlet T0."
    )


def assert_flow_capacity(rows: Sequence[BladeRow], context: str, targets: Optional[Sequence[Optional[float]]] = None) -> None:
    """Raise a ValueError naming every row whose target massflow exceeds its
    sonic (choked) flow capacity at the row's current area/Pt/Tt.

    This is a necessary-only feasibility check: no exit angle can beat
    m~_max, so an infeasible target here can never converge, angle-matching
    mode included. Passing this check does not guarantee convergence - it
    ignores blade throat area, swirl, and spanwise non-uniformity.

    Args:
        rows: Blade rows to check (non Rotor/Stator/IGV rows are skipped).
        context: Short label identifying the caller, used in the error message.
        targets: Optional per-row target massflow overrides, same length as
            `rows`. A `None` entry (or omitting `targets`) uses the row's own
            `total_massflow`.

    Raises:
        ValueError: If one or more rows cannot pass their target massflow.
    """
    failures = []
    for i, row in enumerate(rows):
        if row.row_type not in (RowType.Rotor, RowType.Stator, RowType.IGV):
            continue
        target = targets[i] if targets is not None else None
        cap = row_flow_capacity(row, target)
        if not cap["feasible"]:
            area_eff = cap["area"] * (1.0 - cap["blockage"])
            ratio = cap["min_area"] / area_eff if area_eff > 0 else float("inf")
            failures.append(
                f"  blade row {row.id} ({row.row_type.name}, stage {row.stage_id}): "
                f"required m~ = {cap['m_tilde_required']:.4f} exceeds choked limit "
                f"m~_max = {cap['m_tilde_max']:.4f} for gamma = {cap['gamma']:.3f} "
                f"- no Mach number can satisfy this.\n"
                f"    target massflow : {cap['massflow']:.4f} kg/s  (frame: {cap['frame']})\n"
                f"    P0 = {cap['P0']:.1f} Pa, T0 = {cap['T0']:.1f} K, gamma = {cap['gamma']:.3f}, R = {cap['R']:.1f} J/(kg K)\n"
                f"    annulus area    : {cap['area']:.6f} m^2 (blockage {cap['blockage']:.3f}) -> effective {area_eff:.6f} m^2\n"
                f"    minimum area    : {cap['min_area']:.6f} m^2  ({ratio:.2f}x larger needed)"
            )
    if failures:
        raise ValueError(
            f"{context}: the following blade rows cannot pass their requested massflow "
            "(target massflow exceeds sonic/choked flow capacity):\n" + "\n".join(failures) +
            "\nFix by: enlarging the annulus at these stations, lowering the target massflow, "
            "raising inlet P0, or lowering inlet T0."
        )
