"""Sweep an already-designed, fixed-geometry compressor or turbine across many
operating points (massflow x rpm) to build a performance map/characteristic.

Mutates a persistent `CompressorSpool`/`TurbineSpool` in place rather than
rebuilding per point. Safe because `spool.massflow` and the outlet target are
re-read fresh every `solve()`; `rpm` is not (use `spool.set_rpm`), and stale
`row.power` can degrade convergence if not reset (`flow_math.reset_rotor_power`).

Has no knowledge of loss models or incidence: the massflow/choke axis is
honest immediately (pure continuity); the efficiency/pressure-ratio axis is
only as honest as the loss model plugged into the rows.

Caveat: for `TurbineSpool` in pressure-balance mode, `spool.massflow` is only
a solver seed, not an enforced constraint - the achieved massflow is an
output of the fixed geometry and boundary pressures, so `massflow_points`
often doesn't move it at all (same reason `ShaftMatch.match()` needed two
knobs). `rpm_points` sweeps meaningfully for either machine; `massflow_points`
alone is only meaningful for a `CompressorSpool`.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence, Union

import numpy as np

from .enums import RowType
from .flow_math import reset_rotor_power

__all__ = ["OperatingPoint", "sweep_operating_points"]


@dataclass
class OperatingPoint:
    """One evaluated point of a swept compressor/turbine map."""

    massflow: float
    rpm: float
    corrected_massflow: float   # mdot*sqrt(T0_inlet)/P0_inlet - standard map x-axis
    corrected_speed: float      # rpm*pi/30/sqrt(T0_inlet) - standard speed-line label
    pressure_ratio: float
    efficiency: float
    power: float
    choke_margin_min: float     # min over all rows; <= 0 means something is choked
    converged: bool
    error: Optional[str] = None


def _failed_point(massflow: float, rpm: float, exc: Exception) -> OperatingPoint:
    return OperatingPoint(
        massflow=massflow, rpm=rpm,
        corrected_massflow=float("nan"), corrected_speed=float("nan"),
        pressure_ratio=float("nan"), efficiency=float("nan"), power=float("nan"),
        choke_margin_min=float("nan"),
        converged=False, error=str(exc),
    )


def sweep_operating_points(
    spool,
    massflow_points: Sequence[float],
    rpm_points: Sequence[float],
) -> List[OperatingPoint]:
    """Sweep a fixed-geometry spool across every (rpm, massflow) pair.

    Only `massflow` and `rpm` vary; geometry and loss models stay as built. A
    point that fails to converge (`ValueError`/`RuntimeError`) is recorded
    with `converged=False` and the sweep continues.

    Args:
        spool: An already-constructed `CompressorSpool` or `TurbineSpool`,
            mutated in place. See the module docstring's caveat for
            `TurbineSpool`.
        massflow_points: Massflow values [kg/s] to evaluate at each rpm.
        rpm_points: Shaft speeds [rev/min] to sweep (each becomes one speed
            line on the resulting map).

    Returns:
        One `OperatingPoint` per (rpm, massflow) pair, in nested order:
        all massflow points for the first rpm, then the next rpm, etc.
    """
    results: List[OperatingPoint] = []
    for rpm in rpm_points:
        spool.set_rpm(float(rpm))
        for massflow in massflow_points:
            massflow = float(massflow)
            reset_rotor_power(spool.rows)
            spool.massflow = massflow
            try:
                spool.solve()
            except (ValueError, RuntimeError) as exc:
                results.append(_failed_point(massflow, rpm, exc))
                continue

            inlet = spool.inlet
            T0_inlet = float(np.mean(inlet.T0))
            P0_inlet = float(np.mean(inlet.P0))
            corrected_massflow = massflow * np.sqrt(T0_inlet) / P0_inlet
            corrected_speed = rpm * np.pi / 30.0 / np.sqrt(T0_inlet)

            rows = spool._all_rows()
            margins = [
                float(row.choke_margin_min)
                for row in rows
                if row.row_type in (RowType.Rotor, RowType.Stator, RowType.IGV)
            ]
            choke_margin_min = float(np.min(margins)) if margins else float("nan")

            results.append(OperatingPoint(
                massflow=massflow, rpm=float(rpm),
                corrected_massflow=float(corrected_massflow),
                corrected_speed=float(corrected_speed),
                pressure_ratio=spool.overall_pressure_ratio(),
                efficiency=spool.overall_entropy_efficiency(),
                power=spool.total_power(),
                choke_margin_min=choke_margin_min,
                converged=True,
            ))
    return results
