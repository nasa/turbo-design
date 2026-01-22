"""Convenience helpers for constructing BladeRow objects with minimal boilerplate."""

from __future__ import annotations

from typing import Optional, Sequence, Any

import numpy as np

from .bladerow import BladeRow
from .enums import RowType

__all__ = ["make_blade_row", "make_rotor_row", "make_stator_row"]


def make_blade_row(
    *args: Any,
    row_type: Optional[RowType] = None,
    hub_location: Optional[float] = None,
    shroud_location: Optional[float] = None,
    stage_id: int = 0,
    **attrs: Any,
) -> BladeRow:
    """Generic BladeRow factory that supports legacy calling styles.

    Args mirror the original BladeRow constructor; any extra keyword arguments
    are set as attributes on the created row for convenience.
    """
    # Allow positional legacy (hub_location, row_type, stage_id, shroud_location)
    pos = list(args)
    # Legacy positional usage: (RowType, power=...) or (hub_loc, RowType, ...)
    if pos and isinstance(pos[0], RowType):
        row_type = pos.pop(0)
    if hub_location is None and pos:
        hub_location = pos.pop(0)
    if row_type is None and pos and isinstance(pos[0], RowType):
        row_type = pos.pop(0)
    if shroud_location is None and pos:
        shroud_location = pos.pop(0)

    hub_location = 0.0 if hub_location is None else float(hub_location)
    row_type = row_type if row_type is not None else RowType.Stator

    row = BladeRow(hub_location=hub_location, row_type=row_type, stage_id=stage_id, shroud_location=shroud_location)

    # Apply extra attributes (e.g., loss_function, beta2_metal, etc.)
    for key, val in attrs.items():
        setattr(row, key, val)
    return row

def _maybe_set_pitch(row: BladeRow, pitch_to_chord: Optional[float], solidity: Optional[float]) -> None:
    """Apply pitch/solidity inputs to a blade row if provided."""
    if pitch_to_chord is not None:
        row.pitch_to_chord = pitch_to_chord
    elif solidity is not None and solidity != 0:
        row.pitch_to_chord = 1.0 / solidity


def make_rotor_row(
    hub_location: float,
    metal_exit_angle_deg: Optional[float | Sequence[float]] = None,
    loss_function: Optional[object] = None,
    P0_ratio: float = 1.0,
    pitch_to_chord: Optional[float] = None,
    solidity: Optional[float] = None,
    num_blades: Optional[int] = None,
    axial_chord: Optional[float] = None,
) -> BladeRow:
    """Create a Rotor blade row with common inputs.

    Args:
        hub_location: Streamwise position (0–1) for the row.
        beta2_metal_deg: Exit metal angle(s) in degrees.
        loss_function: Loss model to attach.
        P0_ratio: Total-pressure ratio target across the row.
        pitch_to_chord: Pitch-to-chord ratio (alternative to solidity).
        solidity: Solidity (chord/pitch) if preferred over pitch_to_chord.
        num_blades: Number of blades (used to derive pitch).
        axial_chord: Axial chord length if known.
    """
    row = BladeRow(hub_location=hub_location, row_type=RowType.Rotor)
    row.P0_ratio = P0_ratio
    row.P0_ratio_target = P0_ratio
    _maybe_set_pitch(row, pitch_to_chord, solidity)
    if num_blades is not None:
        row.num_blades = num_blades
    if axial_chord is not None:
        row.axial_chord = axial_chord
    if metal_exit_angle_deg is not None:
        row.metal_exit_angle = np.atleast_1d(metal_exit_angle_deg)
    if loss_function is not None:
        row.loss_function = loss_function
    return row


def make_stator_row(
    hub_location: float,
    metal_exit_angle_deg: Optional[float | Sequence[float]] = None,
    loss_function: Optional[object] = None,
    P0_ratio: float = 1.0,
    pitch_to_chord: Optional[float] = None,
    solidity: Optional[float] = None,
    num_blades: Optional[int] = None,
    axial_chord: Optional[float] = None,
) -> BladeRow:
    """Create a Stator/IGV blade row with common inputs.

    Args:
        hub_location: Streamwise position (0–1) for the row.
        alpha2_metal_deg: Exit metal angle(s) in degrees.
        loss_function: Loss model to attach.
        P0_ratio: Total-pressure ratio target across the row.
        pitch_to_chord: Pitch-to-chord ratio (alternative to solidity).
        solidity: Solidity (chord/pitch) if preferred over pitch_to_chord.
        num_blades: Number of blades (used to derive pitch).
        axial_chord: Axial chord length if known.
    """
    row = BladeRow(hub_location=hub_location, row_type=RowType.Stator)
    row.P0_ratio = P0_ratio
    row.P0_ratio_target = P0_ratio
    _maybe_set_pitch(row, pitch_to_chord, solidity)
    if num_blades is not None:
        row.num_blades = num_blades
    if axial_chord is not None:
        row.axial_chord = axial_chord
    if metal_exit_angle_deg is not None:
        row.metal_exit_angle = np.atleast_1d(metal_exit_angle_deg)
    if loss_function is not None:
        row.loss_function = loss_function
    return row
