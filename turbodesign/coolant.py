from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

from cantera import Solution


@dataclass
class Coolant:
    """Simple container for coolant properties used in spool calculations."""

    fluid: Optional[Solution] = None
    T0: float = 300.0
    P0: float = 101325.0
    massflow_percentage: float = 0.0
    Cp: Optional[float] = None

    def __post_init__(self) -> None:
        # Default Cp to the provided fluid value when available.
        if self.Cp is None:
            self.Cp = float(self.fluid.cp) if self.fluid is not None else 0.0
