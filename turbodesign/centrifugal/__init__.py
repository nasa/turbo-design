"""Centrifugal compressor meanline module -- rothalpy-consistent, radial-safe.

A new module alongside upstream ``turbodesign``, not a patch to it
(``docs/centrifugal/01-design.md`` S1): the axial-oriented core
(``compressor_math.py``, ``compressor_spool.py``, ``radeq.py``,
``flow_math.compute_streamline_areas``, ``deviation/``) is never imported here, so the
axial turbine cases it serves cannot regress.

See ``docs/centrifugal/00-root-cause-analysis.md`` for the bugs this module exists to
avoid, and ``docs/centrifugal/02-tdd-plan.md`` for the slice plan. This is slice 1: an
isentropic impeller only (no slip, no losses, no diffuser).
"""

from .components import Impeller
from .geometry import MeridionalPath, Station
from .losses import (
    EvaluatedLosses,
    ImpellerBladeLoadingCoppage,
    ImpellerChokeAungier,
    ImpellerClearanceJansen,
    ImpellerDiscFrictionDaily,
    ImpellerEntranceDiffusionAungier,
    ImpellerIncidenceConrad,
    ImpellerLeakageAungier,
    ImpellerLossState,
    ImpellerMixingAungier,
    ImpellerRecirculationOh,
    ImpellerSkinFrictionJansen,
    LossModel,
    LossSet,
    NO_LOSSES,
    OhLossSet,
)
from .slip import BusemannSlip, PhysicsError, QiuSlip, SlipModel, StanitzSlip, WiesnerSlip
from .solver import Choked, InletState, OperatingPoint, Stage, StationSeries
from .state import Air, ThermoState

__all__ = [
    "Air",
    "BusemannSlip",
    "Choked",
    "EvaluatedLosses",
    "Impeller",
    "ImpellerBladeLoadingCoppage",
    "ImpellerChokeAungier",
    "ImpellerClearanceJansen",
    "ImpellerDiscFrictionDaily",
    "ImpellerEntranceDiffusionAungier",
    "ImpellerIncidenceConrad",
    "ImpellerLeakageAungier",
    "ImpellerLossState",
    "ImpellerMixingAungier",
    "ImpellerRecirculationOh",
    "ImpellerSkinFrictionJansen",
    "InletState",
    "LossModel",
    "LossSet",
    "MeridionalPath",
    "NO_LOSSES",
    "OhLossSet",
    "OperatingPoint",
    "PhysicsError",
    "QiuSlip",
    "SlipModel",
    "Stage",
    "StanitzSlip",
    "Station",
    "StationSeries",
    "ThermoState",
    "WiesnerSlip",
]
