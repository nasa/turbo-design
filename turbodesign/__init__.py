from .turbine_spool import TurbineSpool
from .stage import Stage
from .enums import LossType, RowType, PassageType, MassflowConstraint
from .inlet import Inlet
from .bladerow import BladeRow
from .coolant import Coolant
from .lossinterp import LossInterp
from .passage import Passage
from .outlet import Outlet
from .deviation import DeviationBaseClass, FixedDeviation

# turbodesign/__init__.py
from importlib import import_module

__all__ = [
    "TurbineSpool", "Stage", "Inlet", "Outlet", "BladeRow", "Coolant",
    "Passage", "RowType", "PassageType", "MassflowConstraint", "LossType",
    "LossInterp", "DeviationBaseClass", "FixedDeviation",
]

_module_map = {
    "TurbineSpool": ("turbodesign.turbine_spool", "TurbineSpool"),
    "Stage": ("turbodesign.stage", "Stage"),
    "Inlet": ("turbodesign.inlet", "Inlet"),
    "Outlet": ("turbodesign.outlet", "Outlet"),
    "BladeRow": ("turbodesign.bladerow", "BladeRow"),
    "Coolant": ("turbodesign.coolant", "Coolant"),
    "Passage": ("turbodesign.passage", "Passage"),
    "LossType": ("turbodesign.enums", "LossType"),
    "RowType": ("turbodesign.enums", "RowType"),
    "PassageType": ("turbodesign.enums", "PassageType"),
    "MassflowConstraint": ("turbodesign.enums", "MassflowConstraint"),
    "LossInterp": ("turbodesign.lossinterp", "LossInterp"),
    "DeviationBaseClass": ("turbodesign.deviation", "DeviationBaseClass"),
    "FixedDeviation": ("turbodesign.deviation", "FixedDeviation"),
}

def __getattr__(name: str):
    try:
        mod_name, attr = _module_map[name]
    except KeyError:
        raise AttributeError(name)
    mod = import_module(mod_name)
    return getattr(mod, attr)
