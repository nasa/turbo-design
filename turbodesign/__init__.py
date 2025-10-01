from .spool import Spool
from .stage import Stage
from .enums import LossType, RowType, PassageType, MassflowConstraint
from .inlet import Inlet
from .bladerow import BladeRow, Coolant
from .lossinterp import LossInterp
from .passage import Passage
from .outlet import Outlet

# turbodesign/__init__.py
from importlib import import_module

__all__ = [
    "Spool", "Stage", "Inlet", "Outlet", "BladeRow", "Coolant",
    "Passage", "RowType", "PassageType", "MassflowConstraint", "LossInterp",
]

_module_map = {
    "Spool": ("turbodesign.spool", "Spool"),
    "Stage": ("turbodesign.stage", "Stage"),
    "Inlet": ("turbodesign.inlet", "Inlet"),
    "Outlet": ("turbodesign.outlet", "Outlet"),
    "BladeRow": ("turbodesign.bladerow", "BladeRow"),
    "Coolant": ("turbodesign.bladerow", "Coolant"),
    "Passage": ("turbodesign.passage", "Passage"),
    "RowType": ("turbodesign.enums", "RowType"),
    "PassageType": ("turbodesign.enums", "PassageType"),
    "MassflowConstraint": ("turbodesign.enums", "MassflowConstraint"),
    "LossInterp": ("turbodesign.lossinterp", "LossInterp"),
}

def __getattr__(name: str):
    try:
        mod_name, attr = _module_map[name]
    except KeyError:
        raise AttributeError(name)
    mod = import_module(mod_name)
    return getattr(mod, attr)
