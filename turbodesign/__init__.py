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
from .agf import AGF_Setup, Inlet_bcs, Outlet_bcs, Settings, Clearance, Domain, read_agf, plot_airfoil_inputs, plot_airfoil_inputs_2D

# turbodesign/__init__.py
from importlib import import_module

__all__ = [
    "TurbineSpool", "Stage", "Inlet", "Outlet", "BladeRow", "Coolant",
    "Passage", "RowType", "PassageType", "MassflowConstraint", "LossType",
    "LossInterp", "DeviationBaseClass", "FixedDeviation",
    "FixedPolytropicEfficiency",
    "AGF_Setup", "Inlet_bcs", "Outlet_bcs", "Settings", "Clearance", "Domain", "read_agf",
    "plot_airfoil_inputs", "plot_airfoil_inputs_2D",
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
    "FixedPolytropicEfficiency": ("turbodesign.loss.fixedpolytropic", "FixedPolytropicEfficiency"),
    "DeviationBaseClass": ("turbodesign.deviation", "DeviationBaseClass"),
    "FixedDeviation": ("turbodesign.deviation", "FixedDeviation"),
    "AGF_Setup": ("turbodesign.agf", "AGF_Setup"),
    "Inlet_bcs": ("turbodesign.agf", "Inlet_bcs"),
    "Outlet_bcs": ("turbodesign.agf", "Outlet_bcs"),
    "Settings": ("turbodesign.agf", "Settings"),
    "Clearance": ("turbodesign.agf", "Clearance"),
    "Domain": ("turbodesign.agf", "Domain"),
    "read_agf": ("turbodesign.agf", "read_agf"),
    "plot_airfoil_inputs": ("turbodesign.agf", "plot_airfoil_inputs"),
    "plot_airfoil_inputs_2D": ("turbodesign.agf", "plot_airfoil_inputs_2D"),
}

def __getattr__(name: str):
    try:
        mod_name, attr = _module_map[name]
    except KeyError:
        raise AttributeError(name)
    mod = import_module(mod_name)
    return getattr(mod, attr)
