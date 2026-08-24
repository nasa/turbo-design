from .turbine_spool import TurbineSpool
from .stage import Stage
from .enums import LossType, RowType, PassageType
from .inlet import Inlet
from .bladerow import BladeRow
from .coolant import Coolant
from .lossinterp import LossInterp
from .passage import Passage
from .outlet import Outlet
from .deviation import DeviationBaseClass, FixedDeviation
from .row_factory import make_blade_row, make_rotor_row, make_stator_row
from .agf import AGF_Setup, Inlet_bcs, Outlet_bcs, Settings, Clearance, Domain, read_agf, plot_airfoil_inputs, plot_airfoil_inputs_2D

# turbodesign/__init__.py
from importlib import import_module

__all__ = [
    "TurbineSpool", "CompressorSpool", "Stage", "Inlet", "Outlet", "BladeRow", "Coolant",
    "Passage", "RowType", "PassageType", "LossType",
    "LossInterp", "DeviationBaseClass", "FixedDeviation",
    "make_blade_row", "make_rotor_row", "make_stator_row",
    "FixedPolytropicEfficiency",
    "FixedPressureLoss",
    "AGF_Setup", "Inlet_bcs", "Outlet_bcs", "Settings", "Clearance", "Domain", "read_agf",
    "plot_airfoil_inputs", "plot_airfoil_inputs_2D",
    "ShaftMatch", "ComponentSizingSpec", "SizingEstimate", "MatchResult", "turbine_massflow",
    "OperatingPoint", "sweep_operating_points",
]

_module_map = {
    "TurbineSpool": ("turbodesign.turbine_spool", "TurbineSpool"),
    "CompressorSpool": ("turbodesign.compressor_spool", "CompressorSpool"),
    "ShaftMatch": ("turbodesign.shaft_match", "ShaftMatch"),
    "ComponentSizingSpec": ("turbodesign.shaft_match", "ComponentSizingSpec"),
    "SizingEstimate": ("turbodesign.shaft_match", "SizingEstimate"),
    "MatchResult": ("turbodesign.shaft_match", "MatchResult"),
    "turbine_massflow": ("turbodesign.shaft_match", "turbine_massflow"),
    "OperatingPoint": ("turbodesign.operating_map", "OperatingPoint"),
    "sweep_operating_points": ("turbodesign.operating_map", "sweep_operating_points"),
    "Stage": ("turbodesign.stage", "Stage"),
    "Inlet": ("turbodesign.inlet", "Inlet"),
    "Outlet": ("turbodesign.outlet", "Outlet"),
    "BladeRow": ("turbodesign.bladerow", "BladeRow"),
    "Coolant": ("turbodesign.coolant", "Coolant"),
    "Passage": ("turbodesign.passage", "Passage"),
    "LossType": ("turbodesign.enums", "LossType"),
    "RowType": ("turbodesign.enums", "RowType"),
    "PassageType": ("turbodesign.enums", "PassageType"),
    "LossInterp": ("turbodesign.lossinterp", "LossInterp"),
    "FixedPolytropicEfficiency": ("turbodesign.loss.fixedpolytropic", "FixedPolytropicEfficiency"),
    "FixedPressureLoss": ("turbodesign.loss.fixedpressureloss", "FixedPressureLoss"),
    "DeviationBaseClass": ("turbodesign.deviation", "DeviationBaseClass"),
    "FixedDeviation": ("turbodesign.deviation", "FixedDeviation"),
    "make_blade_row": ("turbodesign.row_factory", "make_blade_row"),
    "make_rotor_row": ("turbodesign.row_factory", "make_rotor_row"),
    "make_stator_row": ("turbodesign.row_factory", "make_stator_row"),
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
