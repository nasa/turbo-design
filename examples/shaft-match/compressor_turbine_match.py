"""Size and match a turbine to an already-designed compressor on a common shaft.

This is the first example in the repo instantiating both a CompressorSpool and
a TurbineSpool: it takes the Mattingly example-9.1 compressor stage, treats it
as the "known" side of a gas generator, and uses `ShaftMatch` to:

  1. Run a cheap 0-D/1-D sizing pass (shaft power balance + mass-flow
     continuity + the non-dimensional mass flow function) to size a matching
     turbine's annulus.
  2. Scale a hand-built single-stage reference turbine's geometry to hit that
     sizing target, then solve for the (annulus area scale, exit pressure)
     pair that closes both the mass-flow and shaft-power residuals.

The reference turbine's blade angles/loss coefficients are reused as-is
(`ShaftMatch`'s `ScaledReferenceBuilder` approach) - only its annulus radii
and exit static pressure are adjusted by the match. See turbodesign/shaft_match.py
for the two-phase design and its limitations.
"""

import numpy as np

from turbodesign import Inlet, Outlet, Passage, PassageType, TurbineSpool
from turbodesign.compressor_spool import CompressorSpool
from turbodesign.flow_math import radii_for_area
from turbodesign.isentropic import area_for_massflow
from turbodesign.loss.fixedpressureloss import FixedPressureLoss
from turbodesign.row_factory import make_rotor_row, make_stator_row
from turbodesign.shaft_match import ComponentSizingSpec, ShaftMatch


def build_compressor(rpm: float) -> CompressorSpool:
    """Mattingly example 9.1: a single-stage axial compressor stage.

    Station 1 (the compressor face) is sized from the non-dimensional mass
    flow function rather than taken from the textbook table - this is
    "Station 2" in the gas-generator sense (MFP sizes the inlet boundary
    before any blade row has been solved). Stations 2 and 3 (rotor/stator
    exit) stay at Mattingly's prescribed values, since those come from the
    stage's velocity-triangle/loading design, not a freely chosen Mach.
    """
    T01_K = 518.7 / 1.8
    P01_Pa = 14.7 * 6894.76
    massflow_kg_s = 50 * 0.453592
    M1 = 0.7
    rmean = 12 * 0.0254

    area1 = area_for_massflow(massflow_kg_s, P01_Pa, T01_K, M1, gamma=1.4, R=287.0)
    area2 = 179.1 / 39.3701**2
    area3 = 165.3 / 39.3701**2
    rhub1, rshroud1 = radii_for_area(area1, rmean)
    h2 = area2 / (np.pi * 4 * rmean)
    h3 = area3 / (np.pi * 4 * rmean)
    cax = 1 * 0.0254

    xhub_arr = [0, cax, 2 * cax]
    xshroud_arr = [0, cax, 2 * cax]
    rhub_arr = [rhub1, rmean - h2, rmean - h3]
    rshroud_arr = [rshroud1, rmean + h2, rmean + h3]

    passage = Passage(xhub_arr, rhub_arr, xshroud_arr, rshroud_arr, passageType=PassageType.Axial, zero_phi=True)
    inlet = Inlet(hub_location=0)
    inlet.alpha2 = [40]
    inlet.init_total(P01_Pa, T01_K, M=M1)

    rotor = make_rotor_row(hub_location=cax / max(xhub_arr), metal_exit_angle_deg=[-23.87], loss_function=FixedPressureLoss(0))
    stator = make_stator_row(hub_location=2 * cax / max(xhub_arr), metal_exit_angle_deg=[40], loss_function=FixedPressureLoss(0), P0_ratio=1.3)

    outlet = Outlet()
    outlet.init_total(1.3 * P01_Pa, 0.5)

    return CompressorSpool(passage, massflow_kg_s, inlet, outlet, [rotor, stator], num_streamlines=1, rpm=rpm)


def build_reference_turbine(rpm: float) -> TurbineSpool:
    """A hand-picked single-stage turbine whose annulus/exit-pressure ShaftMatch
    will scale to match the compressor - its angles and loss coefficients are
    reused verbatim by the ScaledReferenceBuilder approach."""
    T0, P0, P_exit, rmean = 1200.0, 126680.0, 90000.0, 0.15
    area1, area2, area3 = 0.05, 0.045, 0.05
    h1 = area1 / (np.pi * 4 * rmean)
    h2 = area2 / (np.pi * 4 * rmean)
    h3 = area3 / (np.pi * 4 * rmean)
    cax = 0.02

    xhub_arr = [0, cax, 2 * cax]
    xshroud_arr = [0, cax, 2 * cax]
    rhub_arr = [rmean - h1, rmean - h2, rmean - h3]
    rshroud_arr = [rmean + h1, rmean + h2, rmean + h3]
    passage = Passage(xhub_arr, rhub_arr, xshroud_arr, rshroud_arr, passageType=PassageType.Axial, zero_phi=True)

    inlet = Inlet(hub_location=0, alpha=[0])
    inlet.init_total(P0=P0, T0=T0, M=0.15)
    outlet = Outlet(num_streamlines=1)
    outlet.init_static(P=P_exit, percent_radii=[0.5])

    stator = make_stator_row(hub_location=cax / max(xhub_arr), metal_exit_angle_deg=[65.0], loss_function=FixedPressureLoss(0.03))
    rotor = make_rotor_row(hub_location=2 * cax / max(xhub_arr), metal_exit_angle_deg=[-55.0], loss_function=FixedPressureLoss(0.05))

    turbine = TurbineSpool(passage, 5.0, inlet, outlet, [stator, rotor], num_streamlines=1, rpm=rpm)
    turbine.adjust_streamlines = False
    return turbine


def main() -> None:
    rpm = 1000 * 30 / np.pi  # omega = 1000 rad/s

    compressor = build_compressor(rpm)
    compressor.solve()

    inlet_row = compressor._all_rows()[0]
    print("Station 1 (compressor inlet), sized from the non-dimensional mass flow function:")
    print(f"  area = {float(inlet_row.total_area):.6f} m^2  (at M = {float(np.mean(inlet_row.M)):.2f})")
    print()

    reference_turbine = build_reference_turbine(rpm)

    spec = ComponentSizingSpec(
        inlet_T0=1200.0,          # combustor exit total temperature [K]
        eta_total_guess=0.88,
        exit_static_pressure=101325.0,
        num_stages=1,
        # ngv_choked defaults to True: the NGV throat (station 4) is sized at M=1,
        # the classic choked-metering-orifice assumption.
    )
    shaft = ShaftMatch(compressor, spec, reference=reference_turbine, fuel_air_ratio=0.0, bleed_fraction=0.0)

    sizing = shaft.size()
    print("Station 4 (NGV throat), sizing pass (0-D/1-D, no solver):")
    print(f"  NGV throat Mach          : {sizing.station_M[0]:.2f}  (choked = {spec.ngv_choked})")
    print(f"  target turbine massflow  : {sizing.massflow:.4f} kg/s")
    print(f"  required shaft power     : {sizing.required_power / 1000:.2f} kW")
    print(f"  inlet P0/T0              : {sizing.inlet_P0 / 1000:.1f} kPa / {sizing.inlet_T0:.1f} K")
    print(f"  exit P0/T0/P             : {sizing.exit_P0 / 1000:.1f} kPa / {sizing.exit_T0:.1f} K / {sizing.exit_P / 1000:.1f} kPa")
    print(f"  NGV throat area          : {sizing.station_area[0]:.6f} m^2")
    for w in sizing.warnings:
        print(f"  warning: {w}")

    result = shaft.match(tol_rel=5e-3, max_iter=40)

    print()
    print("Shaft match result:")
    print(f"  converged           : {result.converged} ({result.iterations} iterations)")
    print(f"  compressor massflow : {result.massflow_compressor:.4f} kg/s, power: {result.power_compressor / 1000:.2f} kW")
    print(f"  turbine massflow    : {result.massflow_turbine:.4f} kg/s, power: {result.power_turbine / 1000:.2f} kW")
    print(f"  massflow residual   : {result.massflow_residual:+.4%}")
    print(f"  power residual      : {result.power_residual:+.4%}")
    print(f"  annulus area scale  : {result.area_scale:.3f}x reference")

    shaft.export("compressor_turbine_match", result)
    print()
    print("Wrote compressor_turbine_match_compressor.json, "
          "compressor_turbine_match_turbine.json, compressor_turbine_match_shaft_match.json")


if __name__ == "__main__":
    main()
