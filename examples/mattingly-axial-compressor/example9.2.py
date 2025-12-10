"""
Mattingly Example 9.2 style multi-stage axial compressor sizing.

This uses simple meanline assumptions:
- Constant mean radius/area across stages.
- Rotor raises total pressure (per-stage PR from ψ, η_poly); stator diffuses with no added P0.
- Zero row losses except the polytropic target on the rotor.
"""

import math
from typing import List

import numpy as np

from turbodesign import BladeRow, Inlet, Outlet, Passage, PassageType, RowType
from turbodesign.compressor_spool import CompressorSpool
from turbodesign.loss.fixedpolytropic import FixedPolytropicEfficiency
from turbodesign.loss.fixedpressureloss import FixedPressureLoss


def mass_flow_parameter(M: float, gamma: float = 1.4, R: float = 287.15) -> float:
    expo = -(gamma + 1) / (2 * (gamma - 1))
    return math.sqrt(gamma / R) * M * (1 + (gamma - 1) / 2 * M * M) ** expo


def build_passage(stage_count: int, cax: float, rhub: float, rshroud: float) -> Passage:
    """Create a simple constant-radius passage with 2*stage_count+1 points."""
    npts = 2 * stage_count + 1
    x = np.linspace(0.0, cax * (npts - 1), npts)
    rhub_arr = np.full_like(x, rhub)
    rshroud_arr = np.full_like(x, rshroud)
    return Passage(x, rhub_arr, x, rshroud_arr, passageType=PassageType.Axial)


def main() -> None:
    # Inputs (Mattingly-like)
    T01 = 518.7  # R
    P01 = 14.7  # psia
    massflow = 50.0  # lbm/s
    M1 = 0.7
    alpha1 = 40.0  # deg
    omega = 1000.0  # rad/s
    target_overall_pr = 12.0
    psi = 0.35  # stage loading
    eta_poly = 0.88
    gamma = 1.4
    R = 287.15

    # Derived
    T01_K = T01 / 1.8
    P01_Pa = P01 * 6894.76
    massflow_kg_s = massflow * 0.453592

    mfp = mass_flow_parameter(M1, gamma, R)
    area = massflow_kg_s * math.sqrt(T01_K) / (P01_Pa * math.cos(math.radians(alpha1)) * mfp)

    rmean = 12 * 0.0254
    h = area / (math.pi * 4 * rmean)
    rhub = rmean - h
    rshroud = rmean + h
    cax = 1.0 * 0.0254

    per_stage_pr = (1 + eta_poly * psi * (gamma - 1) / gamma) ** (gamma / (gamma - 1))
    stage_count = math.ceil(math.log(target_overall_pr) / math.log(per_stage_pr))

    passage = build_passage(stage_count, cax, rhub, rshroud)

    inlet = Inlet(hub_location=0)
    inlet.alpha2 = [alpha1]
    inlet.init_total(P01_Pa, T01_K, M=M1)

    rows: List[BladeRow] = []
    for i in range(stage_count):
        rotor_loc = (2 * i + 1) / (2 * stage_count)
        stator_loc = (2 * i + 2) / (2 * stage_count)

        rotor = BladeRow(hub_location=rotor_loc, row_type=RowType.Rotor, stage_id=i)
        rotor.beta2_metal = [30.0]
        rotor.loss_function = FixedPolytropicEfficiency(eta_poly)
        rotor.P0_ratio = per_stage_pr

        stator = BladeRow(hub_location=stator_loc, row_type=RowType.Stator, stage_id=i)
        stator.beta2_metal = [alpha1]
        stator.loss_function = FixedPressureLoss(0.0)
        stator.P0_ratio = 1.0

        rows.extend([rotor, stator])

    outlet = Outlet()
    outlet.init_total(P01_Pa * target_overall_pr, 0.5)

    spool = CompressorSpool(
        passage,
        massflow_kg_s,
        inlet,
        outlet,
        rows,
        rpm=omega * 30 / math.pi,
    )

    spool.solve()
    achieved_pr = spool.overall_pressure_ratio()

    print("Mattingly Example 9.2 (multi-stage meanline)")
    print(f"Target overall PR: {target_overall_pr:0.2f}, per-stage PR (ideal): {per_stage_pr:0.3f}")
    print(f"Stage count used: {stage_count}")
    print(f"Achieved overall PR (inlet/stator-N exit): {achieved_pr:0.3f}")
    print(f"Massflow: {spool.massflow:0.3f} kg/s")


if __name__ == "__main__":
    main()
