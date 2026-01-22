import numpy as np

from turbodesign import Inlet, Outlet, Passage, PassageType
from turbodesign.row_factory import make_rotor_row, make_stator_row
from turbodesign.compressor_spool import CompressorSpool
from turbodesign.loss.fixedpressureloss import FixedPressureLoss
from turbodesign.row_factory import make_rotor_row, make_stator_row


def MFP(M: float, gamma: float = 1.4, R: float = 287.15) -> float:
    """Mass-flow parameter (Mattingly Eq. 9.x)."""
    expo = -(gamma + 1) / (2 * (gamma - 1))
    return np.sqrt(gamma / R) * M * (1 + (gamma - 1) / 2 * M * M) ** expo


def main() -> None:
    # Knowns (Mattingly Example 9.1)
    T01 = 518.7  # R
    P01 = 14.7  # psia
    omega = 1000  # rad/s
    r = 12  # in
    alpha1 = 40  # deg
    massflow = 50  # lbm/s
    M1 = 0.7
    u2_u1 = 1.1

    T01_K = T01 / 1.8
    P01_Pa = P01 * 6894.76
    massflow_kg_s = massflow * 0.453592

    area1 = 207.2 / 39.3701**2 
    area2 = 179.1 / 39.3701**2 
    area3 = 165.3 / 39.3701**2
    rmean = r * 0.0254  # convert inch to meter
    h1 = area1 / (np.pi * 4 * rmean)
    h2 = area2 / (np.pi * 4 * rmean)
    h3 = area3 / (np.pi * 4 * rmean)
    cax = 1 * 0.0254  # Assumed axial chord of 1 inch
    
    xhub_arr = [0, cax, 2 * cax]
    xshroud_arr = [0, cax, 2 * cax]
    # Shift both hub and shroud outward to achieve u2/u1 via mean-radius increase
    
    rhub_arr = [rmean - h1, rmean - h2, rmean - h3]  # Inlet exit, rotor exit, stator exit
    rshroud_arr = [rmean + h1, rmean + h2, rmean + h3]
    
    passage = Passage(xhub_arr, rhub_arr, xshroud_arr, rshroud_arr, passageType=PassageType.Axial, zero_phi=True)
    inlet = Inlet(hub_location=0)
    inlet.alpha2 = [alpha1]
    inlet.init_total(P01_Pa, T01_K, M=M1)

    rotor = make_rotor_row(
        hub_location=cax / max(xhub_arr),
        metal_exit_angle_deg=[-23.87],
        loss_function=FixedPressureLoss(0),
    )

    stator = make_stator_row(
        hub_location=2 * cax / max(xhub_arr),
        metal_exit_angle_deg=[alpha1],
        loss_function=FixedPressureLoss(0),
        P0_ratio=1.3,
    )

    outlet = Outlet()
    outlet.init_total(1.3 * P01_Pa, 0.5)

    spool = CompressorSpool(
        passage,
        massflow_kg_s,
        inlet,
        outlet,
        [rotor, stator],
        rpm=omega * 30 / np.pi,
    )
    spool.solve_balance_pressure()

    print("Mattingly Example 9.1 (single stage)")
    massflow_lbm_s = spool.massflow / 0.453592
    print(f"Massflow: {spool.massflow:0.3f} kg/s ({massflow_lbm_s:0.2f} lbm/s)")
    print(f"Overall total pressure ratio (inlet/stator exit): {spool.overall_pressure_ratio():0.3f}")
    print(f"Rotor exit Mach (meanline): {rotor.M.mean():0.3f}")
    print(f"Rotor exit Relative Mach (meanline): {rotor.M_rel.mean():0.3f}")
    print(f"Stator exit Mach (meanline): {stator.M.mean():0.3f}")


if __name__ == "__main__":
    main()
