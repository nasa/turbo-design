import numpy as np

from turbodesign import Inlet, Outlet, Passage, PassageType
from turbodesign.row_factory import make_rotor_row, make_stator_row
from turbodesign.compressor_spool import CompressorSpool
from turbodesign.loss.fixedpressureloss import FixedPressureLoss


def MFP(M: float, gamma: float = 1.4, R: float = 287.15) -> float:
    """Mass-flow parameter (Mattingly Eq. 9.x)."""
    expo = -(gamma + 1) / (2 * (gamma - 1))
    return np.sqrt(gamma / R) * M * (1 + (gamma - 1) / 2 * M * M) ** expo


def main() -> None:
    # Knowns (Mattingly Example 9.1)
    T01 = 288.16  # K
    P01 = 101.3E3  # pa
    omega = 1000  # rad/s
    rmean = 0.3048  # m
    alpha1 = 40  # deg
    massflow = 22.68  # kg/s
    M1 = 0.7
    u2_u1 = 1.1     # In Mattingly this is the ratio of Vx2 / Vx1
    solidity = 1 # This is actually a calculated quantity based on the hub and shroud curves + number of blades
    phi_rotor = 0.09
    phi_stator = 0.03

    area1 = 0.134 # m^2
    area2 = 0.118 # m^2
    area3 = 0.111 # m^2
      # convert inch to meter
    h1 = area1 / (np.pi * 4 * rmean)
    h2 = area2 / (np.pi * 4 * rmean)
    h3 = area3 / (np.pi * 4 * rmean)
    
    cax = 1 * 0.0254  # Assumed axial chord of 1 inch
    
    xhub_arr = [0, cax, 2 * cax]
    xshroud_arr = [0, cax, 2 * cax]
    # Shift both hub and shroud outward to achieve u2/u1 via mean-radius increase
    number_of_blades = 65
    
    rhub_arr = [rmean - h1, rmean - h2, rmean - h3]  # Inlet exit, rotor exit, stator exit
    rshroud_arr = [rmean + h1, rmean + h2, rmean + h3]
    
    passage = Passage(xhub_arr, rhub_arr, xshroud_arr, rshroud_arr, passageType=PassageType.Axial, zero_phi=True)
    inlet = Inlet(hub_location=0)
    inlet.alpha2 = [alpha1]
    inlet.init_total(P01, T01, M=M1)

    rotor = make_rotor_row(
        hub_location=cax / max(xhub_arr),
        metal_exit_angle_deg=[-23.87],
        loss_function=FixedPressureLoss(phi_rotor),
        solidity=solidity,
        num_blades=int(number_of_blades),
        axial_chord=cax,
    )

    stator = make_stator_row(
        hub_location=2 * cax / max(xhub_arr),
        metal_exit_angle_deg=[alpha1],
        loss_function=FixedPressureLoss(phi_stator),
        P0_ratio=1.3,
        solidity=solidity,
        num_blades=int(number_of_blades),
        axial_chord=cax,
    )

    outlet = Outlet()
    outlet.init_total(1.3 * P01, 0.5)

    spool = CompressorSpool(
        passage,
        massflow,
        inlet,
        outlet,
        [rotor, stator],
        rpm=omega * 30 / np.pi,num_streamlines=1
    )
    spool.solve_balance_pressure()

    print("Mattingly Example 9.1 (single stage)")
    massflow_lbm_s = spool.massflow / 0.453592
    print(f"Massflow: {spool.massflow:0.3f} kg/s ({massflow_lbm_s:0.2f} lbm/s)")
    print(f"Overall total pressure ratio (inlet/stator exit): {spool.overall_pressure_ratio():0.3f}")
    print(f"Rotor exit Mach (meanline): {rotor.M.mean():0.3f}")
    print(f"Rotor exit Relative Mach (meanline): {rotor.M_rel.mean():0.3f}")
    print(f"Stator exit Mach (meanline): {stator.M.mean():0.3f}")
    print(f"Overall polytropic efficiency: {spool.overall_polytropic_efficiency():0.4f}")
    print(f"Rotor flow coefficient (meanline): {rotor.flow_coefficient:0.4f}")
    print(f"Stage loading (meanline): {rotor.stage_loading:0.4f}")
    print(f"Rotor total-to-total efficiency: {rotor.eta_total:0.4f}")
    

if __name__ == "__main__":
    main()
