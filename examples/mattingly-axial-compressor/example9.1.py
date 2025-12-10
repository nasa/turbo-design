import numpy as np

from turbodesign import BladeRow, Inlet, Outlet, Passage, PassageType, RowType
from turbodesign.compressor_spool import CompressorSpool
from turbodesign.loss.fixedpressureloss import FixedPressureLoss


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

    area = massflow_kg_s * np.sqrt(T01_K) / (P01_Pa * np.cos(np.radians(alpha1)) * MFP(M1))

    rmean = r * 0.0254  # convert inch to meter
    h = area / (np.pi * 4 * rmean)
    cax = 1 * 0.0254  # Assumed axial chord of 1 inch
    xhub_arr = [0, cax, 2 * cax]
    xshroud_arr = [0, cax, 2 * cax]
    # Shift both hub and shroud outward to achieve u2/u1 via mean-radius increase
    rmean2 = rmean * u2_u1
    rhub_arr = [rmean - h, rmean2 - h, rmean2 - h]  # Inlet exit, rotor exit, stator exit
    rshroud_arr = [rmean + h, rmean2 + h, rmean2 + h]

    passage = Passage(xhub_arr, rhub_arr, xshroud_arr, rshroud_arr, passageType=PassageType.Axial)
    inlet = Inlet(hub_location=0)
    inlet.alpha2 = [alpha1]
    inlet.init_total(P01_Pa, T01_K, M=M1)

    rotor = BladeRow(hub_location=cax/max(xhub_arr), row_type=RowType.Rotor)
    rotor.beta2_metal = [23.87]
    rotor.loss_function = FixedPressureLoss(0)
    rotor.P0_ratio = 1.3  # target rotor total-pressure ratio

    stator = BladeRow(hub_location=2*cax/max(xhub_arr), row_type=RowType.Stator)
    stator.beta2_metal = [alpha1]
    stator.loss_function = FixedPressureLoss(0)
    stator.P0_ratio = 1.0

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
    spool.massflow_constraint = 
    spool.solve()

    print("Mattingly Example 9.1 (single stage)")
    print(f"Massflow: {spool.massflow:0.3f} kg/s")
    print(f"Overall total pressure ratio (inlet/stator exit): {spool.overall_pressure_ratio():0.3f}")
    print(f"Rotor exit Mach (meanline): {rotor.M.mean():0.3f}")
    print(f"Stator exit Mach (meanline): {stator.M.mean():0.3f}")


if __name__ == "__main__":
    main()
