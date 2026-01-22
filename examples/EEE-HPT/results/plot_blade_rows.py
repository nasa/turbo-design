from typing import List
import matplotlib.pyplot as plt
import numpy.typing as npt
from pathlib import Path

def plot_xz(hub:npt.NDArray,shroud:npt.NDArray,blades:List[npt.NDArray],plot_name:str=''):
    """_summary_

    Args:
        hub (npt.NDArray): _description_
        shroud (npt.NDArray): _description_
        blades (List[npt.NDArray]): _description_
    """
    # Plotting
    plt.figure(figsize=(10, 6))

    # Plot hub
    plt.plot(hub[:, 0], hub[:, 1], label='Hub', color='blue')

    # Plot shroud
    plt.plot(shroud[:, 0], shroud[:, 1], label='Shroud', color='red')

    # Plot blades (X-Z projection)
    for i, blade in enumerate(blades):
        ss = blade[0]
        ps = blade[1]
        for j in range(ss.shape[0]):  # m curves per blade
            x = ss[j, :, 0]
            z = ss[j, :, 2]
            plt.plot(x, z, '.',color='green', alpha=0.6)
            
        for j in range(ps.shape[0]):  # m curves per blade
            x = ps[j, :, 0]
            z = ps[j, :, 2]
            plt.plot(x, z, '.',color='blue', alpha=0.6)

    plt.xlabel('X')
    plt.ylabel('Z')
    plt.title('Blade, Hub, and Shroud in X-Z Plane')
    plt.legend()
    plt.grid(True)
    plt.axis('equal')  # Equal aspect ratio
    script_dir = Path(__file__).resolve().parent
    plt.savefig(str(script_dir / f'XZ - {plot_name}.jpg'),dpi=300)
    plt.show()