"""
EEE 2-Stage HPT Geometry Processing

This script reads IGES geometry files for the EEE 2-Stage High Pressure Turbine
and processes them into pickle files for use with turbo-design.

The processing pipeline:
1. Read hub/shroud IGES files -> save to hub_shroud.pkl
2. Read stator/rotor blade IGES files -> save to stator_rotor.pkl
3. Split blades into suction/pressure sides
4. Resample curves to uniform point counts
5. Fit blade sections to hub/shroud passage
6. Convert units from inches to mm
7. Save processed geometry to eee_hpt_processed.pkl

Usage:
    python eee_hpt_process_geometry.py

Requirements:
    - pyiges[full] (not available on macOS)
    - numpy, scipy, matplotlib
"""

import os
import pickle
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from scipy.interpolate import PchipInterpolator, splev, splprep
from scipy.signal import savgol_filter


def Process_HubShroud_IGES(geometry_dir: Path, output_dir: Path) -> dict:
    """Read hub and shroud curves from IGES files.

    Args:
        geometry_dir: Directory containing case.igs and hub.igs files
        output_dir: Directory to save outputs (csv, pkl, png files)

    Returns:
        Dictionary with 'Hub' and 'Shroud' numpy arrays (units: inches)
    """
    import pyiges

    iges_case = pyiges.read(str(geometry_dir / 'case.igs'))
    iges_hub = pyiges.read(str(geometry_dir / 'hub.igs'))

    curve = iges_case.items[0].to_geomdl()
    curve.delta = 0.001
    case_pts = np.array(curve.evalpts)

    curve = iges_hub.items[0].to_geomdl()
    curve.delta = 0.001
    hub_pts = np.array(curve.evalpts)

    # Save CSV files
    csv_dir = output_dir / 'csv'
    csv_dir.mkdir(exist_ok=True)
    np.savetxt(str(csv_dir / 'shroud.csv'), case_pts, fmt="%f", delimiter=',', header='x,r,theta')
    np.savetxt(str(csv_dir / 'hub.csv'), hub_pts, fmt="%f", delimiter=',', header='x,r,theta')

    # Plot flowpath
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(case_pts[:, 0], case_pts[:, 1], 'b-', linewidth=2, label='Shroud (Casing)')
    ax.plot(hub_pts[:, 0], hub_pts[:, 1], 'r-', linewidth=2, label='Hub')
    ax.fill_between(case_pts[:, 0], hub_pts[:, 1], case_pts[:, 1],
                    alpha=0.1, color='gray', label='Flow passage')
    ax.set_xlabel('Axial Position, x [inches]', fontsize=16)
    ax.set_ylabel('Radial Position, r [inches]', fontsize=16)
    ax.set_title('EEE 2-Stage HPT Flowpath (Meridional View)', fontsize=18, fontweight='bold')
    ax.legend(loc='center right', fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.set_aspect('equal', adjustable='box')
    ax.tick_params(axis='both', labelsize=12)
    fig.tight_layout()
    fig.savefig(str(output_dir / 'flowpath.png'), dpi=300, bbox_inches='tight')
    plt.close(fig)

    # Save pickle file
    hub_shroud = {'Hub': hub_pts, 'Shroud': case_pts}
    with open(str(output_dir / 'hub_shroud.pkl'), 'wb') as f:
        pickle.dump(hub_shroud, f)

    print(f"Hub/Shroud processed: {hub_pts.shape[0]} hub pts, {case_pts.shape[0]} shroud pts")
    return hub_shroud


def Process_StatorRotor_IGES(geometry_dir: Path, output_dir: Path) -> List[List[npt.NDArray]]:
    """Read stator and rotor blade curves from IGES files.

    Note: The IGES file naming is swapped in the original data:
        - hpt_stator1.igs contains rotor1
        - hpt_rotor1.igs contains stator1
        - etc.

    Args:
        geometry_dir: Directory containing blade IGES files
        output_dir: Directory to save outputs

    Returns:
        List of [stator1, rotor1, stator2, rotor2] blade sections
        Each blade is a list of numpy arrays (one per spanwise section)
        Units: inches, coordinates are [x, rtheta, r]
    """
    import pyiges

    # Load IGES files (note: naming is swapped in original data)
    iges_rotor1 = pyiges.read(str(geometry_dir / 'hpt_stator1.igs'))
    iges_stator1 = pyiges.read(str(geometry_dir / 'hpt_rotor1.igs'))
    iges_rotor2 = pyiges.read(str(geometry_dir / 'hpt_stator2.igs'))
    iges_stator2 = pyiges.read(str(geometry_dir / 'hpt_rotor2.igs'))

    curve_delta = 0.001
    csv_dir = output_dir / 'csv'
    csv_dir.mkdir(exist_ok=True)

    def extract_sections(iges_file, name: str, section_range: range) -> List[npt.NDArray]:
        """Extract blade sections from IGES file."""
        sections = []
        for idx, i in enumerate(section_range, start=1):
            curve = iges_file.items[i].to_geomdl()
            curve.delta = curve_delta
            points = np.array(curve.evalpts)
            sections.append(points)
            np.savetxt(str(csv_dir / f'{name}_{idx}.csv'), points,
                       fmt="%f", delimiter=',', header='x,rtheta,r')
        return sections

    # Stage 1
    stator_pts1 = extract_sections(iges_stator1, 'stator1', range(2, 7))
    rotor_pts1 = extract_sections(iges_rotor1, 'rotor1', range(2, 7))

    # Stage 2
    stator_pts2 = extract_sections(iges_stator2, 'stator2', range(2, 6))
    rotor_pts2 = extract_sections(iges_rotor2, 'rotor2', range(2, 7))

    blades = [stator_pts1, rotor_pts1, stator_pts2, rotor_pts2]

    # Save pickle file
    with open(str(output_dir / 'stator_rotor.pkl'), 'wb') as f:
        pickle.dump(blades, f)

    print(f"Blades processed: Stator1={len(stator_pts1)} sections, Rotor1={len(rotor_pts1)} sections, "
          f"Stator2={len(stator_pts2)} sections, Rotor2={len(rotor_pts2)} sections")

    return blades


def resample_curve(curve: npt.NDArray, M: int) -> npt.NDArray:
    """Resample a curve to M points using spline interpolation.

    Args:
        curve: Input curve points (N, 2) or (N, 3)
        M: Number of output points

    Returns:
        Resampled curve with M points
    """
    tck, _ = splprep(curve.T, s=0, per=0)
    u_new = np.linspace(0, 1, M)
    return np.stack(splev(u_new, tck), axis=1)


def split_airfoil_by_angle_distance(coords: npt.NDArray, npts: int = 100,
                                     window_length: int = 11, polyorder: int = 3
                                     ) -> Tuple[npt.NDArray, npt.NDArray]:
    """Split an airfoil curve into suction and pressure sides.

    Uses a scoring method based on angle from centroid and distance to
    identify leading edge (LE) and trailing edge (TE) points, then splits
    the curve accordingly.

    Args:
        coords: Airfoil coordinates (N, 3) as [x, y, z] or [x, rtheta, r]
        npts: Number of points for each resampled side
        window_length: Savitzky-Golay filter window length
        polyorder: Savitzky-Golay filter polynomial order

    Returns:
        Tuple of (suction_side, pressure_side), each (npts, 3)
    """
    coords = np.asarray(coords)

    # Smooth the coordinates
    x = savgol_filter(coords[:, 0], window_length, polyorder, mode='wrap')
    y = savgol_filter(coords[:, 1], window_length, polyorder, mode='wrap')
    z = savgol_filter(coords[:, 2], window_length, polyorder, mode='wrap')
    coords_smooth = np.stack([x, y, z], axis=1)

    # Compute centroid
    centroid = coords_smooth.mean(axis=0)
    cx, cy = centroid[:2]

    xmin, xmax = x.min(), x.max()

    # Precompute normalized distances to centroid
    dx_all = x - cx
    dy_all = y - cy
    dists = np.hypot(dx_all, dy_all)
    max_dist = np.max(dists)
    norm_dists = (dists / max_dist) * np.pi  # match angle scale

    def score_le(pt, norm_dist):
        dx = pt[0] - cx
        dy = pt[1] - cy
        if dx >= 0:
            return -np.inf
        angle = np.arctan2(dy, dx)
        xpos_score = xmax - pt[0]  # bias toward xmin
        return abs(angle) + norm_dist + xpos_score

    def score_te(pt, norm_dist):
        dx = pt[0] - cx
        dy = pt[1] - cy
        if dx <= 0:
            return -np.inf
        angle = np.arctan2(dy, dx)
        xpos_score = pt[0] - xmin  # bias toward xmax
        return abs(angle) + norm_dist + xpos_score

    # Find indices of best scoring LE and TE points
    scores_le = [score_le(pt, nd) for pt, nd in zip(coords_smooth, norm_dists)]
    scores_te = [score_te(pt, nd) for pt, nd in zip(coords_smooth, norm_dists)]
    le_index = np.argmax(scores_le)
    te_index = np.argmax(scores_te)

    # Split the airfoil into two curves
    if le_index < te_index:
        suction_raw = coords_smooth[le_index:te_index + 1]
        pressure_raw = np.vstack([coords_smooth[te_index:], coords_smooth[:le_index + 1]])
    else:
        suction_raw = np.vstack([coords_smooth[le_index:], coords_smooth[:te_index + 1]])
        pressure_raw = coords_smooth[te_index:le_index + 1]

    suction = resample_curve(suction_raw, npts)
    pressure = resample_curve(pressure_raw, npts)

    return suction, pressure


def fit_blade(hub: npt.NDArray, shroud: npt.NDArray,
              blades: List[Tuple[npt.NDArray, npt.NDArray]]) -> None:
    """Fit blade sections to lie within the hub/shroud passage.

    Adjusts the radial (z) coordinate of each blade section so that
    the blade spans from hub to shroud based on its relative position.

    Args:
        hub: Hub curve points (N, 2) as [x, r]
        shroud: Shroud curve points (N, 2) as [x, r]
        blades: List of (suction_side, pressure_side) tuples
                Each side is (nsections, npts, 3)

    Note:
        Modifies blades in place.
    """
    # Create interpolators with slight offset to avoid intersection
    func_hub = PchipInterpolator(hub[:, 0], hub[:, 1] * 0.99)
    func_shroud = PchipInterpolator(shroud[:, 0], shroud[:, 1] * 1.01)

    for blade in blades:
        ss = blade[0]
        ps = blade[1]

        # Make z-coordinate constant across each section
        for i in range(ss.shape[0]):
            ss[i, :, 2] = ss[i, 0, 2]
            ps[i, :, 2] = ss[i, 0, 2]

        # Calculate percent span for suction side
        max_index = np.argmax(ss[:, 0, 2])
        min_index = np.argmin(ss[:, 0, 2])
        height = ss[max_index, :, 2] - ss[min_index, :, 2]
        percent_hub_shroud = (ss[:, :, 2] - ss[min_index, :, 2]) / height

        for i in range(ss.shape[0]):
            for j in range(ss.shape[1]):
                ss[i, j, 2] = ((func_shroud(ss[i, j, 0]) - func_hub(ss[i, j, 0]))
                               * percent_hub_shroud[i, j] + func_hub(ss[i, j, 0]))

        # Calculate percent span for pressure side
        max_index = np.argmax(ps[:, 0, 2])
        min_index = np.argmin(ps[:, 0, 2])
        height = ps[max_index, :, 2] - ps[min_index, :, 2]
        percent_hub_shroud = (ps[:, :, 2] - ps[min_index, :, 2]) / height

        for i in range(ps.shape[0]):
            for j in range(ps.shape[1]):
                ps[i, j, 2] = ((func_shroud(ps[i, j, 0]) - func_hub(ps[i, j, 0]))
                               * percent_hub_shroud[i, j] + func_hub(ps[i, j, 0]))

        # Sort sections from hub to tip
        min_to_max = np.argsort(ss[:, 0, 2])
        blade[0][:] = ss[min_to_max, :, :]
        min_to_max = np.argsort(ps[:, 0, 2])
        blade[1][:] = ps[min_to_max, :, :]


def plot_xz(hub: npt.NDArray, shroud: npt.NDArray,
            blades: List[Tuple[npt.NDArray, npt.NDArray]],
            ax=None, plot_name: str = '') -> None:
    """Plot blades with hub and shroud in the X-Z (meridional) plane.

    Args:
        hub: Hub curve points (N, 2) as [x, r]
        shroud: Shroud curve points (N, 2) as [x, r]
        blades: List of (suction_side, pressure_side) tuples
        ax: Matplotlib axis to plot on (if None, creates new figure)
        plot_name: Optional name suffix for the plot title
    """
    # Plot hub and shroud with fill
    ax.plot(hub[:, 0], hub[:, 1], 'b-', linewidth=2.5, label='Hub')
    ax.plot(shroud[:, 0], shroud[:, 1], 'r-', linewidth=2.5, label='Shroud (Casing)')
    ax.fill_between(shroud[:, 0], hub[:, 1], shroud[:, 1],
                    alpha=0.08, color='blue')

    # Define colors and labels for blade rows
    colors = ['#2ecc71', '#e74c3c', '#9b59b6', '#f39c12']  # Green, Red, Purple, Orange
    labels = ['Stator 1 (Nozzle)', 'Rotor 1', 'Stator 2 (Nozzle)', 'Rotor 2']
    markers = ['o', 's', 'o', 's']  # Circles for stators, squares for rotors

    for i, blade in enumerate(blades):
        ss = blade[0]
        ps = blade[1]
        color = colors[i % len(colors)]
        marker = markers[i % len(markers)]

        for j in range(ss.shape[0]):
            label = labels[i] if j == 0 else None
            # Plot suction side
            ax.plot(ss[j, :, 0], ss[j, :, 2], marker, color=color,
                    markersize=1.5, alpha=0.7, label=label)
            # Plot pressure side
            ax.plot(ps[j, :, 0], ps[j, :, 2], marker, color=color,
                    markersize=1.5, alpha=0.7)

    ax.set_xlabel('Axial Position, x [mm]', fontsize=16)
    ax.set_ylabel('Radial Position, r [mm]', fontsize=16)

    title = 'EEE 2-Stage HPT Meridional View'
    if plot_name:
        title += f' ({plot_name.replace("_", " ").title()})'
    ax.set_title(title, fontsize=18, fontweight='bold')

    ax.legend(loc='center right', fontsize=12, framealpha=0.9)
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.set_aspect('equal', adjustable='box')
    ax.tick_params(axis='both', labelsize=12)

    # Add minor gridlines
    ax.minorticks_on()
    ax.grid(which='minor', linestyle=':', alpha=0.3)


def plot_blade(ss: npt.NDArray, ps: npt.NDArray, output_path: Path, name: str, base_color: str = 'blue') -> None:
    """Plot all spanwise sections of a blade with color gradient indicating percent span.

    Args:
        ss: Suction side points (nsections, npts, 3)
        ps: Pressure side points (nsections, npts, 3)
        output_path: Directory to save the plot
        name: Name for the plot file
        base_color: Base color for the gradient ('blue', 'orange', 'purple', 'green', 'red')
    """
    import matplotlib.colors as mcolors

    fig, ax = plt.subplots(figsize=(12, 8))

    # Define color maps for different blades
    color_maps = {
        'blue': plt.cm.Blues,
        'orange': plt.cm.Oranges,
        'purple': plt.cm.Purples,
        'green': plt.cm.Greens,
        'red': plt.cm.Reds
    }

    cmap = color_maps.get(base_color, plt.cm.Blues)
    n_sections = ss.shape[0]

    # Create color gradient from light to dark (hub to tip)
    colors = [cmap(0.3 + 0.7 * i / (n_sections - 1)) for i in range(n_sections)]

    for section_idx in range(n_sections):
        percent_span = (section_idx / (n_sections - 1)) * 100 if n_sections > 1 else 50
        color = colors[section_idx]

        # Plot suction side (solid line)
        label_ss = f'SS {percent_span:.0f}%'
        ax.plot(ss[section_idx, :, 0], ss[section_idx, :, 1],
               '-', color=color, linewidth=2, alpha=0.9, label=label_ss)

        # Plot pressure side (dashed line)
        label_ps = f'PS {percent_span:.0f}%'
        ax.plot(ps[section_idx, :, 0], ps[section_idx, :, 1],
               '--', color=color, linewidth=2, alpha=0.9, label=label_ps)

    ax.set_xlabel('Axial Position, x [mm]', fontsize=16)
    ax.set_ylabel('Tangential Position, r×θ [mm]', fontsize=16)
    ax.set_title(f'{name} - All Spanwise Sections', fontsize=18, fontweight='bold')
    ax.legend(fontsize=10, ncol=2, loc='best', framealpha=0.9)
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.tick_params(axis='both', labelsize=12)
    ax.axis('scaled')

    fig.tight_layout()
    fig.savefig(str(output_path / f'{name}_all_sections.png'), dpi=300, bbox_inches='tight')
    plt.close()


def plot_all_blades_combined(blades: List[Tuple[npt.NDArray, npt.NDArray]],
                             output_path: Path,
                             blade_names: List[str],
                             blade_colors: List[str]) -> None:
    """Plot all blade rows together on a single figure with color gradients for percent span.

    Args:
        blades: List of (suction_side, pressure_side) tuples for all blade rows
        output_path: Directory to save the plot
        blade_names: Names for each blade row (e.g., ['Stator1', 'Rotor1', ...])
        blade_colors: Color names for each blade row (e.g., ['green', 'red', 'purple', 'orange'])
    """
    import matplotlib.colors as mcolors
    from matplotlib.lines import Line2D

    fig, ax = plt.subplots(figsize=(18, 11))

    # Define color maps for different blades
    color_maps = {
        'blue': plt.cm.Blues,
        'orange': plt.cm.Oranges,
        'purple': plt.cm.Purples,
        'green': plt.cm.Greens,
        'red': plt.cm.Reds
    }

    # Custom legend elements
    legend_elements = []

    for blade_idx, (ss, ps) in enumerate(blades):
        blade_name = blade_names[blade_idx]
        base_color = blade_colors[blade_idx]
        cmap = color_maps.get(base_color, plt.cm.Blues)
        n_sections = ss.shape[0]

        # Create color gradient from dark to light (hub to tip)
        colors = [cmap(0.3 + 0.7 * i / (n_sections - 1)) for i in range(n_sections)]

        # Store hub and tip colors for legend
        hub_color = colors[0]  # Dark = hub (0%)
        tip_color = colors[-1]  # Light = tip (100%)

        for section_idx in range(n_sections):
            color = colors[section_idx]

            # Plot suction side (solid line) - no labels
            ax.plot(ss[section_idx, :, 0], ss[section_idx, :, 1],
                   '-', color=color, linewidth=2, alpha=0.85)

            # Plot pressure side (dashed line) - no labels
            ax.plot(ps[section_idx, :, 0], ps[section_idx, :, 1],
                   '--', color=color, linewidth=2, alpha=0.85)

        # Add custom legend entries for this blade (hub and tip only)
        legend_elements.append(Line2D([0], [0], color=hub_color, linewidth=3,
                                     label=f'{blade_name} Hub (0%)'))
        legend_elements.append(Line2D([0], [0], color=tip_color, linewidth=3,
                                     label=f'{blade_name} Shroud (100%)'))

    ax.set_xlabel('Axial Position, x [mm]', fontsize=22)
    ax.set_ylabel('Tangential Position, r×θ [mm]', fontsize=22)
    ax.set_title('All Blade Rows - Spanwise Sections with Percent Span', fontsize=24, fontweight='bold')

    # Create custom legend with only hub and tip entries
    ax.legend(handles=legend_elements, fontsize=16, ncol=2,
             loc='lower right', framealpha=0.95, borderaxespad=1)

    ax.grid(True, linestyle='--', alpha=0.4)
    ax.tick_params(axis='both', labelsize=16)
    ax.axis('scaled')

    # Add minor gridlines
    ax.minorticks_on()
    ax.grid(which='minor', linestyle=':', alpha=0.2)

    fig.tight_layout()
    fig.savefig(str(output_path / 'all_blades_combined.png'), dpi=300, bbox_inches='tight')
    plt.close()


def process_geometry(script_dir: Path, npts: int = 400) -> dict:
    """Complete geometry processing pipeline.

    Args:
        script_dir: Directory containing geometry/ folder and output location
        npts: Number of points per blade side after resampling

    Returns:
        Dictionary containing:
            - 'hub': Hub curve in mm (N, 2)
            - 'shroud': Shroud curve in mm (N, 2)
            - 'blades': List of (ss, ps) tuples in mm
    """
    geometry_dir = script_dir / 'geometry'

    # Step 1: Read IGES files
    print("Reading IGES files...")
    hub_shroud = Process_HubShroud_IGES(geometry_dir, script_dir)
    blades_raw = Process_StatorRotor_IGES(geometry_dir, script_dir)

    # Step 2: Process blades - split into SS/PS and resample
    print("\nProcessing blade sections...")
    processed_data = []
    for i, blade in enumerate(blades_raw):
        nsections = len(blade)
        n, m = blade[0].shape
        ss = np.zeros(shape=(nsections, npts, m))
        ps = np.zeros(shape=(nsections, npts, m))

        for section_index in range(nsections):
            ss_temp, ps_temp = split_airfoil_by_angle_distance(blade[section_index], npts=npts)
            ss[section_index, :, :] = ss_temp
            ps[section_index, :, :] = ps_temp

        # Sort sections by radial position (hub to tip)
        ss = ss[np.argsort(ss[:, 0, 2]), :, :]
        ps = ps[np.argsort(ps[:, 0, 2]), :, :]

        # Convert from inches to mm
        ss *= 25.4
        ps *= 25.4

        print(f'  Blade-{i} has {nsections} sections')
        processed_data.append((ss, ps))

    # Step 3: Resample and convert hub/shroud
    print("\nProcessing hub/shroud...")
    hub = resample_curve(hub_shroud['Hub'], 200) * 25.4
    shroud = resample_curve(hub_shroud['Shroud'], 200) * 25.4

    # Step 4 & 6: Create comparison plot (before and after fitting)
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 16))

    # Plot before fitting
    plot_xz(hub, shroud, processed_data, ax=ax1, plot_name='Before Fit')

    # Step 5: Fit blades to passage
    print("Fitting blades to passage...")
    fit_blade(hub, shroud, processed_data)

    # Plot after fitting
    plot_xz(hub, shroud, processed_data, ax=ax2, plot_name='After Fit')

    fig.tight_layout()
    fig.savefig(str(script_dir / 'XZ_before_after_fit_comparison.png'),
                dpi=300, bbox_inches='tight')
    plt.close(fig)

    # Step 7: Plot all spanwise sections for each blade with color gradients
    # Colors match the meridional view: Green, Red, Purple, Orange
    labels = ['Stator1', 'Rotor1', 'Stator2', 'Rotor2']
    colors = ['green', 'red', 'purple', 'orange']
    print("Creating blade profile plots with percent span gradients...")
    for i, (ss, ps) in enumerate(processed_data):
        plot_blade(ss, ps, script_dir, labels[i], base_color=colors[i])

    # Step 7b: Plot all blades combined on one figure
    print("Creating combined plot with all blade rows...")
    plot_all_blades_combined(processed_data, script_dir, labels, colors)

    # Step 8: Save processed data
    data = {
        'hub': hub,
        'shroud': shroud,
        'blades': processed_data
    }
    with open(str(script_dir / 'eee_hpt_processed.pkl'), 'wb') as f:
        pickle.dump(data, f)

    print(f"\nProcessed geometry saved to: {script_dir / 'eee_hpt_processed.pkl'}")

    return data


if __name__ == "__main__":
    import platform

    script_dir = Path(__file__).resolve().parent

    if platform.system() == "Darwin":
        print("Warning: pyiges[full] is not available on macOS.")
        print("Please run this script on Linux or Windows.")
        print("Alternatively, use pre-processed pickle files if available.")
    else:
        process_geometry(script_dir, npts=400)
        print("\nGeometry processing complete!")
