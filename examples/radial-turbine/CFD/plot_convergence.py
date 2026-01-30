import matplotlib.pyplot as plt
import os.path as osp
import glob
import re
import numpy as np
import pandas as pd

# Set global plot styling
plt.rcParams.update({
    'font.size': 16,
    'axes.titlesize': 16,
    'axes.labelsize': 16,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 16,
    'figure.titlesize': 16,
    'lines.linewidth': 2,
    'lines.markersize': 4,
    'axes.grid': True,
    'grid.alpha': 0.3,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
    'axes.spines.top': False,
    'axes.spines.right': False,
}) 


def split_floats(text):
    # Regular expression to match floating point numbers (including scientific notation)
    pattern = r'[-+]?[0-9]*\.?[0-9]+(?:[Ee][-+]?[0-9]+)?'

    # Find all matches in the input text
    numbers = re.findall(pattern, text)

    # Convert matches to float
    return [float(num) for num in numbers]


def read_forces(file: str):
    """Read .FORCES file and extract power data (POWERin and ENERGYout)."""
    forces_data = list()
    headers = None

    with open(file, 'r') as f:
        for line in f:
            line_stripped = line.strip()
            if line_stripped.startswith('ITS'):
                # Parse header line
                headers = line_stripped.split()
                break

        # Skip the units line
        next(f, None)

        # Read data lines
        for line in f:
            floats = split_floats(line)
            if len(floats) == len(headers):
                forces_data.append(floats)

    if not forces_data or not headers:
        return None

    forces_data = np.array(forces_data)
    df = pd.DataFrame(forces_data, columns=headers)
    return df


def plot_power(forces_files: list):
    """Plot POWERin from .FORCES files for a single radial turbine."""
    if not forces_files:
        print("No .FORCES files found")
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    for forces_file in forces_files:
        df = read_forces(forces_file)
        if df is None:
            continue

        label = osp.basename(forces_file).replace('.FORCES', '')
        iterations = df['ITS']

        if 'POWERin' in df.columns:
            ax.plot(iterations, -df['POWERin'], label=label)

    ax.set_xlabel('Iterations')
    ax.set_ylabel('Power Out (kW)')
    ax.set_title('Power Out')
    ax.legend(loc='best', framealpha=0.9)

    fig.tight_layout()
    fig.savefig('convergence-power.png', dpi=150, bbox_inches='tight')
    print("Saved convergence-power.png")


def plot_massflow(forces_files: list):
    """Plot FLOWIN and FLOWOT from .FORCES files for a single radial turbine."""
    if not forces_files:
        print("No .FORCES files found")
        return

    fig, ax = plt.subplots(figsize=(10, 6))

    for forces_file in forces_files:
        df = read_forces(forces_file)
        if df is None:
            continue

        label = osp.basename(forces_file).replace('.FORCES', '')
        iterations = df['ITS']

        if 'FLOWIN' in df.columns:
            ax.plot(iterations, df['FLOWIN'], label=f'{label} - In')

        if 'FLOWOT' in df.columns:
            ax.plot(iterations, df['FLOWOT'], linestyle='--', label=f'{label} - Out')

    ax.set_xlabel('Iterations')
    ax.set_ylabel('Mass Flow (kg/s)')
    ax.set_title('Mass Flow')
    ax.legend(loc='best', framealpha=0.9)

    fig.tight_layout()
    fig.savefig('convergence-massflow.png', dpi=150, bbox_inches='tight')
    print("Saved convergence-massflow.png")

def read_convergence(file:str):
    convergence=list()
    with open(file,'r') as f:
        for line in f:
            if line.strip().startswith('ITS'):
                headers = [w for w in line.rstrip("\n").split(' ') if w]
                break
        for line in f:
            floats = split_floats(line)
            if len(floats)==15:
                convergence.append(floats)
    convergence = np.array(convergence)
    df = pd.DataFrame(convergence,columns=headers) # type: ignore
    iterations = df.iloc[:, 0]

    # Plot RHO*U, RHO*V, RHO*W, RHO*E on one figure
    residual_vars = ['RHO*U', 'RHO*V', 'RHO*W', 'RHO*E']
    fig, ax = plt.subplots(figsize=(10, 6))
    for h in residual_vars:
        if h in df.columns:
            ax.plot(iterations, df[h], label=h)
    ax.set_xlabel('Iterations')
    ax.set_ylabel('Residual')
    ax.set_title('Residual Convergence')
    ax.set_yscale('log')
    ax.legend(loc='best', framealpha=0.9)
    fig.tight_layout()
    fig.savefig('convergence-residuals.png', dpi=150, bbox_inches='tight')
    print("Saved convergence-residuals.png")
    plt.close(fig)

    # Determine convergence by looking at the last 1000 iterations, check mean and standard deviation
    last_iterations = int(1000/20)
    rho = df['RHO'].iloc[-last_iterations:]
    rhou = df['RHO*U'].iloc[-last_iterations:]
    rhov = df['RHO*V'].iloc[-last_iterations:]
    rhow = df['RHO*W'].iloc[-last_iterations:]
    rhoe = df['RHO*E'].iloc[-last_iterations:]
    mean_std = np.array([(rho.mean(),rho.std()),
                (rhou.mean(),rhou.std()),
                (rhov.mean(),rhov.std()),
                (rhow.mean(),rhow.std()),
                (rhoe.mean(),rhoe.std())])
    if np.max(mean_std[:,0])<1E-4 and np.max(mean_std[:,1])<1E-4:
        with open('converged.txt', 'w') as f: # Just create the file, do nothing
            pass
    

if __name__=="__main__":
    overall_files = list(glob.glob('*.OVERALL'))
    convergene_files = list(glob.glob('*.CONVERGENCE'))
    forces_files = list(glob.glob('*.FORCES'))

    read_convergence(convergene_files[0])
    plot_power(forces_files)
    plot_massflow(forces_files)
