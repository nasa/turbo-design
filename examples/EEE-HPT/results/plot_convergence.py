import matplotlib.pyplot as plt
import os.path as osp
import glob
import re
import numpy as np
import pandas as pd
from pathlib import Path 


def split_floats(text):
    # Regular expression to match floating point numbers
    pattern = r'[-+]?[0-9]*\.?[0-9]+'
    
    # Find all matches in the input text
    numbers = re.findall(pattern, text)
    
    # Convert matches to float
    return [float(num) for num in numbers]

def read_convergence(file:str):
    script_dir = Path(file).resolve().parent
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
    plot_dir = osp.dirname(file)
    iterations = df.iloc[:, 0]

    # Check the plots
    for h in headers: # type: ignore
        plt.figure(num=1,clear=True,figsize=(10,5))
        plt.plot(iterations, df[h], marker='o', label=h)
        plt.xlabel('Iterations')
        plt.ylabel(f'{h}')
        plt.title(f'{h}')
        plt.grid(True)
        plt.yscale('log')
        h = h.replace('*','-')
        plt.savefig(str(script_dir / f'convergence-{h}.png'),dpi=150)

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
        with open(str(script_dir / 'converged.txt'), 'w') as f: # Just create the file, do nothing
            pass  
    

if __name__=="__main__":
    script_dir = Path(__file__).resolve().parent
    overall_files = list(glob.glob(str(script_dir / '*.OVERALL')))
    convergene_files = list(glob.glob(str(script_dir / '*.CONVERGENCE')))
    read_convergence(convergene_files[0])
