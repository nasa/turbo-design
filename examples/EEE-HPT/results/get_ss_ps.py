'''
    Read blade into suction and pressure sides
'''

import pickle
from typing import Tuple
import numpy as np
import matplotlib.pyplot as plt
import numpy.typing as npt
from pyturbo.helper import pspline, resample_by_curvature, order_points_nearest_neighbor
from scipy.signal import savgol_filter
from scipy.interpolate import splprep, splev
from pathlib import Path

import numpy as np
from scipy.signal import savgol_filter
from scipy.interpolate import splprep, splev

def split_airfoil_by_weighted_angle(coords, npts=100, window_length=11, polyorder=3):
    coords = np.asarray(coords)

    # Smooth coordinates
    x = savgol_filter(coords[:, 0], window_length, polyorder, mode='wrap')
    y = savgol_filter(coords[:, 1], window_length, polyorder, mode='wrap')
    z = savgol_filter(coords[:, 2], window_length, polyorder, mode='wrap')
    coords_smooth = np.stack([x, y, z], axis=1)

    # Centroid
    centroid = coords_smooth.mean(axis=0)
    cx, cy = centroid[:2]

    xmin, xmax = x.min(), x.max()

    def score_le(pt):
        dx = pt[0] - cx
        dy = pt[1] - cy
        if dx >= 0:
            return -np.inf
        angle = np.arctan2(dy, dx)
        dist = np.hypot(dx, dy)
        xpos_score = xmax - pt[0]  # larger near xmin
        return abs(angle) + dist + xpos_score

    def score_te(pt):
        dx = pt[0] - cx
        dy = pt[1] - cy
        if dx <= 0:
            return -np.inf
        angle = np.arctan2(dy, dx)
        dist = np.hypot(dx, dy)
        xpos_score = pt[0] - xmin  # larger near xmax
        return abs(angle) + dist + xpos_score

    # Find indices of best scoring LE and TE points
    scores_le = [score_le(pt) for pt in coords_smooth]
    scores_te = [score_te(pt) for pt in coords_smooth]
    le_index = np.argmax(scores_le)
    te_index = np.argmax(scores_te)

    # Split into suction and pressure sides
    if le_index < te_index:
        suction_raw = coords_smooth[le_index:te_index+1]
        pressure_raw = np.vstack([coords_smooth[te_index:], coords_smooth[:le_index+1]])
    else:
        suction_raw = np.vstack([coords_smooth[le_index:], coords_smooth[:te_index+1]])
        pressure_raw = coords_smooth[te_index:le_index+1]

    # Resample to M points using splines
    def resample_curve(curve, M):
        tck, _ = splprep(curve.T, s=0, per=0)
        u_new = np.linspace(0, 1, M)
        return np.stack(splev(u_new, tck), axis=1) # type: ignore

    ss_unique = suction_raw[np.sort(np.unique(suction_raw, axis=0, return_index=True)[1]),:]
    ps_unique = pressure_raw[np.sort(np.unique(pressure_raw, axis=0, return_index=True)[1]),:]
    
    suction = resample_curve(suction_raw, npts)
    pressure = resample_curve(pressure_raw, npts)
    
    plot_blade(ss_unique,ps_unique,'test_blade')
    plot_blade(suction,pressure,'test_blade2')
    
    return suction, pressure


def split_airfoil_by_angle_distance(coords, npts=100, w_angle=1.0, w_dist=1.0, window_length=11, polyorder=3):
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
        suction_raw = coords_smooth[le_index:te_index+1]
        pressure_raw = np.vstack([coords_smooth[te_index:], coords_smooth[:le_index+1]])
    else:
        suction_raw = np.vstack([coords_smooth[le_index:], coords_smooth[:te_index+1]])
        pressure_raw = coords_smooth[te_index:le_index+1]

    
    
    ss_unique = suction_raw[np.sort(np.unique(suction_raw, axis=0, return_index=True)[1]),:]
    ps_unique = pressure_raw[np.sort(np.unique(pressure_raw, axis=0, return_index=True)[1]),:]
    
    suction = resample_curve(suction_raw, npts)
    pressure = resample_curve(pressure_raw, npts)
    
    plot_blade(ss_unique,ps_unique,'test_blade')
    plot_blade(suction,pressure,'test_blade2')
    
    return suction, pressure


def split_airfoil_smart(coords:npt.NDArray, window_length:int=10, polyorder:int=3,npts:int=100):
    """
    coords: (N, 3) array of [x, y, z] airfoil coordinates
    Returns: suction_side, pressure_side (each (M, 3))
    """
    coords = np.asarray(coords)

    # Step 1: Smooth the coordinates
    x = savgol_filter(coords[:, 0], window_length, polyorder)
    y = savgol_filter(coords[:, 1], window_length, polyorder)
    z = savgol_filter(coords[:, 2], window_length, polyorder)
    coords_smooth = np.stack([x, y, z], axis=1)

    # Step 2: Find leading edge (closest point to centroid)
    centroid = coords_smooth.mean(axis=0)
    dists = np.linalg.norm(coords_smooth - centroid, axis=1)
    le_index = np.argmin(dists)

    # Step 3: Reorder points starting from LE
    coords_smooth = np.roll(coords_smooth, -le_index, axis=0)

    # Step 4: Find trailing edge (opposite side)
    arc_lengths = np.cumsum(np.linalg.norm(np.diff(coords_smooth[:, :2], axis=0), axis=1))
    arc_lengths = np.insert(arc_lengths, 0, 0)
    total_length = arc_lengths[-1]
    te_index = np.argmin(np.abs(arc_lengths - total_length / 2))

    # Step 5: Split into suction and pressure
    ss = coords_smooth[:te_index+1]
    ps = np.vstack([coords_smooth[te_index:], coords_smooth[0:1]])  # close loop

    
    # Step 6: Fit splines to each arc and resample M points
    def resample_curve(curve, M):
        tck, _ = splprep(curve.T, s=0, per=0)
        u_fine = np.linspace(0, 1, M)
        resampled = np.stack(splev(u_fine, tck), axis=1) # type: ignore
        return resampled

    ss_unique = ss[np.sort(np.unique(ss, axis=0, return_index=True)[1]),:]
    ps_unique = ps[np.sort(np.unique(ps, axis=0, return_index=True)[1]),:]
    
    ss = resample_curve(ss_unique,npts)
    ps = resample_curve(ps_unique,npts)
    
    plot_blade(ss_unique,ps_unique,'test_blade')
    plot_blade(ss,ps,'test_blade2')
    return ss, ps

def split_ss_ps(pts:npt.NDArray,npts:int=100) -> Tuple[npt.NDArray,npt.NDArray]:
    """Split the blade points into suction side and pressure side
    Args:
        pts (npt.NDArray): array containing blade points in cartesian coordinates [npts,2]
        npts (int, optional): Number of points to return for suction and pressure sides 
        
    Returns:
        Tuple[npt.NDArray,npt.NDArray]: suction side and pressure side points
    """
    dydx = np.gradient(pts[:,1], pts[:,0])
    
    le_indx1 = -1
    le_indx2 = -1
    te_indx1 = -1 
    te_indx2 = -1 
    
    min_indx = np.argmin(np.abs(pts[:,0] - pts[:,0].min()))
    max_indx = np.argmin(np.abs(pts[:,0] - pts[:,0].max()))

    i = min_indx
    for i in range(0,min_indx):
        if np.sign(dydx[min_indx]) != np.sign(dydx[i]):
            le_indx1 = i
            break
    
    for i in range(min_indx,pts.shape[0]):
        if np.sign(dydx[min_indx]) != np.sign(dydx[i]):
            le_indx2 = i
            break
    
    if le_indx1==-1:
        le_indx = le_indx2
    elif le_indx2 == -1:
        le_indx = le_indx1
    else:
        if pts[le_indx1,0]>pts[le_indx2,0]:
            le_indx = le_indx2
        else:
            le_indx = le_indx1
    
    i = max_indx
    for i in range(max_indx,0,-1):
        if np.sign(dydx[max_indx]) != np.sign(dydx[i]):
            te_indx1 = i
            break
    
    for i in range(max_indx,pts.shape[0]):
        if np.sign(dydx[max_indx]) != np.sign(dydx[i]):
            te_indx2 = i
            break
    
    if te_indx1==-1:
        te_indx = te_indx2
    elif te_indx2 == -1:
        te_indx = te_indx1
    else:
        if pts[te_indx1,0]>pts[te_indx2,0]:
            te_indx = te_indx2
        else:
            te_indx = te_indx1

     # Build SS and PS from Indices
    di = te_indx-le_indx
    new_pts = np.roll(pts,-le_indx,axis=0)
    ss = new_pts[:di+1,:]
    ps = new_pts[di:,:]
    ps = np.vstack([ps,ss[0,:]])
    
    ss_unique = ss[np.sort(np.unique(ss, axis=0, return_index=True)[1]),:]
    ps_unique = ps[np.sort(np.unique(ps, axis=0, return_index=True)[1]),:]
    
    ss_unique = ss_unique[order_points_nearest_neighbor(ss_unique)]
    ps_unique = ps_unique[order_points_nearest_neighbor(ps_unique)]
    ps_unique = np.flipud(ps_unique)
    
    ss = resample_by_curvature(ss_unique,npts,smoothing=0)
    ps = resample_by_curvature(ps_unique,npts,smoothing=0)

    # if bPlot:
    # plot_blade(ss_unique,ps_unique,'test_blade')
    # plot_blade(ss,ps,'test_blade2')
    
    return ss, ps

# Resample to M points each
def resample_curve(curve, M):
    tck, _ = splprep(curve.T, s=0, per=0)
    u_new = np.linspace(0, 1, M)
    return np.stack(splev(u_new, tck), axis=1) # type: ignore
    
def plot_blade(ss:npt.NDArray,ps:npt.NDArray,name:str):
    script_dir = Path(__file__).resolve().parent
    plt.figure(num=0,clear=True,figsize=(10,6))
    plt.plot(ss[:,0],ss[:,1],'-',label='Suction Side')
    plt.plot(ps[:,0],ps[:,1],'-',label='Pressure Side')
    plt.plot(ss[:,0],ss[:,1],'-',label='Suction Side')
    plt.plot(ps[:,0],ps[:,1],'-',label='Pressure Side')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'{name}')
    plt.axis('scaled')
    plt.savefig(str(script_dir / f'{name}.png'),dpi=300)

if __name__ == "__main__":
    script_dir = Path(__file__).resolve().parent
    data = pickle.load(open(str(script_dir / 'stator_rotor.pkl'),'rb'))
    ss,ps = split_ss_ps(data['Stator1'][0])
    print('check')
