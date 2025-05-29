'''
    Read blade into suction and pressure sides 
'''

import pickle
from typing import Tuple
import numpy as np
import matplotlib.pyplot as plt 
import numpy.typing as npt
from pyturbo.helper import pspline, resample_by_curvature, order_points_nearest_neighbor

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
    ps = new_pts[di-1:,:]
    ss_unique = np.unique(ss,axis=0)
    ps_unique = np.unique(ps,axis=0)
    
    ss_unique = ss_unique[order_points_nearest_neighbor(ss_unique)]
    ps_unique = ps_unique[order_points_nearest_neighbor(ps_unique)]
    ss = resample_by_curvature(ss_unique,npts)
    ps = resample_by_curvature(ps_unique,npts)

    # if bPlot:
    # plot_blade(ss_unique,ps_unique,'test_blade')
    # plot_blade(ss,ps,'test_blade2')
    
    return ss, ps

def plot_blade(ss:npt.NDArray,ps:npt.NDArray,name:str):
    plt.figure(num=0,clear=True)
    plt.plot(ss[:,0],ss[:,1],'.',label='Suction Side')
    plt.plot(ps[:,0],ps[:,1],'.',label='Pressure Side')
    plt.plot(ss[:,0],ss[:,1],'.',label='Suction Side')
    plt.plot(ps[:,0],ps[:,1],'.',label='Pressure Side')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'{name}')
    plt.axis('scaled')
    plt.savefig(f'{name}.png',dpi=150)
    
if __name__ == "__main__":
    data = pickle.load(open('stator_rotor.pkl','rb'))
    ss,ps = split_ss_ps(data['Stator1'][0])
    print('check')
