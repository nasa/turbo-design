'''
    Read blade into suction and pressure sides 
'''

import pickle
import numpy as np
import matplotlib.pyplot as plt 

def calculate_slope(pts):
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
    ss = new_pts[:di,:]
    ps = new_pts[di:,:]
    
    plt.plot(ss[:,0],ss[:,1],ps[:,0],ps[:,1])
    plt.axis('scaled')
    plt.show()
    return dydx 


data = pickle.load(open('stator_rotor.pkl','rb'))
dydx = calculate_slope(data['Stator1'][0])

data['Rotor1']
data['Stator2']
data['Rotor2']
print('check')

