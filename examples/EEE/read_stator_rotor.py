import pyiges
from pyiges import examples
import numpy as np
import matplotlib.pyplot as plt 
import pickle
from scipy.interpolate import BSpline, splrep, splev

def Process_HubShroud_IGES():
    iges_case = pyiges.read('case.igs')
    iges_hub = pyiges.read('hub.igs')

    # print an invidiual entity (boring)
    curve = iges_case.items[0].to_geomdl()
    case_pts = np.array(curve.evalpts)
    
    curve = iges_hub.items[0].to_geomdl()
    curve.delta=0.001
    hub_pts = np.array(curve.evalpts)

    np.savetxt('shroud.csv',case_pts,fmt="%f",delimiter=',',header='x,r,theta')
    np.savetxt('hub.csv',hub_pts,fmt="%f",delimiter=',',header='x,r,theta')

    plt.figure(num=0)
    plt.plot(case_pts[:,0],case_pts[:,1])
    plt.plot(hub_pts[:,0],hub_pts[:,1])
    plt.axis('scaled')
    plt.savefig('flowpath.png',transparent=None,dpi=150)


    pickle.dump({'Hub':hub_pts,'Shroud':case_pts},open('hub_shroud.pkl','wb'))


def Process_StatorRotor_IGES():
    # load an example impeller
    iges_rotor1 = pyiges.read('hpt_stator1.igs')
    iges_stator1 = pyiges.read('hpt_rotor1.igs')

    iges_rotor2 = pyiges.read('hpt_stator2.igs')
    iges_stator2 = pyiges.read('hpt_rotor2.igs')
    curve_delta=0.001
    # Stage 1 
    stator_pts1 = list(); indx = 1
    plt.figure(num=2,clear=True)
    for i in range(2,7):
        curve = iges_stator1.items[i].to_geomdl()
        curve.delta=curve_delta
        points = np.array(curve.evalpts); n = points.shape[0]
        ss = points[:n,:]; ps = points[n:,:]
        stator_pts1.append(points)
        np.savetxt(f'csv/stator1_{indx}.csv',stator_pts1[-1],fmt="%f",delimiter=',',header='x,rtheta,r')
        plt.plot(ss[:,0],ss[:,1],'.',label='ss')
        plt.plot(ps[:,0],ps[:,1],'.',label='ps')
        # plt.plot(stator_pts1[-1][:,0],stator_pts1[-1][:,1],'.')
        indx+=1
    plt.axis('scaled')
    plt.title('Stator')
    plt.savefig('Stator1.png',transparent=None,dpi=150)
    
    # print an invidiual entity (boring)
    rotor_pts1 = list(); indx = 1
    plt.figure(num=1,clear=True)
    for i in range(2,7):
        curve = iges_rotor1.items[i].to_geomdl()
        curve.delta=curve_delta
        points = np.array(curve.evalpts)
        rotor_pts1.append(points)
        np.savetxt(f'csv/rotor1_{indx}.csv',rotor_pts1[-1],fmt="%f",delimiter=',',header='x,rtheta,r')
        plt.plot(rotor_pts1[-1][:,0],rotor_pts1[-1][:,1],'.')
        indx+=1
    plt.axis('scaled')
    plt.title('Rotor')
    plt.savefig('Rotor1.png',transparent=None,dpi=150)
    

    # Stage 2
    stator_pts2 = list(); indx = 1 
    plt.figure(num=2,clear=True)
    for i in range(2,6):
        curve = iges_stator2.items[i].to_geomdl()
        curve.delta=curve_delta
        points = np.array(curve.evalpts)
        stator_pts2.append(points)
        np.savetxt(f'csv/stator2_{indx}.csv',stator_pts2[-1],fmt="%f",delimiter=',',header='x,rtheta,r')
        plt.plot(stator_pts2[-1][:,0],stator_pts2[-1][:,1],'.')
        indx+=1
    plt.axis('scaled')
    plt.title('Stator')
    plt.savefig('Stator2.png',transparent=None,dpi=150)
    
    # print an invidiual entity (boring)
    rotor_pts2 = list(); indx = 1
    plt.figure(num=1,clear=True)
    for i in range(2,7):
        curve = iges_rotor2.items[i].to_geomdl()
        curve.delta=curve_delta
        points = np.array(curve.evalpts)
        rotor_pts2.append(points)
        np.savetxt(f'csv/rotor2_{indx}.csv',rotor_pts2[-1],fmt="%f",delimiter=',',header='x,rtheta,r')
        plt.plot(rotor_pts2[-1][:,0],rotor_pts2[-1][:,1],'.')
        indx+=1
    plt.axis('scaled')
    plt.title('Rotor')
    plt.savefig('Rotor2.png',transparent=None,dpi=150)
    
    pickle.dump({
                    'Stator1':stator_pts1,
                    'Rotor1':rotor_pts1,
                    'Stator2':stator_pts2,
                    'Rotor2':rotor_pts2,
                 },open('stator_rotor.pkl','wb'))

def BladeExitLocations():
    data = pickle.load(open('stator_rotor.pkl','rb'))
    data['Stator1']
    data['Rotor1']
    data['Stator2']
    data['Rotor2']
    
if __name__ == "__main__":
    Process_HubShroud_IGES()
    Process_StatorRotor_IGES()