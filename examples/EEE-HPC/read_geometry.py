import numpy as np
import os, pyiges, pickle
import matplotlib.pyplot as plt

def Process_HubShroud_IGES():
    """Exports hub and shroud curves in inches
    """
    iges_case = pyiges.read('geometry/case.igs')
    iges_hub = pyiges.read('geometry/hub.igs')

    # print an invidiual entity (boring)
    curve = iges_case.items[0].to_geomdl()
    case_pts = np.array(curve.evalpts)

    curve = iges_hub.items[0].to_geomdl()
    curve.delta=0.001
    hub_pts = np.array(curve.evalpts)
    os.makedirs('csv',exist_ok=True)

    np.savetxt('csv/shroud.csv',case_pts,fmt="%f",delimiter=',',header='x,r,theta')
    np.savetxt('csv/hub.csv',hub_pts,fmt="%f",delimiter=',',header='x,r,theta')

    plt.figure(num=0,figsize=(10,6))
    plt.plot(case_pts[:,0],case_pts[:,1])
    plt.plot(hub_pts[:,0],hub_pts[:,1])
    plt.axis('scaled')
    plt.savefig('flowpath.png',transparent=None,dpi=300)
    # Hub and Shroud are in inches
    pickle.dump({'Hub':hub_pts,'Shroud':case_pts},open('hub_shroud.pkl','wb'))

def ProcessBlades_IGES():
    """Exports stator and rotor curves in inches
    """
 
    # load an example impeller
    iges_rotor1 = pyiges.read('geometry/compr_rotor1.igs')
    iges_stator1 = pyiges.read('geometry/compr_stator1.igs')

    iges_rotor2 = pyiges.read('geometry/compr_rotor2.igs')
    iges_stator2 = pyiges.read('geometry/compr_stator1.igs')
    
    iges_rotor3 = pyiges.read('geometry/compr_rotor3.igs')
    iges_stator3 = pyiges.read('geometry/compr_stator3.igs')
    
    iges_rotor4 = pyiges.read('geometry/compr_rotor4.igs')
    iges_stator4 = pyiges.read('geometry/compr_stator4.igs')
    
    iges_rotor5 = pyiges.read('geometry/compr_rotor5.igs')
    iges_stator5 = pyiges.read('geometry/compr_stator5.igs')
    
    iges_rotor6 = pyiges.read('geometry/compr_rotor6.igs')
    iges_stator6 = pyiges.read('geometry/compr_stator6.igs')
    
    iges_rotor7 = pyiges.read('geometry/compr_rotor7.igs')
    iges_stator7 = pyiges.read('geometry/compr_stator7.igs')
    
    iges_rotor8 = pyiges.read('geometry/compr_rotor8.igs')
    iges_stator8 = pyiges.read('geometry/compr_stator8.igs')
    
    iges_rotor9 = pyiges.read('geometry/compr_rotor9.igs')
    iges_stator9 = pyiges.read('geometry/compr_stator9.igs')
    
    iges_rotor10 = pyiges.read('geometry/compr_rotor10.igs')
    iges_stator10 = pyiges.read('geometry/compr_stator10.igs')
    
    curve_delta=0.001
    # Stage 1
    os.makedirs('csv',exist_ok=True)

    stator_pts1 = list(); indx = 1
    for i in range(2,7):
        curve = iges_stator1.items[i].to_geomdl()
        curve.delta=curve_delta
        points = np.array(curve.evalpts); n = points.shape[0]
        stator_pts1.append(points)
        os.makedirs('csv', exist_ok=True)
        np.savetxt(f'csv/stator1_{indx}.csv',points,fmt="%f",delimiter=',',header='x,rtheta,r')
        indx+=1

    # print an invidiual entity (boring)
    rotor_pts1 = list(); indx = 1
    for i in range(2,7):
        curve = iges_rotor1.items[i].to_geomdl()
        curve.delta=curve_delta
        points = np.array(curve.evalpts);
        rotor_pts1.append(points)
        np.savetxt(f'csv/rotor1_{indx}.csv',points,fmt="%f",delimiter=',',header='x,rtheta,r')
        indx+=1

    # Stage 2
    stator_pts2 = list(); indx = 1
    for i in range(2,6):
        curve = iges_stator2.items[i].to_geomdl()
        curve.delta=curve_delta
        points = np.array(curve.evalpts)
        stator_pts2.append(points)
        np.savetxt(f'csv/stator2_{indx}.csv',points,fmt="%f",delimiter=',',header='x,rtheta,r')
        indx+=1

    # print an invidiual entity (boring)
    rotor_pts2 = list(); indx = 1
    for i in range(2,7):
        curve = iges_rotor2.items[i].to_geomdl()
        curve.delta=curve_delta
        points = np.array(curve.evalpts)
        rotor_pts2.append(points)
        np.savetxt(f'csv/rotor2_{indx}.csv',points,fmt="%f",delimiter=',',header='x,rtheta,r')
        indx+=1

    pickle.dump([stator_pts1,rotor_pts1,stator_pts2,rotor_pts2],open('stator_rotor.pkl','wb'))

if __name__ == "__main__":
    Process_HubShroud_IGES()
    ProcessBlades_IGES()