from typing import List
import numpy as np
import numpy.typing as npt
import os, pyiges, pickle
import matplotlib.pyplot as plt
from pyiges.geometry import RationalBSplineCurve  # <-- import the class


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
    
    def extract_curves(name:str, iges:pyiges.Iges) -> List[npt.NDArray]:
        curve_delta=0.001
        os.makedirs('csv', exist_ok=True)
        sections = list() # These are the cuts at different radius
        for i in range(len(iges.items)): # these are the cut section
            if isinstance(iges.items[i], RationalBSplineCurve):
                curve = iges.items[i].to_geomdl()
                curve.delta=curve_delta
                points = np.array(curve.evalpts); n = points.shape[0]
                sections.append(points)
                np.savetxt(f'csv/{name}.csv',points,fmt="%f",delimiter=',',header='x,rtheta,r')
        return sections
            
        
    # load an example impeller
    igv = extract_curves("igv", pyiges.read('geometry/compr_igv.igs'))
    rotor1 = extract_curves("rotor1", pyiges.read('geometry/compr_rotor1.igs'))
    stator1 = extract_curves("stator1", pyiges.read('geometry/compr_stator1.igs'))

    rotor2 = extract_curves("rotor2", pyiges.read('geometry/compr_rotor2.igs'))
    stator2 = extract_curves("stator2", pyiges.read('geometry/compr_stator2.igs'))
    
    rotor3 = extract_curves("rotor3", pyiges.read('geometry/compr_rotor3.igs'))
    stator3 = extract_curves("stator3", pyiges.read('geometry/compr_stator3.igs'))
    
    rotor4 = extract_curves("rotor4", pyiges.read('geometry/compr_rotor4.igs'))
    stator4 = extract_curves("stator4", pyiges.read('geometry/compr_stator4.igs'))
    
    rotor5 = extract_curves("rotor5", pyiges.read('geometry/compr_rotor5.igs'))
    stator5 = extract_curves("stator5", pyiges.read('geometry/compr_stator5.igs'))
    
    rotor6 = extract_curves("rotor6", pyiges.read('geometry/compr_rotor6.igs'))
    stator6 = extract_curves("stator6", pyiges.read('geometry/compr_stator6.igs'))
    
    rotor7 = extract_curves("rotor7", pyiges.read('geometry/compr_rotor7.igs'))
    stator7 = extract_curves("stator7", pyiges.read('geometry/compr_stator7.igs'))
    
    rotor8 = extract_curves("rotor8", pyiges.read('geometry/compr_rotor8.igs'))
    stator8 = extract_curves("stator8", pyiges.read('geometry/compr_stator8.igs'))
    
    rotor9 = extract_curves("rotor9", pyiges.read('geometry/compr_rotor9.igs'))
    stator9 = extract_curves("stator9", pyiges.read('geometry/compr_stator9.igs'))
    
    rotor10 = extract_curves("rotor10", pyiges.read('geometry/compr_rotor10.igs'))
    stator10 = extract_curves("stator10", pyiges.read('geometry/compr_stator10.igs'))
    
    pickle.dump([igv,rotor1,stator1,rotor2,stator2,
                 rotor3,stator3,rotor4,stator4,
                 rotor5,stator5,rotor6,stator6,
                 rotor7,stator7,rotor8,stator8,
                 rotor9,stator9,rotor10,stator10],open('rotor_stator.pkl','wb'))

if __name__ == "__main__":
    Process_HubShroud_IGES()
    ProcessBlades_IGES()