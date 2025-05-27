import pyiges
from pyiges import examples
import numpy as np
import matplotlib.pyplot as plt 
import pickle
from scipy.interpolate import BSpline, splrep, splev
from agf import Inlet_bcs, Outlet_bcs, Settings, AGF_Setup, Clearance
import subprocess


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
        stator_pts1.append({'ss':ss,'ps':ps})
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
        points = np.array(curve.evalpts); n = points.shape[0]
        ss = points[:n,:]; ps = points[n:,:]
        rotor_pts1.append({'ss':ss,'ps':ps})
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
        points = np.array(curve.evalpts); n = points.shape[0]
        ss = points[:n,:]; ps = points[n:,:]
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
        points = np.array(curve.evalpts); n = points.shape[0]
        ss = points[:n,:]; ps = points[n:,:]
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
    
    blades = pickle.load(open('stator_rotor.pkl','rb'))
    hub_shroud = pickle.load(open('hub_shroud.pkl','rb'))

    blades['Stator1']
    blades['Rotor1']
    blades['Stator2']
    blades['Rotor2']
    nblades = [46,76,48,70] # Vanes, Rotors, Vanes, Rotors
    
    T0 = 1588       # K
    P = 100         # kPa
    P0 = 4.933 * P  # kPa 
    CorrectedSpeed = 33.19 # rad/(sec * sqrt(K))
    RPM = CorrectedSpeed * np.sqrt(T0) * 30/np.pi
    
    inlet = Inlet_bcs(ptin=P0,ttin=1588,machin=0.05,alpin=0,phiin=0,pspan=50)
    outlet = Outlet_bcs(rpm=RPM,gamma=1.4,
                        psout=P,twall=0,molwt=28.96)
    
    settings = Settings(ifang=10)
    settings.nblades = nblades[0]
    
    clearance = Clearance(tlecl=0.005,tmccl=0.005,ttecl=0.005,hlecl=0,hmccl=0,htecl=0)
    
    agf = AGF_Setup(template_file='template.agf')
    agf.add_inlet(inlet=inlet)
    agf.add_outlet(outlet=outlet)
    agf.add_clearance(clearance=clearance)
    agf.add_settings(settings=settings)
    
    # cen.plot()
    agf.add_passage(hub=hub_shroud['hub'], shroud=hub_shroud['shroud'])
    # agf.add_blade(blades['Stator1'][])
    agf.build(output_filename='radial_turbine.agf')
    
    # Run Wand
    result = subprocess.run('./runwand.sh', shell=True, capture_output=True, text=True)

    # Check for negative cells 
    
    # Print the output
    with open("wand.stdout", "w") as file:
        file.write(result.stdout)
    # Print any errors
    if "NO NEGATIVE VOLUME in result.stdout":
        result = subprocess.run('./runleo.sh', shell=True, capture_output=True, text=True)
        with open("leo.stdout", "w") as file:
            file.write(result.stdout)
        if "A valid ADS license could not be acquired." not in "leo.stdout":
            # Plot convergence
            from plot_convergence import read_convergence
            import glob
            overall_files = list(glob.glob('*.OVERALL'))
            convergene_files = list(glob.glob('*.CONVERGENCE'))
            read_convergence(convergene_files[0])
            import post_process