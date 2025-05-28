import pyiges
import numpy as np
import matplotlib.pyplot as plt 
import pickle
from agf import Inlet_bcs, Outlet_bcs, Settings, AGF_Setup, Clearance
import subprocess
import platform
import os
from get_ss_ps import split_ss_ps

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

def BladeExitLocations():
    data = pickle.load(open('stator_rotor.pkl','rb'))

    
if __name__ == "__main__":
    if platform.system() != "Darwin": # pyiges[full] does not work on MacOS        
        Process_HubShroud_IGES()
        Process_StatorRotor_IGES()
    
    blades = pickle.load(open('stator_rotor.pkl','rb'))
    processed_data = []
    for blade in blades:
        section_data = list() 
        for section in blade:
            ss, ps = split_ss_ps(section)
            section_data.append({'ss':ss,'ps':ps})
        processed_data.append(section_data)

    hub_shroud = pickle.load(open('hub_shroud.pkl','rb'))

    
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
