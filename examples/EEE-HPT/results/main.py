from typing import List, Tuple
import numpy as np
import numpy.typing as npt 
import matplotlib.pyplot as plt 
from turbodesign.agf import Inlet_bcs, Outlet_bcs, Settings, AGF_Setup, Clearance
from turbodesign.row_factory import make_rotor_row, make_stator_row
import pickle, os, subprocess, platform, pyiges
from get_ss_ps import split_ss_ps, split_airfoil_smart,split_airfoil_by_angle_distance, resample_curve
from plot_blade_rows import plot_xz
from scipy.interpolate import PchipInterpolator
from mpl_toolkits.mplot3d import Axes3D

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

    plt.figure(num=0,figsize=(10,6))
    plt.plot(case_pts[:,0],case_pts[:,1])
    plt.plot(hub_pts[:,0],hub_pts[:,1])
    plt.axis('scaled')
    plt.savefig('flowpath.png',transparent=None,dpi=300)
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

def fit_blade(hub:npt.NDArray,shroud:npt.NDArray,blades:List[Tuple[npt.NDArray,npt.NDArray]]):
    func_hub = PchipInterpolator(hub[:,0],hub[:,1]*0.99)
    func_shroud = PchipInterpolator(shroud[:,0],shroud[:,1]*1.01)
    
    for s,blade in enumerate(blades):
        ss = blade[0]
        ps = blade[1]
        for i in range(ss.shape[0]):
            ss[i,:,2] = ss[i,0,2]
            ps[i,:,2] = ss[i,0,2]
        
        max_index = np.argmax(ss[:,0,2])
        min_index = np.argmin(ss[:,0,2])
        height = ss[max_index,:,2] - ss[min_index,:,2]
        percent_hub_shroud = (ss[:,:,2]-ss[min_index,:,2])/height # Location of the profiles as a percentage of hub and shroud 
        for i in range(ss.shape[0]):
            for j in range(ss.shape[1]):
                ss[i,j,2] = (func_shroud(ss[i,j,0]) -  func_hub(ss[i,j,0])) * percent_hub_shroud[i,j] + func_hub(ss[i,j,0])
        
        max_index = np.argmax(ps[:,0,2])
        min_index = np.argmin(ps[:,0,2])
        height = ps[max_index,:,2] - ps[min_index,:,2]
        percent_hub_shroud = (ps[:,:,2]-ps[min_index,:,2])/height # Location of the profiles as a percentage of hub and shroud 
        for i in range(ps.shape[0]):
            for j in range(ps.shape[1]):
                ps[i,j,2] = (func_shroud(ps[i,j,0]) -  func_hub(ps[i,j,0])) * percent_hub_shroud[i,j] + func_hub(ps[i,j,0])
        
        min_to_max = np.argsort(ss[:,0,2])
        ss = ss[min_to_max,:,:]
        min_to_max = np.argsort(ps[:,0,2])
        ps = ps[min_to_max,:,:]
        
        # fig = plt.figure(num=1,clear=True,dpi=200,figsize=(10,6))
        # ax = fig.add_subplot(111, projection='3d')

        # for i in range(ss.shape[0]):
        #     ax.plot(ss[i,:,0],ss[i,:,1],ss[i,:,2],'.')
        #     ax.plot(ps[i,:,0],ps[i,:,1],ps[i,:,2],'.')
        #     ax.axis('equal')
        #     plt.savefig(f'fit_blade_{s}_section_{i}.jpg',dpi=150)
        
        
if __name__ == "__main__":
    if platform.system() != "Darwin": # pyiges[full] does not work on MacOS        
        Process_HubShroud_IGES()
        Process_StatorRotor_IGES()
    
    blades = pickle.load(open('stator_rotor.pkl','rb'))
    processed_data = []
    npts = 400
    for i,blade in enumerate(blades):
        nsections = len(blade)
        n,m = blade[0].shape
        ss = np.zeros(shape=(nsections,npts,m))
        ps = np.zeros(shape=(nsections,npts,m))
        for section_index in range(nsections):
            ss_temp, ps_temp = split_airfoil_by_angle_distance(blade[section_index],npts=npts)
            # plot_blade(ss,ps,'test')
            ss[section_index,:,:] = ss_temp
            ps[section_index,:,:] = ps_temp
        ss = ss[np.argsort(ss[:,0,2]),:,:]
        ps = ps[np.argsort(ps[:,0,2]),:,:]
        ss *= 25.4 # convert to mm
        ps *= 25.4 # convert to mm
        print(f'Blade-{i} has {nsections} sections')
        processed_data.append((ss, ps))

    hub_shroud = pickle.load(open('hub_shroud.pkl','rb'))
    hub = resample_curve(hub_shroud['Hub'],200) * 25.4
    shroud = resample_curve(hub_shroud['Shroud'],200) * 25.4
    
    # Lets verify by plotting
    plot_xz(hub,shroud,processed_data,'before fit')
    fit_blade(hub,shroud,processed_data)
    plot_xz(hub,shroud,processed_data, 'after fit')
    
    # Output the AGF Files for ADS 
    nblades = [46,76,48,70] # Vanes, Rotors, Vanes, Rotors
    PT_Ratio = [2.25,2.11]  # P01/P03,P03/P05 These are stage total pressure ratio with 01 being total pressure of stator inlet 
    P0 = 1257514       # Pa
    T0 = 1587       # K
    P = 230300        # Pa
    
    CorrectedSpeed = 33.19 # rad/(sec * sqrt(K))
    RPM = 12400 # CorrectedSpeed * np.sqrt(T0) * 30/np.pi
    RPM = [0,RPM,0,RPM]     # Stator - Rotor - Stator - Rotor
    inlet = Inlet_bcs(ptin=P0,ttin=1587,machin=0.2,alpin=0,phiin=0,pspan=50)
    clearance = Clearance(tlecl=0.000,tmccl=0.000,ttecl=0.000,hlecl=0,hmccl=0,htecl=0)
    
    for i in range(len(processed_data)):
        ss = processed_data[i][0]; ps = processed_data[i][1]
        settings = Settings(ifang=0)
        settings.nblades = nblades[i]
        settings.ity = 5 # x,y,z
        settings.lete = 10 
        outlet = Outlet_bcs(rpm=RPM[i],gamma=1.3,
                        psout=P,twall=0,molwt=28.96)
        agf = AGF_Setup(template_file='template.agf',name=f"EEE-HPT-{i}")
        agf.add_inlet(inlet=inlet)
        agf.add_outlet(outlet=outlet)
        agf.add_clearance(clearance=clearance)
        agf.add_settings(settings=settings)
 
        # cen.plot()
        agf.add_passage(hub=hub, shroud=shroud)
        agf.add_blade(ss,ps)
        agf.build(output_filename=f'blade{i}.agf')
    
    # Run Wand
    result = subprocess.run('./runwand.sh', shell=True, capture_output=True, text=True)
    print(result.stdout)
    # # Check for negative cells 
    
    # Print the output
    with open("wand.stdout", "w") as file:
        file.write(result.stdout)
    # Print any errors
    if "NO NEGATIVE VOLUME in result.stdout":
        result = subprocess.run('./runleo.sh', shell=True, capture_output=True, text=True)
        print(result.stdout)
        with open("leo.stdout", "w") as file:
            file.write(result.stdout)
            
        if "A valid ADS license could not be acquired." not in "leo.stdout":
            # Plot convergence
            from plot_convergence import read_convergence
            import glob
            overall_files = list(glob.glob('*.OVERALL'))
            convergene_files = list(glob.glob('*.CONVERGENCE'))
            read_convergence(convergene_files[0])
