from typing import List, Tuple
import numpy as np
import numpy.typing as npt
import os, pyiges, pickle
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
from pyiges.geometry import RationalBSplineCurve  # <-- import the class
from scipy.interpolate import PchipInterpolator
from copy import deepcopy

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
                points = np.array(curve.evalpts)
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

def flip_blade(blade:List[npt.NDArray]):
    for section in blade:
        xmin, xmax = section[:,0].min(), section[:,0].max()
        section[:,0] = xmin + xmax - section[:,0]

def random_colors(n, seed=None):
    """
    Return n random RGB tuples in [0,1].
    Use a seed for reproducibility.
    """
    rng = np.random.default_rng(seed)
    return [tuple(rng.random(3)) for _ in range(n)]

def distinct_hsv_colors(n, seed=None, shuffle=True):
    """
    Return n visually distinct colors by sampling evenly around HSV hue.
    Optionally shuffle to reduce adjacent similarity.
    """
    hues = np.linspace(0, 1, n, endpoint=False)
    sat  = np.full(n, 0.65)   # medium saturation
    val  = np.full(n, 0.95)   # bright
    hsv  = np.stack([hues, sat, val], axis=1)
    rgb  = mcolors.hsv_to_rgb(hsv)
    idx  = np.arange(n)
    if shuffle:
        rng = np.random.default_rng(seed)
        rng.shuffle(idx)
    return [tuple(rgb[i]) for i in idx]

def tab_palette(n):
    """
    Pull colors from Matplotlib's qualitative palettes (loops if n is large).
    """
    cmap = plt.get_cmap('tab20')  # or 'tab10'
    return [cmap(i % cmap.N) for i in range(n)]

from scipy.signal import savgol_filter
from scipy.interpolate import splprep, splev

def resample_curve(curve, M):
    tck, _ = splprep(curve.T, s=0, per=0)
    u_new = np.linspace(0, 1, M)
    return np.stack(splev(u_new, tck), axis=1) # type: ignore

def split_airfoil_by_angle_distance(coords, npts=100, window_length=11, polyorder=3):
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

    # plot_blade(ss_unique,ps_unique,'test_blade')
    # plot_blade(suction,pressure,'test_blade2')

    return suction, pressure

def plot_blade(ss:npt.NDArray,ps:npt.NDArray,name:str):
    plt.figure(num=0,clear=True,figsize=(10,6))
    plt.plot(ss[:,0],ss[:,1],'-',label='Suction Side')
    plt.plot(ps[:,0],ps[:,1],'-',label='Pressure Side')
    plt.plot(ss[:,0],ss[:,1],'-',label='Suction Side')
    plt.plot(ps[:,0],ps[:,1],'-',label='Pressure Side')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(f'{name}')
    plt.axis('scaled')
    plt.savefig(f'{name}.png',dpi=300)
    
    
if __name__ == "__main__":
    # Process_HubShroud_IGES()
    # ProcessBlades_IGES()
    
    rotor_stator = []
    with open('rotor_stator.pkl','rb') as fp:
        rotor_stator = pickle.load(fp)
    with open('hub_shroud.pkl','rb') as fp:
        hub_shroud = pickle.load(fp)

    processed_data = []
    npts = 400
    for i,blade in enumerate(rotor_stator):
        nsections = len(blade)
        n,m = blade[0].shape
        ss = np.zeros(shape=(nsections,npts,m))
        ps = np.zeros(shape=(nsections,npts,m))
        for section_index in range(nsections):
            ss_temp, ps_temp = split_airfoil_by_angle_distance(blade[section_index],npts=npts,window_length=5)
            # plot_blade(ss_temp,ps_temp,'test')
            ss[section_index,:,:] = ss_temp
            ps[section_index,:,:] = ps_temp
            
        ss = ss[np.argsort(ss[:,0,2]),:,:]
        ps = ps[np.argsort(ps[:,0,2]),:,:]
        ss *= 25.4 # convert to mm
        ps *= 25.4 # convert to mm
        print(f'Blade-{i} has {nsections} sections')
        processed_data.append((ss, ps))

    hub_shroud = pickle.load(open('hub_shroud.pkl','rb'))
    hub = resample_curve(hub_shroud['Hub'],200) * 25.4          # Convert in to mm
    shroud = resample_curve(hub_shroud['Shroud'],200) * 25.4
    
    labels = ['IGV',
              'R1','S1','R2','S2',
              'R3','S3','R4','S4',
              'R5','S5','R6','S6',
              'R7','S7','R8','S8',
              'R9','S9','R10','S10']
    colors = distinct_hsv_colors(len(labels), seed=42)
    
    plt.figure(num=1)
    for i, sections in enumerate(processed_data):
        first = True
        ss,ps = sections[0], sections[1]
        for j in range(ss.shape[0]):
            coords = np.vstack([ss[j,:,:],ps[j,1:,:]])
            plt.plot(
                coords[:,0], coords[:,1],
                color=colors[i],
                label=labels[i] if first else "_nolegend_"
            )
            first = False

    # nice legend without duplicates
    plt.legend(loc="lower left", mode="expand", ncol=4)
    plt.xlabel('x - Axial')
    plt.ylabel('rtheta')
    plt.axis('equal')
    plt.show()
    
    plt.figure(num=2)
    plt.plot(hub[:,0],hub[:,1])
    plt.plot(shroud[:,0],shroud[:,1])
    for i, sections in enumerate(processed_data):
        first = True
        ss,ps = sections[0], sections[1]
        for j in range(ss.shape[0]):
            coords = np.vstack([ss[j,:,:],ps[j,1:,:]])
            plt.plot(
                coords[:,0], coords[:,2],
                color=colors[i],
                label=labels[i] if first else "_nolegend_"
            )
            first = False

    # nice legend without duplicates
    plt.legend(loc="lower left", mode="expand", ncol=4)
    plt.xlabel('x - Axial')
    plt.ylabel('rtheta')
    plt.axis('equal')
    plt.title('IGES Data position in the passage')
    plt.show()

    # Lets fit the blade
    fit_blade(hub,shroud,processed_data)
    
    # Re-Plot 
    plt.figure(num=3)
    plt.plot(hub[:,0],hub[:,1])
    plt.plot(shroud[:,0],shroud[:,1])
    for i, sections in enumerate(processed_data):
        first = True
        ss,ps = sections[0], sections[1]
        for j in range(ss.shape[0]):
            coords = np.vstack([ss[j,:,:],ps[j,1:,:]])
            plt.plot(
                coords[:,0], coords[:,2],
                color=colors[i],
                label=labels[i] if first else "_nolegend_"
            )
            first = False

    # nice legend without duplicates
    plt.legend(loc="lower left", mode="expand", ncol=4)
    plt.xlabel('x - Axial')
    plt.ylabel('rtheta')
    plt.axis('equal')
    plt.title('Blades fitted to the passage')
    plt.show()
    
    data = {}
    data['hub'] = hub
    data['shroud'] = shroud
    data['blades'] = processed_data
    pickle.dump(data,open('e3_hpc_processed.pkl','wb'))