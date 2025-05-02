import json
import numpy as np 
from pyturbo.helper import bezier, arc, line2D, xr_to_mprime
from endwall import build_endwalls 

def read_input(input_filename):
    x = dict()
    with open(input_filename, "r") as f: 
        for line in f:
            split_val = line.split('=')
            if len(split_val)==2: # x1 = 2 # Grab the 2
                x[split_val[0].strip()] = float(split_val[1])
    return x
 
def print_output():
    import json 
    output = dict()
    with open('output.json') as o:
        data = json.load(o)
        output['M0'] = np.mean(data['blade_rows'][0]['M'])
        output['Ps0'] = np.mean(data['blade_rows'][0]['P'])
        output['Ts0'] = np.mean(data['blade_rows'][0]['T'])
        output['P0'] = np.mean(data['blade_rows'][0]['P0'])
        output['T0'] = np.mean(data['blade_rows'][0]['T0'])
        
        output['P01'] = np.mean(data['blade_rows'][1]['P0'])
        output['T01'] = np.mean(data['blade_rows'][1]['T0'])
        output['P1'] = np.mean(data['blade_rows'][1]['P'])
        output['T1'] = np.mean(data['blade_rows'][1]['T'])
        output['M1'] = np.mean(data['blade_rows'][1]['M'])
        output['M1_rel'] = np.mean(data['blade_rows'][1]['M_rel'])
        output['Vm1'] = np.mean(data['blade_rows'][1]['Vm'])
        output['W1'] = np.mean(data['blade_rows'][1]['W'])
        output['W1t'] = np.mean(data['blade_rows'][1]['Wt'])
        output['alpha1'] = np.mean(data['blade_rows'][1]['alpha1'])
        output['alpha2'] = np.mean(data['blade_rows'][1]['alpha2'])
        output['Cp1'] = np.mean(data['blade_rows'][1]['Cp'])
        output['gamma1'] = np.mean(data['blade_rows'][1]['gamma'])
        output['Vt1'] = np.mean(data['blade_rows'][1]['Vt'])
        
        output['P02'] = np.mean(data['blade_rows'][2]['P0'])
        output['T02'] = np.mean(data['blade_rows'][2]['T0'])
        output['P2'] = np.mean(data['blade_rows'][2]['P'])
        output['T2'] = np.mean(data['blade_rows'][2]['T'])
        output['M2'] = np.mean(data['blade_rows'][2]['M'])
        output['M2_rel'] = np.mean(data['blade_rows'][2]['M_rel'])
        output['Vm2'] = np.mean(data['blade_rows'][2]['Vm'])
        output['W2'] = np.mean(data['blade_rows'][2]['W'])
        output['W2t'] = np.mean(data['blade_rows'][2]['Wt'])
        output['beta1'] = np.mean(data['blade_rows'][2]['beta1'])
        output['beta2'] = np.mean(data['blade_rows'][2]['beta2'])
        output['Cp2'] = np.mean(data['blade_rows'][2]['Cp'])
        output['gamma2'] = np.mean(data['blade_rows'][2]['gamma'])
        output['Vt2'] = np.mean(data['blade_rows'][2]['Vt'])
        
        output['power'] = data['blade_rows'][2]['Power']
        output['eta_total'] = data['blade_rows'][2]['eta_total']
        output['eta_static'] = data['blade_rows'][2]['eta_static']
        output['euler_power'] = data['blade_rows'][2]['euler_power']
        output['degree_of_reaction'] = data['blade_rows'][2]['rp']
        output['P0_P'] = np.mean(data['blade_rows'][2]['P0_P'])
        output['Calculated_Massflow_Stator'] = data['blade_rows'][1]['calculated_massflow']
        output['Calculated_Massflow_Rotor'] = data['blade_rows'][2]['calculated_massflow']
        output['Calculated_Massflow_diff'] = abs(data['blade_rows'][1]['calculated_massflow'] - data['blade_rows'][2]['calculated_massflow'])
        output['dVt'] = output['Vt1'] - output['Vt2']
    
        output['massflow_deviation'] = abs(output['Calculated_Massflow_Rotor'] - 1)
        output['power_deviation'] = abs(output['power']-80E3)
        output['etal_total_inv'] = 1/output['eta_total']
        
        with open("output.txt", "w") as f:        
            for k,v in output.items():
                f.write(f'{k} = {v:0.6f}\n')
 
if __name__ == '__main__':
    x = read_input("input.dat")
    # Call Rosebrock test function 
    from single_execution import radial_turbine
    x = read_input('input.dat')
    
    #%% Define the Passage 
    # Shroud is defined using a thickness offset from the hub to construct a spline
    hub_inlet, hub, hub_outlet, shroud_inlet, shroud, shroud_outlet = build_endwalls(radius=x['radius'],
                                                                                hub_outlet_radius_scale=x['hub_outlet_radius_scale'],
                                                                                shroud_inlet_offset=x['shroud_inlet_offset'],shroud_outlet_offset_ratio=x['shroud_outlet_offset_ratio'],rhub_out=x['rhub_out'],
                                                                                x_stretch_factor=x['x_stretch_factor'],
                                                                                inlet_ext_percent=0.3,
                                                                                outlet_ext_percent=0.4)

    hub_all = np.vstack([hub_inlet[:-1,:],hub,hub_outlet[1:,:]])
    shroud_all = np.vstack([shroud_inlet[:-1,:],shroud,shroud_outlet[1:,:]])
    
    import matplotlib.pyplot as plt 
    plt.figure(num=1,clear=True)
    plt.plot(hub_inlet[:,0],hub_inlet[:,1],'k',linewidth=1.5,label='hub-inlet')
    plt.plot(hub[:,0],hub[:,1],'b',linewidth=1.5,label='hub')
    plt.plot(hub_outlet[:,0],hub_outlet[:,1],'m',linewidth=1.5,label='hub-outlet')
    
    plt.plot(shroud_inlet[:,0],shroud_inlet[:,1],'k',linewidth=1.5,label='shroud-inlet')
    plt.plot(shroud[:,0],shroud[:,1],'b',linewidth=1.5,label='shroud')
    plt.plot(shroud_outlet[:,0],shroud_outlet[:,1],'m',linewidth=1.5,label='shroud-outlet')
    
    plt.xlabel('x-axial')
    plt.ylabel('radius')
    plt.axis('equal')
    plt.legend()
    plt.savefig('design_passage.jpg',dpi=150)
    plt.close('all')
    total_hub_arc_len = xr_to_mprime(hub_all)[1][-1]
    inlet_arc_len = xr_to_mprime(hub_inlet)[1][-1]
    rotor_arc_len = xr_to_mprime(hub)[1][-1]
    
    blade_position = (0,(inlet_arc_len+rotor_arc_len)/total_hub_arc_len)
    P0 = x['Inlet_P0']
    P = x['Outlet_P']
    T0 = x['Inlet_T0']
    P0_P = P0/P
    RPM = -50000
    radial_turbine(P0=P0,P0_P=P0_P,
                   T0=T0,Design_RPM=RPM,
                   alpha2=x['alpha2'],beta3=x['beta3'],
                   massflow=x['massflow'],hub=hub_all,shroud=shroud_all,blade_position=blade_position)
    
    print_output()