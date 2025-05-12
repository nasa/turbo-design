from typing import Dict
import tecplot as tp
from tecplot.exception import *
from tecplot.constant import *
import pandas as pd
import math, json
import matplotlib.pyplot as plt 
import numpy as np 

R = 287.15
gamma = 1.35
Cp = gamma*R/(gamma-1)
n_blades = 5 
rpm = 50000

def set_fluid():
    tp.macro.execute_extended_command(command_processor_id='CFDAnalyzer4',
    command=f"SetFluidProperties Incompressible='F' Density=1 SpecificHeat={Cp} UseSpecificHeatVar='F' SpecificHeatVar=1 GasConstant={R} UseGasConstantVar='F' GasConstantVar=1 Gamma={gamma} UseGammaVar='F' GammaVar=1 Viscosity=1 UseViscosityVar='F' ViscosityVar=1 Conductivity=1 UseConductivityVar='F' ConductivityVar=1")
    
    # Set the field variables 
    tp.macro.execute_extended_command(command_processor_id='CFDAnalyzer4',
    command="SetFieldVariables ConvectionVarsAreMomentum='F' UVarNum=6 VVarNum=7 WVarNum=8 ID1='Pressure' Variable1=31 ID2='Temperature' Variable2=5")

def apply_equations():
    tp.data.operate.execute_equation(equation='{gam}='+f'{gamma}',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{Cp}='+f'{Cp}',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{R}='+f'{R}',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{RPM}=' + f'-{rpm}',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{Xm}={X}/1000',
    ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{Ym}={Y}/1000',
    ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{Zm}={Z}/1000',
    ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{P}={StaticPressure[kPa]}*1000',
    ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{T}={StaticTemperature[K]}',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{T0R}={RelativeTotalTemperature[K]}',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{P0R}={P}*({T0R}/{StaticTemperature[K]})**({gam}/({gam}-1))',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{radius}=sqrt({Ym}**2+{Zm}**2)',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{U}={RPM}*pi/30*{radius}',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{W}=sqrt({RelativeVelocity-X[m/s]}**2+{RelativeVelocity-Y[m/s]}**2+{RelativeVelocity-Z[m/s]}**2)',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{Vr}=(1/{radius})*({Zm}*{RelativeVelocity-Z[m/s]} + {Ym}*{RelativeVelocity-Y[m/s]})',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{Vm}=sqrt({RelativeVelocity-X[m/s]}**2 + {Vr}**2)',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{theta}=atan2({Zm},{Ym})',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{Wt}=(1/{radius})*({Ym}*{RelativeVelocity-Z[m/s]}-{Zm}*{RelativeVelocity-Y[m/s]})',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{beta} = atan2({Wt},{Vm})*180/pi',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{Vt} = {Wt} + {U}',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{V}=sqrt({Vt}**2+{Vm}**2)',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{W}=sqrt({Wt}**2+{Vm}**2)',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{a}=sqrt({gam}*{R}*{T})',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{M_rel}={W}/{a}',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{M}={V}/{a}',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{M_rel}={W}/{a}',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{T0}={T} * (1 + ({gam}-1)/2*{M}**2)',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{M_isen}=sqrt(2/({gam}-1)*({T0}/{T}-1))',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{alpha} = atan2({Vt},{Vm})*180/pi',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{T0_version2}={T}+{V}**2/(2*{Cp})',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{P0}={P}*({T0}/{StaticTemperature[K]})**({gam}/({gam}-1))',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{rho}={P}/({R}*{T})',
        ignore_divide_by_zero=True)
    tp.data.operate.execute_equation(equation='{rhoVm}={rho}*{Vm}',
        ignore_divide_by_zero=True)

def ConvertFrame(dataset:tp.data.dataset.Dataset,cylindrical:bool=False):
    if cylindrical:
        tp.active_frame().plot().axes.x_axis.variable_index=dataset.variable_names.index('Xm')
        tp.active_frame().plot().axes.y_axis.variable_index=dataset.variable_names.index('radius')
        tp.active_frame().plot().axes.z_axis.variable_index=dataset.variable_names.index('theta')
    else:
        tp.active_frame().plot().axes.x_axis.variable_index=dataset.variable_names.index('Xm')
        tp.active_frame().plot().axes.y_axis.variable_index=dataset.variable_names.index('Ym')
        tp.active_frame().plot().axes.z_axis.variable_index=dataset.variable_names.index('Zm')
    
def ExtractSlice(dataset:tp.data.dataset.Dataset, x:float=None,r:float=None):
    """_summary_

    Args:
        dataset (tp.data.dataset.Dataset): _description_
        x (float, optional): slice in constant x. Defaults to None.
        r (float, optional): slice in constant r. Defaults to None.

    Returns:
        zone number: _description_
    """
    tp.active_frame().plot().rgb_coloring.red_variable_index=3
    tp.active_frame().plot().rgb_coloring.green_variable_index=8
    tp.active_frame().plot().rgb_coloring.blue_variable_index=3
    tp.active_frame().plot().contour(0).variable_index=3
    tp.active_frame().plot().contour(1).variable_index=4
    tp.active_frame().plot().contour(2).variable_index=5
    tp.active_frame().plot().contour(3).variable_index=6
    tp.active_frame().plot().contour(4).variable_index=7
    tp.active_frame().plot().contour(5).variable_index=8
    tp.active_frame().plot().contour(6).variable_index=9
    tp.active_frame().plot().contour(7).variable_index=10
    tp.active_frame().plot(PlotType.Cartesian3D).show_slices=True
    if x:
        tp.active_frame().plot().slice(0).orientation=SliceSurface.XPlanes
        tp.active_frame().plot().slice(0).origin.x=x
        y = tp.active_frame().plot().slice(0).origin.y
        z = tp.active_frame().plot().slice(0).origin.z
    else:
        tp.active_frame().plot().slice(0).orientation=SliceSurface.YPlanes
        x = tp.active_frame().plot().slice(0).origin.x
        tp.active_frame().plot().slice(0).origin.y=r
        z = tp.active_frame().plot().slice(0).origin.z
        
    tp.active_frame().plot().slices(0).extract(transient_mode=TransientOperationMode.AllSolutionTimes)
    # tp.data.extract.extract_slice(
    #     origin=(x,y,z),
    #     normal=normal, 
    #     source=tp.constant.SliceSource.SurfaceZones, 
    #     dataset=dataset)
    return dataset.num_zones

def GetValuesFromSlice(zone:int,scalar_var:int=0,AverageType:str='MassFlowWeightedAverage',WeightedAvgVar:int=49) -> float:
    """Gets a value from slice 

    Args:
        zone (int): slice surface zone index 
        scalar_var (int): Tecplot variable index 
        AverageType (str): Selects the type of averaging to use. Defaults to MassFlowRate. Options: MassFlowRate, MassFlowWeightedAverage

    Returns:
        float: Value given from a partcular scalar_var 
    """
    frame = tp.active_frame()
    frame.plot()
    tp.macro.execute_extended_command(command_processor_id='CFDAnalyzer4',
        command=f"Integrate [{zone}] VariableOption='{AverageType}' XOrigin=0 YOrigin=0 ZOrigin=0 ScalarVar={scalar_var}" + f" Absolute='F' ExcludeBlanked='F' XVariable={WeightedAvgVar} YVariable=2 ZVariable=3" + " IntegrateOver='Cells' IntegrateBy='Zones' IRange={MIN =1 MAX = 0 SKIP = 1} JRange={MIN =1 MAX = 0 SKIP = 1} KRange={MIN =1 MAX = 0 SKIP = 1} PlotResults='F' PlotAs='Result' TimeMin=291.5947 TimeMax=291.5947")
    total = float(frame.aux_data['CFDA.INTEGRATION_TOTAL'])
    return total

def get_massflow(zone:int):
    frame = tp.active_frame()
    frame.plot()
    tp.macro.execute_extended_command(command_processor_id='CFDAnalyzer4', command=f"Integrate [{zone}] VariableOption='MassFlowRate' XOrigin=0 YOrigin=0 ZOrigin=0 ScalarVar=31 Absolute='F' ExcludeBlanked='F' XVariable=1 YVariable=2 ZVariable=3 IntegrateOver='Cells' IntegrateBy='Zones' " + "IRange={MIN =1 MAX = 0 SKIP = 1} JRange={MIN =1 MAX = 0 SKIP = 1} KRange={MIN =1 MAX = 0 SKIP = 1} PlotResults='F' PlotAs='Result' TimeMin=0 TimeMax=0")
    return float(frame.aux_data['CFDA.INTEGRATION_TOTAL'])
    
def get_key_index(d, key): # Get index of a key from a dictionary 
    try: 
        return list(d.keys()).index(key)
    except ValueError:
        return -1  # Return -1 if the key is not found

def calculate_properties(station01:Dict[str,float],station02:Dict[str,float],IsRotor:bool=True):
    P0_P = station01['P0']/station02['P']
    T3_is = station01['T0'] * (1/P0_P)**((gamma-1)/gamma)
   
    a = math.sqrt(gamma*R*station02['T'])
    T03_is = T3_is * (1+(gamma-1)/2*station02['M']**2) # This is probably V not Vm
    station02['total-total_power'] = (station01['T0'] - station02['T0'])*station02['massflow']
    station02['total-total_efficiency'] = (station01['T0'] - station02['T0'])/(station01['T0'] - T03_is)
    station02['total-static_efficiency'] = (station01['T0'] - station02['T0'])/(station01['T0'] - T3_is)
    station02['torque_efficiency'] = station02['total_power_torque']/(n_blades*station02['massflow']*(station01['T0'] - T03_is))
    station02['Total-Total_Power_kW_per_blade'] = station02['massflow'] * Cp * (station01['T0'] - station02['T0']) / 1000
    station02['Total-Total_Power_kW'] = n_blades*station02['Total-Total_Power_kW_per_blade']
    station02['Euler_Power_kW_per_blade'] = station02['massflow'] * (station01['U']*station01['Vt'] - station02['U']*station02['Vt'])/1000
    station02['Euler_Power_kW'] = n_blades*station02['Euler_Power_kW_per_blade']
    station02['P0_loss'] = (station01['P0R'] - station02['P0R'])/(station01['P0R'] - station02['P'])
    station02['T03_is'] = T03_is
    station02['T3_is'] = T3_is
    
def read_forces(file_path:str):
    table = []
    with open(file_path, 'r') as file:
        [file.readline() for _ in range(2)] # Skip first two blank lines
        header = file.readline().strip().split()
        file.readline() # Skip the units 
        for line in file:
            # Strip the line and split by whitespace (handles multiple spaces)
            row = line.strip().split()
            row = [float(r) for r in row]
            table.append(row)
    df = pd.DataFrame(table,columns=header)
    
    return df

def plot_velocity_triangles(station01:Dict[str,float],station02:Dict[str,float]):
    """Plots the velocity triangles for each blade row
        Made for turbines 
    """
    prop = dict(arrowstyle="-|>,head_width=0.4,head_length=0.8",
        shrinkA=0,shrinkB=0)
    
    x_start = 0
    y_max = 0; y_min = 0
    plt.figure(num=1,clear=True)
    for station in [station01,station02]:
        Vt = station['Vt']
        Wt = station['Wt']
        U = station['U']

        y_max = (Vt if Vt>y_max else y_max)
        y_max = (Wt if Wt>y_max else y_max)
        
        y_min = (Vt if Vt<y_min else y_min)
        y_min = (Wt if Wt<y_min else y_min)
        
    for station,row_name in zip([station01,station02],['stator','rotor']):
        x_end = x_start + station['Vm']
        dx = x_end - x_start
            
        Vt = station['Vt']
        Wt = station['Wt']
        U = station['U']
            
        plt.annotate("", xy=(x_end,Vt), xytext=(x_start,0), arrowprops=prop)        # V
        plt.text((x_start+x_end)/2,Vt/2*1.1,"V",fontdict={"fontsize":"xx-large"})
        
        plt.annotate("", xy=(x_end,Wt), xytext=(x_start,0), arrowprops=prop)        # W
        plt.text((x_start+x_end)/2,Wt/2*1.1,"W",fontdict={"fontsize":"xx-large"})
        
        if (abs(station['Vt']) > abs(station['Wt'])):
            # Shift Vt to right just a bit so we can see it
            plt.annotate("", xy=(x_end,Wt), xytext=(x_end,0), arrowprops=prop) # Wt
            plt.text(x_end+dx*0.1,Wt/2,"Wt",fontdict={"fontsize":"xx-large"})

            plt.annotate("", xy=(x_end,U+Wt), xytext=(x_end,Wt), arrowprops=prop) # U
            plt.text(x_end+dx*0.1,(Wt+U)/2,"U",fontdict={"fontsize":"xx-large"})
        else:
            # Shift Wt to right just a bit so we can see it
            plt.annotate("", xy=(x_end,Vt), xytext=(x_end,0), arrowprops=prop) # Vt
            plt.text(x_end+dx*0.1,Vt/2,"Vt",fontdict={"fontsize":"xx-large"})

            plt.annotate("", xy=(x_end,Wt+U), xytext=(x_end,Wt), arrowprops=prop) # U
            plt.text(x_end+dx*0.1,Wt+U/2,"U",fontdict={"fontsize":"xx-large"})

        if np.sign(Vt)>0:
            y = y_min/2
        else:
            y = y_max/2
        plt.text((x_start+x_end)/2,y,row_name,fontdict={"fontsize":"xx-large"})
        x_start += station['Vm']
        plt.axis([0,x_end+dx, y_min, y_max])
    plt.ylabel("Tangental Velocity [m/s]")
    plt.xlabel("Vm [m/s]")
    plt.title(f"Velocity Triangles")
    plt.savefig(f"velocity_triangles.png",transparent=False,dpi=150)
    

# Read .FORCES File
df = read_forces('CFD/radial-turbine.FORCES')
Torque = df['TQ-PF'].iloc[-1] + df['TQ-VF'].iloc[-1]
Power_torque_def = Torque * rpm*math.pi/30
Total_power_torque = Power_torque_def * n_blades # kW

tp.session.connect()
tp.active_page().name='Untitled'
tp.add_page()
tp.new_layout()

dataset = tp.data.load_tecplot_szl('CFD/radial-turbine.szplt')
tp.macro.execute_command('$!RedrawAll')
tp.active_frame().plot().fieldmaps(0,1,2,3,4).surfaces.surfaces_to_plot=SurfacesToPlot.BoundaryFaces
tp.active_frame().plot().show_mesh=True
tp.macro.execute_command('$!RedrawAll')
tp.active_frame().plot().show_shade=False
tp.macro.execute_command('$!RedrawAll')

apply_equations()
set_fluid()

variables_of_interest = ['P','T','P0','T0','T0R','P0R','rho','rhoVm']
variables_of_interest.extend(['U','W','V','Wt','Vt','Vm','U','Vr','M','M_rel'])
variables_of_interest.extend(['beta','alpha'])

ConvertFrame(dataset,cylindrical=False)
outlet_zone_id = ExtractSlice(dataset,x=11/1000)
outlet_data = {}
for v in variables_of_interest:
    outlet_data[v] = GetValuesFromSlice(outlet_zone_id,dataset.variable_names.index(v)+1)
outlet_data['massflow'] = get_massflow(outlet_zone_id)
outlet_data['total_massflow'] = get_massflow(outlet_zone_id)*n_blades
outlet_data['total_power_torque'] = Total_power_torque
outlet_data['torque_per_blade'] = Torque        # kN*m
outlet_data['torque_total'] = Torque*n_blades   # kN*m

ConvertFrame(dataset,cylindrical=True)
inlet_zone_id = ExtractSlice(dataset,r=0.055)
ConvertFrame(dataset,cylindrical=False)
inlet_data = {}
for v in variables_of_interest:
    inlet_data[v] = GetValuesFromSlice(inlet_zone_id,dataset.variable_names.index(v)+1)# ,'WeightedAverage',dataset.variable_names.index('rhoVm')+1)
inlet_data['massflow'] = get_massflow(inlet_zone_id)
inlet_data['total_massflow'] = get_massflow(inlet_zone_id)*n_blades

# Calculations
calculate_properties(inlet_data,outlet_data)

# Example dictionary
data = {
    'inlet':inlet_data,
    'outlet':outlet_data
}
plot_velocity_triangles(inlet_data,outlet_data)

# Write to JSON file
with open('results.json', 'w') as json_file:
    json.dump(data, json_file, indent=4)