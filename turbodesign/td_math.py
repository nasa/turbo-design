from typing import List, Optional, Tuple
import numpy as np
import math
import numpy.typing as npt
from .bladerow import BladeRow, compute_gas_constants
from .enums import RowType, LossType
from scipy.integrate import trapezoid
from .passage import Passage
from .isentropic import IsenP

def T0_coolant_weighted_average(row:BladeRow) -> float:
    """Calculate the new weighted Total Temperature array considering coolant

    Args:
        coolant (Coolant): Coolant
        massflow (np.ndarray): massflow mainstream

    Returns:
        float: Total Temperature drop
    """
    
    massflow = row.massflow
    total_massflow_no_coolant = row.total_massflow_no_coolant
    Cp = row.Cp
    
    Cpc = row.coolant.Cp
    T0c = row.coolant.T0
    massflow_coolant = row.coolant.massflow_percentage*total_massflow_no_coolant*row.massflow[1:]/row.massflow[-1] 
    if massflow_coolant.mean()>0:
        if row.row_type == RowType.Stator:
            T0= row.T0
            dT0 = T0.copy() * 0 
            T0_new = (massflow[1:]*Cp*T0[1:] + massflow_coolant*Cpc*T0c) \
                        /(massflow[1:]*Cp + massflow_coolant*Cpc)
            dT0[1:] = T0_new - row.T0[1:]
            dT0[0] = dT0[1]
        else:
            T0R = row.T0R
            T0R_new = T0R.copy()
            Cp = row.Cp
            T0R_new[1:] = (massflow[1:]*Cp*T0R[1:] + massflow_coolant*Cpc*T0c) \
                        /(massflow[1:]*Cp + massflow_coolant*Cpc)
            T0R_new[0] = T0R_new[1]
            
            T = T0R_new - row.W**2/(2*Cp)   # Dont change the velocity triangle but adjust the static temperature 
            T0_new = T+row.V**2/(2*Cp)      # Use new static temperature to calculate the total temperature 
            dT0 = T0_new - row.T0
        return dT0
    else:
        return row.T0*0

def compute_massflow(row:BladeRow) -> None:
    """Populates row.massflow and row.calculated_massflow 

    Calculated_massflow is massflow[-1]

    Args:
        row (BladeRow): current blade row. All quantities are at exit
        upstream (BladeRow): upstream blade row. All quantities are at exit
    """
    massflow_fraction =  np.linspace(0,1,len(row.percent_hub_shroud))
    massflow = row.percent_hub_shroud*0
    total_area = 0 
    for j in range(1,len(row.percent_hub_shroud)):
        Vm = (row.Vm[j]+row.Vm[j-1])/2
        rho = (row.rho[j]+row.rho[j-1])/2
        if np.abs((row.x[j]-row.x[j-1]))<1E-5: # Axial Machines
            total_area += np.pi*(row.r[j]**2-row.r[j-1]**2)
            massflow[j] = Vm * rho * np.pi* (row.r[j]**2-row.r[j-1]**2) + massflow[j-1]
        else:   # Radial Machines
            dx = row.x[j]-row.x[j-1]
            S = (row.r[j]-row.r[j-1])
            C = np.sqrt(1+((row.r[j]-row.r[j-1])/dx)**2)
            area = 2*np.pi*C*(S/2*dx**2+row.r[j-1]*dx)
            total_area += area
            massflow[j] = Vm * rho *area + massflow[j-1]
    
    row.total_massflow_no_coolant = massflow[-1]
    if row.coolant != None:
        massflow += massflow_fraction*row.coolant.massflow_percentage*row.total_massflow_no_coolant    # Take into account the coolant massflow
    row.massflow = massflow
    row.calculated_massflow = massflow[-1]
    row.total_massflow = massflow[-1]
    row.area = total_area

def compute_reynolds(rows:List[BladeRow],passage:Passage):
    """Calculates the Reynolds Number 

    Args:
        rows (List[BladeRow]): Blade row to calculate the Reynolds number
        passage (Passage): Passage 
    """
    
    for i in range(1,len(rows)):
        row = rows[i]
        xr = passage.get_xr_slice(0.5,(rows[i-1].location,row.percent_hub))
        dx = np.diff(xr[:,0])
        dr = np.diff(xr[:,1])
        c = np.sum(np.sqrt(dx**2+dr**2))
        mp = [2/(xr[i,1]+xr[i-1,1])*np.sqrt(dr[i-1]**2 + dx[i-1]**2) for i in range(1,len(xr[:,1]))]
        mp = np.hstack([[0],np.cumsum(mp)])
        
        if row.row_type == RowType.Rotor:
            V = row.W.mean()
        else:
            V = row.V.mean()
        rho = row.rho.mean()
        mu = row.mu
        row.Reynolds = c*V*rho/mu
        row.mprime = mp
        row.axial_chord = max(c,1E-12) # Axial chord
        # row.num_blades = int(2*np.pi*row.r.mean() / row.pitch_to_chord * row.axial_chord)

def compute_power(row:BladeRow,upstream:BladeRow) -> None:
    """Calculates the power

    Args:
        row (BladeRow): _description_
        upstream (BladeRow): _description_
    """
    if row.row_type == RowType.Stator:
        row.power = 0
        row.eta_static = 0
        row.eta_total = 0
        row.stage_loading = 0
        row.euler_power = 0
        row.T_is = 0 * row.T0 
        row.T0_is = 0 * row.T0 # Make it an array
    else:
        P0_P = (upstream.P0/row.P).mean()
        row.T_is = upstream.T0 * (1/P0_P)**((row.gamma-1)/row.gamma)
        a = np.sqrt(row.gamma*row.R*row.T_is)
        row.T0_is = row.T_is * (1+(row.gamma-1)/2*(row.V/a)**2)
        
        row.power = row.massflow[-1] * (row.Cp * (upstream.T0 - row.T0)).mean()
        # row.power = sum(v * w for v, w in zip(row.power[1:], np.diff(row.massflow))) # Massflow weighted average 
        row.eta_static = row.power/ (row.massflow[-1]*row.Cp*(upstream.T0.mean()-row.T_is.mean()))
        row.eta_total = (upstream.T0.mean() - row.T0.mean()) / (upstream.T0.mean() - row.T0_is.mean())
        row.stage_loading = row.Cp*(upstream.T0.mean() - row.T0.mean())/row.U.mean()**2
        row.euler_power = row.massflow[-1]* (upstream.U*upstream.Vt - row.U*row.Vt).mean()
    
def compute_quantities(row:BladeRow,upstream:BladeRow):
    """Calculation of all quantites after radial equilibrium has been solved assuming we know the static pressure at the exit.

    Note:
        Radial Equilibrium gives P0, T0, Vm. This code assumes the loss either enthalpy or pressure loss has already been calculated 

        compute_velocity has been called so we know W, Wt, V, Vt, U, M, M_rel

        Static Pressure and Temperature should come from Total Temperature and Pressure + Velocity 
    Args:
        row (BladeRow): current blade row. All quantities are at exit
        upstream (BladeRow): upstream blade row. All quantities are at exit
    """

    if row.row_type == RowType.Rotor:
        Cp_avg = (row.Cp+upstream.Cp)/2
        # Factor any coolant added and changes in streamline radius
        row.T0R = upstream.T0R - T0_coolant_weighted_average(row) # - (upstream.U**2-row.U**2)/(2*Cp_avg) 
        row.P = upstream.P0_stator_inlet/row.P0_P
        
        if row.loss_function.loss_type == LossType.Pressure: 
            # This affects the velocity triangle
            row.P0R = upstream.P0R - row.Yp*(upstream.P0R-row.P)
            row.T = (row.P/row.P0R)**((row.gamma-1)/row.gamma) * row.T0R
            row.T0 = (1+(row.gamma-1)/2 * row.M**2) * row.T
            row.power_distribution = row.massflow * row.Cp * (upstream.T0 - row.T0)
            row.power = np.trapezoid(row.power_distribution,row.r-row.r[0])
            row.power_mean = row.massflow[-1] * row.Cp * (upstream.T0.mean()-row.T0.mean())

        elif row.loss_function.loss_type == LossType.Enthalpy:
            ' For Enthalpy related loss, assume the static quantities do not change '
            row.T = (row.P/row.P0R)**((row.gamma-1)/row.gamma) * row.T0R
            row.T0 = (1+(row.gamma-1)/2 * row.M**2) * row.T

            def calculate_power(T0:npt.NDArray):
                row.power_distribution = row.massflow * row.Cp * (upstream.T0 - T0)
                row.power = np.trapezoid(row.power_distribution,row.r-row.r[0])
                row.power_mean = row.massflow[-1] * row.Cp * (upstream.T0.mean() - T0.mean())

                # Factor in T0R_drop. Convert T0R drop to absolute terms
                T_drop = (upstream.T0R - row.T0R) - row.W**2/(2*row.Cp) # row.T0R contains the drop
                T0_drop = T_drop*(1+(row.gamma-1)/2*row.M**2)

                T0 = upstream.T0 - row.power/row.eta_total/(row.total_massflow*row.Cp) + T0_drop
                return T0

            T0 = row.T0.copy()
            for _ in range(5):  # interate on the convergence of T0_drop
                T0 = calculate_power(T0)

            row.T0 = T0 
            row.T = row.T0-row.W**2/(2*row.Cp)
            row.P0R = row.P * (row.T0R/row.T)**((row.gamma)/(row.gamma-1))
            row.P0 = row.P * (row.T0/row.T)**((row.gamma)/(row.gamma-1))
            
    elif row.row_type == RowType.Stator:
        ' For the stator we already assume the upstream P0 already applied '
        if row.loss_function == LossType.Pressure:
            row.P0 = upstream.P0 - row.Yp*(upstream.P0-row.P)
        else:
            row.P0 = upstream.P0
        row.T0 = upstream.T0 - T0_coolant_weighted_average(row)
        row.T = row.T0 / (1+(row.gamma-1)/2*row.M**2)
        row.P = row.P0 * (row.T/row.T0)**((row.gamma)/(row.gamma-1))
        row.T0R = row.T + row.W**2 / (2*row.Cp)
        row.P0R = row.P*(row.T0R/row.T)**((row.gamma)/(row.gamma-1))
   
def stator_calc(row:BladeRow,upstream:BladeRow,downstream:Optional[BladeRow]=None,calculate_vm:bool=True):
    """Given P0, T0, P, alpha2 of stator calculate all other quantities

    Usage:
        Set row.P0 = upstream.P0 - any pressure loss
        row.T0 = upstream.T0 - any cooling
        row.P = row.rp*(row.P0 - rotor.P) + rotor.P 
        Set alpha2 
        
    Args:
        row (BladeRow): Stator Row
        upstream (BladeRow): Stator or Rotor Row 
        downstream (BladeRow): Stator or Rotor Row. Defaults to None
         
    """
    ## degree of reaction (rp) is assumed 
    # downstream.P = upstream.P0 * 1/downstream.P0_P 
    # if downstream is not None:
    #     # Use the upstream P0 value then later factor in the loss
    #     row.P = downstream.rp*(upstream.P0 - downstream.P) + downstream.P
    # else:
    #     row.P = upstream.P    
    
    # Static Pressure is assumed 
    row.P0 = upstream.P0 - row.Yp*(upstream.P0-row.P)
    
    if downstream is not None:
        row.P0_P = float((row.P0/downstream.P).mean())
        row.rp = ((row.P-downstream.P)/(upstream.P0-downstream.P)).mean()
        
    if calculate_vm:
        row.M = ((row.P0/row.P)**((row.gamma-1)/row.gamma) - 1) * 2/(row.gamma-1)
        row.M = np.sqrt(row.M)
        T0_T = (1+(row.gamma-1)/2 * row.M**2)
        row.T0 = upstream.T0 - T0_coolant_weighted_average(row)
        row.T = row.T0/T0_T
        row.V = row.M*np.sqrt(row.gamma*row.R*row.T)
        row.Vm = row.V*np.cos(row.alpha2)
        row.Vx = row.Vm*np.cos(row.phi)
        row.Vr = row.Vm*np.sin(row.phi)
        row.Vt = row.Vm*np.tan(row.alpha2)
    else: # We know Vm, P0, T0, P
        row.Vx = row.Vm*np.cos(row.phi)
        row.Vr = row.Vm*np.sin(row.phi)
        row.Vt = row.Vm*np.tan(row.alpha2)
        row.V = np.sqrt(row.Vx**2 + row.Vr**2 + row.Vt**2)
        
        row.T = row.P/(row.R*row.rho)   # We know P, this is a guess
        row.M = row.V/np.sqrt(row.gamma*row.R*row.T)
        
    if upstream.row_type == RowType.Rotor:
        row.alpha1 = upstream.alpha2 # Upstream rotor absolute frame flow angle
    row.beta1 = upstream.beta2
    row.rho = row.P/(row.R*row.T)
    row.U = row.omega*row.r
    row.Wt = row.Vt-row.U
    row.P0_stator_inlet = upstream.P0

def rotor_calc(row:BladeRow,upstream:BladeRow,calculate_vm:bool=True):
    """Calculates quantities given beta2 

    Args:
        row (BladeRow): Rotor Row
        upstream (BladeRow): Stator Row or Rotor Row
    """
    row.P0_stator_inlet = upstream.P0_stator_inlet
    ## P0_P is assumed 
    # row.P = row.P0_stator_inlet*1/row.P0_P
    
    # Static Pressure is assumed
    row.P0_P = (row.P0_stator_inlet/row.P).mean()
    upstream_radius = upstream.r
    row.U = row.omega*row.r
    # Upstream Relative Frame Calculations 
    upstream.U = upstream.rpm*np.pi/30 * upstream_radius # rad/s 
    upstream.Wt = upstream.Vt - upstream.U
    upstream.W = np.sqrt(upstream.Vx**2 + upstream.Wt**2 + upstream.Vr**2)
    upstream.beta2 = np.arctan2(upstream.Wt,upstream.Vm)
    upstream.T0R = upstream.T+upstream.W**2/(2*upstream.Cp)
    upstream.P0R = upstream.P * (upstream.T0R/upstream.T)**((upstream.gamma)/(upstream.gamma-1))      
    upstream.M_rel = upstream.W/np.sqrt(upstream.gamma*upstream.R*upstream.T)
    
    upstream_rothalpy = upstream.T0R*upstream.Cp - 0.5*upstream.U**2 # H01R - 1/2 U1^2 
    if np.any(upstream_rothalpy < 0):
        print('U is too high, reduce RPM or radius')
    # Rotor Exit Calculations
    row.beta1 = upstream.beta2
    #row.Yp # Evaluated earlier 
    row.P0R = upstream.P0R - row.Yp*(upstream.P0R-row.P)
    
    # Total Relative Temperature stays constant through the rotor. Adjust for change in radius from rotor inlet to exit
    row.T0R = upstream.T0R # (upstream_rothalpy + 0.5*row.U**2)/row.Cp # - T0_coolant_weighted_average(row) 
    P0R_P = row.P0R / row.P
    T0R_T = P0R_P**((row.gamma-1)/row.gamma)
    row.T = (row.T0R/T0R_T)     # Exit static temperature
    if calculate_vm:    # Calculates the T0 at the exit
        row.W = np.sqrt(2*row.Cp*(row.T0R-row.T)) #! nan popups here a lot for radial machines 
        if np.isnan(np.sum(row.W)):
            # Need to adjust T
            raise ValueError(f'nan detected: check flow path. Turbine inlet cut should be horizontal')
        row.Vr = row.W*np.sin(row.phi)
        row.Vm = row.W*np.cos(row.beta2)
        row.Wt = row.W*np.sin(row.beta2)
        row.Vx = row.Vm*np.cos(row.phi)
        row.Vt = row.Wt + row.U 
        row.V = np.sqrt(row.Vr**2+row.Vt**2+row.Vx**2)
        row.M = row.V/np.sqrt(row.gamma*row.R*row.T)
        row.Vm = np.sqrt(row.Vx**2+row.Vr**2)
        row.T0 = row.T + row.V**2/(2*row.Cp)
        row.P0 = row.P*(row.T0/row.T)**(row.gamma/(row.gamma-1))
        row.alpha2 = np.arctan2(row.Vt,row.Vm)
    else: # We know Vm, P0, T0
        row.Vr = row.Vm*np.sin(row.phi)
        row.Vx = row.Vm*np.cos(row.phi)
        
        row.W = np.sqrt(2*row.Cp*(row.T0R-row.T))
        row.Wt = row.W*np.sin(row.beta2)
        row.U = row.omega * row.r 
        row.Vt = row.Wt+row.U
        
        row.alpha2 = np.arctan2(row.Vt,row.Vm)
        row.V = np.sqrt(row.Vm**2*(1+np.tan(row.alpha2)**2))
        
        row.M = row.V/np.sqrt(row.gamma*row.R*row.T)
        T0_T = (1+(row.gamma-1)/2 * row.M**2)
        row.P0 = row.P * T0_T**(row.gamma/(row.gamma-1))
    
    row.M_rel = row.W/np.sqrt(row.gamma*row.R*row.T)
    row.T0 = row.T+row.V**2/(2*row.Cp)

def inlet_calc(row:BladeRow):
    """Calculates the conditions for the Inlet 

    Args:
        row (BladeRow): _description_
    """
    
    area = row.Vm.copy()*0
    # Estimate the density
    row.T = row.T0
    row.P = row.P0
    row.rho = row.P/(row.T*row.R)
    total_area = 0
    for iter in range(5): # Lets converge the Mach and Total and Static pressures
        for j in range(1,len(row.percent_hub_shroud)):
            rho = row.rho[j]
            tube_massflow = row.massflow[j]-row.massflow[j-1]
            if np.abs((row.x[j]-row.x[j-1]))<1E-6: # Axial Machines  
                total_area += np.pi*(row.r[j]**2-row.r[j-1]**2)
                row.Vm[j] = tube_massflow/(rho*np.pi*(row.r[j]**2-row.r[j-1]**2))
            else:   # Radial Machines
                dx = row.x[j]-row.x[j-1]
                S = (row.r[j]-row.r[j-1])
                C = np.sqrt(1+((row.r[j]-row.r[j-1])/dx)**2)
                area[j] = 2*np.pi*C*(S/2*dx**2+row.r[j-1]*dx)
                total_area += area[j]
                row.Vm[j] = tube_massflow/(rho*area[j])
        avg_mach = np.mean(row.M)
        row.Vm[0] = 1/(len(row.Vm)-1)*row.Vm[1:].sum() # Initialize the value at the hub to not upset the mean
        row.Vr = row.Vm*np.sin(row.phi)
        row.Vt = row.Vm*np.tan(row.alpha2)
        row.V = np.sqrt(row.Vt**2+row.Vm**2)
        # Fine tune the Temperature and Pressure and density
        row.M = row.V/np.sqrt(row.gamma*row.R*row.T)
        row.T = row.T0 * 1/(1+(row.gamma-1)/2*row.M**2)
        row.P = row.P0 * (row.T/row.T0)**(row.gamma/(row.gamma-1))
        compute_gas_constants(row)
        
    if np.mean(row.M)>0.5:
        raise ValueError(f"High inlet mach can lead to errors iter:{iter} Mach:{avg_mach}")
    
    if np.mean(row.M)<0.01:
        print(f"Unusually slow flow:{iter} Mach:{avg_mach}")