from typing import Optional
from scipy.interpolate import interp1d,PchipInterpolator
from scipy.integrate import solve_ivp
import numpy as np  
from .bladerow import BladeRow
from .enums import RowType
import math 

def radeq(row:BladeRow,upstream:BladeRow,downstream:Optional[BladeRow]=None) -> BladeRow:
    """Solves the radial equilibrium equation for axial machines and returns the convergence. 

    Note:
        This function will give you T0, P0, Vm as a function of the radius. 
        
    Args:
        row (BladeRow): Current row
        upstream (BladeRow): Previous row
        downstream (BladeRow): Next row
        
    Returns:
        BladeRow: current row with T0, P0, and Vm calculated
    """
    row_radius = row.r     # Use these for gradient 
    up_radius = upstream.r
        
    def ode_radeq_streamtube(r:np.ndarray,y:np.ndarray):
        """Solves the radial equilibrium equation for a streamtube 

        Args:
            r (np.ndarray): radius not as a percent
            y (np.ndarray): Array containing [P0,Vt,VtU,T0]
            
        """
        P0 = y[0]
        T0 = y[1]
        Vm = y[2]
        r = row.r.mean()+r
        if r>row_radius[-1]:
            return [0,0,0]
        elif r<row_radius[0]:
            return [0,0,0]

        Cp = row.Cp
        phi = interp1d(row_radius, row.phi)(r)
        alpha = interp1d(row_radius, row.alpha2)(r)
        T = interp1d(row_radius, row.T)(r)
        P = interp1d(row_radius, row.P)(r)
        rm = interp1d(row_radius, row.rm)(r)
        rho = row.rho.mean()
        
        if (row.row_type == RowType.Rotor):
            omega = row.rpm*np.pi/30
            U = omega*r
        else:
            omega = 0
            U = 0
        gamma = row.gamma

        # Solve the Radial Equlibrium 
        Vt = Vm*np.tan(alpha)
        Vr = Vm*np.sin(phi)
        # Estimations 
        dVm_dr = float(interp1d(row_radius, np.gradient(row.Vm, row_radius))(r))        
        up_Vm = interp1d(row_radius, upstream.Vm)(r)
        
        if downstream:
            if downstream.row_type == RowType.Outlet:
                down_Vm = Vm
            else:
                down_Vm = interp1d(row_radius, downstream.Vm)(r)
        else:
            down_Vm = Vm
        up_m = interp1d(row_radius, upstream.m)(r)
        
        # Get a rough guess of dVm/dm
        if downstream!=None:
            down_m = interp1d(row_radius, downstream.m)(r)
            row_m = interp1d(row_radius, row.m)(r)
            if down_m != row_m:
                func_Vm_m = PchipInterpolator([up_m, row_m, down_m],[up_Vm, Vm, down_Vm])
            else:
                func_Vm_m = PchipInterpolator([up_m, row_m],[up_Vm, Vm])    
        else:
            func_Vm_m = PchipInterpolator([up_m, row_m],[up_Vm, Vm])     # type: ignore
        dVm_dm = func_Vm_m.derivative()(row_m) # type: ignore
        
        # Upstream 
        dT_dr = float(interp1d(row_radius, np.gradient(row.T,row_radius))(r))
        dP_dr = float(interp1d(row_radius, np.gradient(row.P,row_radius))(r))
        P = float(interp1d(row_radius, np.gradient(row.P,row_radius))(r))
        dT0_dr = dT_dr + Vm/Cp * (1 + np.tan(alpha)**2)*dVm_dr 
        dP0_dr = dP_dr * (T0/T)**(gamma/(gamma-1)) + P*gamma/(gamma-1) * (T0/T)**(1/(gamma-1)) * (T*dT0_dr-T0*dT_dr)/T**2
        # dT0_dr = float(interp1d(row.percent_hub_shroud, np.gradient(row.T0,row_radius))((r-row_radius[0])/(row_radius[-1]-row_radius[0]))) 
        # dP0_dr = float(interp1d(row.percent_hub_shroud, np.gradient(row.P0,row_radius))((r-row_radius[0])/(row_radius[-1]-row_radius[0]))) 
        
        C = (1 + np.tan(alpha)**2) * Vm**2/(2*Cp*T0)
        if (C>1):
            raise Exception(f"Invalid value of C {C:0.2f} which causes Vm to be nan.\nChange reduce alpha/beta for {row.row_type} {row.id}")
        B = (1-C)**(gamma/(gamma-1))
        A = -P0 * gamma/(gamma-1) * (1-C)**(1/(gamma-1)) * (1 + np.tan(alpha)**2)/(2*Cp)
        
        eqn15_rhs = Vt**2/r - Vm**2/rm*np.cos(phi) - Vr*dVm_dm # right hand side of equation 15
        eqn15_rhs_simple = Vt**2/r # right hand side of equation 15 simplified for axial machines
        
        epsilon = 1e-10  # or another small threshold
        if abs(rm) > epsilon:
            dVm_dr = T0/(2*Vm*A) * (rho*eqn15_rhs - B*dP0_dr) + Vm/(2*T0) * dT0_dr # Eqn 21
        else:
            dVm_dr = T0/(2*Vm*A) * (rho*eqn15_rhs_simple - B*dP0_dr) + Vm/(2*T0) * dT0_dr  # Eqn 21, simple 
        ydot = np.array([dP0_dr,dT0_dr,dVm_dr])

        return ydot

    T0 = row.T0
    
    P0 = row.P0
    Vm = row.Vm

    # Estimate the Vt based on a given turning angle 
    mean_radius = row_radius.mean()
    tip_radius = row_radius.max()
    hub_radius = row_radius.min()

    T0m = interp1d(row.percent_hub_shroud,T0)(0.5); 
    P0m = interp1d(row.percent_hub_shroud,P0)(0.5); Vmm = interp1d(row.percent_hub_shroud,Vm)(0.5)
    # We are solving for the values of these quantities at row exit
    ics = np.array([P0m,T0m,Vmm])

    # hub_to_tip = np.linspace(hub_radius,tip_radius)
    # res = odeint(ode_radeq_streamtube, ics, hub_to_tip, tfirst=True)
    
    # P0_new = interp1d(hub_to_tip,res[:,0])(row_radius)
    # T0_new = interp1d(hub_to_tip,res[:,1])(row_radius)
    # Vm_new = interp1d(hub_to_tip,res[:,2])(row_radius)
    r_eval = row.r - mean_radius 
    # mean_radius_to_tip = np.linspace(0,tip_radius-mean_radius,len(row_radius)*5)
    res1 = solve_ivp(ode_radeq_streamtube, t_span =[0, tip_radius-mean_radius], y0 = ics,
                     t_eval=np.linspace(0,tip_radius-mean_radius,len(row_radius)*2))
    
    # mean_radius_to_hub = np.linspace(0,hub_radius-mean_radius,len(row_radius)*5)
    res2 = solve_ivp(ode_radeq_streamtube, t_span = [0,hub_radius-mean_radius], y0 = ics,
                     t_eval=np.linspace(0,hub_radius-mean_radius,len(row_radius)*2))
    
    mid_to_tip_vals = res1.y.transpose()
    mid_to_tip_r = res1.t + mean_radius
    mid_to_hub_vals = res2.y.transpose()
    mid_to_hub_r = res2.t + mean_radius
    mid_to_hub_vals = np.flipud(mid_to_hub_vals)
    mid_to_hub_r = np.flipud(mid_to_hub_r)
    
    hub_to_tip_vals = np.concatenate([mid_to_hub_vals[:-1,:],mid_to_tip_vals])
    
    r = np.concatenate([mid_to_hub_r[:-1], mid_to_tip_r])
    
    P0_new = interp1d(r,hub_to_tip_vals[:,0])(row_radius)
    T0_new = interp1d(r,hub_to_tip_vals[:,1])(row_radius)
    Vm_new = interp1d(r,hub_to_tip_vals[:,2])(row_radius)
    
    row.P0 = P0_new
    row.T0 = T0_new
    row.Vm = Vm_new
    return row

