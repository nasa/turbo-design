from scipy.interpolate import interp1d
from scipy.integrate import odeint
import numpy as np 
from .bladerow import BladeRow
from .enums import RowType


def radeq(row:BladeRow,upstream:BladeRow) -> BladeRow:
    """Solves the radial equilibrium equation for axial machines and returns the convergence. 

    Note:
        This function will give you T0, P0, Vm as a function of the radius. 
        
    Args:
        row (BladeRow): Current row
        upstream (BladeRow): Previous row

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
        if r>row_radius[-1]:
            return [0,0,0]
        elif r<row_radius[0]:
            return [0,0,0]

        Cp = row.Cp
        phi = interp1d(row_radius, row.phi)(r)
        alpha = interp1d(row_radius, row.alpha2)(r)
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
        Vt = Vm*np.cos(phi)*np.tan(alpha)
        Vr = Vm*np.sin(phi)
        # Estimations 
        dVm_dr = float(interp1d(row_radius, np.gradient(row.Vm, row_radius))(r))
        dVt_dr = dVm_dr*np.cos(phi)*np.tan(alpha)
        dVr_dr = 0 #dVm_dr*np.sin(phi)

        # Upstream 
        dT0up_dr = float(interp1d(upstream.percent_hub_shroud, np.gradient(upstream.T0,up_radius))((r-row_radius[0])/(row_radius[-1]-row_radius[0]))) # use percentage to get the T0 upstream value
        dP0up_dr = float(interp1d(upstream.percent_hub_shroud, np.gradient(upstream.P0,up_radius))((r-row_radius[0])/(row_radius[-1]-row_radius[0]))) # use percentage to get the T0 upstream value
        
        dP0_dr = float(interp1d(row.percent_hub_shroud, np.gradient(row.P0,row_radius))((r-row_radius[0])/(row_radius[-1]-row_radius[0]))) # use percentage to get 

        U_up = float(interp1d(upstream.percent_hub_shroud,up_radius*omega)((r-row_radius[0])/(row_radius[-1]-row_radius[0]))) # use percentage to get the T0 upstream value
        dVtup_dr = float(interp1d(upstream.percent_hub_shroud,np.gradient(upstream.Vt,up_radius))((r-row_radius[0])/(row_radius[-1]-row_radius[0])))
        Vtup = float(interp1d(upstream.percent_hub_shroud,upstream.Vt)((r-row_radius[0])/(row_radius[-1]-row_radius[0])))
        
        
        dT0_dr = dT0up_dr - 1/row.Cp*(U_up*dVtup_dr + Vtup*omega - (U*dVt_dr+Vt*omega)) # Eqn 8
        # if row.loss_function.LossType == LossType.Pressure: # type: ignore
        #     dP0_dr = dP0up_dr-row.Yp*(dP0up_dr - dP_dr)     # Eqn 9

        C = Vm**2*(1+np.cos(phi)**2 * np.tan(alpha)**2)/(2*Cp*T0)
        B = (1-C)**(gamma/(gamma-1))      
        A = P0 * gamma/(gamma-1) * (1-C)**(1/(gamma-1))
        dVm_dr = 1/(2*Vm*A) * (rho*(Vt/r - Vm**2/rm * np.cos(phi) - Vr*dVr_dr) - dP0_dr*B) + 1/(2*T0) *dT0_dr  # Eqn 6
        
        ydot = np.array([dP0_dr,dT0_dr,dVm_dr])

        return ydot

    T0 = row.T0
    
    P0 = row.P0
    Vm = row.Vm

    # Estimate the Vt based on a given turning angle 
    mean_radius = row_radius.mean()
    tip_radius = row_radius[-1]
    hub_radius = row_radius[0]

    T0m = interp1d(row.percent_hub_shroud,T0)(0.5); 
    P0m = interp1d(row.percent_hub_shroud,P0)(0.5); Vmm = interp1d(row.percent_hub_shroud,Vm)(0.5)
    # We are solving for the values of these quantities at row exit
    ics = np.array([P0m,T0m,Vmm])

    rm_to_tip = np.linspace(mean_radius,tip_radius)
    res1 = odeint(ode_radeq_streamtube, ics, rm_to_tip, tfirst=True)
    
    rm_to_hub = np.flip(np.linspace(hub_radius,mean_radius))
    res2 = odeint(ode_radeq_streamtube, ics, rm_to_hub, tfirst=True)
    
    res2 = np.flipud(res2)
    res = np.concatenate([res2[:-1,:],res1])
    r = np.concatenate([np.flip(rm_to_hub)[:-1], rm_to_tip])
    
    P0_new = interp1d(r,res[:,0])(row_radius)
    T0_new = interp1d(r,res[:,1])(row_radius)
    Vm_new = interp1d(r,res[:,2])(row_radius)
    
    row.P0 = P0_new
    row.T0 = T0_new
    row.Vm = Vm_new
    if row.row_type == RowType.Rotor:
        # U(VT1-VT2) = Power/massflow; VT2 = VT1 - Power/massflow 
        row.Vt = upstream.Vt-row.power/(row.total_massflow*row.U)
        row.alpha2 = np.arctan2(row.Vt,row.Vx)
    elif row.row_type == RowType.Stator:
        row.Vt = row.Vm*np.cos(row.phi)*np.tan(row.alpha2)
    row.Vr = row.Vm*np.sin(row.phi)
    row.Vx = row.Vm*np.cos(row.phi)
    
    return row


def radeq_normalized(row:BladeRow,upstream:BladeRow) -> BladeRow:
    """Solves the radial equilibrium equation for axial or radial machines and returns the convergence. 

    Args:
        row (BladeRow): Current row
        upstream (BladeRow): Previous row

    Returns:
        BladeRow: current row with Vt, T0, P0, and Vm calculated
    """
    _,row_radius = row.streamline.get_point(row.percent_hub_shroud)     # Use these for gradient 
    _,up_radius = upstream.streamline.get_point(upstream.percent_hub_shroud)

    def ode_radeq_streamtube(t:np.ndarray,y:np.ndarray):
        """Solves the radial equilibrium equation for a streamtube 

        Args:
            t (np.ndarray): percent from hub to shroud
            y (np.ndarray): Array containing [P0,Vt,VtU,T0]
            
        """
        P0 = y[0]
        T0 = y[1]
        Vm = y[2]
        if t>1:
            return [0,0,0]
        elif t<0:
            return [0,0,0]

        _,r = row.streamline.get_point()
        Cp = row.Cp
        # Interpolate angle of inclination (phi), exit flow angle (alpha), radius of curvature (rm) at a particular percentage from hub to shroud
        phi = interp1d(row.percent_hub_shroud, row.phi)(t)
        alpha = interp1d(row.percent_hub_shroud, row.alpha2)(t)
        rm = interp1d(row.percent_hub_shroud,row.rm)(t)
        rho = interp1d(row.percent_hub_shroud,row.rho)(t)
    
        
        if (row.row_type == RowType.Rotor):
            omega = row.rpm*np.pi/30
            U = omega*r
        else:
            omega = 0
            U = 0
        gamma = row.gamma

        # Solve the Radial Equlibrium 
        Vt = Vm*np.cos(phi)*np.tan(alpha)
        Vr = float(interp1d(row.percent_hub_shroud, row.Vr)(t))
        # Estimations: need the radius of the streamline to compute gradients  
        dVm_dr = interp1d(row_radius,np.gradient(row.Vm, row_radius))(r)
        dVt_dr = dVm_dr*np.cos(phi)*np.tan(alpha)
        
        dVm_dm = interp1d(row_radius, np.gradient(Vm, r)) + interp1d(x, np.gradient(Vm, r))
        
        # Upstream: We interpolate the gradient based on the percentage from hub to shroud
        dT0up_dr = interp1d(upstream.percent_hub_shroud,
                            np.gradient(upstream.T0,up_radius))(t)
        dP0up_dr = interp1d(upstream.percent_hub_shroud,
                            np.gradient(upstream.P0,up_radius))(t)
        dP0_dr = interp1d(upstream.percent_hub_shroud,
                          np.gradient(upstream.P0,up_radius))(t)

        U_up = interp1d(upstream.percent_hub_shroud,up_radius*omega)(t) # use percentage to get the T0 upstream value
        dVtup_dr = interp1d(upstream.percent_hub_shroud,np.gradient(upstream.Vt,up_radius))(t)
        Vtup = interp1d(upstream.percent_hub_shroud,upstream.Vt)(t)
        
        dP_dr = interp1d(row.percent_hub_shroud,np.gradient(row.P,row_radius))(t)
        dT0_dr = dT0up_dr - 1/row.Cp*(U_up*dVtup_dr + Vtup*omega - (U*dVt_dr+Vt*omega)) # Eqn 8
        # if row.loss_function.LossType == LossType.Pressure: # type: ignore
        #     dP0_dr = dP0up_dr-row.Yp*(dP0up_dr - dP_dr)     # Eqn 9

        C = Vm**2*(1+np.cos(phi)**2 * np.tan(alpha)**2)/(2*Cp*T0)
        B = (1-C)**(gamma/(gamma-1))      
        A = P0 * gamma/(gamma-1) * (1-C)**(1/(gamma-1))
        dVm_dr = 1/(2*Vm*A) * (rho*(Vt/r - Vm**2/rm * np.cos(phi)-Vr*dVm_dm) - dP0_dr*B) + 1/(2*T0) *dT0_dr # Eqn 6
        
        ydot = np.array([dP0_dr,dT0_dr,dVm_dr])

        return ydot

    T0 = row.T0
    P0 = row.P0
    Vm = row.Vm

    # Estimate the Vt based on a given turning angle 
    _, mean_radius = row.streamline.get_point(0.5)
    _, tip_radius = row.streamline.get_point(1)
    _, hub_radius = row.streamline.get_point(0)

    T0m = interp1d(row.percent_hub_shroud,T0)(0.5); 
    P0m = interp1d(row.percent_hub_shroud,P0)(0.5); 
    Vmm = interp1d(row.percent_hub_shroud,Vm)(0.5)
    # We are solving for the values of these quantities at row exit
    ics = np.array([P0m,T0m,Vmm])

    mid_to_tip = np.linspace(0,1)
    res1 = odeint(ode_radeq_streamtube, ics, mid_to_tip, tfirst=True) # Results
    
    mid_to_hub = np.flip(np.linspace(hub_radius,mean_radius))
    res2 = odeint(ode_radeq_streamtube, ics, mid_to_hub, tfirst=True) # Results
    
    res2 = np.flipud(res2)
    res = np.concatenate([res2[:-1,:],res1])    # Combine the results
    t = np.concatenate([np.flip(mid_to_hub)[:-1], mid_to_tip])
    
    P0_new = interp1d(t,res[:,0])(row.percent_hub_shroud)
    T0_new = interp1d(t,res[:,1])(row.percent_hub_shroud)
    Vm_new = interp1d(t,res[:,2])(row.percent_hub_shroud)
    
    row.P0 = P0_new
    row.T0 = T0_new
    row.Vm = Vm_new
    if row.row_type == RowType.Rotor:
        # U(VT1-VT2) = Power/massflow; VT2 = VT1 - Power/massflow 
        row.Vt = upstream.Vt-row.power/(row.total_massflow*row.U)
        row.alpha2 = np.arctan2(row.Vt,row.Vx)
    elif row.row_type == RowType.Stator:
        row.Vt = row.Vm*np.cos(row.phi)*np.tan(row.alpha2)
    row.Vr = row.Vm*np.sin(row.phi)
    row.Vx = row.Vm*np.cos(row.phi)
    
    return row
  

        
