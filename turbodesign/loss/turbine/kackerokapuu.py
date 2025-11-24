import pickle, os
from typing import Dict
from ...bladerow import BladeRow, sutherland
from ...lossinterp import LossInterp
from ...enums import RowType, LossType
import numpy as np
import pathlib
from ..losstype import LossBaseClass
import requests


def _mean_value(value):
    """Return the mean of an array-like as a Python float."""
    return float(np.asarray(value).mean())


def _mean_radians(value):
    """Alias for readability when working with angular quantities already in radians."""
    return _mean_value(value)


class KackerOkapuu(LossBaseClass):

    def __init__(self):
        """KackerOkapuu model is an improvement to the Ainley Mathieson model. 
        
        Limitations:
            - Doesn't factor incidence loss 
            - For steam turbines and impulse turbines
            
        Reference:
            Kacker, S. C., and U. Okapuu. "A mean line prediction method for axial flow turbine efficiency." (1982): 111-119.
            
        """
        super().__init__(LossType.Pressure)
        path = pathlib.Path(os.path.join(os.environ['TD3_HOME'],"kackerokapuu"+".pkl"))
        
        if not path.exists():
            url = "https://github.com/nasa/turbo-design/raw/main/references/Turbines/KackerOkapuu/kackerokapuu.pkl"
            response = requests.get(url, stream=True)
            with open(path.absolute(), mode="wb") as file:
                for chunk in response.iter_content(chunk_size=10 * 1024):
                    file.write(chunk)
        
        with open(path.absolute(),'rb') as f:
            self.data = pickle.load(f) # type: ignore
    
        
    
    def __call__(self,row:BladeRow, upstream:BladeRow) -> float:
        """Kacker Okapuu is an updated version of Ainley Mathieson and Dunham Came. This tool uses the pressure loss definition. 

        Note: 
            All equation numbers are from the Kacker Okapuu paper
        
        Reference:
            Kacker, S. C., and U. Okapuu. "A mean line prediction method for axial flow turbine efficiency." (1982): 111-119.
        
        Args:
            upstream (BladeRow): Upstream blade row
            row (BladeRow): downstream blade row

        Returns:
            float: Pressure Loss 
        """
        # Get the Inlet incoming mach number relative to the blade
        if upstream.row_type == RowType.Stator:
            M1 = _mean_value(upstream.M)
        else:
            M1 = _mean_value(upstream.M_rel) 
        
        c = row.chord
        b = row.axial_chord
        if row.row_type == RowType.Stator:
            alpha1_rad = _mean_radians(row.alpha1)
            beta1_rad = _mean_radians(row.beta1_metal)
            alpha2_rad = _mean_radians(row.alpha2)
            M2 = _mean_value(row.M)
            h = 0
            Rec = _mean_value(row.V*row.rho*row.chord / sutherland(row.T))
        else:
            h = row.tip_clearance * (row.r[-1]-row.r[0])
            alpha1_rad = _mean_radians(row.beta1)
            beta1_rad = _mean_radians(row.beta1_metal)
            alpha2_rad = _mean_radians(row.beta2)
            M2 = _mean_value(row.M_rel)
            Rec = _mean_value(row.W*row.rho*row.chord / sutherland(row.T))

        alpha1_deg = float(np.degrees(alpha1_rad))
        alpha2_deg = float(np.degrees(alpha2_rad))
        beta1_deg = float(np.degrees(beta1_rad))
            
        Yp_beta0 = self.data['Fig01_beta0'](float(row.pitch_to_chord), float(alpha2_deg))
        Yp_beta1_alpha2 = self.data['Fig02'](float(row.pitch_to_chord), float(alpha2_deg))
        t_max_c = self.data['Fig04'](float(np.abs(beta1_deg)+np.abs(alpha2_deg)))
        
        ratio = beta1_rad / alpha2_rad
        Yp_amdc = (Yp_beta0 + np.abs(ratio) * ratio * (Yp_beta1_alpha2-Yp_beta0)) * ((t_max_c)/0.2)**(ratio) # Eqn 2, AMDC = Ainley Mathieson Dunham Came
        
        # Shock Loss
        if M1>=0.4: # You'll have imaginary numbers if M1<0.4
            dP_q1_hub = 0.75*(M1-0.4)**1.75 # Eqn 4, this is at the hub
            dP_q1_shock = row.r[-1]/row.r[0] * dP_q1_hub # Eqn 5
            Y_shock = dP_q1_shock * upstream.P/row.P * (1-(1+(upstream.gamma-1)/2*M1**2))/(1-(1+(row.gamma-1)/2*M2**2)) # Eqn 6
            Y_shock = _mean_value(Y_shock)
        else:
            Y_shock = 0
        
        if M2 <= 0.2:
            K1 = 1
        else:
            K1 = 1-1.25*(M2-0.2)
        K2 = (M1/M2)**2 
        Kp = 1-K2*(1-K1)
        
        if M2>1:
            CFM = 1+60*(M2-1)**2    # Eqn 9 
        else:
            CFM = 0
        
        Yp = 0.914 * (2/3*Yp_amdc*Kp + Y_shock) # Eqn 8 Subsonic regime 
        if M2>1:
            Yp = Yp*CFM
        
        f_ar = (1-0.25*np.sqrt(2-h/c)) / (h/c) if h/c<=2 else 1/(h/c)
        alpham = np.arctan(0.5*(np.tan(alpha1_rad) - np.tan(alpha2_rad)))
        Cl_sc = 2*(np.tan(alpha1_rad)+np.tan(alpha2_rad))*np.cos(alpham)
        Ys_amdc = 0.0334 *f_ar *np.cos(alpha2_rad)/np.cos(beta1_rad) * (Cl_sc)**2 * np.cos(alpha2_rad)**2 / np.cos(alpham)**3
        # Secondary Loss
        if h>0: # h is calculated from tip clearance. When h is 0 there is no tip clearance  
            K3 = 1/(h/(b))**2       # Fig 13, it's actually bx in the picture which is the axial chord; h is 0 this causes nan
            Ks = 1-K3*(1-Kp)        # Eqn 15
            Ys = 1.2*Ys_amdc*Ks     # Eqn 16
        else:
            K3 = 0
            Ks = 0 
            Ys = 0 
        
        # Trailing Edge
        if np.abs(alpha1_deg-alpha2_deg)<5:
            delta_phi2 = self.data['Fig14_Impulse'](float(row.te_pitch*row.pitch / row.throat))
        else:
            delta_phi2 = self.data['Fig14_Axial_Entry'](float(row.te_pitch*row.pitch / row.throat))
        
        Ytet = (1-(row.gamma-1)/2 * M2**2 * (1/(1-delta_phi2)-1)) **(-row.gamma/(row.gamma-1)) - 1 # Equation 18
        Ytet = Ytet / (1-(1+(row.gamma-1)/2*M2**2)**(-row.gamma/(row.gamma-1)))
        
        # Tip Clearance
        if h > 0:
            kprime = row.tip_clearance/(3)**0.42 # Number of seals 
            Ytc = 0.37*c/h * (kprime/c)**0.78 * Cl_sc**2 * np.cos(alpha2_rad)**2 / np.cos(alpham)**3
        else:
            Ytc = 0 
            
        if Rec <= 2E5:
            f_re = (Rec/2E5)**-0.4 
        elif Rec<1E6:
            f_re = 1 
        else:
            f_re = (Rec/1E6)**-0.2
        
        Yt = Yp*f_re + Ys + Ytet + Ytc 
        return Yt
        
