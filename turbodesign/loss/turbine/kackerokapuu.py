import pickle, os
from typing import Dict
from ...bladerow import BladeRow, sutherland
from ...lossinterp import LossInterp
from ...enums import RowType, LossType
import numpy as np
import numpy.typing as npt
import pathlib
from ..losstype import LossBaseClass
import requests


def _mean_value(value):
    """Return the mean of an array-like as a Python float."""
    return float(np.asarray(value).mean())

class KackerOkapuu(LossBaseClass):
    UseCFM:bool = False
    
    def __init__(self,UseCFM:bool=False):
        """KackerOkapuu model is an improvement to the Ainley Mathieson model. 
        
        Limitations:
            - Doesn't factor incidence loss 
            - For steam turbines and impulse turbines
        
        Args:
            UseCFM (bool): Factor in supersonic drag rise. Authors state in AMDC loss that this is not accurate. It's a multplier to pressure loss. 
            
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
        self.UseCFM = UseCFM
        
    
    def __call__(self,row:BladeRow, upstream:BladeRow) -> npt.NDArray:
        """Kacker Okapuu is an updated version of Ainley Mathieson and Dunham Came. This tool uses the pressure loss definition. 

        Note: 
            All equation numbers are from the Kacker Okapuu paper
        
        Reference:
            Kacker, S. C., and U. Okapuu. "A mean line prediction method for axial flow turbine efficiency." (1982): 111-119.
        
        Args:
            row (BladeRow): Blade row being evaluated.
            upstream (BladeRow): Upstream blade row providing inlet conditions.

        Returns:
            numpy.ndarray: Pressure loss coefficient array matching ``row.r``.
        """
        # Get the Inlet incoming mach number relative to the blade
        c = row.chord
        b = row.axial_chord
        if row.row_type == RowType.Stator:
            beta1_rad = np.abs(_mean_value(np.radians(row.beta1_metal))) # Metal angle from fig 3
            alpha1_rad = np.abs(_mean_value(row.alpha1)) # Flow angle
            alpha2_rad = np.abs(_mean_value(row.alpha2)) # Flow angle at exit which is metal angle
            beta2_rad = alpha2_rad                  
            alpham_rad = (alpha1_rad + alpha2_rad)*0.5
            M1 = _mean_value(upstream.M)
            M2 = _mean_value(row.M)
            h = 0
            Rec = _mean_value(row.V*row.rho*row.chord / sutherland(row.T))
            
            beta1_deg = np.abs(np.degrees(beta1_rad))
            alpha2_deg = np.abs(np.degrees(alpha2_rad))
            Yp_beta0 = self.data['Fig01_beta0'](float(row.pitch_to_chord), alpha2_deg)  # when beta1 = 0 
            Yp_beta1_alpha2 = self.data['Fig02'](float(row.pitch_to_chord), alpha2_deg) # When beta1 = alpha2
            t_max_c = self.data['Fig04'](beta1_deg+alpha2_deg)
        else:
            h = row.tip_clearance * (row.r[-1]-row.r[0])
            alpha1_rad = np.abs(_mean_value(row.beta1))
            beta1_rad = np.abs(_mean_value(np.radians(row.beta1_metal))) # metal angles are stored as degrees 
            beta2_rad = np.abs(_mean_value(np.radians(row.beta2_metal)))
            alpha2_rad = np.abs(_mean_value(row.beta2))
            alpham_rad = (beta1_rad + beta2_rad)*0.5
            M1 = _mean_value(upstream.M_rel) 
            M2 = _mean_value(row.M_rel)
            Rec = _mean_value(row.W*row.rho*row.chord / sutherland(row.T))

            beta1_deg = np.abs(np.degrees(beta1_rad))
            alpha2_deg = np.abs(np.degrees(beta2_rad))
            Yp_beta0 = self.data['Fig01_beta0'](float(row.pitch_to_chord), alpha2_deg)  # when beta1 = 0 
            Yp_beta1_alpha2 = self.data['Fig02'](float(row.pitch_to_chord), alpha2_deg) # When beta1 = alpha2
            t_max_c = self.data['Fig04'](beta1_deg+alpha2_deg)
        
        ratio = beta1_rad / alpha2_rad
        Yp_amdc = (Yp_beta0 + np.abs(ratio) * ratio * (Yp_beta1_alpha2-Yp_beta0)) * ((t_max_c)/0.2)**(ratio) # Eqn 2, AMDC = Ainley Mathieson Dunham Came
        
        # Shock Loss
        if M1>=0.4: # You'll have imaginary numbers if M1<0.4
            dP_q1_hub = 0.75*(M1-0.4)**1.75 # Eqn 4, this is at the hub
            dP_q1_shock = row.r[-1]/row.r[0] * dP_q1_hub # Eqn 5
            # Eqn 6: convert the inlet-q shock loss to the exit dynamic-head basis using the
            # isentropic dynamic-head fraction q/P0 = 1 - (1 + (gamma-1)/2 M^2)^(-gamma/(gamma-1)).
            # The (-gamma/(gamma-1)) exponent was dropped in the original port and is restored
            # here (cf. the trailing-edge denominator below, which uses it correctly). The
            # spurious static-pressure ratio P1/P2 is removed: the q/P0 terms already carry the
            # per-row total-pressure normalization, so the conversion is the q1/q2 ratio.
            q1_frac = 1-(1+(upstream.gamma-1)/2*M1**2)**(-upstream.gamma/(upstream.gamma-1))
            q2_frac = 1-(1+(row.gamma-1)/2*M2**2)**(-row.gamma/(row.gamma-1))
            Y_shock = dP_q1_shock * q1_frac/q2_frac # Eqn 6
            Y_shock = _mean_value(Y_shock)
        else:
            Y_shock = 0
        
        if M2 <= 0.2:
            K1 = 1
        else:
            K1 = 1-1.25*(M2-0.2)
        K2 = (M1/M2)**2 
        Kp = 1-K2*(1-K1)
        
        if (M2>1) and (self.UseCFM is True):
            CFM = 1+60*(M2-1)**2    # Eqn 9 
        else:
            CFM = 1
        
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
        if np.abs(beta1_deg-np.degrees(beta2_rad))<5: # impulse turbine the inlet and exit angles are the same
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
        return Yt+row.r*0
        
