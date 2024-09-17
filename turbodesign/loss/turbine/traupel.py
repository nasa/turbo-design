import pickle, os
from typing import Dict
from ...bladerow import BladeRow, sutherland
from ...lossinterp import LossInterp
from ...enums import RowType, LossType
import numpy as np
import pathlib
from ..losstype import LossBaseClass
import requests

class Traupel(LossBaseClass):
    def __init__(self):
        super().__init__(LossType.Enthalpy)
        path = pathlib.Path(os.path.join(os.environ['TD3_HOME'],"traupel"+".pkl"))

        if not path.exists():
            url = "https://github.com/nasa/turbo-design/raw/main/references/Turbines/Traupel/traupel.pkl"
            response = requests.get(url, stream=True)
            with open(path.absolute(), mode="wb") as file:
                for chunk in response.iter_content(chunk_size=10 * 1024):
                    file.write(chunk)   
        
        with open(path.absolute(),'rb') as f:
            self.data = pickle.load(f) # type: ignore
        
    def __call__(self,row:BladeRow, upstream:BladeRow) -> float:
        """Enthalpy loss is computed for the entire stage. 

        Args:
            upstream (BladeRow): Stator Row
            row (BladeRow): Rotor Row

        Returns:
            float: Efficiency
        """
        
        alpha1 = 90-np.degrees(upstream.alpha1.mean())
        alpha2 = 90-np.degrees(upstream.alpha2.mean())
        beta2 = 90 - np.degrees(row.beta1.mean())
        beta3 = 90 - np.degrees(row.beta2.mean())
            
        g = upstream.pitch # G is the pitch 
        h_stator = upstream.r[-1] - upstream.r[0]
        h_rotor = row.r[-1] - row.r[0]

        if row.row_type == RowType.Rotor:
            turning = np.abs(np.degrees(upstream.beta2-row.beta2).mean())
            F = self.data['Fig06']((upstream.W/row.W).mean(),turning) # Inlet velocity
        else:
            turning = np.abs(np.degrees(upstream.alpha2-row.alpha2).mean())
            F = self.data['Fig06']((upstream.V/row.V).mean(),turning) # Inlet velocity

        H = self.data['Fig07'](alpha1-beta2,alpha2-beta3)
        
        zeta_s = F*g/h_stator  # (h1-h1s)/(0.5*c1s**2) # no idea what h1s or h2s is
        zeta_r = F*g/h_rotor # (h2-h2s)/(0.5*w2s**2)
        x_p_stator = self.data['Fig01'](alpha1,alpha2) # not sure if this is the right figure
        x_p_rotor = self.data['Fig01'](beta2,beta3) # not sure if this is the right figure
        zeta_p_stator = self.data['Fig02'](alpha1,alpha2)
        x_m_stator = self.data['Fig03_0'](upstream.M)
        zeta_p_rotor = self.data['Fig02'](beta2,beta3)
        x_m_rotor = self.data['Fig03_0'](row.M_rel)
        
        
        e_te = upstream.te_pitch * g
        o = upstream.throat 
        ssen_alpha2 = e_te/o # Thickness of Trailing edge divide by throat 
        ssen_beta2 = row.te_pitch*g / row.throat
        
        x_delta_stator = self.data['Fig05'](ssen_alpha2,alpha2)
        zeta_delta_stator = self.data['Fig04'](ssen_alpha2,alpha2)
        x_delta_rotor = self.data['Fig05'](ssen_beta2,beta3)
        zeta_delta_rotor = self.data['Fig04'](ssen_beta2,beta3)
        
        Dm = 2* (upstream.r[-1] + upstream.r[0])/2 # Is this the mean diameter? I dont know
        zeta_f = 0.5 * (h_stator/Dm)**2
        
        zeta_pr_stator = zeta_p_stator * x_p_stator * x_m_stator * x_delta_stator + zeta_delta_stator + zeta_f
        
        Dm = 2* (row.r[-1] + row.r[0])/2 # Is this the mean diameter? I dont know
        zeta_f = 0.5 * (h_rotor/Dm)**2
        
        zeta_pr_rotor = zeta_p_rotor * x_p_rotor * x_m_rotor * x_delta_rotor + zeta_delta_rotor + zeta_f
        
        if row.row_type == RowType.Stator:
            zeta_cl = 0 
        else:
            zeta_cl = self.data['Fig08'](row.tip_clearance) # For simplicity assume unshrouded blade 
            
        zeta_z = 0 # Do not factor this in, a bit complicated
        # 1 - (internal) - (external) 
        zeta_v = 0 
        zeta_off = 0 
        eta_stator = 1- (zeta_pr_stator + zeta_s + 0 + zeta_z) - (zeta_r+zeta_v) - zeta_off # Presentation slide 9
        eta_rotor = 1 - (zeta_pr_rotor + zeta_r + zeta_cl + zeta_z) - (zeta_r+zeta_v) - zeta_off
        return eta_stator+eta_rotor
        
        
        
        