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
        
    def __call__(self,row:BladeRow, upstream:BladeRow) -> npt.NDArray:
        """Compute Traupel stage enthalpy efficiency from an upstream/downstream pair.

        Args:
            row (BladeRow): Blade row being evaluated (stator or rotor).
            upstream (BladeRow): Upstream blade row supplying inlet conditions.

        Returns:
            numpy.ndarray: Spanwise efficiency array matching ``row.r``.
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
            F = self.data['Fig06'](float((upstream.W/row.W).mean()), float(turning)) # Inlet velocity
        else:
            turning = np.abs(np.degrees(upstream.alpha2-row.alpha2).mean())
            F = self.data['Fig06'](float((upstream.V/row.V).mean()), float(turning)) # Inlet velocity

        H = self.data['Fig07'](float(alpha1-beta2), float(alpha2-beta3))
        
        zeta_s = F*g/h_stator  # Stator loss factor scaled by pitch-to-span
        zeta_r = F*g/h_rotor   # Rotor loss factor scaled by pitch-to-span
        x_p_stator = self.data['Fig01'](float(alpha1), float(alpha2))
        x_p_rotor = self.data['Fig01'](float(beta2), float(beta3))
        zeta_p_stator = self.data['Fig02'](float(alpha1), float(alpha2))
        x_m_stator = self.data['Fig03_0'](float(np.mean(upstream.M)))
        zeta_p_rotor = self.data['Fig02'](float(beta2), float(beta3))
        x_m_rotor = self.data['Fig03_0'](float(np.mean(row.M_rel)))
        
        
        e_te = upstream.te_pitch * g
        o = upstream.throat 
        ssen_alpha2 = e_te/o # Thickness of Trailing edge divide by throat 
        ssen_beta2 = row.te_pitch*g / row.throat
        
        x_delta_stator = self.data['Fig05'](float(ssen_alpha2), float(alpha2))
        zeta_delta_stator = self.data['Fig04'](float(ssen_alpha2), float(alpha2))
        x_delta_rotor = self.data['Fig05'](float(ssen_beta2), float(beta3))
        zeta_delta_rotor = self.data['Fig04'](float(ssen_beta2), float(beta3))
        
        Dm = 2* (upstream.r[-1] + upstream.r[0])/2  # Mean diameter used for annulus friction
        zeta_f = 0.5 * (h_stator/Dm)**2
        
        zeta_pr_stator = zeta_p_stator * x_p_stator * x_m_stator * x_delta_stator + zeta_delta_stator + zeta_f
        
        Dm = 2* (row.r[-1] + row.r[0])/2  # Mean diameter used for annulus friction
        zeta_f = 0.5 * (h_rotor/Dm)**2
        
        zeta_pr_rotor = zeta_p_rotor * x_p_rotor * x_m_rotor * x_delta_rotor + zeta_delta_rotor + zeta_f
        
        if row.row_type == RowType.Stator:
            zeta_cl = 0 
        else:
            zeta_cl = self.data['Fig08'](float(row.tip_clearance))  # Clearance loss for unshrouded blades
            
        zeta_z = 0  # Disk friction loss not modeled
        # 1 - (internal) - (external)
        zeta_v = 0 
        zeta_off = 0 
        eta_stator = 1- (zeta_pr_stator + zeta_s + 0 + zeta_z) - (zeta_r+zeta_v) - zeta_off  # Per Traupel formulation
        eta_rotor = 1 - (zeta_pr_rotor + zeta_r + zeta_cl + zeta_z) - (zeta_r+zeta_v) - zeta_off
        return (eta_stator+eta_rotor) + row.r*0
        
        
        
        
