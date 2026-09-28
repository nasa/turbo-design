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
        """Compute Traupel stage total-to-total efficiency from a stator/rotor pair.

        Experimental: the digitised figures and the kinetic-energy weighting of the row losses are
        not validated against Traupel (1977).

        The stage loss is assembled at the rotor, as CraigCox does, because the turbine solver
        applies the whole stage loss to the rotor Yp and treats stators as lossless.

        Note:
            TurbineSpool's LossType.Enthalpy path (turbine_spool.py, balance_loop) does not
            currently work for any enthalpy model: it converts this array with float(), and
            its search for Yp does not depend on Yp. Calling this model directly is fine.

        Args:
            row (BladeRow): Blade row being evaluated. Stators return zeros.
            upstream (BladeRow): Upstream (stator) row supplying inlet conditions.

        Returns:
            numpy.ndarray: Spanwise efficiency array matching ``row.r``.
        """
        if row.row_type != RowType.Rotor:
            return 0.0 * row.r

        alpha1 = 90-np.degrees(upstream.alpha1.mean())
        alpha2 = 90-np.degrees(upstream.alpha2.mean())
        beta2 = 90 - np.degrees(row.beta1.mean())
        beta3 = 90 - np.degrees(row.beta2.mean())
            
        g = upstream.pitch # Stator pitch
        g_rotor = row.pitch
        h_stator = upstream.r[-1] - upstream.r[0]
        h_rotor = row.r[-1] - row.r[0]

        turning = np.abs(np.degrees(upstream.beta2-row.beta2).mean())
        F = self.data['Fig06'](float((upstream.W/row.W).mean()), float(turning)) # Inlet velocity

        # Fig07 (divergence factor) is not used: whether it is additive or a multiplier is unverified.
        zeta_s = F*g/h_stator  # Stator loss factor scaled by pitch-to-span
        zeta_r = F*g_rotor/h_rotor   # Rotor loss factor scaled by pitch-to-span
        x_p_stator = self.data['Fig01'](float(alpha1), float(alpha2))
        x_p_rotor = self.data['Fig01'](float(beta2), float(beta3))
        zeta_p_stator = self.data['Fig02'](float(alpha1), float(alpha2))
        x_m_stator = self.data['Fig03_0'](float(np.mean(upstream.M)))
        zeta_p_rotor = self.data['Fig02'](float(beta2), float(beta3))
        x_m_rotor = self.data['Fig03_0'](float(np.mean(row.M_rel)))
        
        
        e_te = upstream.te_pitch * g
        o = upstream.throat 
        ssen_alpha2 = e_te/o # Thickness of Trailing edge divide by throat 
        ssen_beta2 = row.te_pitch*g_rotor / row.throat
        
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
        
        zeta_cl = self.data['Fig08'](float(row.tip_clearance))  # Clearance loss for unshrouded blades

        zeta_z = 0  # Disk friction loss not modeled
        zeta_v = 0  # Ventilation loss not modeled
        zeta_off = 0  # Leaving loss not modeled

        # Each row carries its own secondary loss; clearance applies to the rotor only.
        # Row losses are summed unweighted; Traupel weights them by each row's exit kinetic energy.
        zeta_stator = zeta_pr_stator + zeta_s
        zeta_rotor = zeta_pr_rotor + zeta_r + zeta_cl
        eta_stage = 1.0 - (zeta_stator + zeta_rotor + zeta_z + zeta_v + zeta_off)
        return eta_stage + row.r*0
