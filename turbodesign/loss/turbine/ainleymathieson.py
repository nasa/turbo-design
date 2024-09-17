import pickle, os
from typing import Dict
from ...bladerow import BladeRow, sutherland
from ...lossinterp import LossInterp
from ...enums import RowType, LossType
import numpy as np
import pathlib
from ..losstype import LossBaseClass
import requests 

class AinleyMathieson(LossBaseClass):
        
    def __init__(self):
        """Ainley Mathieson loss correlation used for loss estimation for Gas and Steam Turbines. If you are designing an impulse turbine, this could be a good model to start with.
        Use only with axial turbines. Loss correlations used with Reynolds numbers 1E5 to 3E5.
        
        Limitations:
            - Doesn't factor incidence Loss. 
            - Steamturbines
            - Impulse Turbines beta1=alpha2
            - Reynolds Range 
            
        Reference:
            Ainley, D. G., and Mathieson, G. C. R., "A Method of Performance Estimation for Axial-Flow Turbines," British ARC, R & M 2974, 1951.

        """
        super().__init__(LossType.Pressure)
        path = pathlib.Path(os.path.join(os.environ['TD3_HOME'],"ainleymathieson"+".pkl"))

        if not path.exists():
            url = "https://github.com/nasa/turbo-design/raw/main/references/Turbines/AinleyMathieson/ainleymathieson.pkl"
            response = requests.get(url, stream=True)
            with open(path.absolute(), mode="wb") as file:
                for chunk in response.iter_content(chunk_size=10 * 1024):
                    file.write(chunk)
        
        with open(path.absolute(),'rb') as f:
            self.data = pickle.load(f) # type: ignore
            
    def __call__(self,row:BladeRow, upstream:BladeRow) -> float:
        """Ainley Mathieson predicts the pressure loss of a turbine nozzle or rotor. Since these correlations are from Cascade experiments, the user should be familar with the reynolds and mach number requirements for each equation and figure. Using something outside the bounds can give inaccuare approximations of loss. Additionally these correlations were done on unoptimized blades so efficiencies maybe lower than what's attainable. Massflow can also be affected because exit P0 is affected. 
        
        This code will attempt use the correct equations and warn the user if mach number is out of range.
        
        Note:
            alpha: gas flow angle relative to axial direction
            beta: blade angle relative to axial direction 
            
        Args:
            upstream (BladeRow): Upstream blade row
            row (BladeRow): downstream blade row

        Returns:
            float: Pressure Loss at zero incidence 
        """
        Kp = 12300 # Ft / (pdl C); 1 pdl = 0.138254954376 N :(
        At = row.throat
        s = row.pitch
        Area_upstream = np.pi*(upstream.r[-1]**2-upstream.r[0]**2)
        Area_downstream = np.pi*(row.r[-1]**2-row.r[0]**2)
        e = row.r.mean() # mean radius 
        h = (row.r[-1]-row.r[0])
        
        s_c = row.pitch/row.chord
        
        if row.row_type == RowType.Stator:
            k = 0
            alpha1 = -np.abs(row.alpha1)
            beta1 = -np.abs(np.radians(row.beta1_metal))
            alpha2 = np.abs(row.alpha2)
            if row.M<0.5:
                alpha2 = self.data['Fig05'](np.degrees(np.cos(row.throat/row.pitch)))
            elif row.M<0.95:
                X = 0.7
                alpha2 = np.arctan( 
                                   ((1-X*k/h)*(np.cos(beta1)/np.cos(alpha2)))*np.tan(alpha2)
                                   + X*k/h*np.cos(beta1)/np.cos(alpha2)*np.tan(beta1)
                                   )    # Eqn 4
            elif row.M<=1 and row.M>0.5:
                if Area_upstream/Area_downstream <1.02 and Area_upstream/Area_downstream>0.98: # Not flared
                    alpha2 = -np.arccos(At/Area_downstream) # Eqn 2
                else: # Flared so use equation 3
                    At = At/s * (5*Area_downstream + Area_upstream)/6 # Eqn 3
            B = 0.25 # Eqn 6 
        else:   # Rotor
            k = row.tip_clearance
            alpha1 = np.abs(row.beta1)
            beta1 = np.abs(np.radians(row.beta1_metal))
            alpha2 = -np.abs(row.beta2)
            if row.M_rel<0.5:
                alpha2 = self.data['Fig05'](np.degrees(np.cos(row.throat/row.pitch)))
            elif row.M_rel<0.95:
                X = 1.35 # Shrouded Blade
                alpha2 = np.arctan( 
                                   ((1-X*k/h)*(np.cos(beta1)/np.cos(alpha2)))*np.tan(alpha2)
                                   + X*k/h*np.cos(beta1)/np.cos(alpha2)*np.tan(beta1)
                                   )    # Eqn 4
            elif row.M<=1.05 and row.M>=0.95:
                alpha2 = -np.arccos(At/Area_downstream) # Eqn 2
                
            B = 0.5 # Eqn 6 
            
        
        t_c = 0.2; # Impulse turbine 0.15 < t/c < 0.25    
        
        ID = (upstream.r[-1]+row.r[-1])/2; OD = (upstream.r[-1]+row.r[-1])/2
        A1 = np.pi*(upstream.r[-1]**2 - upstream.r[0]**2)*np.cos(beta1)
        A2 = np.pi*(row.r[-1]**2 - row.r[0]**2)*np.cos(alpha2)
        lam = self.data['Fig08']((A2/A1)**2/(1+ID/OD)) # Eqn but using Figure 8
        alpha_m = np.arctan((np.tan(alpha1) - np.tan(alpha2))/2)
        Cl_s_c = 2*(np.tan(alpha1)-np.tan(alpha2))*np.cos(alpha_m)
        # calculated but not used 
        Y_secondary_clearance = (lam + B * k/h) * (Cl_s_c)**2 * (np.cos(alpha2)/np.cos(alpha_m)) # Eqn 6
        
        
        Yp_beta0 = self.data['Fig04a'](s_c, np.degrees(alpha2))
        Yp_beta1_eq_alpha2 = self.data['Fig04b'](s_c, np.degrees(alpha2))
        
        Yp_i0 = (Yp_beta0 + (beta1/alpha2)**2 *(Yp_beta1_eq_alpha2 - Yp_beta0)) *(t_c/0.2)**(-beta1/alpha2) # Fig 4 and Eqn 5
        
        Yt = Yp_i0 + Y_secondary_clearance 
        return Yt   # Profile loss at zero incidence 
        