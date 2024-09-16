import math
from ...lossinterp import LossInterp
from ...enums import RowType, LossType
from typing import Any, Callable, Dict, List, Union 
from ...bladerow import BladeRow
from scipy.stats import linregress
import numpy as np 
from ..losstype import LossBaseClass


class TD2(LossBaseClass):

    
    def __init__(self):
        super().__init__(LossType.Pressure)

        
    @property
    def LossType(self):
        return self._loss_type
    
    def __call__(self,row:BladeRow, upstream:BladeRow) -> float:
        """TD2-2_Manual equations 12a and 12b. This is the loss equation used inside original TD2.
        #! Wrong very wrong from the book definition
        Total Pressure Loss (Y) = (P_in - P_ex) / (P_ex - P_ex,static)
        Where P_in and P_ex are total quantities 
        
        Authors Comments:
            Given the assumptions, the use of this loss correlation should be limited to starting a solution. 

        Assumptions:
            1. Rotor and stator loss coefficients were equal when their relative design requirements are identical 
            2. Stage reaction at meanline is 50%
            3. Axial velocity was constant through the stage
            4. Stator exit Mach number is 0.8 

        Args:
            beta_in (float): Inlet flow angle in degrees. 
            beta_ex (float): Exit flow angle in degrees.
            V_ratio (float): Ratio if inlet velocity magnitude with relative exit velocity magnitude

        Returns:
            float: Total Pressure loss coefficient
        """
        beta_in = row.beta1.mean() # absolute flow angle entering the blade row
        beta_ex = row.beta2.mean()  # flow angle leaving the blade row
        if row.row_type == RowType.Stator:
            V_ratio = row.V.mean() / upstream.V.mean()    # Vin/Vex Equation 12a and 12b. Relative to Stator
        elif row.row_type == RowType.Rotor:
            V_ratio = row.W.mean() /upstream.W.mean()     # Vin/Vex Equation 12a and 12b. Relative to Rotor 

        a = [0.055, 0.15, 0.6, 0.6, 0.8, 0.03, 0.157255, 3.6] # Coefficients from pw3hp1.json
        # mu should be in 
        if V_ratio < a[2]: 
            Y = abs(math.tan(beta_in) - math.tan(beta_ex)) / (a[3] + a[4]*math.cos(beta_ex)) * (a[5] + a[6]*V_ratio**a[7])
        else:
            Y = abs(math.tan(beta_in) - math.tan(beta_ex)) / (a[3] + a[4]*math.cos(beta_ex)) * (a[0] + a[1]*(V_ratio - a[2]))

        return Y
    
class TD2_Reynolds_Correction(LossBaseClass):
    
    def __init__(self):
        super().__init__(LossType.Pressure)

        self.TD2 = TD2()
    
    @property
    def LossType(self):
        return self._loss_type
    
    def __call__(self,upstream:BladeRow, row:BladeRow) -> float:
        """TD2_Pressure_Loss with Reynolds Correction Factor. This is from NASA SP-290 (Vol.1, p.62)
    
        The correlations come from td2-2.f Line 2771 
        WYECOR(I)=WYECOR(I)*(0.35+0.65*18.21)/(0.35+0.65*(FLWP/VISC/RST(MEAN))**(0.2))
        
        I have assumed that 18.21 is some reference reynolds number. There is very little documentation on what these numbers are. The entire fortran code is frustrating and probably should never have happened in the first place. 

        Args:
            row (BladeRow): Current Blade Row
        
        Returns:
            float: Total Pressure loss coefficient
        """
        Y = self.TD2(upstream,row)
        A = 0.35
        B = 0.65
        Y = Y * (A+B*18.21)/(0.35+0.65*row.massflow/(row.mu* row.r.mean()))
        row.Yp = Y
        return Y
    
# def Soderberg(upstream:BladeRow,row:BladeRow) -> float: 
#     """Soderberg Loss for axial machines. Takes into account the aspect ratio 

#     Args:
#         upstream (BladeRow): _description_
#         row (BladeRow): _description_

#     Returns:
#         float: Enthalpy Loss Coefficients 
#     """
#     H = row.r.max()-row.r.min()
#     l = row.x[-1]-row.x[0]
#     xi = 0.04+0.06*((row.beta2-row.beta1)/100)**2 
#     if row.row_type == RowType.Stator:
#         mu = sutherland(row.T.mean())
#         Re = row.rho*row.V*(l)/mu
#         enthalpy_loss = (10E5/Re)**0.25 * ((1+xi)*(0.993+0.075*l/H)-1)
#     else:
#         mu = sutherland(row.T.mean())
#         Re = row.rho*row.W*(row.x[-1]-row.x[0])/mu
#         enthalpy_loss = (10E5/Re)**0.25 * ((1+xi)*(0.975+0.075*l/H)-1)
#     return enthalpy_loss

# def AinleyMathieson(upstream:BladeRow,row:BladeRow) -> float:
#     """Ainley Mathieson equation for computing total pressure loss (Yp)

#     Args:
#         upstream (BladeRow): _description_
#         row (BladeRow): _description_
#     """
#     #! Need to extract data from plots 
#     pass 

# def AinleyMathiesonUpdated(upstream:BladeRow,row:BladeRow) -> float:
#     """Updated version of Ainley Mathieson
#     https://www.mdpi.com/2504-186X/7/2/14
#     These derivations are valid for steam turbines.

#     Note:
#         According to the authors, Pressure Loss divided in to 3 components: Profile Loss, secondary loss, and tip/shroud clearance loss

#     Args:
#         upstream (BladeRow): _description_
#         row (BladeRow): _description_

#     Returns:
#         float: _description_
#     """
    
#     pass 
