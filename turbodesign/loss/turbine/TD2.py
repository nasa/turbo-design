import math
from ...lossinterp import LossInterp
from ...enums import RowType, LossType
from typing import Any, Callable, Dict, List, Union 
from ...bladerow import BladeRow
from scipy.stats import linregress
import numpy as np
import numpy.typing as npt
from ..losstype import LossBaseClass


class TD2(LossBaseClass):

    
    def __init__(self):
        super().__init__(LossType.Pressure)

        
    @property
    def LossType(self):
        return self._loss_type
    
    def __call__(self,row:BladeRow, upstream:BladeRow) -> npt.NDArray:
        """TD2-2 manual equations 12a/12b for total pressure loss coefficient.

        The implementation mirrors the original TD2 code path, which differs from
        the textbook definition but preserves legacy behavior. Use primarily for
        initial estimates.

        Assumptions:
            1. Rotor and stator loss coefficients are equal when design requirements match.
            2. Stage reaction at meanline is 50%.
            3. Axial velocity is constant through the stage.
            4. Stator exit Mach number is 0.8.

        Args:
            row (BladeRow): Blade row being evaluated.
            upstream (BladeRow): Upstream blade row supplying inlet conditions.

        Returns:
            numpy.ndarray: Total pressure loss coefficient repeated across ``row.r``.
        """
        beta_in = row.beta1.mean()  # Inlet flow angle at the mean radius
        beta_ex = row.beta2.mean()  # Exit flow angle at the mean radius
        if row.row_type == RowType.Stator:
            V_ratio = row.V.mean() / upstream.V.mean()    # Vin/Vex Equation 12a and 12b. Relative to Stator
        elif row.row_type == RowType.Rotor:
            V_ratio = row.W.mean() /upstream.W.mean()     # Vin/Vex Equation 12a and 12b. Relative to Rotor 

        a = [0.055, 0.15, 0.6, 0.6, 0.8, 0.03, 0.157255, 3.6] # Coefficients from pw3hp1.json
        if V_ratio < a[2]: 
            Y = abs(math.tan(beta_in) - math.tan(beta_ex)) / (a[3] + a[4]*math.cos(beta_ex)) * (a[5] + a[6]*V_ratio**a[7])
        else:
            Y = abs(math.tan(beta_in) - math.tan(beta_ex)) / (a[3] + a[4]*math.cos(beta_ex)) * (a[0] + a[1]*(V_ratio - a[2]))

        return Y+row.r*0
    
class TD2_Reynolds_Correction(LossBaseClass):
    
    def __init__(self):
        super().__init__(LossType.Pressure)

        self.TD2 = TD2()
    
    @property
    def LossType(self):
        return self._loss_type
    
    def __call__(self,row:BladeRow, upstream:BladeRow) -> npt.NDArray:
        """Apply TD2 Reynolds correction (NASA SP-290 Vol.1, p.62).

        The correction follows td2-2.f line 2771:
        WYECOR = WYECOR*(0.35+0.65*18.21)/(0.35+0.65*(FLWP/VISC/RST(MEAN))**0.2)

        Args:
            row (BladeRow): Blade row receiving the correction.
            upstream (BladeRow): Upstream blade row supplying inlet conditions.

        Returns:
            numpy.ndarray: Reynolds-corrected total pressure loss coefficient.
        """
        Y = self.TD2(row, upstream)
        # Reynolds group massflow/(mu * r_mean) is dimensionless and corresponds to
        # FLWP/VISC/RST(MEAN) in td2-2.f. The legacy 0.2 turbulent exponent was dropped
        # in the original Python port and is restored here.
        Re_group = row.massflow / (row.mu * row.r.mean())
        Y = Y * (0.35 + 0.65*18.21) / (0.35 + 0.65*Re_group**0.2)
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
