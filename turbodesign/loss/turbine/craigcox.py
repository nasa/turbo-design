import pickle, os
from typing import Dict
import numpy.typing as npt
from ...bladerow import BladeRow, sutherland
from ...lossinterp import LossInterp
from ...enums import RowType, LossType
import numpy as np
import pathlib
from ..losstype import LossBaseClass
import requests

class CraigCox(LossBaseClass):
    
    def __init__(self):
        """Craig and Cox used to estimate loss for subsonic higher mach number flows in axial turbines. Assumes suction side has low convex surface close to throat. 

        Reference:
            Craig, H. R. M., and H. J. A. Cox. "Performance estimation of axial flow turbines." Proceedings of the Institution of Mechanical Engineers 185.1 (1970): 407-424.
            
        """
        super().__init__(LossType.Enthalpy)
        path = pathlib.Path(os.path.join(os.environ['TD3_HOME'],"craigcox"+".pkl"))
        
        if not path.exists():
            # Download data from Github 
            url = "https://github.com/nasa/turbo-design/raw/main/references/Turbines/CraigCox/craigcox.pkl"
            response = requests.get(url, stream=True)
            with open(path.absolute(), mode="wb") as file:
                for chunk in response.iter_content(chunk_size=10 * 1024):
                    file.write(chunk)
        
        with open(path.absolute(),'rb') as f:
            self.data = pickle.load(f) # type: ignore
        
        self.C = 1/(200*32.2*778.16) # https://www.sciencedirect.com/science/article/pii/S2666202721000574 
    
    
    def __call__(self,row:BladeRow, upstream:BladeRow) -> npt.NDArray:
        """Craig and Cox uses the enthalpy definition of loss to calculate the loss of a turbine stage. 
        
        Note: 
            Losses are organized as Group 1 which include profile losses and secondary flows. 
            Group 2 losses include rotor tip leakage, balance losses, guide gland losses, lacing wire losses.
            
            All equation numbers are from the craig cox paper
            Craig, H. R. M., and H. J. A. Cox. "Performance estimation of axial flow turbines." Proceedings of the Institution of Mechanical Engineers 185.1 (1970): 407-424.
        
        Equations:
            Eta_t = (Work done - Group 2 losses) / (Work done + Group 1 Losses) 
            
            Group1 Losses = (Xp + Xs + Xa)*C1**2/(200gJ) + (Xp + Xs + Xa)*W2**2/(200gJ)
            Group1 Losses = Stator Component + Rotor Component where C1 and W2 are exit velocities for stator and rotor 
            
            i+i_stall = (i+i_stall)_basic + (\delta i + stall)_sb + (\delta i + stall)_cb 

            (i+i_stall)_basic from Figure 11
            (delta incidence + stall)_sb + (delta incidence + stall)_cb from Figure 12

            Profile Loss
            Xp = x_pb N_pr N_pi N_pt + (\delta x_p)_t + (\delta x_p)_s/e + (\delta x_p)_m

            - x_pb from Figure 5 but use Figure 4 to calculate Fl. Fl*x*s/b is the x axis for Figure 5 
            - N_pr from figure 3 
            - N_pi from Figure 10
            - N_pt from Figure 6
            - (\delta x_p)_t 
            - (\delta x_p)_s/e from Figure 9
            - (\delta x_p)_m from Figure 8

        Args:
            row (BladeRow): Downstream blade row being evaluated.
            upstream (BladeRow): Upstream blade row providing inlet conditions.

        Returns:
            numpy.ndarray | int: Stage efficiency; returns 0 for stators, spanwise array for rotors.
        """
        if row.row_type == RowType.Stator:
            # Craig-Cox assembles a combined stage efficiency at the rotor row; the stator
            # branch contributes nothing on its own. Return a spanwise zero array (not a bare
            # int 0) so the return type matches the rotor branch and downstream array ops.
            return 0.0 * row.r
        else:
            V_inlet = upstream.V.mean()
            V = row.W.mean()
        
        def lookup_constants(currentRow:BladeRow):
            # Contraction Ratio find s/b and use Fig07
            s_b = currentRow.pitch/currentRow.camber
            if currentRow.row_type == RowType.Stator:
                V = currentRow.V.mean()
                # See weird velocity triangle in Figure 1
                inlet_flow_angle = 90-np.abs(np.degrees(currentRow.alpha1).mean())
                outlet_flow_angle =90-np.abs(np.degrees(currentRow.alpha2).mean())
                incidence_angle = np.degrees(currentRow.alpha1 - currentRow.beta1_metal).mean()
                M_out = currentRow.M.mean()
            else:
                V = currentRow.W.mean()
                inlet_flow_angle = 90-np.abs(np.degrees(currentRow.beta1).mean())
                outlet_flow_angle =90-np.abs(np.degrees(currentRow.beta2).mean())
                incidence_angle = np.degrees(currentRow.beta1 - currentRow.beta1_metal).mean()
                M_out = currentRow.M_rel.mean()
                
            Re =  V * currentRow.rho*currentRow.chord / sutherland(currentRow.T.mean())
            Re = Re.mean()
            
            if currentRow.beta1_fixed:
                blade_inlet_angle = 90-np.abs(np.degrees(currentRow.beta1_metal.mean())) # Alpha
            else:
                if currentRow.row_type == RowType.Stator:
                    blade_inlet_angle = 90-np.abs(np.degrees(currentRow.alpha1.mean())) # Alpha1
                else:
                    blade_inlet_angle = 90-np.abs(np.degrees(currentRow.beta1.mean())) # beta1
                
            te = currentRow.te_pitch * currentRow.pitch
            e_s = 0.3 # Pitch to back radius ratio, assumed. lower value = less loss
            if currentRow.beta1_fixed:
                imin = 90-currentRow.beta1_metal.mean()
            else:
                imin = 90-currentRow.beta1.mean() # Incidence required for minimum loss
            
            asin_os = np.degrees(np.arcsin(currentRow.throat/currentRow.pitch))
            
            N_pr = self.data['Fig03'](float(Re), 0.05) # use a good finish for the geometry
            if (inlet_flow_angle-imin < 10):
                Fl = 13
            else:
                Fl = self.data['Fig04'](float(outlet_flow_angle), float(inlet_flow_angle-imin))
            
            x = 1-np.sin(np.radians(outlet_flow_angle))/np.sin(np.radians(inlet_flow_angle))
            contraction_ratio = self.data['Fig07'](float(x), float(s_b)) # contraction ratio

            # Fig 5's ordinate is the loss *parameter* x_pb*(s/b)*sin(beta), not x_pb itself
            # (axis label: "Basic Profile Loss Parameter Xp(s/b)sin(beta)"). Recover x_pb by
            # dividing out (s/b)*sin(beta), using the low-speed exit flow angle for beta per the
            # figure's note. Eqn 10 then consumes x_pb directly.
            xpb_s_b_sinB = self.data['Fig05'](float(Fl*s_b), float(contraction_ratio))
            X_pb = xpb_s_b_sinB / (s_b * np.sin(np.radians(outlet_flow_angle)))
            delta_X_pt = self.data['Fig06_delta_Xpt'](float(currentRow.te_pitch))
            N_pt = self.data['Fig06_Npt'](float(currentRow.te_pitch), float(outlet_flow_angle))
            delta_Xpm = self.data['Fig08'](float(M_out), float(np.degrees(np.arcsin((currentRow.throat+te)/currentRow.pitch))))
            delta_Xp_se = self.data['Fig09'](float(e_s), float(M_out)) 
            
            Fi = self.data['Fig15'](float(blade_inlet_angle), float(s_b))
            # Incidence Effects 
            if currentRow.beta1_fixed:
                if incidence_angle>0: # Positive incidence
                    stall_incidence_angle = self.data['Fig11'](float(currentRow.beta1.mean()), float(asin_os))
                    
                    incidence_ratio = (incidence_angle - imin)/(stall_incidence_angle-imin)

                    i_plus_istall_sb = self.data['Fig12_sb'](float(s_b), float(asin_os))
                    i_plus_istall_cor = self.data['Fig12_cr'](float(contraction_ratio), float(asin_os))
                    
                    if blade_inlet_angle<=90:
                        i_plus_istall_basic = self.data['Fig11'](float(inlet_flow_angle), float(asin_os))
                        i_plus_istall = i_plus_istall_basic + i_plus_istall_sb + i_plus_istall_cor # Eqn 5 
                    else:
                        i_plus_istall_basic = self.data["Fig14_i+istall"](float(blade_inlet_angle), float(asin_os))
                        i_plus_istall = i_plus_istall_basic + (1-(blade_inlet_angle-90)/(90-asin_os))*(i_plus_istall_sb + i_plus_istall_cor) # Eqn 7 
                else:
                    i_minus_istall_sb = self.data['Fig13'](float(s_b), float(asin_os))

                    if blade_inlet_angle<=90: 
                        i_minus_istall_basic = self.data['Fig13_alpha1'](float(s_b), float(asin_os))
                        i_minus_istall = i_minus_istall_basic + i_minus_istall_sb # Eqn 6
                    else:
                        i_minus_istall_basic = self.data["Fig14_i-istall"](float(blade_inlet_angle), float(asin_os))
                        i_minus_istall = i_minus_istall_basic + (1-(blade_inlet_angle - 90)/(90-asin_os)) * i_minus_istall_sb   # Eqn 8 
                        
                imin = (i_plus_istall + Fi * (i_minus_istall))/(1+Fi) # type: ignore # Eqn 9
                N_pi = self.data['Fig10'](float(imin), float(incidence_ratio))
            else:
                N_pi = 1 # No effect
            
            Xp = X_pb*N_pr*N_pi*N_pt + delta_X_pt + delta_Xp_se + delta_Xpm # Eqn 10
            
            # Secondary Loss 
            Ns_hb = self.data['Fig17'](float(1/currentRow.aspect_ratio))
            x_sb = self.data['Fig18'](float((V_inlet/V)**2), float(s_b*Fl))
                
            # Craig & Cox apply the same Reynolds/surface-finish factor to the secondary loss
            # as to the profile loss (Fig 3, N_pr); reuse it here as N_sr (Eqn 11) rather than 1.
            Nsr = N_pr
            Xs = Nsr*Ns_hb*x_sb 
            # Annulus Loss Factor
            Xa = 0            
            return Xp, Xs, Xa 
        

        Xp1,Xs1,Xa1 = lookup_constants(upstream)
        Xp2,Xs2,Xa2 = lookup_constants(row)
        Group1_Loss = (Xp1 + Xs1 + Xa1) * V_inlet**2 *3.28**2 * self.C + (Xp2 + Xs2 + Xa2*V_inlet**2/V**2) *  V**2 *3.28**2 * self.C # type: ignore
                # Eqn 4, convert V from m^2/s^2to ft^2/s^2
                
        # According to Equation 3, Group 1 loss is an enthalpy loss Cp*T0 Btu/lb. Need to convert to Pressure Loss
        
        # Btu/lbf to J/Kg
       
        T0_Loss = Group1_Loss * 2326 / row.Cp # Eqn 3 in Kelvin
        T0_T = (row.T0/row.T).mean() 
        T02 = row.T0.mean()-T0_Loss # P02 Changes 
        
        # According to Equation 3, Group 1 loss is an enthalpy loss Cp*T0. Need to convert to Pressure Loss
        eta_total = (upstream.T0.mean() - row.T0.mean())/(upstream.T0.mean()-(row.T0.mean()-T0_Loss))
        return eta_total + row.r*0
        
