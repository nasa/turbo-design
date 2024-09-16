# Read in the excel data
# Create bspline surface with points
import sys, os, pathlib
sys.path.insert(0,"../../../")
from turbodesign import LossInterp
import pickle

if __name__=="__main__": 
    data = {
            'Fig03':LossInterp("Fig03_profile_loss_ratio.csv",
                                        xlabel="Reynolds Number",
                                        ylabel="Profile Loss Ratio",
                                        clabel="Finish", logx10=True),
            
            'Fig04':LossInterp("Fig04_Lift_Parameter_FL.csv",
                                        xlabel="Outlet Flow Angle [Deg]",
                                        ylabel="Lift Parameter [Fl]",
                                        clabel="Fluid Inlet Angle at minimum loss"),
            
            'Fig05':LossInterp("Fig05_base_profile_loss.csv",
                                        xlabel="Modified Lift Coefficient",
                                        ylabel="Basic Profile Loss Parameter Xp(s/b)sin(beta)",
                                        clabel="Contraction Ratio"),
            
            'Fig06_delta_Xpt':LossInterp("Fig06_TE_Thickness_Loss_Increment.csv",
                                        xlabel="te/s",
                                        ylabel="Loss Increment Xpt"),
            
            'Fig06_Npt':LossInterp("Fig06_TE_Thickness_Profile_Loss.csv",
                               xlabel="te/s",
                               ylabel="Profile Loss Ratio Npt"),
            
            'Fig07':LossInterp("Fig07_contraction_ratio.csv",
                        xlabel='1-sin(alpha_out)/sin(alpha_in)',
                        ylabel='Contraction Ratio',
                        clabel='s/b'),
            
            'Fig08':LossInterp("Fig08_Profile_Loss_Increment.csv",
                        xlabel='Outlet Isentropic Mach Number',
                        ylabel='Profile Loss Increment',
                        clabel='arcsin((o+te)/s)'),
            
            'Fig09':LossInterp("Fig09_Profile_Loss_Increment.csv",
                        xlabel='Blade Pitch to Back Radius Ratio',
                        ylabel='Profile Loss Increment xp_se',
                        clabel='outlet isentropic mach'),
            
            'Fig10':LossInterp("Fig10_incidence_loss.csv",
                        xlabel='Incidence Ratio (i-imin)/(istall-imin)',
                        ylabel='Profile Loss Ratio Npi'),
            
            'Fig11':LossInterp("Fig11_Basic_Positive_stalling_incidence.csv",
                        xlabel='Blade Inlet Angle',
                        ylabel='Basic Positive Stalling Incidence (i+stall)',
                        clabel='asin(o/s)'),
            
            'Fig12_sb':LossInterp("Fig12_Incidence_Correction_i+istall.csv",
                        xlabel='Pitch to backbone length ratio s/b',
                        ylabel='incidence correction i-istall_sb',
                        clabel='asin(o/s)'),
            
            'Fig12_cr':LossInterp("Fig12_Incidence_Correction_i+istall_cr.csv",
                        xlabel='Blade Contraction Ratio',
                        ylabel='Incidence Correction i+istall',
                        clabel='asin(o/s)'),
            
            'Fig13':LossInterp("Fig13_Incidence_Correction_i-istall.csv",
                        xlabel='Pitch to backbone length ratio s/b',
                        ylabel='Incidence correction i-istall_sb',
                        clabel='asin(o/s)'),
            
            'Fig13_alpha1':LossInterp("Fig13_Negative_Stalling_Incidence_i-istall.csv",
                        xlabel='Blade Inlet Angle',
                        ylabel='Basic Negative Stalling Incidence (i-istall)',
                        clabel='asin(o/s)'),
            
            'Fig14_i-istall':LossInterp("Fig14_Basic_stalling_incidences_blade_angle_gt_90_i-istall.csv",
                        xlabel='Blade Inlet Angle',
                        ylabel='i-istall',
                        clabel='asin(o/s)'),
            
            'Fig14_i+istall':LossInterp("Fig14_Basic_stalling_incidences_blade_angle_gt_90_i+istall.csv",
                        xlabel='Blade Inlet Angle',
                        ylabel='i+istall',
                        clabel='asin(o/s)'),
            
            'Fig15':LossInterp("Fig15_Minimum_Loss_Incidence-Range_Ratio_Fi.csv",
                        xlabel='Blade Inlet Angle',
                        ylabel='incidence range ratio Fi',
                        clabel='s/b'),
            
            'Fig17':LossInterp("Fig17_Secondary_Loss_Ns_hb.csv",
                        xlabel='Inv Aspect Ratio',
                        ylabel='Secondary Loss Ratio'),
            
            'Fig18':LossInterp("Fig18_base_secondary_factor.csv",
                        xlabel='Vinlet^2/Voutlet^2',
                        ylabel='Basic Secondary Loss Factor',
                        clabel='Fl*s/b'),
            
            }
    default_home = os.path.join(os.path.expanduser("~"), ".cache")
    os.environ['TD3_HOME'] = os.path.join(default_home,'TD3_LossModels')
    os.makedirs(os.environ['TD3_HOME'],exist_ok=True)
    path = pathlib.Path(os.path.join(os.environ['TD3_HOME'],"craigcox"+".pkl"))

    with open(path.absolute(),'wb') as f:
        pickle.dump(data,f)
    
    import shutil
    shutil.copyfile(path.absolute(),os.path.join(os.getcwd(),path.name))