# Read in the excel data
# Create bspline surface with points
import sys, os, pathlib
sys.path.insert(0,"../../../")
from td3 import LossInterp
import pickle

if __name__=="__main__": 
    default_home = os.path.join(os.path.expanduser("~"), ".cache")
    os.environ['TD3_HOME'] = os.path.join(default_home,'TD3_LossModels')
    os.makedirs(os.environ['TD3_HOME'],exist_ok=True)
    path = pathlib.Path(os.path.join(os.environ['TD3_HOME'],"ainleymathieson"+".pkl"))

    with open(path.absolute(),'wb') as f:
        pickle.dump({
            'Fig04a':LossInterp("Fig04a.csv",
                                        xlabel="Pitch/Chord",
                                        ylabel="Profile Loss Coefficient Yp, Nozzle/Stator Beta1=0",
                                        clabel="Outlet Gas Angle"),
            
            'Fig04b':LossInterp("Fig04b.csv",
                                        xlabel="Pitch/Chord",
                                        ylabel="Profile Loss Coefficient Yp, Impulse",
                                        clabel="Outlet Values of Exit Gas Angle, alpha2"),
            
            'Fig05':LossInterp("Fig05.csv",
                                        xlabel="arccos(o/s)",
                                        ylabel="Exlt Flow Angle Alpha2*"),
            
            'Fig06':LossInterp("Fig06.csv",
                                        xlabel="i/istall",
                                        ylabel="Yp/Yp(i=0)"),
            
            'Fig07_top':LossInterp("Fig07a_top.csv",
                               xlabel="pitch/chord",
                               ylabel="is-is(s/c=0.75)"),
            
            'Fig07_mid':LossInterp("Fig07a_mid.csv",
                        xlabel='pitch/chord',
                        ylabel='alpha2/alpha2(s_c=0.75)'),
            
            'Fig07_bot':LossInterp("Fig07b_bottom.csv",
                        xlabel='beta1/alpha2',
                        ylabel='stalling incidence s_c=0.75',
                        clabel='alpha2'),
            
            'Fig08':LossInterp("Fig08.csv",
                        xlabel='(A2/A1)**2 / (1+ID/OD)',
                        ylabel='Lambda'),
            
            'Fig09':LossInterp("Fig09.csv",
                        xlabel='te/s',
                        ylabel='Yt/Yt(te/s=0.02)'),
            
            'Fig10':LossInterp("Fig10.csv",
                        xlabel='Y',
                        ylabel='Non-dimensional massflow',
                        clabel='gamma'),
            
            'Fig11':LossInterp("Fig11.csv",
                        xlabel='Y',
                        ylabel='(W/P)_Q=Qmax',
                        clabel='gamma'),
            
            'Fig12':LossInterp("Fig12.csv",
                        xlabel='Q/Q_max',
                        ylabel='W/W_Q=Qmax'),
            
            'Fig14':LossInterp("Fig14.csv",
                        xlabel='Incidence i',
                        ylabel='Yt'),
            
            'Fig14_Nozzle_Inlet_Angle':LossInterp("Fig14_Top_Nozzle.csv",
                        xlabel='Blade Outlet Mach',
                        ylabel='Nozzle Inlet Angle'),
            
            'Fig14_Nozzle_Outlet_Angle':LossInterp("Fig14_Nozzle_Outlet.csv",
                        xlabel='Blade Outlet Mach',
                        ylabel='Nozzle Outlet Angle'),
            
            'Fig14_Rotor_Outlet_Angle':LossInterp("Fig14_Top_Rotor.csv",
                        xlabel='Blade Outlet Mach',
                        ylabel='Rotor Outlet Angle'),
            
            'Fig14_Rotor_Inlet_Angle':LossInterp("Fig14_Top.csv",
                        xlabel='Blade Outlet Mach',
                        ylabel='Rotor Inlet Angle'),
            
            'Fig15':LossInterp("Fig15.csv",
                        xlabel='Overall Pressure Ratio P3/Pi',
                        ylabel='Non-dimensional Massflow'),
            
            },f)
    # Fig04ab 
    # Fig04a = LossInterp("Fig04a.csv",
    #         xlabel="Pitch/Chord",
    #         ylabel="Profile Loss Coefficient Yp, Nozzle/Stator Beta1=0",
    #         clabel="Outlet Gas Angle"),
    # Fig04a.plot()
    
    # Fig04b = LossInterp("Fig04b.csv",
    #                     xlabel="Pitch/Chord",
    #                     ylabel="Profile Loss Coefficient Yp, Impulse",
    #                     clabel="Outlet Values of Exit Gas Angle, alpha2"),
    # Fig04b.plot()
    # Fig08
    Fig08 = LossInterp("Fig08.csv",
                        xlabel='(A2/A1)**2 / (1+ID/OD)',
                        ylabel='Lambda')
    Fig08.plot()
    print('check')