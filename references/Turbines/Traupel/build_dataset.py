# Read in the excel data
# Create bspline surface with points
import sys, os, pathlib
sys.path.insert(0,"../../../")
from turbodesign import LossInterp
import pickle

if __name__=="__main__": 
    default_home = os.path.join(os.path.expanduser("~"), ".cache")
    os.environ['TD3_HOME'] = os.path.join(default_home,'TD3_LossModels')
    os.makedirs(os.environ['TD3_HOME'],exist_ok=True)
    path = pathlib.Path(os.path.join(os.environ['TD3_HOME'],"traupel"+".pkl"))

    with open(path.absolute(),'wb') as f:
        pickle.dump({
            'Fig01':LossInterp("Fig01_ζ0.csv",
                                        xlabel="alpha1,beta2",
                                        ylabel="Zeta0",
                                        clabel="alpha2,beta3"),
            'Fig02':LossInterp("Fig02_ZetaP0.csv",
                                        xlabel="ZetaP0",
                                        ylabel="alpha1,beta2",
                                        clabel="alpha2,beta3"),
            'Fig03_0':LossInterp("Fig03_0.csv",
                                        xlabel="M",
                                        ylabel="Xm"),
            'Fig03_1':LossInterp("Fig03_1.csv",
                                        xlabel="M",
                                        ylabel="Xm"),
            'Fig04':LossInterp("Fig04.csv",
                                xlabel="delta alpha-tau/g",
                                ylabel="Zeta_Delta",
                                clabel="alpha2,beta3"),
            'Fig05':LossInterp("Fig05.csv",
                               xlabel="delta alpha",
                               ylabel="X_delta",
                               clabel='alpha2,beta3'),
            'Fig06':LossInterp("Fig06.csv",
                        xlabel='c1/c2,w1/w2',
                        ylabel='F',
                        clabel='Delta C, Delta W'),
            'Fig07':LossInterp("Fig07.csv",
                        xlabel='alpha1-beta2',
                        ylabel='H',
                        clabel='alpha2-beta3'),
            'Fig08':LossInterp("Fig08.csv",
                        xlabel='Clearance',
                        ylabel='Xi_Cl'),
            'Fig09':LossInterp("Fig09.csv",
                        xlabel='alpha1-beta2',
                        ylabel='G',
                        clabel='alpha2-beta3'),
            'Fig10':LossInterp("Fig10.csv",
                        xlabel='delta/c',
                        ylabel='K',
                        clabel='0-stator, 1-rotor'),
            'Fig11':LossInterp("Fig11.csv",
                        xlabel='mu',
                        ylabel='Zeta_cl',
                        clabel='0-stator, 1-rotor'),
            'Fig12':LossInterp("Fig12.csv",
                        ylabel='Zeta_z',
                        xlabel='f'),
            'Fig16':LossInterp("Fig16.csv",
                        xlabel='Z',
                        ylabel='(alpha1-alpha1,opt),(beta2-beta2,opt)',
                        clabel='0=a,1=b'),
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
                        xlabel='Clearance',
                        ylabel='Xi_Cl')
    
    import shutil
    shutil.copyfile(path.absolute(),os.path.join(os.getcwd(),path.name))