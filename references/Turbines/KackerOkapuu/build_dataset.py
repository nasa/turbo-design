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
    path = pathlib.Path(os.path.join(os.environ['TD3_HOME'],"kackerokapuu"+".pkl"))

    with open(path.absolute(),'wb') as f:
        pickle.dump({
            'Fig01_beta0':LossInterp("Fig01_beta0.csv",
                                        xlabel="Pitch/Chord",
                                        ylabel="Yp",
                                        clabel="alpha2"),
            
            'Fig02':LossInterp("Fig02_beta=alpha2.csv",
                                        xlabel="Pitch/Chord",
                                        ylabel="Yp",
                                        clabel="alpha2"),
            
            'Fig04':LossInterp("Fig04_tmax_c.csv",
                                        xlabel="beta1+beta2",
                                        ylabel="tmax/c"),
            
            'Fig05':LossInterp("Fig05_stagger.csv",
                                        xlabel="beta1",
                                        ylabel="stagger",
                                        clabel="beta2"),
            
            'Fig08_K1':LossInterp("Fig08_K1.csv",
                                        xlabel="M2",
                                        ylabel="K1"),
            
            'Fig09_K2':LossInterp("Fig09_K2.csv",
                                        xlabel="M1/M2",
                                        ylabel="K2"),
            
            'Fig14_Axial_Entry':LossInterp("Fig14_axial_entry_nozzle.csv",
                                        xlabel="t/o",
                                        ylabel="delta phi tet 2"),
            
            'Fig14_Impulse':LossInterp("Fig14_impulse_blading.csv",
                                        xlabel="t/o",
                                        ylabel="delta phi tet 2")
            
            },f)
        
        Fig01 = LossInterp("Fig01_beta0.csv",
                                        xlabel="Pitch/Chord",
                                        ylabel="Yp",
                                        clabel="alpha2")
        Fig01.plot()
    
    import shutil
    shutil.copyfile(path.absolute(),os.path.join(os.getcwd(),path.name))