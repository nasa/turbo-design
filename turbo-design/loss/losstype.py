from typing import Any, Dict
from ..lossinterp import LossInterp
import os 
from ..enums import LossType

class LossBaseClass:  
    data: Dict[str,LossInterp]
    _loss_type:LossType
    
    def __init__(self,lossType:LossType):
        # Make the environment directory 
        default_home = os.path.join(os.path.expanduser("~"), ".cache")
        os.environ['TD3_HOME'] = os.path.join(default_home,'TD3_LossModels')
        os.makedirs(os.environ['TD3_HOME'],exist_ok=True)
        
        self._loss_type = lossType

    def __call__(self, row:Any, upstream:Any):
        raise Exception("Method needs to be overridden")

    
    
    @property
    def loss_type(self):
        return self._loss_type