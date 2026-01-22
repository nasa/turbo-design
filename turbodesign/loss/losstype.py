from abc import ABC, abstractmethod
from typing import Any, Dict, Iterable, Tuple
from ..lossinterp import LossInterp
import numpy.typing as npt
import os 
from ..enums import LossType

class LossBaseClass(ABC):  
    data: Dict[str,LossInterp]
    _loss_type:LossType
    
    def __init__(self,lossType:LossType):
        # Make the environment directory 
        default_home = os.path.join(os.path.expanduser("~"), ".cache")
        os.environ['TD3_HOME'] = os.path.join(default_home,'TD3_LossModels')
        os.makedirs(os.environ['TD3_HOME'],exist_ok=True)
        
        self._loss_type = lossType

    @abstractmethod
    def __call__(self, row:Any, upstream:Any) -> npt.NDArray:
        """Evaluate the loss for the supplied blade row."""
        raise NotImplementedError
    
    @property
    def loss_type(self):
        return self._loss_type
