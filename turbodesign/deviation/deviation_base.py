from abc import ABC, abstractmethod
from typing import Any, Dict, Iterable, Tuple
import numpy.typing as npt
import os 

class DeviationBaseClass(ABC):  
    """All deviation models should inherit from this abstract base class 
    """
    def __init__(self):
        # Make the environment directory 
        default_home = os.path.join(os.path.expanduser("~"), ".cache")
        os.environ['TD3_HOME'] = os.path.join(default_home,'TD3_LossModels')
        os.makedirs(os.environ['TD3_HOME'],exist_ok=True)
        
    @abstractmethod
    def __call__(self, row:Any, upstream:Any) -> npt.NDArray:
        """Evaluate the loss for the supplied blade row."""
        raise NotImplementedError
    

