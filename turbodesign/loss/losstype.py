from abc import ABC, abstractmethod
from typing import Any, Dict, Iterable, Tuple
from ..lossinterp import LossInterp
import numpy.typing as npt
import os 
from ..enums import LossType

class LossBaseClass(ABC):  
    data: Dict[str,LossInterp]
    _loss_type:LossType
    
    def __init__(self,lossType:LossType, is_parasitic:bool=False):
        # Make the environment directory
        default_home = os.path.join(os.path.expanduser("~"), ".cache")
        os.environ['TD3_HOME'] = os.path.join(default_home,'TD3_LossModels')
        os.makedirs(os.environ['TD3_HOME'],exist_ok=True)

        self._loss_type = lossType
        self._is_parasitic = is_parasitic

    @abstractmethod
    def __call__(self, row:Any, upstream:Any) -> npt.NDArray:
        """Evaluate the loss for the supplied blade row."""
        raise NotImplementedError

    @property
    def loss_type(self):
        return self._loss_type

    @property
    def is_parasitic(self) -> bool:
        """Whether this model carries parasitic work, and therefore cannot be
        expressed as a pressure-loss coefficient at all.

        Internal loss (incidence, friction, clearance, mixing, shock, blade
        loading, diffusion): destroys total pressure, does not change the
        work. It can be expressed as a pressure-loss coefficient (Yp).

        Parasitic loss (disc friction, recirculation, leakage, or an
        aggregate that includes any of these): adds work (raises T0) and
        destroys no total pressure. It belongs in the denominator of
        efficiency, eta = (dh_Euler - dh_internal) / (dh_Euler +
        dh_parasitic). It CANNOT be expressed as a Yp -- there is no pressure
        term to attribute it to.

        Defaults to False (internal); subclasses that carry parasitic work
        must pass ``is_parasitic=True`` to this constructor.
        """
        return self._is_parasitic
