from abc import ABC, abstractmethod
from typing import Any, Dict, Iterable, Tuple
from ..lossinterp import LossInterp
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
    def __call__(self, row:Any, upstream:Any) -> float:
        """Evaluate the loss for the supplied blade row."""
        raise NotImplementedError

    
    
    @property
    def loss_type(self):
        return self._loss_type


class CompositeLossModel(LossBaseClass):
    """Combines multiple loss models of the same type."""

    def __init__(self, models: Iterable[LossBaseClass]):
        models_tuple: Tuple[LossBaseClass, ...] = tuple(models)
        if not models_tuple:
            raise ValueError("CompositeLossModel requires at least one loss model.")

        loss_type = models_tuple[0].loss_type
        for model in models_tuple[1:]:
            if model.loss_type != loss_type:
                raise ValueError("All loss models must share the same LossType.")

        super().__init__(loss_type)
        self._models: Tuple[LossBaseClass, ...] = models_tuple

    def __call__(self, row: Any, upstream: Any) -> float:
        total_loss = 0.0
        for model in self._models:
            total_loss += float(model(row, upstream))
        return total_loss

    @property
    def models(self) -> Tuple[LossBaseClass, ...]:
        return self._models
