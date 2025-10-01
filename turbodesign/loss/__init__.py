from .losstype import LossBaseClass, LossType

__all__ = ["LossBaseClass", "FixedPressureLoss"]

def __getattr__(name: str):
    if name == "FixedPressureLoss":
        from .fixedpressureloss import FixedPressureLoss
        return FixedPressureLoss
    raise AttributeError(name)
