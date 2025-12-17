from typing import List
import numpy as np 


def convert_to_list(x:float) -> List[float]:
    return x if isinstance(x, list) else [x]

def convert_to_ndarray(t) -> np.ndarray:
    """Converts a scalar or list to a numpy array
    Args:
        t (float,list): [description]
    Returns:
        np.ndarray: variable as an array
    """
    if hasattr(t, "default_factory"):
        try:
            t = t.default_factory()
        except Exception:
            t = np.array([0.0])
    if type(t) is not np.ndarray and type(t) is not list: # Scalar
        t = np.array([t],dtype=float)
    elif (type(t) is list):
        t = np.array(t,dtype=float)
    return t


def safe_interpolate(values, src_r, dst_r, default: float = 0.0, radians: bool = False, interp_func=None):
    """Safely convert, default, and interpolate quantities onto target radii."""
    from .bladerow import interpolate_quantities  # local import to avoid cycles

    arr = convert_to_ndarray(values)
    if hasattr(arr, "default_factory"):
        arr = arr.default_factory()
    if arr.size == 0:
        arr = np.array([default], dtype=float)
    src = convert_to_ndarray(src_r)
    dst = convert_to_ndarray(dst_r)
    if src.size == 0:
        src = np.linspace(0, 1, len(arr))
    if arr.size > 1 and src.size != arr.size:
        src = np.linspace(0, 1, len(arr))
    if arr.size == 1:
        arr = arr[0] * np.ones_like(dst)
    else:
        f = interp_func if interp_func is not None else interpolate_quantities
        arr = f(arr, src, dst)
    if radians:
        arr = np.radians(arr)
    return arr
