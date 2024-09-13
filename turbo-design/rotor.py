from dataclasses import dataclass, field
from .enums import RowType
from .bladerow import BladeRow
from .arrayfuncs import convert_to_ndarray
import numpy as np 

class Rotor(BladeRow):
    """Station defined at rotor exit

    Inherits:
        (BladeRow): Defines the properties of the blade row
    """

    beta2:np.ndarray = field(default_factory=lambda: np.array([72]))    # Degrees   
    M_rel:np.ndarray = field(default_factory=lambda: np.array([1.2]))   # Relative mach number
    power:float = 0         # Watts

    def __init__(self,power:float=50000):
        """Initialize the Rotor

        Args:
            coolant_P0 (float, optional): Coolant massflow used to cool the rotor internal cooling, tip leakage, hub disk leakage. Defaults to 0.
            coolant_T0 (float, optional): _description_. Defaults to 500.
            power (float, optional): _description_. Defaults to 50000.
        """
        self.row_type = RowType.RotorExit
        self.power = power
        super().__init__(self)

    def initialize_rotor_exit(self,beta:float,M_exit:float):
        """Uses the beta and exit mach number to predict a value for Vm

        Args:
            beta (float): exit relative flow angle
            M_exit (float): exit relative mach number
        """
        self.beta2 = convert_to_ndarray(beta)
        self.M_rel = convert_to_ndarray(M_exit)
