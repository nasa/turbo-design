from dataclasses import dataclass, field
from cantera import Solution

@dataclass
class Coolant:
    fluid:Solution = field(default = Solution('air.yaml')) # cantera solution
    T0:float = field(default=900)                               # Kelvin
    P0:float = field(default=50*101325)                         # Pascal
    massflow_percentage:float = field(default=0.03)     # Fraction of total massflow going through compressor 
    