from enum import Enum

class LossType(Enum):
    Pressure = 1
    Enthalpy = 2
    Entropy = 3

class RowType(Enum):
    Stator = 1
    Rotor = 2 
    CounterRotating = 2
    Inlet = 3
    Outlet = 4

class MassflowConstraint(Enum):
    MatchMassFlow:int = 1 # Changes the exit angles to match the massflow
    BalanceMassFlow:int = 2 # Keeps the exit angle but balances the massflow between the stages as best it can. This will affect the static pressure at the stage exit

class PowerType(Enum):
    """The code for BladeRow will assume a PowerType automatically depending on what you specify. If you specify the blade row to have P0_P which is the stator inlet total pressure to rotor exit static pressure then that will be used to calculate all the quantities.
    
    Or if you specify the power for the blade row, the Total Temperature at Rotor Exit will change. 
    
    P0_P: Assumes a Total-Static Pressure Ratio. Using this method will not allow you to match a particular massflow 
    
    T0: Total Temperature at exit is specified based on power required
    
    """
    # 
    P0_P:int = 1 
    T0:int = 2 


class PassageType(Enum):
    Centrifugal:int=0
    Axial:int=1