from dataclasses import dataclass,field

@dataclass
class Stage :
    blade_rows: list = None
    cp: float = None
    streamlines:list = field(default_factory=list)
