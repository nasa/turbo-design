from typing import List
from cantera import composite
from cantera.composite import Solution
from .bladerow import BladeRow,RowType
from .spool import Spool
import json

class CompressorSpool(Spool):
    def __init__(self,rhub:List[float],xhub:List[float],
                 rshroud:List[float],xshroud:List[float],
                 massflow:float,rows=List[BladeRow],
                 num_streamlines:int=3,
                 fluid:Solution=Solution('air.yaml'),rpm:float=-1):
        """Initializes a Compressor Spool

        Args:
            rhub (List[float]): hub radius. Units meters
            xhub (List[float]): x location of each hub radius station. Units meters
            rshroud (List[float]): shroud radius. Units meters
            xshroud (List[float]): x location of each shroud radius. Units meters
            massflow (float): massflow at spool inlet 
            rows (List[BladeRow], optional): List of blade rows. Defaults to List[BladeRow].
            num_streamlines (int, optional): number of streamlines. Defaults to 3.
            gas (ct.Solution, optional): cantera gas solution. Defaults to ct.Solution('air.yaml').
            rpm (float, optional): RPM for the entire spool Optional, you can also set rpm of the blade rows individually. Defaults to -1.
            power_target (float, optional): Sets a target power in kW
        """
        super().__init__(rhub, xhub, rshroud, xshroud, massflow, rows,num_streamlines, fluid, rpm)
        
        pass
    
    def export_properties(self,filename:str="compressor_spool.json"):
        """Export the spool object to json 

        Args:
            filename (str, optional): name of export file. Defaults to "spool.json".
        """
        blade_rows = list()
        total_total_efficiency = list() 

        for indx,row in enumerate(self.blade_rows):
            blade_rows.append(row.to_dict()) # Appending data 
            if row.row_type == RowType.Stator:
                pass
        
        data = {            
            "blade_rows": blade_rows,
            "massflow":self.massflow,
            "rpm":self.rpm,
            "r_streamline":self.r_streamline.tolist(),
            "x_streamline":self.x_streamline.tolist(),
            "rhub":self.rhub.tolist(),
            "rshroud":self.rshroud.tolist(),
            "xhub":self.xhub.tolist(),
            "xshroud":self.xshroud.tolist(),
            "num_streamlines":self.num_streamlines,
        }
        # Dump all the Python objects into a single JSON file.
        with open(filename, "w") as f:
            json.dump(data, f, indent=4)