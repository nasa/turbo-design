# type: ignore[arg-type, reportUnknownArgumentType]
from dataclasses import field
import json
from typing import Dict, List, Union
import matplotlib.pyplot as plt
from .bladerow import BladeRow
import numpy as np
import numpy.typing as npt
from .enums import RowType, MassflowConstraint
from cantera.composite import Solution
from .loss.turbine import TD2
from pyturbo.helper import line2D
from scipy.interpolate import interp1d
from .passage import Passage
from .inlet import Inlet
from .outlet import Outlet 

class Spool:
    blade_rows:List[BladeRow]=[]
    massflow:float=0
    rpm:float = 0
    
    passage:Passage
    t_streamline: npt.NDArray = field(default_factory=lambda: np.zeros((10,)))
    num_streamlines:int=0
    
    # Fluid entering the spool
    _fluid: Solution = Solution("air.yaml")

    massflow_constraint: MassflowConstraint
    
    # Inlet Conditions 
    def __init__(self,passage:Passage,
                 massflow:float,rows=List[BladeRow],
                 num_streamlines:int=3,
                 fluid:Solution=Solution('air.yaml'),
                 rpm:float=-1,
                 massflow_constraint:MassflowConstraint=MassflowConstraint.MatchMassFlow):
        """Initializes a Spool


        Note:
            When it comes to massflow, the exit angle of the blade rows along with the upstream total pressure is used to set the massflow through the spool. If turning is too high then massflow is limited and cannot. 
            This code gives you the option of setting a massflow and varying the exit angles to match that particular massflow. (MatchMassflow)
            Or keeping the exit angle and modifying the inlet condition to set a specific massflow through the stage. (BalanceMassflow)
            
        Args:
            passage (Passage): Passage defining hub and shroud
            massflow (float): massflow at spool inlet 
            rows (List[BladeRow], optional): List of blade rows. Defaults to List[BladeRow].
            num_streamlines (int, optional): number of streamlines. Defaults to 3.
            gas (ct.Solution, optional): cantera gas solution. Defaults to ct.Solution('air.yaml').
            rpm (float, optional): RPM for the entire spool Optional, you can also set rpm of the blade rows individually. Defaults to -1.
            massflow_constraint (MassflowConstraint, optional): MatchMassflow - Matches the massflow defined in the spool. BalanceMassflow - Balances the massflow between BladeRows, matches the lowest massflow.
        """
        self.passage = passage
        self.fluid = fluid
        self.massflow = massflow
        self.num_streamlines = num_streamlines
        self.blade_rows = rows
        self.massflow_constraint = massflow_constraint
        
        self.rpm = rpm
        
        for i in range(len(self.blade_rows)):
            '''Defining RPM in the bladerows
                There's 2 ways of setting the RPM. 
                1. Conventional: Stator-Rotor-Stator-Rotor-etc. Set the RPM equally across all
                2. Counter Rotation: We either use the RPM already persecribed for each blade row.
            '''    
            if (type(self.blade_rows[i]) != Inlet) and (type(self.blade_rows[i]) != Outlet):
                self.blade_rows[i].fluid = self.fluid
                self.blade_rows[i].rpm = rpm
                self.blade_rows[i].axial_chord = self.blade_rows[i].axial_location * self.passage.hub_length
            

    @property
    def fluid(self):
        return self._fluid

    @fluid.setter
    def fluid(self,newFluid:Solution):
        """Change the type of gas used in the spool

        Args:
            new_gas (ct.Solution.Solution): New gas mixture
        """
        self._fluid = newFluid


    def set_blade_row_rpm(self,index:int,rpm:float):
        """sets the rpm of a particular blade row"""    
        self.blade_rows[index].rpm = rpm
        
    def set_blade_row_type(self,blade_row_index:int,rowType:RowType):
        """Sets a blade row to be a different type.

        Args:
            blade_row_index (int): Index of blade row. Keep in mind Blade row 0 is inlet, -1 = outlet
        """
        self.blade_rows[blade_row_index].row_type = rowType

    def set_blade_row_exit_angles(self,radius:Dict[int,List[float]],beta:Dict[int,List[float]],IsSupersonic:bool=False):
        """Set the intended exit flow angles for the spool. 

            Useful if you already have a geometry that you want to simulate
            
            You can set the values for each blade row or simply the inlet and exit
 
        Args:
            radius (np.ndarray): Matrix containing [[r1,r2,r3],[r1,r2,r3]] blade_rows x radial locations
            beta (Dict[int,List[float]]): Metal Angle of the geometry, it is assumed that the fluid will leave at this angle.
            IsSupersonic(bool): if solution is supersonic at spool exit

        Example: 
            # Station 0 has angles of [73.1,74.2,75.4] 
            # Station 5 has angles of [69,69,69] 
            beta =  {0: [73.1,74.2,75.4], 5: [69,69,69]}

        """
        for k,v in radius.items():                                   
            self.blade_rows[k].radii_geom = v 
        for k,v in beta.items():                                   
            self.blade_rows[k].beta_geom = v 
            self.blade_rows[k].beta_fixed = True
        for br in self.blade_rows:
            if not IsSupersonic:
                br.solution_type = SolutionType.subsonic
            else:
                br.solution_type = SolutionType.supersonic
        
    
    def initialize_streamlines(self):
        """Initialize streamline and compute the curvature. 
            This function is called once
        """
             
        for row in self.blade_rows: 
            row.phi = np.zeros((self.num_streamlines,))
            row.rm = np.zeros((self.num_streamlines,))
            row.r = np.zeros((self.num_streamlines,))
            
            t_radial = np.linspace(0,1,self.num_streamlines)
            self.calculate_streamline_curvature(row,t_radial)
                
            # Set the loss function if it's not set
            if (row.loss_function == None):  
                row.loss_function = TD2()
            
    def calculate_streamline_curvature(self,row:BladeRow,t_radial:Union[List[float],npt.NDArray]):
        """Called to calculate new streamline curvature

        Args:
            row (BladeRow): blade row
            t_radial (Union[List[float],npt.NDArray]): normalized streamline radial locations
        """
        for i,tr in enumerate(t_radial):
            t_streamline, x_streamline, r_streamline = self.passage.get_streamline(tr)                
            phi, rm, r = self.passage.streamline_curvature(x_streamline,r_streamline)
            row.phi[i] = float(interp1d(t_streamline,phi)(row.axial_location))
            row.rm[i] = float(interp1d(t_streamline,rm)(row.axial_location))
            row.r[i] = float(interp1d(t_streamline,r)(row.axial_location))
          
    def solve(self):
        raise NotImplementedError('Solve is not implemented')

    
    def plot(self):
        """Plots the hub, shroud, and all the streamlines
        """
        plt.figure(num=1,clear=True,dpi=150,figsize=(15,10))
        plt.plot(self.passage.xhub_pts,self.passage.rhub_pts,label='hub',linestyle='solid',linewidth=2,color='black')
        plt.plot(self.passage.xshroud_pts,self.passage.rshroud_pts,label='shroud',linestyle='solid',linewidth=2,color='black')
        
        hub_length = np.sum(np.sqrt(np.diff(self.passage.xhub_pts)**2 + np.diff(self.passage.rhub_pts)**2))
        # Populate the streamlines 
        x_streamline = np.zeros((self.num_streamlines,len(self.blade_rows)))
        r_streamline = np.zeros((self.num_streamlines,len(self.blade_rows)))
        for i in range(len(self.blade_rows)):
            x_streamline[:,i] = self.blade_rows[i].x
            r_streamline[:,i] = self.blade_rows[i].r
        
        for i in range(1,len(self.blade_rows)-1): # plot dashed lines for each steamline
            plt.plot(x_streamline[:,i],r_streamline[:,i],'--b',linewidth=1.5)
        
        for i in range(len(self.blade_rows)):
            
            row = self.blade_rows[i]
            plt.plot(row.x,row.r,linestyle='dashed',linewidth=1.5,color='blue',alpha=0.4)
            plt.plot(x_streamline[:,i],r_streamline[:,i],'or')
            
            if i == 0: 
                pass # Inlet 
            else:  # i>0
                upstream = self.blade_rows[i-1]
                if upstream.row_type== RowType.Inlet:
                    cut_line1,_,_ = self.passage.get_cutting_line((row.axial_location*hub_length +(0.5*row.blade_to_blade_gap*row.axial_chord) - row.axial_chord)/hub_length)
                else:
                    cut_line1,_,_ = self.passage.get_cutting_line((upstream.axial_location*hub_length)/hub_length)
                cut_line2,_,_ = self.passage.get_cutting_line((row.axial_location*hub_length-(0.5*row.blade_to_blade_gap*row.axial_chord))/hub_length)
                
            if self.blade_rows[i].row_type == RowType.Stator:
                x1,r1 = cut_line1.get_point(np.linspace(0,1,10))
                plt.plot(x1,r1,'m')
                x2,r2 = cut_line2.get_point(np.linspace(0,1,10))
                plt.plot(x2,r2,'m')
                x_text = (x1+x2)/2; r_text = (r1+r2)/2
                plt.text(x_text.mean(),r_text.mean(),"Stator",fontdict={"fontsize":"xx-large"})
            elif self.blade_rows[i].row_type == RowType.Rotor:
                x1,r1 = cut_line1.get_point(np.linspace(0,1,10))
                plt.plot(x1,r1,color='brown')
                x2,r2 = cut_line2.get_point(np.linspace(0,1,10))
                plt.plot(x2,r2,color='brown')
                x_text = (x1+x2)/2; r_text = (r1+r2)/2
                plt.text(x_text.mean(),r_text.mean(),"Rotor",fontdict={"fontsize":"xx-large"})
        plt.axis('scaled')
        plt.savefig(f"Meridional.png",transparent=False,dpi=150)
        plt.show()
        
    def plot_velocity_triangles(self):
        """Plots the velocity triangles for each blade row
            Made for turbines 
        """
        prop = dict(arrowstyle="-|>,head_width=0.4,head_length=0.8",
            shrinkA=0,shrinkB=0)
        
        
        for j in range(self.num_streamlines):
            x_start = 0
            y_max = 0; y_min = 0
            plt.figure(num=1,clear=True)
            for i in range(1,len(self.blade_rows)):
                row = self.blade_rows[i]
                x_end = x_start + row.Vx.mean()
                dx = x_end - x_start
                
                Vt = row.Vt[j]
                Wt = row.Wt[j]
                U = row.U[j]

                y_max = (Vt if Vt>y_max else y_max)
                y_max = (Wt if Wt>y_max else y_max)
                
                y_min = (Vt if Vt<y_min else y_min)
                y_min = (Wt if Wt<y_min else y_min)
                
                # V 
                plt.annotate("", xy=(x_end,Vt), xytext=(x_start,0), arrowprops=prop)
                plt.text((x_start+x_end)/2,Vt/2*1.1,"V",fontdict={"fontsize":"xx-large"})
                
                # W
                plt.annotate("", xy=(x_end,Wt), xytext=(x_start,0), arrowprops=prop)
                plt.text((x_start+x_end)/2,Wt/2*1.1,"W",fontdict={"fontsize":"xx-large"})
                
                if (abs(row.Vt[j]) > abs(row.Wt[j])):
                    # Shift Vt to right just a bit so we can see it
                    plt.annotate("", xy=(x_end,Wt), xytext=(x_end,0), arrowprops=prop) # Wt
                    plt.text(x_end+dx*0.1,Wt/2,"Wt",fontdict={"fontsize":"xx-large"})

                    plt.annotate("", xy=(x_end,U+Wt), xytext=(x_end,Wt), arrowprops=prop) # U
                    plt.text(x_end+dx*0.1,(Wt+U)/2,"U",fontdict={"fontsize":"xx-large"})
                else:
                    # Shift Wt to right just a bit so we can see it
                    plt.annotate("", xy=(x_end,Vt), xytext=(x_end,0), arrowprops=prop) # Vt
                    plt.text(x_end+dx*0.1,Vt/2,"Vt",fontdict={"fontsize":"xx-large"})

                    plt.annotate("", xy=(x_end,Wt+U), xytext=(x_end,Wt), arrowprops=prop) # U
                    plt.text(x_end+dx*0.1,Wt+U/2,"U",fontdict={"fontsize":"xx-large"})

                plt.text((x_start+x_end)/2,-y_max*0.95,row.row_type.name,fontdict={"fontsize":"xx-large"})
                x_start += row.Vx[j]
                plt.axis([0,x_end+dx, y_min, y_max])
            plt.title(f"Velocity Triangles for Streamline {j}")
            plt.savefig(f"streamline_{j:04d}.png",transparent=False,dpi=150)
    


   
    
def JSON_to_Spool(filename:str="spool.json"):
    """Reads a JSON file with the properties 
    #! Need to instantiate and return spool. Fluid isn't loaded fyi 

    Args:
        filename (str, optional): json filename with properties . Defaults to "spool.json".
    """
    with open(filename, "r") as f:
        data = json.load(f)
