from typing import List
from cantera.composite import Solution
from .bladerow import BladeRow, interpolate_streamline_radii
from .enums import RowType, MassflowConstraint, LossType, PassageType
from .spool import Spool
import json
from .passage import Passage
from scipy.interpolate import interp1d
import numpy as np
import numpy.typing as npt
from .td_math import inlet_calc,rotor_calc, stator_calc, compute_massflow, compute_power, compute_gas_constants
from .solve_radeq import adjust_streamlines, radeq
from scipy.optimize import minimize_scalar, minimize, fmin_slsqp


class TurbineSpool(Spool):
    def __init__(self,passage:Passage,
                 massflow:float,rows:List[BladeRow],
                 num_streamlines:int=3,
                 fluid:Solution=Solution('air.yaml'),
                 rpm:float=-1,
                massflow_constraint:MassflowConstraint=MassflowConstraint.MatchMassFlow):
        """Initializes a Turbine Spool

        Args:
            passage (Passage): Passage defining hub and shroud
            massflow (float): massflow at spool inlet 
            rows (List[BladeRow], optional): List of blade rows. Defaults to List[BladeRow].
            num_streamlines (int, optional): number of streamlines. Defaults to 3.
            fluid (ct.Solution, optional): cantera gas solution. Defaults to ct.Solution('air.yaml').
            rpm (float, optional): RPM for the entire spool Optional, you can also set rpm of the blade rows individually. Defaults to -1.
            massflow_constraint (MassflowConstraint, optional): MatchMassflow - Matches the massflow defined in the spool. BalanceMassflow - Balances the massflow between BladeRows, matches the lowest massflow.

        """
        super().__init__(passage, massflow, rows,num_streamlines, fluid, rpm)
        self.massflow_constraint = massflow_constraint
        pass

    def initialize_quantities(self):
        """Initializes the massflow throughout the rows 
        """
        # Massflow from inlet already defined
        
        # Inlet
        W0 = self.massflow
        inlet = self.blade_rows[0]
        if self.blade_rows[0].row_type == RowType.Inlet:
            self.blade_rows[0].massflow = np.linspace(0,1,self.num_streamlines)*W0
            self.blade_rows[0].total_massflow = W0 
            self.blade_rows[0].total_massflow_no_coolant = W0 
            
            interpolate_streamline_radii(self.blade_rows[0],self.passage)

            # Set the gas to total values for now             
            self.blade_rows[0].fluid.TP = self.blade_rows[0].T0.mean(), self.blade_rows[0].P0.mean()
            self.blade_rows[0].Cp = self.blade_rows[0].fluid.cp
            self.blade_rows[0].Cv = self.blade_rows[0].fluid.cv
            self.blade_rows[0].R = self.blade_rows[0].Cp-self.blade_rows[0].Cv
            self.blade_rows[0].gamma = self.blade_rows[0].Cp/self.blade_rows[0].Cv
            self.blade_rows[0].rho[:] = self.blade_rows[0].fluid.density
            inlet_calc(self.blade_rows[0])
        
        for row in self.blade_rows:
            interpolate_streamline_radii(row,self.passage)

        outlet = self.blade_rows[-1]
        for j in range(self.num_streamlines):
            P0 = inlet.get_total_pressure(inlet.percent_hub_shroud[j])
            percents = np.zeros(shape=(len(self.blade_rows)-2)) + 0.3
            percents[-1] = 1
            Ps_range = outlet_pressure(percents=percents,inletP0=inlet.P0[j],outletP=outlet.P[j])
            for i in range(1,len(self.blade_rows)-1):
                self.blade_rows[i].P[j] = Ps_range[i-1]
            
        # Pass T0 and P0 to all the other blade_rows
        for i in range(1,len(self.blade_rows)-1):
            upstream = self.blade_rows[i-1] # Inlet conditions solved before this step
            if i+1<len(self.blade_rows):
                downstream = self.blade_rows[i+1]
            else:
                downstream = None 
            
            row = self.blade_rows[i]
            if (row.coolant is not None):
                T0c = self.blade_rows[i].coolant.fluid.T
                P0c = self.blade_rows[i].coolant.fluid.P
                W0c = self.blade_rows[i].coolant.massflow_percentage * self.massflow
                Cpc = self.blade_rows[i].coolant.fluid.cp
            else:
                T0c = 100
                P0c = 0
                W0c = 0
                Cpc = 0
            
            T0 = upstream.T0
            P0 = upstream.P0
            Cp = upstream.Cp
            
            T0 = (W0*Cp*T0 + W0c*Cpc*T0c)/(Cpc * W0c + Cp*W0)
            P0 = (W0*Cp*P0 + W0c*Cpc*P0c)/(Cpc * W0c + Cp*W0)   
            Cp = (W0*Cp + W0c*Cpc)/(W0c + W0)                   # Weighted 
            
            if row.row_type == RowType.Stator:
                T0 = upstream.T0
            else:
                T0 = upstream.T0 - row.power / (Cp*(W0 + W0c))
            
            W0 += W0c
            row.T0 = T0
            row.P0 = P0
            row.Cp = Cp
            row.total_massflow = W0
            row.massflow = np.linspace(0,1,self.num_streamlines)*row.total_massflow
            # Pass Quantities: rho, P0, T0
            
            row.rho = upstream.rho
            row.gamma = upstream.gamma
            row.R = upstream.R
            
            if row.row_type == RowType.Stator:
                stator_calc(row,upstream,downstream)
            elif row.row_type == RowType.Rotor:
                rotor_calc(row,upstream)
                compute_massflow(row)
                compute_power(row,upstream)        
    
    def solve(self):
        """
            Solve for the exit flow angles to match the massflow distribution at the stage exit
        """
        self.initialize_streamlines()
        self.initialize_quantities()
        
        if self.massflow_constraint ==MassflowConstraint.MatchMassFlow:
            self.__match_massflow()
        elif self.massflow_constraint == MassflowConstraint.BalanceMassFlow:
            self.__balance_massflow()
    
    
    def __match_massflow(self):
        """ Matches the massflow between streamtubes by changing exit angles. Doesn't use radial equilibrium.
        """
        for _ in range(3):
            # Step 1: Solve a blade row for exit angle to maintain massflow
            for i in range(len(self.blade_rows)):
                row = self.blade_rows[i]
                # Upstream Row
                if i == 0:
                    upstream = self.blade_rows[i]
                else:
                    upstream = self.blade_rows[i-1]
                if i<len(self.blade_rows)-1:
                    downstream = self.blade_rows[i+1]
                else:
                    downstream = None
                
                if row.row_type == RowType.Stator:
                    bounds = [0,80]
                elif row.row_type == RowType.Rotor:
                    bounds = [-80,0]
                if row.row_type != RowType.Inlet:
                    for j in range(1,self.num_streamlines):
                        res = minimize_scalar(massflow_loss_function, bounds=bounds,args=(j,row,upstream,downstream),tol=1E-2)
                        if row.row_type == RowType.Rotor:
                            row.beta2[j] = np.radians(res.x)
                             # Initialize the value at the hub to not upset the mean
                            row.beta2[0] = 1/(len(row.beta2)-1)*row.beta2[1:].sum()
                        elif row.row_type == RowType.Stator:
                            row.alpha2[j] = np.radians(res.x)
                            row.alpha2[0] = 1/(len(row.alpha2)-1)*row.alpha2[1:].sum()
                    upstream = compute_gas_constants(upstream)
                    row = compute_gas_constants(row)
                    
                    
            # Step 3: Adjust streamlines to evenly divide massflow
            adjust_streamlines(self.blade_rows,self.passage)
            
    def __balance_massflow(self):
        """ Balances the massflow between rows. Use radial equilibrium.

            Types of stages:
            1. Stator - Rotor         | Stator - Rotor
            2. Rotor                  | Stator - Rotor | Stator - Rotor
            3. Stator - Rotor         | CounterRotating | Stator - Rotor
            4. Rotor-Counter Rotating | Stator - Rotor
            5. Counter Rotating - Rotor | Stator - Rotor 
            
            Steps:
                1. Split the blade rows into stages stator-rotor pairs or rotor rotor or rotor 
                2. Change degree of reaction to match the total massflow
                3. Adjust the streamlines for each blade row to balance the massflow
        """
        
        # Balance the massflow between Stages
        def balance_massflows(x0:List[float],blade_rows:List[List[BladeRow]],P0:npt.NDArray,P:npt.NDArray):
            total_massflow = list(); massflow_stage = list()
            stage_ids = list(set([row.stage_id for row in self.blade_rows if row.stage_id>=0])); s = 0 
            sign = 1 
            for j in range(self.num_streamlines):                
                Ps_range = outlet_pressure(x0,P0[j],P[j])
                for i in range(1,len(blade_rows)-1):
                    blade_rows[i].P[j] = Ps_range[i-1]
            blade_rows[-1].P = P
            calculate_massflows(blade_rows,True)
            for row in blade_rows[1:]:
                total_massflow.append(row.total_massflow_no_coolant)
            for s in stage_ids:
                for row in blade_rows:
                    if row.stage_id == s and row.row_type == RowType.Rotor:
                        massflow_stage.append(sign*row.total_massflow_no_coolant)
                        sign*=-1
            if len(stage_ids) % 2 == 1:
                massflow_stage.append(massflow_stage[-1]*sign)
            print(x0)
            return np.std(total_massflow)*2 # + abs(sum(massflow_stage))  # Equation 28
        
        # Break apart the rows to stages
        outlet_P=list(); outlet_P_guess = list()
        
        for i in range(1,len(self.blade_rows)-2):
            outlet_P.append(self.blade_rows[i].inlet_to_outlet_pratio)
            outlet_P_guess.append(np.mean(self.blade_rows[i].inlet_to_outlet_pratio))
        
        if len(outlet_P) == 1:
            res1 = minimize_scalar(fun=balance_massflows,args=(self.blade_rows[:-1],self.blade_rows[0].P0,self.blade_rows[-1].P), 
                        bounds=outlet_P[0],tol=0.001,options={'disp': True})
            x = res1.x
        else:
            x = fmin_slsqp(func=balance_massflows,args=(self.blade_rows[:-1],self.blade_rows[0].P0,self.blade_rows[-1].P), 
                        bounds=outlet_P, x0=outlet_P_guess,epsilon=0.001,iter=100) # ,tol=0.001,options={'disp': True})
        
        # Adjust the inlet: Set the massflow
        self.blade_rows[0].massflow = np.linspace(0,1,self.num_streamlines)*self.blade_rows[1].total_massflow_no_coolant
        inlet_calc(self.blade_rows[0]) # adjust the inlet to match massflow 
        for _ in range(3):
            adjust_streamlines(self.blade_rows[:-1],self.passage)
            self.blade_rows[-1].transfer_quantities(self.blade_rows[-2])
            self.blade_rows[-1].P = self.blade_rows[-1].get_static_pressure(self.blade_rows[-1].percent_hub_shroud)
            err = balance_massflows(x,self.blade_rows[:-1],self.blade_rows[0].P0,self.blade_rows[-1].P)
        if err>5E-2:
            print(f"Massflow is not convergenced error:{err}")
        else: 
            print(f"Massflow converged to less than 0.05kg/s error:{err}")
    
    def export_properties(self,filename:str="turbine_spool.json"):
        """Export the spool object to json 

        Args:
            filename (str, optional): name of export file. Defaults to "spool.json".
        """
        blade_rows = list()
        degree_of_reaction = list() 
        total_total_efficiency = list() 
        total_static_efficiency = list()
        stage_loading = list() 
        euler_power = list() 
        
        x_streamline = np.zeros((self.num_streamlines,len(self.blade_rows)))
        r_streamline = np.zeros((self.num_streamlines,len(self.blade_rows)))
        massflow = list()
        for indx,row in enumerate(self.blade_rows):
            blade_rows.append(row.to_dict()) # Appending data 
            if row.row_type == RowType.Rotor:
                # Calculation for these are specific to Turbines 
                degree_of_reaction.append(((self.blade_rows[indx-1].P- row.P)/(self.blade_rows[indx-2].P-row.P)).mean())
                
                total_total_efficiency.append(row.eta_total)
                total_static_efficiency.append(row.eta_static)

                stage_loading.append(row.stage_loading)
                euler_power.append(row.euler_power)
            
            if row.row_type!=RowType.Inlet and row.row_type!=RowType.Outlet:
                massflow.append(row.massflow[-1])
            
            for j,p in enumerate(row.percent_hub_shroud):
                t,x,r = self.passage.get_streamline(p)
                x_streamline[j,indx] = float(interp1d(t,x)(row.axial_location))
                r_streamline[j,indx] = float(interp1d(t,r)(row.axial_location))
                
        data = {            
            "blade_rows": blade_rows,
            "massflow":np.mean(massflow),
            "rpm":self.rpm,
            "r_streamline":r_streamline.tolist(),
            "x_streamline":x_streamline.tolist(),
            "rhub":self.passage.rhub_pts.tolist(),
            "rshroud":self.passage.rshroud_pts.tolist(),
            "xhub":self.passage.xhub_pts.tolist(),
            "xshroud":self.passage.xshroud_pts.tolist(),
            "num_streamlines":self.num_streamlines,
            "total-total_efficiency":total_total_efficiency,
            "total-static_efficiency":total_static_efficiency,
            "stage_loading":stage_loading,
            "degree_of_reaction":degree_of_reaction
        }
        # Dump all the Python objects into a single JSON file.
        class NumpyEncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, np.ndarray):
                    return obj.tolist()
                return super().default(obj)

        with open(filename, "w") as f:
            json.dump(data, f, indent=4,cls=NumpyEncoder)


def calculate_massflows(blade_rows:List[BladeRow],calculate_vm:bool=False):
    """Calculates the massflow 

    Args:
        blade_rows (List[BladeRow]): _description_
        passage (Passage): _description_
        calculate_vm (bool, optional): _description_. Defaults to False.
    """
    for p in range(3):
        for i in range(1,len(blade_rows)):
            row = blade_rows[i]
            # Upstream Row
            if i == 0:
                upstream = blade_rows[i]
            else:
                upstream = blade_rows[i-1]
            if i<len(blade_rows)-1:
                downstream = blade_rows[i+1]
            else:
                downstream = None 
            
                # Pressure loss = shift in entropy which affects the total pressure of the row
            if row.row_type == RowType.Inlet:
                row.Yp = 0
            else:
                if row.loss_function.loss_type == LossType.Pressure:
                    row.Yp = row.loss_function(row,upstream)
                    for _ in range(2):
                        if row.row_type == RowType.Rotor:
                            rotor_calc(row,upstream,calculate_vm=True)
                            # Finds Equilibrium between Vm, P0, T0
                            row = radeq(row,upstream) 
                            row = compute_gas_constants(row)
                            rotor_calc(row,upstream,calculate_vm=False)
                        elif row.row_type == RowType.Stator:
                            stator_calc(row,upstream,downstream,calculate_vm=True)
                            # Finds Equilibrium between Vm, P0, T0
                            row = radeq(row,upstream)
                            row = compute_gas_constants(row)
                            stator_calc(row,upstream,downstream,calculate_vm=False)
                        row = compute_gas_constants(row)
                        compute_massflow(row)
                        compute_power(row,upstream)
                        
                elif row.loss_function.loss_type == LossType.Enthalpy: 
                    if row.row_type == RowType.Rotor:
                        row.Yp = 0
                        rotor_calc(row,upstream,calculate_vm=calculate_vm)
                        eta_total = float(row.loss_function(row,upstream))
                        def find_yp(Yp,row,upstream):
                            row.Yp = Yp
                            rotor_calc(row,upstream,calculate_vm=True)
                            row = radeq(row,upstream)
                            row = compute_gas_constants(row)
                            rotor_calc(row,upstream,calculate_vm=False)
                            return abs(row.eta_total - eta_total)
                        
                        res = minimize_scalar(find_yp,bounds=[0,0.6],args=(row,upstream))
                        row.Yp = res.x
                    elif row.row_type == RowType.Stator:
                        row.Yp = 0
                        stator_calc(row,upstream,downstream,calculate_vm=True)
                        row = radeq(row,upstream) 
                        row = compute_gas_constants(row)
                        stator_calc(row,upstream,downstream,calculate_vm=False)
                    row = compute_gas_constants(row)
                    compute_massflow(row)
                    compute_power(row,upstream)

def massflow_loss_function(exit_angle:float,index:int,row:BladeRow,upstream:BladeRow,downstream:BladeRow=None):
    """Finds the blade exit angles that balance the massflow throughout the stage 

    Args:
        exit_angle (float): exit flow angle of the rotor row 
        index (int): streamline index for the current row 
        row (BladeRow): current blade row
        upstream (BladeRow): upstream blade row
        downstream (BladeRow): downstream blade row 

    Returns:
        float: massflow loss 
    """
    # Pressure loss = shift in entropy which affects the total pressure of the row
    if row.row_type == RowType.Inlet:
        row.Yp = 0
    else:
        if row.loss_function.loss_type == LossType.Pressure:
            row.Yp = row.loss_function(row,upstream)
            if row.row_type == RowType.Rotor:
                row.beta2[index] = np.radians(exit_angle)
                rotor_calc(row,upstream)
            elif row.row_type == RowType.Stator:
                row.alpha2[index] = np.radians(exit_angle)
                stator_calc(row,upstream,downstream)
            upstream = compute_gas_constants(upstream)
            row = compute_gas_constants(row)
        elif row.loss_function.loss_type == LossType.Enthalpy:
            # Search for pressure loss that results in the correct total temperature drop
            if row.row_type == RowType.Rotor:
                row.Yp = 0
                row.beta2[index] = np.radians(exit_angle)
                rotor_calc(row,upstream)
                T0_drop = row.loss_function(row,upstream)
                T0_target = row.T0.mean()-T0_drop
                def find_yp(Yp):
                    row.Yp = Yp
                    rotor_calc(row,upstream)
                    upstream = compute_gas_constants(upstream)
                    row = compute_gas_constants(row)
                    return abs(row.T0.mean() - T0_target)
                res = minimize_scalar(find_yp,bounds=[0,0.6])
                row.Yp = res.x
            elif row.row_type == RowType.Stator:
                row.Yp = 0
                row.alpha2[index] = np.radians(exit_angle)
                stator_calc(row,upstream,downstream)
                upstream = compute_gas_constants(upstream)
                row = compute_gas_constants(row)
        

    # if use_radeq:
    #     row = radeq(row,upstream) # Finds Equilibrium between Vm, P0, T0
    
    compute_massflow(row)
    compute_power(row,upstream)
    
    if row.row_type!=RowType.Inlet:
        if row.row_type == RowType.Rotor:
            T3_is = upstream.T0 * (1/row.P0_P)**((row.gamma-1)/row.gamma)
        else:
            T3_is = upstream.T0 * (row.P0/row.P)**((row.gamma-1)/row.gamma)
        a = np.sqrt(row.gamma*row.R*T3_is)
        T03_is = T3_is * (1+(row.gamma-1)/2*(row.Vm/a)**2)
        row.eta_total = (upstream.T0.mean() - row.T0.mean())/(upstream.T0.mean()-T03_is.mean())
        
    return np.abs(row.total_massflow*index/(len(row.massflow)-1) - row.massflow[index])
    
def outlet_pressure(percents:List[float],inletP0:float,outletP:float) -> npt.NDArray:
    """Given a list of percents from 0 to 1 for each row, output each row's outlet static pressure

    Args:
        percents (List[float]): List of floats as percents [[0 to 1],[0 to 1]]
        inletP0 (float): Inlet Total Pressure 
        outletP (float): Outlet Static Pressure

    Returns:
        npt.NDArray: Array of static pressures
    """
    maxP = inletP0
    minP = outletP
    if isinstance(percents, float):
        Ps = [percents*(minP-maxP)+maxP]
    else:
        Ps = np.zeros(shape=(len(percents),1)); i = 0
        for p in percents:
            Ps[i] = p*(minP-maxP)+maxP
            maxP = Ps[i]
            i+=1
    return Ps