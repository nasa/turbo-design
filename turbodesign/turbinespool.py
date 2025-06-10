from typing import List, Optional
from cantera.composite import Solution
from .bladerow import BladeRow, interpolate_streamline_radii
from .enums import RowType, MassflowConstraint, LossType, PassageType
from .spool import Spool
import json, copy
from .passage import Passage
from scipy.interpolate import interp1d
import numpy as np
import numpy.typing as npt
from .td_math import inlet_calc,rotor_calc, stator_calc, compute_massflow, compute_power, compute_gas_constants, compute_reynolds
from .solve_radeq import adjust_streamlines, radeq
from scipy.optimize import minimize_scalar, differential_evolution, fmin_slsqp
from .inlet import Inlet
from .outlet import Outlet
from pyturbo.helper import convert_to_ndarray

class TurbineSpool(Spool):
    
    def __init__(self,passage:Passage,
                 massflow:float,rows:List[BladeRow],
                 num_streamlines:int=3,
                 fluid:Optional[Solution]=Solution('air.yaml'),
                 rpm:float=-1,
                 massflow_constraint:MassflowConstraint=MassflowConstraint.MatchMassFlow): # type: ignore
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
        if self.fluid:
            inlet.initialize_fluid(self.fluid) # type: ignore
        else:
            inlet.initialize_fluid(R=self.blade_rows[1].R, # type: ignore
                                    gamma=self.blade_rows[1].gamma,
                                    Cp=self.blade_rows[1].Cp)
        
        inlet.total_massflow = W0 
        inlet.total_massflow_no_coolant = W0
        inlet.massflow = np.linspace(0,1,self.num_streamlines)*W0
        
        inlet.initialize_inputs(self.num_streamlines)
        inlet.initialize_velocity(self.passage,self.num_streamlines)     # type: ignore
        interpolate_streamline_radii(inlet,self.passage,self.num_streamlines)

        compute_gas_constants(inlet,self.fluid)
        inlet_calc(inlet)
        
        for row in self.blade_rows:
            interpolate_streamline_radii(row,self.passage,self.num_streamlines)

        outlet = self.blade_rows[-1]
        for j in range(self.num_streamlines):
            P0 = inlet.get_total_pressure(inlet.percent_hub_shroud[j]) # type: ignore
            percents = np.zeros(shape=(len(self.blade_rows)-2)) + 0.3
            percents[-1] = 1
            Ps_range = outlet_pressure(percents=percents,inletP0=inlet.P0[j],outletP=outlet.P[j]) # type: ignore
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
                T0c = self.blade_rows[i].coolant.T0
                P0c = self.blade_rows[i].coolant.P0
                W0c = self.blade_rows[i].coolant.massflow_percentage * self.massflow
                Cpc = self.blade_rows[i].coolant.Cp
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
                stator_calc(row,upstream,downstream) # type: ignore
                compute_massflow(row)
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
            self.__match_massflow()  # Matches massflow by changing turning angle
        elif self.massflow_constraint == MassflowConstraint.BalanceMassFlow:
            self.__balance_massflow()   # Balances massflow by changing row exit static pressure

    
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
                else:# row.row_type == RowType.Rotor:
                    bounds = [-80,0]
                if row.row_type != RowType.Inlet:
                    for j in range(1,self.num_streamlines):
                        res = minimize_scalar(massflow_loss_function, bounds=bounds,args=(j,row,upstream,downstream),tol=1E-3) 
                        if row.row_type == RowType.Rotor:
                            row.beta2[j] = np.radians(res.x)
                             # Initialize the value at the hub to not upset the mean
                            row.beta2[0] = 1/(len(row.beta2)-1)*row.beta2[1:].sum()
                        elif row.row_type == RowType.Stator:
                            row.alpha2[j] = np.radians(res.x)
                            row.alpha2[0] = 1/(len(row.alpha2)-1)*row.alpha2[1:].sum()
                    compute_gas_constants(upstream,self.fluid)
                    compute_gas_constants(row,self.fluid)
                
                    
            # Step 3: Adjust streamlines to evenly divide massflow
            adjust_streamlines(self.blade_rows,self.passage)
        compute_reynolds(self.blade_rows,self.passage)
    
    @staticmethod # Private static method
    def __massflow_std__(blade_rows:List[BladeRow]):
        """Returns the standard deviation of the massflow

        Args:
            blade_rows (List[BladeRow]): List of blade rows

        Returns:
            _type_: _description_
        """
        total_massflow = list(); s = 0; massflow_stage = list()
        stage_ids = list(set([row.stage_id for row in blade_rows if row.stage_id>=0]))
        for row in blade_rows: # Ignore inlet and outlet
            total_massflow.append(row.total_massflow_no_coolant)
            sign = 1
            for s in stage_ids:
                for row in blade_rows:
                    if row.stage_id == s and row.row_type == RowType.Rotor:
                        massflow_stage.append(sign*row.total_massflow_no_coolant)
                        sign*=-1
            if len(stage_ids) % 2 == 1:
                massflow_stage.append(massflow_stage[-1]*sign)
        deviation = np.std(total_massflow)*2
        if deviation>1.0:
            print("high massflow deviation detected")
        return np.std(total_massflow)*2 # + abs(sum(massflow_stage))  # Equation 28  
    
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
        def balance_massflows(x0:List[float],blade_rows:List[BladeRow],P0:npt.NDArray,P:npt.NDArray,balance_mean_pressure:bool=True):
            """Balance Massflows. 
            
            Steps:
                1. Balance the mean static pressure in between the blade rows. X0 = [0.2,0.5,...] size = num_rows
                    P = outlet static pressure 
                    
                2. Keep the mean 

            Args:
                x0 (List[float]): Percentage of P0 exiting each row
                blade_rows (List[List[BladeRow]]): _description_
                P0 (npt.NDArray): _description_
                P (npt.NDArray): (1) Outlet Static Pressure. (2) 
                balance_mean_pressure (bool, optional): _description_. Defaults to True.

            Returns:
                _type_: _description_
            """
            # blade_rows_backup = copy.deepcopy(blade_rows)
            # try:
            if balance_mean_pressure:
                for j in range(self.num_streamlines):
                    Ps = outlet_pressure(x0,P0[j],P[j])
                    for i in range(1,len(blade_rows)-2):
                        blade_rows[i].P[j] = float(Ps[i-1]) # type: ignore
                blade_rows[-2].P = P # type: ignore
            else:
                for i in range(1,len(blade_rows)-1):
                    for j in range(self.num_streamlines):
                        blade_rows[i].P[j] = P[j]*x0[(i-1)*self.num_streamlines+j]    # type: ignore # x0 size = num_streamlines -1 
            # try:  
            calculate_massflows(blade_rows,True,self.fluid)
            print(x0)
            return self.__massflow_std__(blade_rows[1:-1]) # do not consider inlet and outlet
            # except Exception as e:
            #     print(e)
            # finally:
            #     blade_rows = blade_rows_backup
            #     return np.inf # Return a high error
            

        # Break apart the rows to stages
        outlet_P=list(); outlet_P_guess = list() # Outlet P is the bounds, outlet_p_guess is the guessed values 
        
        for i in range(1,len(self.blade_rows)-2):
            outlet_P.append(self.blade_rows[i].inlet_to_outlet_pratio)
            outlet_P_guess.append(np.mean(self.blade_rows[i].inlet_to_outlet_pratio))
        
        print(f"Looping to converge massflow")
        past_err = -100; loop_iter = 0; err = 0.001
        while (np.abs((err-past_err)/err)>0.05) and loop_iter<10:
            if len(outlet_P) == 1:
                # x = balance_massflows(0.22896832148169688,self.blade_rows,self.blade_rows[0].P0,self.blade_rows[-1].P)
                res = minimize_scalar(fun=balance_massflows,args=(self.blade_rows,self.blade_rows[0].P0,self.blade_rows[-1].P),bounds=outlet_P[0],tol=0.001,options={'disp': True},method='bounded')
                x = res.x
                print(x)
            else:
                x = fmin_slsqp(func=balance_massflows,args=(self.blade_rows,self.blade_rows[0].P0,self.blade_rows[-1].P), 
                            bounds=outlet_P, x0=outlet_P_guess,epsilon=0.001,iter=100) # ,tol=0.001,options={'disp': True})
                outlet_P_guess = x 
        
            # Adjust the inlet: Set the massflow
            self.blade_rows[0].massflow = np.linspace(0,1,self.num_streamlines)*self.blade_rows[1].total_massflow_no_coolant
            self.blade_rows[0].total_massflow_no_coolant = self.blade_rows[1].total_massflow_no_coolant
            self.blade_rows[0].total_massflow = self.blade_rows[1].total_massflow_no_coolant
            self.blade_rows[0].calculated_massflow = self.blade_rows[0].total_massflow_no_coolant
            inlet_calc(self.blade_rows[0]) # adjust the inlet to match massflow 
        
            if self.adjust_streamlines:
                adjust_streamlines(self.blade_rows[:-1],self.passage)
                
            self.blade_rows[-1].transfer_quantities(self.blade_rows[-2])    # This would be the outlet
            self.blade_rows[-1].P = self.blade_rows[-1].get_static_pressure(self.blade_rows[-1].percent_hub_shroud)
            
            past_err = err
            err = self.__massflow_std__(self.blade_rows)
            loop_iter += 1 
            print(f"Loop {loop_iter} massflow convergenced error:{err}")
        
        # calculate Reynolds number
        compute_reynolds(self.blade_rows,self.passage)
        
    
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
        enthalpy_power = list()
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
                enthalpy_power.append(row.power)
            if row.row_type!=RowType.Inlet and row.row_type!=RowType.Outlet:
                massflow.append(row.massflow[-1])
            
            for j,p in enumerate(row.percent_hub_shroud):
                t,x,r = self.passage.get_streamline(p)
                x_streamline[j,indx] = float(interp1d(t,x)(row.percent_hub))
                r_streamline[j,indx] = float(interp1d(t,r)(row.percent_hub))
        
        Pratio_Total_Total = np.mean(self.blade_rows[0].P0 / self.blade_rows[-2].P0)
        Pratio_Total_Static = np.mean(self.blade_rows[0].P0 / self.blade_rows[-2].P)
        FlowFunction = np.mean(massflow)*np.sqrt(self.blade_rows[0].T0)*self.blade_rows[0].P0/1000 # kg sqrt(K)/(sec kPa)
        CorrectedSpeed = self.rpm * np.pi/30 / np.sqrt(self.blade_rows[0].T0.mean())   # rad/s * 1/sqrt(K)
        EnergyFunction = (self.blade_rows[0].T0 - self.blade_rows[-2].T0) * 0.5* (self.blade_rows[0].Cp + self.blade_rows[-2].Cp) / self.blade_rows[0].T0 # J/(KgK)
        EnergyFunction = np.mean(EnergyFunction)
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
            "euler_power": euler_power,
            "enthalpy_power":enthalpy_power,
            "total-total_efficiency":total_total_efficiency,
            "total-static_efficiency":total_static_efficiency,
            "stage_loading":stage_loading,
            "degree_of_reaction":degree_of_reaction,
            "Pratio_Total_Total":Pratio_Total_Total,
            "Pratio_Total_Static":Pratio_Total_Static,
            "FlowFunction":FlowFunction,
            "CorrectedSpeed":CorrectedSpeed,
            "EnergyFunction":EnergyFunction
        }
        # Dump all the Python objects into a single JSON file.
        class NumpyEncoder(json.JSONEncoder):
            def default(self, obj):
                if isinstance(obj, np.ndarray):
                    return obj.tolist()
                return super().default(obj)

        with open(filename, "w") as f:
            json.dump(data, f, indent=4,cls=NumpyEncoder)


def calculate_massflows(blade_rows:List[BladeRow],calculate_vm:bool=False,fluid:Optional[Solution]=None):
    """Calculates the massflow 

    Args:
        blade_rows (List[BladeRow]): _description_
        passage (Passage): _description_
        calculate_vm (bool, optional): _description_. Defaults to False.
    """
    for i in range(1,len(blade_rows)-1):
        row = blade_rows[i]
        # Upstream Row
        if i == 0:
            upstream = blade_rows[i]
        else:
            upstream = blade_rows[i-1]
        if i<len(blade_rows)-1:
            downstream = blade_rows[i+1]
        
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
                        row = radeq(row,upstream,downstream) 
                        compute_gas_constants(row,fluid)
                        rotor_calc(row,upstream,calculate_vm=False)
                    elif row.row_type == RowType.Stator:
                        stator_calc(row,upstream,downstream,calculate_vm=True)
                        # Finds Equilibrium between Vm, P0, T0
                        row = radeq(row,upstream,downstream)
                        compute_gas_constants(row,fluid)
                        stator_calc(row,upstream,downstream,calculate_vm=False)
                    compute_gas_constants(row,fluid)
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
                        compute_gas_constants(row,fluid)
                        rotor_calc(row,upstream,calculate_vm=False)
                        return abs(row.eta_total - eta_total)
                    
                    res = minimize_scalar(find_yp,bounds=[0,0.6],args=(row,upstream))
                    row.Yp = res.x
                elif row.row_type == RowType.Stator:
                    row.Yp = 0
                    stator_calc(row,upstream,downstream,calculate_vm=True)
                    row = radeq(row,upstream) 
                    row = compute_gas_constants(row,fluid)
                    stator_calc(row,upstream,downstream,calculate_vm=False)
                row = compute_gas_constants(row,fluid)
                compute_massflow(row)
                compute_power(row,upstream)
    
def massflow_loss_function(exit_angle:float,index:int,row:BladeRow,upstream:BladeRow,downstream:BladeRow=None,fluid:Solution=None):
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
            upstream = compute_gas_constants(upstream,fluid)
            row = compute_gas_constants(row,fluid)
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
                    upstream = compute_gas_constants(upstream,fluid)
                    row = compute_gas_constants(row,fluid)
                    return abs(row.T0.mean() - T0_target)
                res = minimize_scalar(find_yp,bounds=[0,0.6])
                row.Yp = res.x
            elif row.row_type == RowType.Stator:
                row.Yp = 0
                row.alpha2[index] = np.radians(exit_angle)
                stator_calc(row,upstream,downstream)
                upstream = compute_gas_constants(upstream,fluid)
                row = compute_gas_constants(row,fluid)
        

    # if use_radeq:
    #     row = radeq(row,upstream) # Finds Equilibrium between Vm, P0, T0
    
    compute_massflow(row)
    compute_power(row,upstream)
    
    if row.row_type!=RowType.Inlet:
        # if row.row_type == RowType.Rotor:
        T3_is = upstream.T0 * (1/row.P0_P)**((row.gamma-1)/row.gamma)
        # else:
        #     T3_is = upstream.T0 * (row.P0/row.P)**((row.gamma-1)/row.gamma)
        a = np.sqrt(row.gamma*row.R*T3_is)
        T03_is = T3_is * (1+(row.gamma-1)/2*(row.V/a)**2)
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
    percents = convert_to_ndarray(percents)
    Ps = np.zeros((len(percents),))
    for i in range(len(percents)):
        Ps[i] = float(interp1d((0,1),(inletP0,outletP))(percents[i]))
        inletP0 = Ps[i] 
    return Ps
    