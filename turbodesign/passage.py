from typing import List, Tuple
import numpy as np
import numpy.typing as npt
from scipy.interpolate import interp1d
from pyturbo.helper import line2D
from .enums import PassageType
from scipy.optimize import minimize_scalar
from findiff import FinDiff
from pyturbo.helper import convert_to_ndarray 
import matplotlib.pyplot as plt 

class Passage:
    xhub:interp1d
    rhub:interp1d
    xhub_pts:npt.NDArray
    rhub_pts:npt.NDArray
    
    xshroud:interp1d
    rshroud:interp1d
    xshroud_pts:npt.NDArray
    rshroud_pts:npt.NDArray
    
    n:int
    passageType:PassageType

    x_streamlines:npt.NDArray
    r_streamlines:npt.NDArray
    
    def __init__(self,xhub:List[float],rhub:List[float],
                 xshroud:List[float],rshroud:List[float],
                 passageType:PassageType=PassageType.Axial):
        """_summary_

        Args:
            xhub (List[float]): xhub coordinates
            rhub (List[float]): rhub coordinates 
            xshroud (List[float]): xshroud coordinates 
            rshroud (List[float]): rshroud coordinates 
            passageType (PassageType, optional): Axial or Centrifugal. Defaults to PassageType.Axial.
        """
        assert len(xhub) == len(xshroud), "xHub and xShroud should be the same length"
        assert len(rhub) == len(rshroud), "rHub and rShroud should be the same length"

        self.n = len(xhub)
        t_streamline = np.linspace(0,1,len(xhub))
        self.xhub = interp1d(t_streamline,xhub)
        self.rhub = interp1d(t_streamline,rhub)
        self.xshroud = interp1d(t_streamline,xshroud)
        self.rshroud = interp1d(t_streamline,rshroud)
        
        self.xhub_pts = convert_to_ndarray(xhub)
        self.rhub_pts = convert_to_ndarray(rhub)
        self.xshroud_pts = convert_to_ndarray(xshroud)
        self.rshroud_pts = convert_to_ndarray(rshroud)
        
        self.passageType = passageType
        
    def get_streamline(self,t_radial:float) -> Tuple[npt.NDArray,npt.NDArray, npt.NDArray]:
        """Gets the streamline at a certain percent radius 

        Args:
            t_radial (float): percent between hub and shroud

        Returns:
            Tuple containing:
            
                t_streamline (npt.NDArray): Non dimensional length of the streamline
                x_streamline (npt.NDArray): x-coordinate along the streamline
                r_streamline (npt.NDArray): r-coordinate along the streamline
        """
        t_streamline = np.linspace(0,1,self.n)
        # first streamline is at the hub
        r_streamline = t_streamline.copy()*0
        x_streamline = t_streamline.copy()*0
        for i,t in enumerate(t_streamline):
            xhub = self.xhub(t)
            rhub = self.rhub(t)
            xshroud = self.xshroud(t)
            rshroud = self.rshroud(t)
            x_streamline[i],r_streamline[i] = line2D([xhub,rhub],[xshroud,rshroud]).get_point(t_radial)
        
        return t_streamline,x_streamline,r_streamline

    @staticmethod
    def streamline_curvature(x_streamline:npt.NDArray,r_streamline:npt.NDArray) -> Tuple[npt.NDArray,npt.NDArray,npt.NDArray]:
        """Hub and casing values of streamline angles of inclination and curvature 

            x_streamline[axial,radial]
            r_streamline[axial,radial]
            
        Args:
            x_streamline (np.ndarray): Axial position as a matrix with shape [number of stations, number of x-positions]
            r_streamline (np.ndarray): Annulus Radii of streamlines arranged with shape [number of stations, number of x-positions]

        Returns:
            Tuple: containing

                *phi* (np.ndarray): Array containing angles of inclination for each radi at each station. Rows = radius, columns = station
                *rm* (np.ndarray): Array containing the curvature for each station and streamline
                *r* (np.ndarray): Annulus radius

        References:
            https://stackoverflow.com/questions/28269379/curve-curvature-in-numpy

        """
        phi = np.zeros(shape=x_streamline.shape)
        r  = np.zeros(shape=x_streamline.shape)
        radius_curvature = np.zeros(shape=x_streamline.shape)
        # Have to make sure there isn't a divide by zero which could happen if there is a vertical line somewhere
        indices = np.where(np.abs(np.diff(x_streamline))>np.finfo(float).eps)[0]
    
        d_dx = FinDiff(0,x_streamline[indices[0]:indices[-1]],1)
        d2_dx2 = FinDiff(0,x_streamline[indices[0]:indices[-1]],2)
        dr_dx = d_dx(r_streamline[indices[0]:indices[-1]])
        d2r_dx2 = d2_dx2(r_streamline[indices[0]:indices[-1]])    
            
        radius_curvature[indices[0]:indices[-1]] = np.power((1+np.power(dr_dx,2)),1.5)
        radius_curvature[indices[0]:indices[-1]] = np.divide(radius_curvature[indices[0]:indices[-1]], np.abs(d2r_dx2))
        radius_curvature = np.nan_to_num(radius_curvature,nan=0)
        
        def vertical_line_phi(start:int,end:int):
            # Lets get the phi and slope for the vertical parts 
            for i in range(start,end):     
                dx = x_streamline[i] - x_streamline[i-1]
                dr = r_streamline[i] - r_streamline[i-1]
                if (dr < 0) & (np.abs(dx) < np.finfo(float).eps):
                    phi[i-1] = -np.pi/2
                elif (dr > 0) & (np.abs(dx) < np.finfo(float).eps):
                    phi[i-1] = np.pi/2
                radius_curvature[i-1] = 1000000 # Initialize to high number, used in radeq. 
                
        vertical_line_phi(1,indices[0])
        vertical_line_phi(indices[1],len(x_streamline))
        
        rm = radius_curvature     # https://www.cuemath.com/radius-of-curvature-formula/ should be 1/curvature
        phi[indices[0]:indices[-1]] = np.arctan(dr_dx)
        r = r_streamline
            
        return phi, rm, r
        
    def get_cutting_line(self, t_hub:float) -> line2D:
        """Gets the cutting line perpendicular to hub and shroud 

        Args:
            t_hub (float): percentage along the axial direction 

        Returns:
            (Tuple) containing:
        
                cut (line2D): line from hub to shroud
                t_hub (float): t corresponding to xhub location
                t_shroud (float): t corresponding to intersection of bisector of hub 
                                
        """
        xhub = self.xhub(t_hub)
        rhub = self.rhub(t_hub)
        
        if t_hub>0 and t_hub<1:
            dx = self.xhub(t_hub+0.0001) - self.xhub(t_hub-0.0001) 
            dr = self.rhub(t_hub+0.0001) - self.rhub(t_hub-0.0001)
        elif t_hub>0:
            dx = self.xhub(t_hub) - self.xhub(t_hub-0.0001) 
            dr = self.rhub(t_hub) - self.rhub(t_hub-0.0001)
        elif t_hub<1:
            dx = self.xhub(t_hub+0.0001) - self.xhub(t_hub)
            dr = self.rhub(t_hub+0.0001) - self.rhub(t_hub)
        
        if self.passageType == PassageType.Centrifugal:
            if np.abs(dr)>1e-6:
                # Draw a line perpendicular to the hub. 
                # Find the intersection point to the shroud. 
                h = -dx/dr # Slope of perpendicular line
            
                f = lambda t: h*(self.xshroud(t) - xhub)+rhub # line from hub to shroud 
                fun = lambda t: np.abs(f(t)-self.rshroud(t)) # find where it intersects
                res = minimize_scalar(fun,bounds=[0,1],tol=1E-3) 
                t_shroud = res.x
            else:
                t_shroud = t_hub # Vertical line 
        else:
            t_shroud = t_hub
        
        xshroud = self.xshroud(t_shroud)
        rshroud = self.rshroud(t_shroud)
        return line2D([xhub,rhub],[xshroud,rshroud]), t_hub, t_shroud
    
    def get_xr_slice(self,t_span:float,axial_location:Tuple[float,float]):
        """Returns the xr coordinates of a streamline, a line that is parallel to both hub and shroud
            
        Args:
            t_span (float): _description_
            axial_location (float): _description_

        Returns:
            np.NDArray: _description_
        """
        t_hub = np.linspace(axial_location[0],axial_location[1],100)
        
        shroud_pts_cyl = np.vstack([self.xshroud(t_hub),self.rshroud(t_hub)]).transpose()
        hub_pts_cyl = np.vstack([self.xhub(t_hub),self.rhub(t_hub)]).transpose()
        n = len(t_hub)
            
        xr = np.zeros((n,2))
        for j in range(n):
            l = line2D(hub_pts_cyl[j,:],shroud_pts_cyl[j,:])
            xr[j,0],xr[j,1] = l.get_point(t_span)
        return xr
        
    @property
    def hub_length(self):
        """returns the computed length of the hub 
        Returns:
            _type_: _description_
        """
        return np.sum(np.sqrt(np.diff(self.xhub_pts)**2 + np.diff(self.rhub_pts)**2))
    
    def plot_cuts(self,percent_axial:List[float]=[]):
        """_summary_

        Args:
            percent_axial (List[float], optional): _description_. Defaults to [].
        """
        
        plt.figure(num=1,clear=True,dpi=150,figsize=(15,10))
        plt.plot(self.xhub_pts,self.rhub_pts,label='hub',linestyle='solid',linewidth=2,color='black')
        plt.plot(self.xshroud_pts,self.rshroud_pts,label='shroud',linestyle='solid',linewidth=2,color='black')
        for p in percent_axial:
            cut,_,_ = self.get_cutting_line(p)
            x,r = cut.get_point(np.linspace(0,1,10))
            plt.plot(x,r,label=f'{p}',linestyle='dashed')
        
        
        plt.ylim([-self.rshroud_pts.max()*0.1, self.rshroud_pts.max()])
        plt.legend()
        plt.axis('scaled')
        plt.show()
        
