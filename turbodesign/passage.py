from typing import List, Optional, Tuple, Union
import numpy as np
import numpy.typing as npt
from scipy.interpolate import PchipInterpolator, interp1d
from pyturbo.helper import line2D
from .enums import PassageType
from scipy.optimize import minimize_scalar
from findiff import FinDiff
from pyturbo.helper import convert_to_ndarray, xr_to_mprime

import matplotlib.pyplot as plt 

class Passage:
    xhub:PchipInterpolator
    rhub:PchipInterpolator
    xhub_pts:npt.NDArray
    rhub_pts:npt.NDArray
    
    xshroud:PchipInterpolator
    rshroud:PchipInterpolator
    xshroud_pts:npt.NDArray
    rshroud_pts:npt.NDArray
    
    n:int
    passageType:PassageType

    x_streamlines:npt.NDArray
    r_streamlines:npt.NDArray
    hub_arc_len:float 
    
    def __init__(self,xhub:Union[npt.NDArray,List[float]],rhub:Union[npt.NDArray,List[float]],
                 xshroud:Union[npt.NDArray,List[float]],rshroud:Union[npt.NDArray,List[float]],
                 passageType:PassageType=PassageType.Axial, zero_phi: bool = False):
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
        self.zero_phi = zero_phi

        hub_arc_len = xr_to_mprime(np.vstack([xhub,rhub]).transpose())[1]
        self.hub_arc_len = hub_arc_len[-1]
        
        self.xhub = PchipInterpolator(hub_arc_len/hub_arc_len[-1],xhub)         # Get the xhub,rhub in terms of the hub arc len
        self.rhub = PchipInterpolator(hub_arc_len/hub_arc_len[-1],rhub)
        self.xshroud = PchipInterpolator(hub_arc_len/hub_arc_len[-1],xshroud)
        self.rshroud = PchipInterpolator(hub_arc_len/hub_arc_len[-1],rshroud)
        
        if len(xhub) < 10:
            self.n = 10
        else:
            self.n = len(xhub)
        
        self.xhub_pts = convert_to_ndarray(xhub) # type: ignore
        self.rhub_pts = convert_to_ndarray(rhub) # type: ignore
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
            xhub = float(self.xhub(t))
            rhub = float(self.rhub(t))
            xshroud = float(self.xshroud(t))
            rshroud = float(self.rshroud(t))
            x_streamline[i] ,r_streamline[i] = line2D((xhub,rhub),(xshroud,rshroud)).get_point(t_radial)
        return t_streamline,x_streamline,r_streamline

    def streamline_curvature(self, x_streamline:npt.NDArray,r_streamline:npt.NDArray) -> Tuple[npt.NDArray,npt.NDArray,npt.NDArray]:
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
        if self.zero_phi:
            return phi, radius_curvature, r_streamline
        # Have to make sure there isn't a divide by zero which could happen if there is a vertical line somewhere
        indices = np.where(np.abs(np.diff(x_streamline))>np.finfo(float).eps)[0]
    
        d_dx = FinDiff(0,x_streamline[indices[0]:indices[-1]],1)
        d2_dx2 = FinDiff(0,x_streamline[indices[0]:indices[-1]],2)
        
        dr_dx = d_dx(r_streamline[indices[0]:indices[-1]]) # type: ignore
        d2r_dx2 = d2_dx2(r_streamline[indices[0]:indices[-1]])     # type: ignore
            
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
        
    def get_area(self,t_hub:float) -> float:
        """Get Area

        Args:
            t_hub (float): Percent arc length along the hub 

        Returns:
            float: Area
        """
        n = 100
        line = self.get_cutting_line(t_hub)[0]
        x,r = line.get_point(np.linspace(0,1,n))
        total_area = 0 
        for j in range(1,n):
            if np.abs((x[-1]-x[0]))<1E-12: # Axial Machines
                total_area += np.pi*(r[j]**2-r[j-1]**2)
            else:   # Radial Machines
                dx = x[j]-x[j-1]
                S = (r[j]-r[j-1])
                C = np.sqrt(1+((r[j]-r[j-1])/dx)**2)
                area = 2*np.pi*C*(S/2*dx**2+r[j-1]*dx)
                total_area += area
        return total_area
        
    def get_cutting_line(self, t_hub:float,t_shroud:Optional[float]=None) -> Tuple[line2D,float,float]:
        """Gets the cutting line perpendicular to hub and shroud 

        Args:
            t_hub (float): percentage along the axial direction 

        Returns:
            (Tuple) containing:
        
                cut (line2D): line from hub to shroud
                t_hub (float): Percentage along hub arc length
                t_shroud (Optional[float]): t corresponding to intersection of bisector of hub. Defaults to None
                                
        """
        xhub = float(self.xhub(t_hub))
        rhub = float(self.rhub(t_hub))
        if t_shroud is None:            
            if t_hub>0 and t_hub<1:
                dx = self.xhub(t_hub+0.0001) - self.xhub(t_hub-0.0001) 
                dr = self.rhub(t_hub+0.0001) - self.rhub(t_hub-0.0001)
            elif t_hub>0:
                dx = self.xhub(t_hub) - self.xhub(t_hub-0.0001) 
                dr = self.rhub(t_hub) - self.rhub(t_hub-0.0001)
            else: # t_hub<1:
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
                    t_shroud = res.x # type: ignore
                else:
                    t_shroud = t_hub # Vertical line 
            else:
                t_shroud = t_hub
                
        xshroud = float(self.xshroud(t_shroud))
        rshroud = float(self.rshroud(t_shroud))
        
        return line2D((xhub,rhub),(xshroud,rshroud)), t_hub, t_shroud # type: ignore
    
    def get_xr_slice(self, t_span: float, percent_hub: Tuple[float, float], 
                     percent_shroud: Optional[Tuple[float, float]] = None, resolution: int = 100) -> npt.NDArray[np.float64]:
        """
        Return the (x, r) coordinates of a *straight* streamline segment that
        connects corresponding hub and shroud points, sampled uniformly along
        each surface between the given percent limits.

        The point returned on each connecting line is at parametric position
        `t_span` in [0, 1], where 0 = hub point and 1 = shroud point.

        Args:
            t_span: Interpolation parameter along each hub→shroud connector (0..1).
            percent_hub: (start, end) fractional arc-length positions along the hub (0..1).
            percent_shroud: Optional (start, end) along the shroud (0..1). If None,
                the shroud uses the same normalized range as `percent_hub`.
            resolution: Number of sample points along the streamwise direction.

        Returns:
            (resolution, 2) array of [x, r] coordinates.
        """
        # ---- validation
        if not (0.0 <= t_span <= 1.0):
            raise ValueError("t_span must be in [0, 1].")
        if resolution < 2:
            raise ValueError("resolution must be >= 2.")
        if not (0.0 <= percent_hub[0] <= 1.0 and 0.0 <= percent_hub[1] <= 1.0):
            raise ValueError("percent_hub values must be in [0, 1].")
        if percent_shroud is not None and not (
            0.0 <= percent_shroud[0] <= 1.0 and 0.0 <= percent_shroud[1] <= 1.0
        ):
            raise ValueError("percent_shroud values must be in [0, 1].")

        # ---- parameterize along hub and shroud (use each surface's own length!)
        t_hub = np.linspace(percent_hub[0], percent_hub[1], resolution) * self.hub_length
        if percent_shroud is None:
            t_shroud = np.linspace(percent_hub[0], percent_hub[1], resolution) * self.shroud_length
        else:
            t_shroud = np.linspace(percent_shroud[0], percent_shroud[1], resolution) * self.shroud_length

        # ---- sample hub & shroud curves (x, r)
        hub_pts = np.column_stack([self.xhub(t_hub), self.rhub(t_hub)])          # (N, 2)
        shroud_pts = np.column_stack([self.xshroud(t_shroud), self.rshroud(t_shroud)])  # (N, 2)

        # ---- vectorized interpolation along each connector: hub + t*(shroud - hub)
        xr = hub_pts + (shroud_pts - hub_pts) * float(t_span)  # (N, 2)

        return xr.astype(np.float64, copy=False)

    
    def get_m(self,t_span:float,resolution:int=100) -> npt.NDArray:
        """Meridional cooridnates

        Args:
            t_span (float): _description_
            resolution (int, optional): _description_. Defaults to 100.

        Returns:
            npt.NDArray: _description_
        """
        xr = self.get_xr_slice(t_span=t_span,percent_hub=(0,1),resolution=resolution)
        dx = np.diff(xr[:,0])
        dr = np.diff(xr[:,1])
        m = np.concat([[0],np.cumsum(np.sqrt(dx**2 + dr**2))])
        return m
    
    def get_dm(self,t_span:float,location:float,resolution:int=1000) -> float:
        """return the derivative in the meridional direction at a particular point

        Args:
            t_span (float): percent span 
            location (float): hub location of the blade
            resolution (int, optional): number of points to represent the hub curve. Defaults to 1000.

        Returns:
            (float) : returns the derivative 
        """
        m = self.get_m(t_span,resolution)
        return PchipInterpolator(np.linspace(0,1,resolution),np.diff(m))(location) # type: ignore
    
    @property
    def hub_length(self):
        """returns the computed length of the hub 
        Returns:
            _type_: _description_
        """
        return np.sum(np.sqrt(np.diff(self.xhub_pts)**2 + np.diff(self.rhub_pts)**2))
    
    @property
    def shroud_length(self):
        """returns the computed length of the shroud 
        Returns:
            _type_: _description_
        """
        return np.sum(np.sqrt(np.diff(self.xshroud_pts)**2 + np.diff(self.rshroud_pts)**2))
    
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
        
