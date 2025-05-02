from pyturbo.helper import arc, line2D, bezier
import numpy as np 

def build_endwalls(radius:float,
                   hub_outlet_radius_scale:float=1.1,
                   shroud_inlet_offset:float=0.1,
                   shroud_outlet_offset_ratio:float=1.2,
                   rhub_out:float=0.01,
                   x_stretch_factor:float=1,
                   inlet_ext_percent:float=0.1,
                   outlet_ext_percent=0.2):
    """Build endwalls

    Args:
        radius (float): radius of the hub 
        hub_outlet_radius_scale (float, optional): Amount to grow the radius. Defaults to 1.1.
        shroud_inlet_offset (float, optional): Initial hub to shroud offset. Defaults to 0.1.
        shroud_outlet_offset_ratio (float, optional): Final hub to shroud offset. Defaults to 1.2.
        rhub_out (float, optional): how much to shift the hub radius by above the rotation axis. Defaults to 0.01.
        x_stretch_factor (float, optional): _description_. Defaults to 1.
        inlet_ext_percent (float, optional): _description_. Defaults to 0.1.
        outlet_ext_percent (float, optional): _description_. Defaults to 0.2.

    Returns:
        _type_: _description_
    """
    # Build the hub then offset the shroud 
    b = bezier([0,0.1,0.8,1],[1, 1, hub_outlet_radius_scale, hub_outlet_radius_scale])
    _,ratio = b.get_point(np.linspace(0,1,1000))
    
    hub = arc(xc=0,yc=0,radius=radius*ratio,alpha_start=180,alpha_stop=270)
    [xhub,rhub] = hub.get_point(np.linspace(0,1,1000))
    xhub*=x_stretch_factor

    hub = np.vstack([xhub,rhub]).transpose()
    hub[:,1] += -hub[:,1].min() + rhub_out
    xhub += -xhub.min()
    
    # Build the shroud by offsetting from the hub 
    offset_bezier = bezier([0,0.1,0.8,1],[shroud_inlet_offset, shroud_inlet_offset, shroud_outlet_offset_ratio*shroud_inlet_offset, shroud_outlet_offset_ratio*shroud_inlet_offset])
    
    offset = offset_bezier.get_point(np.linspace(0,1,hub.shape[0]))[1]
    shroud = offset_curve(hub[:,0],hub[:,1],offset)
    
    hub_inlet_ext = np.array([[hub[0,0],hub[0,1]+radius*inlet_ext_percent],
                                [hub[0,0],hub[0,1]]])
    hub_inlet_ext_pts = np.array(line2D(hub_inlet_ext[0,:],hub_inlet_ext[1,:]).get_point(np.linspace(0,1,500))).transpose()
    
    shroud_inlet_ext = np.array([[shroud[0,0],shroud[0,1]+radius*inlet_ext_percent],
                                    [shroud[0,0],shroud[0,1]]])
    shroud_inlet_ext_pts = np.array(line2D(shroud_inlet_ext[0,:],shroud_inlet_ext[1,:]).get_point(np.linspace(0,1,500))).transpose()
    
    # Extend Outlet
    hub_outlet_ext = np.array([[hub[-1,0],hub[-1,1]],
                                [hub[-1,0]+radius*outlet_ext_percent,hub[-1,1]]])
    hub_outlet_ext_pts = np.array(line2D(hub_outlet_ext[0,:],hub_outlet_ext[1,:]).get_point(np.linspace(0,1,500))).transpose()
    
    shroud_outlet_ext = np.array([[shroud[-1,0],shroud[-1,1]],
                                    [shroud[-1,0]+radius*outlet_ext_percent,shroud[-1,1]]])
    shroud_outlet_ext_pts = np.array(line2D(shroud_outlet_ext[0,:],shroud_outlet_ext[1,:]).get_point(np.linspace(0,1,500))).transpose()
    
    hub_new = np.vstack([hub_inlet_ext_pts, hub[1:-1,:], hub_outlet_ext_pts])
    shroud_new = np.vstack([shroud_inlet_ext_pts, shroud[1:-1,:], shroud_outlet_ext_pts])

  
    return hub_inlet_ext_pts,hub,hub_outlet_ext_pts,shroud_inlet_ext_pts,shroud,shroud_outlet_ext_pts

def compute_normals(x, y):
    # Compute first derivatives
    dx = np.gradient(x)
    dy = np.gradient(y)
    
    # Compute normal vectors (perpendicular to tangent)
    length = np.hypot(dx, dy)
    nx = -dy / length
    ny = dx / length
    
    return nx, ny

def offset_curve(x, y, offset_distance):
    nx, ny = compute_normals(x, y)

    # Offset points along the normal direction
    x_offset = x + offset_distance * nx
    y_offset = y + offset_distance * ny
    
    return np.vstack([x_offset, y_offset]).transpose()
    