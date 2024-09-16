from typing import List, Tuple, Union
import pandas as pd 
import numpy as np 
import numpy.typing as npt
from scipy.interpolate import bisplrep, bisplev, interp1d, LSQBivariateSpline
import matplotlib.pyplot as plt 


class LossInterp:
    """Interpret the loss of an XY Graph with an additional 3rd variable 
    """
    df: pd.DataFrame
    x:np.ndarray
    y:np.ndarray
    c:np.ndarray    

    x_max_c:np.ndarray
    x_min_c:np.ndarray
    
    c_max:float
    c_min:float

    fxc_max: interp1d
    fxc_min: interp1d

    weights: np.ndarray
    xlabel:str
    ylabel:str
    clabel:str
    _name:str 

    is_xy:bool
    func: interp1d
    
    xlogscale:bool
    
    @property
    def name(self):
        return self._name

    @name.setter
    def name(self,val:str):
        self._name = val 

    def __init__(self,csv_filename:str,xlabel:str="",ylabel:str="",clabel:str="",logx10:bool=False):
        """Initialize the loss interpolation with data

        Note:
            The data contains x,y, and c. Y is predicted from x and c
        Args:
            csv_filename (str): csv data with 3 columns, x,y,c.
            xlabel (str, optional): xlabel. Defaults to "".
            ylabel (str, optional): ylabel. Defaults to "".
            clabel (str, optional): clabel. Defaults to "".
            logx10 (bool, optional): apply log base 10 to x axis. Defaults to False
        """
        self.df = pd.read_csv(csv_filename)
        self.df = self.df.sort_values(self.df.columns[0],ascending=True)
        self.name = csv_filename
        data = self.df.to_numpy()
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.logX10 = logx10
        if data.shape[1] == 3:    
            self.x = data[:,0]
            self.y = data[:,1]
            self.c = data[:,2]
            
            self.clabel = clabel
            
            # Find the minimum x value for every item in row
            self.x_max_c = np.zeros(shape=(len(self.df.iloc[:,2].unique()),2))
            self.x_min_c = np.zeros(shape=(len(self.df.iloc[:,2].unique()),2))

            i = 0 
            for inlet_flow in self.df.iloc[:,2].unique():
                df_slice = self.df[self.df.iloc[:,2] == inlet_flow]
                df_slice = df_slice.sort_values(df_slice.columns[0],ascending=True)
                self.x_max_c[i,0] = inlet_flow
                self.x_max_c[i,1] = df_slice.iloc[:,0].max()
                
                self.x_min_c[i,0] = inlet_flow
                self.x_min_c[i,1] = df_slice.iloc[:,0].min()
                i+=1 
            self.fxc_max = interp1d(self.x_max_c[:,0],self.x_max_c[:,1]) # xmax as a function of c
            self.fxc_min = interp1d(self.x_min_c[:,0],self.x_min_c[:,1]) # xmin as a function of c

            # self.weights = bisplrep(self.x,self.c,self.y,kx=5, ky=5)
            if self.logX10:
                self.func = LSQBivariateSpline(np.log10(self.x),self.c,self.y,np.log10(np.unique(self.x[:4])),np.unique(self.c[:4]))
            else:
                self.func = LSQBivariateSpline(self.x,self.c,self.y,np.unique(self.x[:3]),np.unique(self.c[:3]))
            self.is_xy = False # 3 columns
            self.c_max = self.df.to_numpy()[:,-1].max()
            self.c_min = self.df.to_numpy()[:,-1].min()
        else:
            self.is_xy = True
            self.x = data[:,0]
            self.y = data[:,1]
            if self.logX10:
                self.func = interp1d(np.log10(self.x),self.y)
            else:
                self.func = interp1d(self.x,self.y)

    def __call__(self,x:Union[npt.NDArray,float],c:Union[npt.NDArray,float,None]=None) -> float:
        """Pass in an array of x and c values

        Args:
            x (Union[npt.NDArray,float]): value on x axis 
            c (Union[npt.NDArray,float,None]): Third axis value

        Returns:
            float: y
        """
            
        if self.is_xy or c==None:
            if self.logX10:
                x = np.log10(x)
            if isinstance(x, np.ndarray):
                return self.func(x)
            elif isinstance(x, float):
                return float(self.func(x))
        else:
            if isinstance(x, np.ndarray):
                xmax = self.fxc_max(c)
                xmin = self.fxc_min(c)
                x[x>xmax] = xmax
                x[x<xmin] = xmin
                if self.logX10:
                    y = self.func(np.log10(x),c)
                else:
                    y = self.func(x,c)
            elif isinstance(x, float):
                if (c>self.c_max):
                    c = self.c_max
                elif c<self.c_min:
                    c = self.c_min
                xmax = float(self.fxc_max(c))
                xmin = float(self.fxc_min(c))
                x = xmin if x<xmin else x
                x = xmax if x>xmax else x
                # y[j] = bisplev(x[j],cc,self.weights)
                if self.logX10:
                    y = self.func(np.log10(x),c)[0][0] # type: ignore
                else:
                    y = self.func(x,c)[0][0] # type: ignore
        return y
    
    def plot(self):
        """Plot the data with predicted values
        """
        plt.figure(num=1,figsize=(10,6)) 
        if not self.is_xy:
            c_unique = np.unique(self.c)
            graycolors = plt.cm.gray(np.linspace(0,1,len(c_unique)))
            coolcolors = plt.cm.cool(np.linspace(0,1,len(c_unique)))

        
            # Plot the actual data 
            i = 0
            for c in c_unique:
                df2 = self.df[self.df.iloc[:,2] == c]
                x = df2.iloc[:,0].to_numpy()
                plt.plot(x,df2.iloc[:,1],'.-',label=f'{c}',color=graycolors[i])
                y = self(x,c)
                plt.plot(x,y,'o-',label=f'predicted-{c}',color=coolcolors[i])
                i+=1
        else:
            x = self.df.iloc[:,0].to_numpy()            
            plt.plot(x,self.df.iloc[:,1],'ro',linewidth=2,label='actual')
            y = self(x)
            plt.plot(x,y,'b-',label='predicted')
            
        plt.ylabel(self.ylabel)
        plt.xlabel(self.xlabel)
        plt.title(self.name)
        plt.legend()
        plt.show()
