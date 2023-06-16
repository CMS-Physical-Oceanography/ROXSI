import numpy as np
from scipy.optimize import curve_fit

def get_cell_roughness(bath,
                       bath_X,
                       cell_stats=False):
    """
    This function calculates seafloor roughness within a grid cell using 
    high resolution bathymetry (bath) as z0 = std(z_mean - z'). z_mean is a plane 
    given by z_mean = a1*x + a2*y + x and is fit to the hi resolution bathymetry.
    z' are the vertical variations about z_mean and is given by z' = bath - z_mean.
    Seafloor roughness is given by z0 = std(z').
    
    INPUTS:
        bath       - High resolution bathemetry soundings within the grid cell.
        bath_X     - Coordinates of bath where x = bath_X[:,0] and y = bath_X[:,1]
        cell_stats - Defult False. If True the number of high res soundings, the
                     mean slope of the fit and the mean depth of the fit are 
                     returned.
    Outputs:
        rgh        - seafloor roughness.
        mean_h     - Mean depth of z_mean.
        mean_slope - magnitude of the slope of z_mean
        N          - Number of high res soundings in the cell."""
    
        
    def plane(x,a1,a2,c):
        return a1*X[:,0] + a2*X[:,1] +  c
    
    non_nans = ~np.isnan(bath)

    pop,cov = curve_fit(plane,bath_X[non_nans],bath[non_nans])
    
    fit = plane(bath_X,pop[0],pop[1],pop[2])
    rgh = np.std(bath[non_nans]-fit)
    mean_h = np.nanmean(fit)
    mean_slope = np.sqrt(pop[0]**2 + pop[1]**2)
    N = len(bath)
    
    if cell_stats == False:
        return rgh
    else:
        return rgh,mean_h,mean_slope,N
    
