import numpy as np
from .cell_roughness import get_cell_roughness
from .grab_data import get_data,get_cell_data




def grid_roughness(grid_lat,
                   grid_lon,
                   bath,
                   bath_lat,
                   bath_lon,
                   point_thresh,
                   parallel = True):
    """
    This function calculates roughness for the 100m and 500m grids.
    
    Inputs:
        grid_lat (2d array) - lat coordinates of grid.
        grid_lon (2d array) - lon coordinates of grid.
        bath     (1d array) - bathymetry soundings.
        bath_lat (1d array) - latitude of soundings.
        bath_lon (1d array) - longitude of soundings.
        point_thresh        - minimum number of soundings within a cell
                              for values to be assigned.
        parallel (default = True) - True if verticies of grid cells run
                                    parallel with lines of latitude and 
                                    longitude. 
    Outputs:
        roughness_grid (2d array) - gridded cell roughness as std(z')
        mean_h         (2d array) - mean grid cell depth.
        n              (2d array) - number of points within grid cell."""
    
    roughness_grid = np.zeros((grid_lat.shape[0]-1,grid_lat.shape[1]-1))
    mean_h = np.zeros_like(roughness_grid)
    n = np.zeros_like(roughness_grid)

    for row in range(roughness_grid.shape[0]):
        for col in range(roughness_grid.shape[1]):
   
            if parallel == True:
                max_lat = grid_lat[row,col]
                min_lat = grid_lat[row+1,col]
                max_lon = grid_lon[row,col+1]
                min_lon = grid_lon[row,col]
                
                cell_h,cell_lats,cell_lons = get_data(bath,
                                                      bath_lat,
                                                      bath_lon,
                                                      min_lat,
                                                      max_lat,
                                                      min_lon,
                                                      max_lon)
            elif parallel == False:
                
                cell_h,cell_X = get_cell_data(grid_lat,
                                              grid_lon,
                                              row,
                                              col,
                                              bath,
                                              np.column_stack((bath_lon,bath_lat)))
                cell_lats = cell_X[:,1]
                cell_lons = cell_X[:,0]
            else:
                raise ValueError('parallel must be boolian value')
            
                
            if len(cell_h[~np.isnan(cell_h)]) < point_thresh:   
                
                roughness_grid[row,col],mean_h[row,col],n[row,col] = [np.nan for i in range(3)]
            else:
                roughness_grid[row,col],mean_h[row,col],n[row,col] = get_cell_roughness(cell_h,cell_lats,cell_lons)

            
    return roughness_grid,mean_h,n




def get_slopes(h,res=2):
    """
    This function calculates the magnitude 
    of the slope between bathymetry points h
    spaced res units apart.
    
    Inputs:
        h - 2d array of bathymetry soundings.
        res - horizontal spacing of soundings
    Outputs:
        abs_slope - magnitude of the slope between
                    soundings."""
    
    east_slopes = (h[:,2:] - h[:,:-2])/2
    north_slopes = (h[2:,:] - h[:-2,:])/2
    
    abs_slope_ = np.sqrt(north_slopes[:,1:-1]**2+east_slopes[1:-1,:]**2)
    abs_slope = np.full((h.shape[0],h.shape[1]),np.nan)
    abs_slope[1:-1,1:-1] = abs_slope_
    
    return abs_slope


def grid_slopes(grid_lat,
                grid_lon,
                slopes,
                slopes_X):
    
    """
    This function calculated the standerd deviation and mean of the 
    z' slopes within grid cells.
    
    Used for the 100m and 500m grid.
    
    Inputs:
        - grid_lat (2d array) - latitude points of the grid.
        - grid_lon (2d array) - longitude points of the grid.
        - slopes   (1d array) - slopes between hi res bathymetry points.
        - slopes_X (1d array) - coordinates of slopes where slopes_X[:,0] = x
                                and slopes_X[:,1] = y.
    Outputs:
        - zprime_slopes - mean slope between z' points in each cell.
        - zprime_slope_stds - standerd deviation of slopes between z'
                              points in each cell."""
    
    zprime_slopes = np.zeros((grid_lat.shape[0]-1,grid_lat.shape[1]-1))
    zprime_slope_stds = np.zeros_like(zprime_slopes)
    
    for row in range(zprime_slopes.shape[0]):
        for col in range(zprime_slopes.shape[1]):
            
            max_lat = grid_lat[row,col]
            min_lat = grid_lat[row+1,col]
       
            max_lon = grid_lon[row,col+1]
            min_lon = grid_lon[row,col]
   

            cell_slopes,cell_X = get_data(slopes,
                                          slopes_X,
                                          min_lat,
                                          max_lat,
                                          min_lon,
                                          max_lon)
 
            if len(cell_slopes[~np.isnan(cell_slopes)]) < 50:   
                
                zprime_slopes[row,col],zprime_slope_stds[row,col]  = [np.nan for i in range(2)]
            else:
                zprime_slopes[row,col],zprime_slope_stds[row,col] = np.nanmean(cell_slopes),np.nanstd(cell_slopes)

                
    return zprime_slopes,zprime_slope_stds