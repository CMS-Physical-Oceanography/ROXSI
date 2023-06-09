import numpy as np
from cell_roughness import get_cell_roughness

"""
This script contains funcitons used to calculate
sigma z' for all model grids.

for the 100m and 500m grids:
    get_data()
    grid_roughness()
    grid_slopes()
    
for the 20m grid:
    get_cell_data()
    grid20_roughness()
    grid20_slopes()
    
for all grids:
    get_slopes()
"""

def get_data(bath,
             bath_X,
             min_y,
             max_y,
             min_x,
             max_x):
    
    """
    This funtion retrieves bathymetry data within the bounds of min_lat, max_lat, 
    and min_lon, max_lon.
    
    Used for the 100m and 500m grid.
    
    INPUTS:
        -bath     - bathymetry soundings.
        -bath_X  -  Coordinates of soundings where bath_X[:,0]= x and bath_X[:,1]=y
        -min_lat  - minimum latitude.
        -max_lat  - maximum latitude. 
        -min_lon,max_lon - see above.
    OUTPUTS:
        Bathymetry soudnings and thier coordinates within the specified grid cell."""
    
    good_idxs = np.where(  (bath_X[:,1] > min_y) & 
                           (bath_X[:,1] < max_y) & 
                           (bath_X[:,0] > min_x) & 
                           (bath_X[:,0] < max_x) )
    

    return bath[good_idxs],bath_X[good_idxs]




def get_cell_data(grid_lat,grid_lon,row,col,bath,bath_X):
    
    """
    Assumes grid[row,col] is the top corner of the cell.
    
    This function retrievs data from a grid cell that has 
    verticies which do not run parallel to lines of latitude 
    and longitude. This is done by taking the 2d cross product 
    of the vectors between the verticies of the grid cell specified
    by row, col and the vectors from each bathymetry point to the 
    vertex.
    
    Used for the 20m grid.
    
    Inputs:
        grid_lat - (2d array) latitide points of the grid.
        grid_lon - longitude points of the grid
        row      - row index of the cell where data is needed.
        col      - column index of the cell.
        bath     - All hi resolution bathymetry
        bath_X   - Coordinates of bath where bath_X[:,0] = x and
                   bath_X[:,1] = y.
    Outputs:
        Bathymetry soundings withing the grid cell and thier 
        coordinates."""

    
    def cross(A,B):
        return A[0]*B[:,1] - A[1]*B[:,0]
    
    tc = np.array([grid_lon[row,col],grid_lat[row,col]]) # top corner
    lc = np.array([grid_lon[row+1,col],grid_lat[row+1,col]]) # left corner
    bc = np.array([grid_lon[row+1,col+1],grid_lat[row+1,col+1]]) # bottom corner
    rc = np.array([grid_lon[row,col+1],grid_lat[row,col+1]]) # right corner
    
    # vectors between verticies
    tclc = tc-lc
    lcbc = lc-bc
    bcrc = bc-rc
    rctc = rc-tc
    
    tcpv,lcpv,bcpv,rcpv = [np.zeros((len(bath),2)) for i in range(4)]
    
    
    for i in range(2): # get vectors between verticies and bath points
        tcpv[:,i] = bath_X[:,i]-tc[i]
        lcpv[:,i] = bath_X[:,i]-lc[i]
        bcpv[:,i] = bath_X[:,i]-bc[i]
        rcpv[:,i] = bath_X[:,i]-rc[i]

    # take cross product
    tclc_check = cross(tclc,tcpv)
    lcbc_check = cross(lcbc,lcpv)
    bcrc_check = cross(bcrc,bcpv)
    rctc_check = cross(rctc,rcpv)
    
    # get indicies of bathymetry data within the grid cell
    in_bounds = np.where((tclc_check < 0) & 
                         (lcbc_check < 0) & 
                         (bcrc_check < 0) & 
                         (rctc_check < 0))
    
    return bath[in_bounds],bath_X[in_bounds,:][0]



def grid_roughness(grid_lat,
                   grid_lon,
                   bath,
                   bath_X):
    """
    This function calculates roughness for the 100m and 500m grids.
    
    Inputs:
        grid_lat (2d array) - lat coordinates of grid.
        grid_lon (2d array) - lon coordinates of grid.
        bath     (1d array) - bathymetry soundings.
        bath_X   (shape = (len(bath),2) - coordinates of soundings where 4
                 bath_X[:,0] = x and bath_X[:,1] = y.
    Outputs:
        roughness_grid (2d array) - gridded cell roughness as std(z')"""
    
    roughness_grid = np.zeros((grid_lat.shape[0]-1,grid_lat.shape[1]-1))

    for row in range(roughness_grid.shape[0]):
        print(row)
        for col in range(roughness_grid.shape[1]):
            
            max_lat = grid_lat[row,col]
            min_lat = grid_lat[row+1,col]
       
            max_lon = grid_lon[row,col+1]
            min_lon = grid_lon[row,col]
   

            cell_h,cell_X = get_data(bath,
                                            bath_X,
                                            min_lat,
                                            max_lat,
                                            min_lon,
                                            max_lon)
 
            if len(cell_h[~np.isnan(cell_h)]) < 50:   
                
                roughness_grid[row,col] = np.nan
            else:
                roughness_grid[row,col] = get_cell_roughness(cell_h,cell_lats,cell_lons)

            
    return roughness_grid


def grid20_roughness(grid_lat,grid_lon,bath,bath_X):
    
    """
    This function calculates roughness for the 20m grid as 
    std(z')
    
    Inputs and outputs are the same as grid_roughness() above."""
    
    
    roughness = np.zeros((grid_lat.shape[0]-1,grid_lat.shape[1]-1))


    for row in range(roughness.shape[0]):
        print(row/roughness.shape[0])
        for col in range(roughness.shape[1]):
            
            cell_h,cell_X = get_cell_data(grid_lat,grid_lon,row,col,bath,bath_X)
            
            if len(cell_h) <5:      
                roughness[row,col] = np.nan          
            else:
                roughness[row,col],rgh_square[row,col],mean_h[row,col],mean_b[row,col] = get_cell_roughness(cell_h,cell_X)

            
    return roughness


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
    
    Used for the 100m and 50m grid.
    
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
        print(row)
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



def grid20_slopes(grid_lat,
                  grid_lon,
                  slopes,
                  slopes_X):
    
    """
    Same as grid_slopes but for the 20m grid."""
    
    zprime_slopes = np.zeros((grid_lat.shape[0]-1,grid_lat.shape[1]-1))
    zprime_slope_stds = np.zeros_like(zprime_slopes)
    
    for row in range(zprime_slopes.shape[0]):
        print(row/zprime_slopes.shape[0])
        for col in range(zprime_slopes.shape[1]):
               

            cell_slopes,slopes_X = get_cell_data(grid_lat,
                                                 grid_lon,
                                                 row,
                                                 col,
                                                 slopes,
                                                 slopes_X)
 
            if len(cell_slopes[~np.isnan(cell_slopes)]) < 5:   
                
                zprime_slopes[row,col],zprime_slope_stds[row,col]  = [np.nan for i in range(2)]
            else:
                zprime_slopes[row,col],zprime_slope_stds[row,col] = np.nanmean(cell_slopes),np.nanstd(cell_slopes)

                if np.nanstd(cell_slopes) == 0:
                    print(len(cell_slopes[~np.isnan(cell_slopes)]),row,col,cell_slopes[~np.isnan(cell_slopes)])
            
    return zprime_slopes,zprime_slope_stds