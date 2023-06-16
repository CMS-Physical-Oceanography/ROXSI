import numpy as np


def make_finer_grid(grid_lat,grid_lon,fact):
    """
    This function inputs a grid and returns a
    different grid over the same area that is fact
    times finer than the input.
    
    For example if the inputted grid has a 100m 
    resolution and fact=2 the output will be 50m'
    resolution.
    Inputs:
        -grid_lat (2d array) : lat of input grid
        -grid_lon (2d array) : lon of input grid, must be the
                               the same shape as grid_lat
        -fact (integer ) : factor by which the resolution is lowered
    Outputs:
        -new_lat (2d array) : lower resolution grid lat.
        -new_lon (2d array) : lower resolution grid lon."""
        
    new_lat = np.full((grid_lat.shape[0]*fact,grid_lat.shape[1]*fact),np.nan)
    new_lon = np.full(new_lat.shape,np.nan)
    
    row_idxs = np.linspace(0,new_lat.shape[0]-1,new_lat.shape[0],dtype=int)
    col_idxs = np.linspace(0,new_lat.shape[1]-1,new_lat.shape[1],dtype=int)
    
    new_lat[row_idxs%fact==0,0] = grid_lat[:,0]    
    new_lon[0,col_idxs%fact==0] = grid_lon[0,:]
    
    lat_non_nans = row_idxs[~np.isnan(new_lat[:,0])]
    lon_non_nans = col_idxs[~np.isnan(new_lon[0,:])]

    lat_vals = np.interp(row_idxs,lat_non_nans,new_lat[lat_non_nans,0])
    lon_vals = np.interp(col_idxs,lon_non_nans,new_lon[0,lon_non_nans])


    for col in range(new_lat.shape[1]):
        new_lat[:,col] = lat_vals
    for row in range(new_lat.shape[0]):
        new_lon[row,:] = lon_vals
    
    return new_lat,new_lon


def make_coarser_grid(grid_lat,grid_lon,fact):
    """
    This function inputs a grid and returns a different grid
    over the save area that is fact times coarser than the
    input.
    
    For example if the inputted grid has 100m resolution and 
    fact =3 the output will be 200m resolution.
    Inputs:
        -grid_lat (2d array) : lat of input grid
        -grid_lon (2d array) : lon of input grid, must be the
                               the same shape as grid_lat.
        -fact (integer ) : factor by which the resolution is lowered
    Outputs:
        -new_lat (2d array) : lower resolution grid lat.
        -new_lon (2d array) : lower resolution grid lon."""
    
    idxs = np.linspace(0,grid_lat.shape[0]-1,grid_lat.shape[0],dtype=int)
    row_idxs = idxs[idxs%fact==0]
    idxs = np.linspace(0,grid_lat.shape[1]-1,grid_lat.shape[1],dtype=int)
    col_idxs = idxs[idxs%fact==0]
    
    new_lat_ = grid_lat[:,col_idxs]
    new_lon_ = grid_lon[:,col_idxs]

    new_lat = new_lat_[row_idxs,:]
    new_lon = new_lon_[row_idxs,:]
    
    return new_lat,new_lon