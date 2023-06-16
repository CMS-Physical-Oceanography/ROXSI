import numpy as np

def get_data(bath,
             bath_lat,
             bath_lon,
             min_lat,
             max_lat,
             min_lon,
             max_lon):
    
    """
    This funtion retrieves bathymetry data within the bounds of min_lat, max_lat, 
    and min_lon, max_lon.
    
    Used for grid cells with verticies running parallel with lines of latitude and
    longitude.
    
    INPUTS:
        -bath     - bathymetry soundings.
        -bath_lat,
         bath_lon -  Coordinates of soundings.
        -min_lat  - minimum latitude.
        -max_lat  - maximum latitude. 
        -min_lon,max_lon - see above.
    OUTPUTS:
        Bathymetry soundings and thier coordinates within the specified grid cell."""
    
    good_idxs = np.where(  (bath_lat > min_lat) & 
                           (bath_lat < max_lat) & 
                           (bath_lon > min_lon) & 
                           (bath_lon < max_lon) )
    

    return bath[good_idxs],bath_lat[good_idxs],bath_lon[good_idxs]


def get_cell_data(grid_lat,grid_lon,row,col,bath,bath_X):
    
    """
    Assumes grid[row,col] is the top corner of the cell.
    
    This function retrievs data from a grid cell that has 
    verticies which do not run parallel to lines of latitude 
    and longitude. This is done by taking the 2d cross product 
    of the vectors between the verticies of the grid cell specified
    by row, col and the vectors from each bathymetry point to the 
    vertex.
    
    Used for the 20m diamond shaped model grid.
    
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