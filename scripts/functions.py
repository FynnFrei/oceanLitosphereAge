import numpy as np
from scipy.io import netcdf_file
import os


# the theory of cooling oceanic litosphere in geo3 script1 P. 124

# TODO plot actual depth over expected depth for all grid points


# the not ideal way to downsample data. Lot of processing.
def compress_grid(grid, block_size):
    compressed_grid = np.zeros((grid.shape[0] // block_size, grid.shape[1] // block_size))
    compressed_lat = np.linspace(-90, 90, grid.shape[0] // block_size)
    compressed_lon = np.linspace(-180, 180, grid.shape[1] // block_size)
    print(f'Shape of compressed grid: {np.shape(compressed_grid)}')
    for i in range(0, grid.shape[0] - 1, block_size):  # calculate blocks
        for j in range(0, grid.shape[1] - 1, block_size):
            # fill compressed_grid with values
            compressed_grid[i // block_size, j // block_size] = grid[i:i + block_size, j:j + block_size].mean()

    return compressed_lon, compressed_lat, compressed_grid


# Does not work properly right now
def downsample_grid(grid, block_size):
    y, x = grid.shape
    new_y = (y // block_size) * block_size
    new_x = (x // block_size) * block_size
    trimmed_grid = grid[:new_y, :new_x]  # trimmed grid to match the divisible size
    downsampled_lat = np.linspace(-90, 90, new_y // block_size)
    downsampled_lon = np.linspace(-180, 180, new_x // block_size)

    # Reshape into 4 dimensions-> each new pixel contains a 2D block -> 4D
    reshaped = trimmed_grid.reshape(new_x // block_size, block_size, new_y // block_size, block_size)
    downsampled_grid = reshaped.mean(axis=(1, 3))  # Average across block rows and block columns
    print(np.shape(downsampled_grid))

    return downsampled_lat, downsampled_lon, downsampled_grid


def write_to_netcdf(output_filename, lon, lat, grid):
    """
    Write lon, lat, grid into a new .nc file.

    Parameters:
    - output_filename: str, name of the output NetCDF file
    - new_b_lon: np.ndarray, longitudes
    - new_b_lat: np.ndarray, latitudes
    - new_b_grid: np.ndarray, data grid
    """
    # Ensure the 'grids' subfolder exists
    output_folder = "grids"
    os.makedirs(output_folder, exist_ok=True)

    # Combine folder and filename
    full_output_path = os.path.join(output_folder, output_filename)

    # Open the NetCDF file for writing
    with netcdf_file(full_output_path, 'w') as nc_file:
        # Define dimensions
        nc_file.createDimension('lon', lon.size)
        nc_file.createDimension('lat', lat.size)

        # Create variables
        lon_var = nc_file.createVariable('lon', 'f8', ('lon',))
        lat_var = nc_file.createVariable('lat', 'f8', ('lat',))
        grid_var = nc_file.createVariable('elevation', 'f8', ('lat', 'lon'))

        # Write data to variables
        lon_var[:] = lon
        lat_var[:] = lat
        grid_var[:, :] = grid

        # Add metadata (optional)
        lon_var.units = 'degrees_east'
        lat_var.units = 'degrees_north'
        grid_var.units = 'elevation'

    print(f"NetCDF file '{full_output_path}' created successfully.")


def age_to_depth(grid):
    '''m = -0.3
    depth_grid = m * np.sqrt(1000000 * grid)'''

    alpha = 3 * 10 ** -5  # coefficient of thermal expansion in 1/K
    T_m = 1570  # mantle temperature in K (e.g. Vl.1, P.115)
    rho_m = 3300  # initial density of the hot mantle at T_m in Kg/m^3 (example value)
    rho_w = 1000  # density of water in kg/m^3
    kappa = 10 ** -6  # thermal diffusivity in m^2/s (data from Vl.1, P. 115)
    t = 3600 * 24 * 365 * 1000000 * grid  # time in s

    d = -(2 * alpha * rho_m * T_m / (rho_m - rho_w)) * np.sqrt((kappa / np.pi) * t)
    return d


def exclude_continental_plate(grid1, grid2):
    # replaces all grid points where continental plate exists by nan -> oceanic plate remains
    grid2_ocean = np.where(np.isnan(grid1), np.nan, grid2)
    return grid2_ocean


def shift_longitude(lon, grid, new_center):
    # shift center to wherever you want
    '''new_lon = np.where(lon < 0, lon + 360, lon)
    sorted_idx = np.argsort(new_lon)
    new_lon = new_lon[sorted_idx]
    new_grid = grid[:, sorted_idx]'''
    cut_deg = new_center
    if 0 <= new_center <= 180:
        cut_deg -= 180
    elif -180 <= new_center < 0:
        cut_deg += 180
    else:
        raise ValueError("ERROR: Choose value between -180 and 180 for new_center!")
    # TODO fix the index problem
    right_lon = lon[(lon >= cut_deg) | (lon < -180+cut_deg)]
    left_lon = lon[(lon < cut_deg) | (lon > 180+cut_deg)]
    new_lon = np.concatenate((right_lon, left_lon))

    right_idx = np.where((lon >= cut_deg) | (lon < -180+cut_deg))[0]
    left_idx = np.where((lon < cut_deg) | (lon > 180+cut_deg))[0]
    new_idx = np.concatenate((right_idx, left_idx))

    new_grid = grid[:, new_idx]
    '''
    # Calculate cell edges for pcolormesh
    new_lon_edges = np.concatenate(([new_lon[0] - 0.5 * (new_lon[1] - new_lon[0])],
                                    0.5 * (new_lon[:-1] + new_lon[1:]),
                                    [new_lon[-1] + 0.5 * (new_lon[-1] - new_lon[-2])]))'''
    return new_lon, new_grid


def trim_grid(lat, lon, grid, trim_lat, trim_lon):
    # TODO automatically trim to same size
    # cut off the last rows and columns for equal shape
    new_lat = lon[trim_lat:]
    new_lon = lat[:-trim_lon]
    new_grid = grid[trim_lat:, :-trim_lon]  # slice a_grid to same shape as b_grid
    print(f"New shape: {np.shape(new_grid)}")
    return new_lat, new_lon, new_grid
