import numpy as np
from scipy.io import netcdf_file
import os


# the theory of cooling oceanic litosphere in geo3 script1 P. 124


# the not ideal way to downsample data. Lot of processing.
def compress_grid(grid, block_size):
    compressed_grid = np.zeros((grid.shape[0] // block_size, grid.shape[1] // block_size))
    compressed_lat = np.linspace(-90, 90, grid.shape[0] // block_size)
    compressed_lon = np.linspace(-180, 180, grid.shape[1] // block_size)
    print(f'Target shape: {np.shape(compressed_grid)}')
    for i in range(0, grid.shape[0] - 1, block_size):  # calculate blocks
        # print(f'column{i} of {grid.shape[0]-1}', end='\r')
        for j in range(0, grid.shape[1] - 1, block_size):
            # fill compressed_grid with values
            compressed_grid[i // block_size, j // block_size] = grid[i:i + block_size, j:j + block_size].mean()
    print(f'grid downsampled successfully')

    return compressed_lat, compressed_lon, compressed_grid


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


def write_to_netcdf(output_filename, lat, lon, grid):
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
    print('Writing .nc file...')
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


def get_land_grid(b_grid):
    land_grid = np.where(b_grid >= 0, 1, np.nan)
    return land_grid


def calc_ridge_depth(a_grid, b_grid, max_age):
    # calculates the mean depth for a given max crust age
    a_grid_ridge = np.where(a_grid > max_age, np.nan, a_grid)   # grid only contains young lithosphere -> ridge
    ridge_depth = np.where(a_grid_ridge > 0, b_grid, np.nan)    # gives depht to ridge grid points
    avr_ridge_depth = np.nanmean(ridge_depth)   # calculate mean depth
    num_of_points = np.sum(a_grid_ridge > 0)
    print(avr_ridge_depth, num_of_points)

    return avr_ridge_depth


def age_to_depth(grid, d_0):

    alpha = 2 * 10 ** -5  # coefficient of thermal expansion in 1/K
    T_m = 1570  # mantle temperature in K (e.g. Vl.1, P.115)
    rho_m = 3300  # initial density of the hot mantle at T_m in Kg/m^3 (example value)
    rho_w = 1030  # density of water in kg/m^3
    kappa = 10 ** -6  # thermal diffusivity in m^2/s (data from Vl.1, P. 115)
    t = 3600 * 24 * 365 * 1000000 * grid  # time in s
    # d_0 is average ocean ridge depth
    # d is not absolute depth but depth change compared to ocean ridge

    d = -(2 * alpha * rho_m * T_m / (rho_m - rho_w)) * np.sqrt((kappa / np.pi) * t) + d_0
    return d


def exclude_continental_plate(grid1, grid2):
    # replaces all grid points where continental plate exists by nan -> oceanic plate remains
    grid2_ocean = np.where(np.isnan(grid1), np.nan, grid2)
    return grid2_ocean


def shift_longitude(lon, grid, new_center):
    shift = - (len(lon)/360) * new_center
    shifted_lon = np.roll(lon, shift)
    shifted_grid = np.roll(grid, shift, axis=1)
    return shifted_lon, shifted_grid


def trim_grid(lat, lon, grid, ref_grid):
    # cut off the last rows and columns for equal shape
    trim_lat = len(grid[0]) - len(ref_grid[0])
    trim_lon = len(grid[1]) - len(ref_grid[1])
    new_lat = lon[trim_lat:]
    new_lon = lat[:-trim_lon]
    new_grid = grid[trim_lat:, :-trim_lon]  # slice a_grid to same shape as b_grid
    print(f"New shape: {np.shape(new_grid)}")
    return new_lat, new_lon, new_grid


def compare_grid_points(grid1, grid2):
    grid1_flat = np.reshape(grid1, -1)
    grid2_flat = np.reshape(grid2, -1)
    mask_nan = ~np.isnan(grid1_flat) & ~np.isnan(grid2_flat)
    grid1_flat_clean = grid1_flat[mask_nan]
    grid2_flat_clean = grid2_flat[mask_nan]
    diff = grid2_flat_clean - grid1_flat_clean
    return diff, grid1_flat_clean, grid2_flat_clean

