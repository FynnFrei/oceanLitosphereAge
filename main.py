from scipy.io import netcdf_file
import numpy as np
from matplotlib import pyplot as plt

path = "grids/"
age_file = netcdf_file(path + 'age.2020.1.GTS2012.6m.grd', 'r', mmap=False)
elevation_file = netcdf_file(path + 'GRIDONE_2D_2008_1m.nc', mmap=False)


def compress_grid(grid, block_size):
    compressed_grid = np.zeros((grid.shape[0]//block_size, grid.shape[1]//block_size))
    compressed_lat = np.linspace(-90, 90, grid.shape[0]//block_size)
    compressed_lon = np.linspace(-180, 180, grid.shape[1]//block_size)
    print(np.shape(compressed_grid))
    for i in range(0, grid.shape[0]-1, block_size):     # calculate blocks
        for j in range(0, grid.shape[1]-1, block_size):
            compressed_grid[i//block_size, j//block_size] = (grid[i, j] + grid[i, j + 1] + grid[i + 1, j] + grid[i + 1, j + 1]) / 4

    return compressed_lon, compressed_lat, compressed_grid


# dataset a -> age_file
a_lon_var = age_file.variables['lon']
a_lat_var = age_file.variables['lat']
a_grid_var = age_file.variables['z']
a_lon = a_lon_var[:]
a_lat = a_lat_var[:]
a_grid = a_grid_var[:]      # age in Ma

# dataset e -> elevation_file
e_lon_var = elevation_file.variables['lon']
e_lat_var = elevation_file.variables['lat']
e_grid_var = elevation_file.variables['elevation']
e_lon = e_lon_var[:]
e_lat = e_lat_var[:]
e_grid = e_grid_var[:]      # depth in m

#print(age_file.variables.keys())
#print(elevation_file.variables.keys())
age_file.close()
elevation_file.close()

new_e_lon, new_e_lat, new_e_grid = compress_grid(e_grid, 12)
new_a_lon, new_a_lat, new_a_grid = compress_grid(a_grid, 1)

#print(np.shape(new_e_grid), np.shape(new_a_grid))

# create subplots:
fig, axes = plt.subplots(2, 1, figsize=(12, 12))

axes[0].pcolormesh(new_e_lon, new_e_lat, new_e_grid)
axes[0]. grid()

axes[1].pcolormesh(new_a_lon, new_a_lat, new_a_grid)
axes[1]. grid()

plt.show()

# this is test commit.