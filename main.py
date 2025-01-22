from scipy.io import netcdf_file
import numpy as np
from matplotlib import pyplot as plt
from scripts import functions

path = "grids/"
age_file = netcdf_file(path + 'age.2020.1.GTS2012.6m.grd', 'r', mmap=False)
bathymetry_file = netcdf_file(path + 'GRIDONE_2D_2008_1m.nc', mmap=False)
# print(age_file.variables.keys())
# print(bathymetry_file.variables.keys())

# dataset a -> age_file
a_lon_var = age_file.variables['lon']
a_lat_var = age_file.variables['lat']
a_grid_var = age_file.variables['z']
a_lon = a_lon_var[:]        # reads longitude values into list
a_lat = a_lat_var[:]
a_grid = a_grid_var[:]      # age in Ma

# dataset b -> bathymetry_file
b_lon_var = bathymetry_file.variables['lon']
b_lat_var = bathymetry_file.variables['lat']
b_grid_var = bathymetry_file.variables['elevation']
b_lon = b_lon_var[:]
b_lat = b_lat_var[:]
b_grid = b_grid_var[:]      # depth in m

age_file.close()
bathymetry_file.close()

# downsample grids:
new_b_lon, new_b_lat, new_b_grid = functions.compress_grid(b_grid, 20)
new_a_lon, new_a_lat, new_a_grid = functions.compress_grid(a_grid, 1)

# create subplots:
fig, axes = plt.subplots(2, 1, figsize=(12, 12))

axes[0].pcolormesh(new_b_lon, new_b_lat, new_b_grid)
axes[0]. grid()

axes[1].pcolormesh(new_a_lon, new_a_lat, new_a_grid)
axes[1]. grid()

plt.show()
