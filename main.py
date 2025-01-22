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
new_b_lon, new_b_lat, new_b_grid = functions.compress_grid(b_grid, 6)
#new_a_lon, new_a_lat, new_a_grid = functions.compress_grid(a_grid, 1)
d_grid = functions.age_to_depth(a_grid)

d_grid_trim = d_grid[1:, :-1]   # slice d_grid to same shape as b_grid
depth_anomaly = new_b_grid - d_grid_trim  # array with the depth difference between bathymetry grid and calculated depth

# create subplots:
fig, axes = plt.subplots(2, 2, figsize=(16, 9))

mesh1 = axes[0][0].pcolormesh(a_lon, a_lat, a_grid)
axes[0][0].set_title("Age")
fig.colorbar(mesh1, ax=axes[0][0], label="Age")
axes[0][0]. grid()

mesh2 = axes[1][0].pcolormesh(a_lon, a_lat, d_grid)
axes[1][0].set_title("Calculated depth")
fig.colorbar(mesh2, ax=axes[1][0], label="Depth")
axes[1][0]. grid()

mesh3 = axes[0][1].pcolormesh(new_b_lon, new_b_lat, new_b_grid)
axes[0][1].set_title("Bathymetric map")
fig.colorbar(mesh3, ax=axes[0][1], label="Depth")
axes[0][1]. grid()

mesh4 = axes[1][1].pcolormesh(new_b_lon, new_b_lat, depth_anomaly)
axes[1][1].set_title("Depth anomaly")
fig.colorbar(mesh4, ax=axes[1][1], label="Δ Depth")
axes[1][1]. grid()

plt.show()
