from scipy.io import netcdf_file
import numpy as np
from matplotlib import pyplot as plt
from scripts import functions
#from scripts.functions import write_to_netcdf
from matplotlib.colors import Normalize, ListedColormap

path = "grids/"
age_file = netcdf_file(path + 'age.2020.1.GTS2012.6m.grd', 'r', mmap=False)
bathymetry_file = netcdf_file(path + 'GRIDONE_2D_2008_6m.nc', mmap=False)
#print(age_file.variables.keys())
#print(bathymetry_file.variables.keys())

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


# adjust grids
b_lon_shift, b_grid_shift = functions.shift_longitude(b_lon, b_grid, -180)  # shift lon to new center
a_lon_shift, a_grid_shift = functions.shift_longitude(a_lon, a_grid, -180)
a_lat_trim, a_lon_trim_shift, a_grid_trim_shift = functions.trim_grid(a_lat, a_lon_shift, a_grid_shift, 1, 1)   # trim to same shape
b_grid = functions.exclude_continental_plate(a_grid_trim_shift, b_grid_shift)     # make continents NaN

# calculate difference:
d_grid = functions.age_to_depth(a_grid_trim_shift)
depth_anomaly = b_grid_shift - d_grid  # array with the depth difference between bathymetry grid and calculated depth

# create subplots:
fig, ax = plt.subplots(2, 2, figsize=(20, 9))

cmap1 = plt.cm.viridis
cmap2 = plt.cm.turbo_r
cmap3 = plt.cm.seismic
continent_color = (0.3, 0.3, 0.3)
cmap1.set_bad(color=continent_color)
cmap2.set_bad(color=continent_color)
cmap3.set_bad(color=continent_color)

mesh1 = ax[0][0].pcolormesh(a_lon_shift, a_lat, a_grid_shift, cmap=cmap2, vmax=200)
ax[0][0].set_title("Age in Ma")
ax[0][0].set_aspect(1)
fig.colorbar(mesh1, ax=ax[0][0], label="Age")
ax[0][0].grid()

mesh2 = ax[0][1].pcolormesh(b_lon_shift, b_lat, b_grid, cmap=cmap1, vmin=-7000, vmax=0)
ax[0][1].set_title("Bathymetric map")
ax[0][1].set_aspect(1)
fig.colorbar(mesh2, ax=ax[0][1], label="Depth")
ax[0][1].grid()

mesh3 = ax[1][0].pcolormesh(a_lon_shift, a_lat, d_grid, cmap=cmap1, vmin=-7000, vmax=0)
ax[1][0].set_title("Calculated depth")
ax[1][0].set_aspect(1)
fig.colorbar(mesh3, ax=ax[1][0], label="Depth")
ax[1][0].grid()

mesh4 = ax[1][1].pcolormesh(b_lon_shift, b_lat, depth_anomaly, cmap=cmap3, vmin=-3000, vmax=3000)
ax[1][1].set_title("Depth anomaly")
ax[1][1].set_aspect(1)
fig.colorbar(mesh4, ax=ax[1][1], label="Δ Depth")
ax[1][1].grid()

plt.show()


'''fig2 = plt.figure(figsize=(16, 9))
plt.pcolormesh(b_lon, b_lat, b_grid)
plt.title("Bathymetric map")
plt.colorbar(label="Depth")
plt.grid()
plt.show()'''
