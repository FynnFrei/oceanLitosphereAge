from scipy.io import netcdf_file
from scripts import functions
from scripts.functions import depth_to_age
import numpy as np

path = "grids/"
accuracy = 2  # how many minutes° are one grid point (2 or 6)
map_center = 200    # central longitude of plotted map
max_ridge_age = 0.01    # in Ma, for calculating avr ridge depth
age_file = netcdf_file(path + f'age.2020.1.GTS2012.{accuracy}m.grd', 'r', mmap=False)
bathymetry_file = netcdf_file(path + f'GRIDONE_2D_2008_{accuracy}m.nc', mmap=False)
# print(age_file.variables.keys())

# dataset a -> age_file
a_lon_var = age_file.variables['lon']
a_lat_var = age_file.variables['lat']
a_grid_var = age_file.variables['z']
a_lon = a_lon_var[:]  # reads longitude values into list
a_lat = a_lat_var[:]
a_grid = a_grid_var[:]  # age in Ma
print(f'grid point: {a_grid[1000][3600]}')
a = np.array([
    [10, 20, 30, 1],
    [40, 50, 60, 2],
    [70, 80, 90, 3]
])
print(len(a[:]))
y_index = - np.repeat(np.arange(len(a[:])), len(a[0])) + 1
a_flat = np.reshape(a, -1)
a_y = np.stack((a_flat, y_index), axis=1)
print(a_y)
print(np.cos(90/180*np.pi))

# dataset b -> bathymetry_file
b_lon_var = bathymetry_file.variables['lon']
b_lat_var = bathymetry_file.variables['lat']
b_grid_var = bathymetry_file.variables['elevation']
b_lon = b_lon_var[:]
b_lat = b_lat_var[:]
b_grid = b_grid_var[:]  # depth in m

age_file.close()
bathymetry_file.close()

# new_lat, new_lon, new_grid = functions.compress_grid(b_grid, 2)
#write_to_netcdf("GRIDONE_2D_2008_2m.nc", new_lat, new_lon, new_grid)

b_lon_shift, b_grid_shift = functions.shift_longitude(b_lon, b_grid, map_center)  # shift lon to new center
a_lon_shift, a_grid_shift = functions.shift_longitude(a_lon, a_grid, map_center)  # shift lon to new center

# adjust grids
a_lat_trim, a_lon_trim_shift, a_grid_trim_shift = functions.trim_grid_to_ref(a_lat, a_lon_shift, a_grid_shift,
                                                                      b_grid)  # trim to same shape
b_grid_shift_ocean = functions.exclude_continental_plate(a_grid_trim_shift, b_grid_shift)  # make continents NaN

# calculate difference:
avr_ridge_depth = functions.calc_ridge_depth(a_grid_trim_shift, b_grid_shift_ocean, max_ridge_age)
d_grid = functions.age_to_depth(a_grid_trim_shift, avr_ridge_depth)     # dataset d -> calculated depth
depth_anomaly = abs(d_grid) - abs(b_grid_shift_ocean)  # array shows how much calculated depth deviates from actual depth
depth_anomaly_abs = abs(depth_anomaly)
land_grid = functions.get_land_grid(b_grid_shift)


'''print(f'{depth_to_age(-8000, avr_ridge_depth)} Ma old')   # to get age from depth (doublecheck)'''




# Glossar---------------------------------------------------------------------------------------------------------------
# b_            bathymetric data
# a_            age data
# d_            calculated depth data
# shift         new longitude center (shifted x-axis)
# trim          cut edges of grid to same shape of all grids
# ocean         only ocean lithosphere has values - continents are NaN
# clean         without NaN
# window        only a window of the whole world map (e.g. only northern hemisphere, western part)