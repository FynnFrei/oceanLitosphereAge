from scipy.io import netcdf_file
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize, ListedColormap
from scripts import functions

path = "grids/"
age_file = netcdf_file(path + 'age.2020.1.GTS2012.2m.grd', 'r', mmap=False)
bathymetry_file = netcdf_file(path + 'GRIDONE_2D_2008_2m.nc', mmap=False)
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

# new_lat, new_lon, new_grid = functions.compress_grid(b_grid, 2)
#write_to_netcdf("GRIDONE_2D_2008_2m.nc", new_lat, new_lon, new_grid)

b_lon_shift, b_grid_shift = functions.shift_longitude(b_lon, b_grid, 200)  # shift lon to new center
a_lon_shift, a_grid_shift = functions.shift_longitude(a_lon, a_grid, 200)  # shift lon to new center

# adjust grids
a_lat_trim, a_lon_trim_shift, a_grid_trim_shift = functions.trim_grid(a_lat, a_lon_shift, a_grid_shift, b_grid)   # trim to same shape
b_grid_shift = functions.exclude_continental_plate(a_grid_trim_shift, b_grid_shift)     # make continents NaN

# calculate difference:
d_grid = functions.age_to_depth(a_grid_trim_shift)
depth_anomaly = d_grid - b_grid_shift  # array shows how much calculated depth deviates from actual depth

# create subplots:

cmap1 = plt.cm.viridis
cmap2 = plt.cm.turbo_r
cmap3 = plt.cm.seismic
continent_color = (0.5, 0.5, 0.5)
cmap1.set_bad(color=continent_color)
cmap2.set_bad(color=continent_color)
cmap3.set_bad(color=continent_color)

'''
fig1, ax = plt.subplots(2, 2, figsize=(18, 8))

#mesh1 = ax[0][0].pcolormesh(a_lon_shift, a_lat, a_grid_shift, cmap=cmap2, vmax=200)
mesh1 = ax[0][0].imshow(
    a_grid_shift,
    extent=[a_lon_shift.min(), a_lon_shift.max(), a_lat.min(), a_lat.max()],
    origin="lower",
    cmap=cmap2,
    vmax=200,
    aspect=1
)
ax[0][0].set_title("Age in Ma")
ax[0][0].set_aspect(1)
fig1.colorbar(mesh1, ax=ax[0][0], label="Age")
ax[0][0].grid()

mesh2 = ax[0][1].imshow(
    b_grid_shift,
    extent=[b_lon_shift.min(), b_lon_shift.max(), b_lat.min(), b_lat.max()],
    origin="lower",
    cmap=cmap1, vmin=-7000, vmax=0, aspect=1
)
ax[0][1].set_title("Bathymetric map")
ax[0][1].set_aspect(1)
fig1.colorbar(mesh2, ax=ax[0][1], label="Depth")
ax[0][1].grid()

mesh3 = ax[1][0].imshow(
    d_grid,
    extent=[b_lon_shift.min(), b_lon_shift.max(), b_lat.min(), b_lat.max()],
    origin="lower",
    cmap=cmap1, vmin=-7000, vmax=0, aspect=1
)
ax[1][0].set_title("Calculated depth")
ax[1][0].set_aspect(1)
fig1.colorbar(mesh3, ax=ax[1][0], label="Depth")
ax[1][0].grid()

mesh4 = ax[1][1].imshow(
    depth_anomaly,
    extent=[b_lon_shift.min(), b_lon_shift.max(), b_lat.min(), b_lat.max()],
    origin="lower",
    cmap=cmap3, vmin=-3000, vmax=3000, aspect=1
)
ax[1][1].set_title("Depth anomaly")
ax[1][1].set_aspect(1)
fig1.colorbar(mesh4, ax=ax[1][1], label="Δ Depth")
ax[1][1].grid()

plt.show()
'''


fig2 = plt.figure(figsize=(16, 8))
plt.imshow(
    depth_anomaly,
    extent=[b_lon_shift.min(), b_lon_shift.max(), b_lat.min(), b_lat.max()],
    origin="lower",
    cmap=cmap3, vmin=-3000, vmax=3000, aspect=1
)
plt.title("Depth Anomaly (comparing calculated with real depth)")
plt.colorbar(label="Δ Depth")
plt.grid()
plt.show()

# functions.plot_scatter(d_grid, b_grid_shift)
'''
grid1_flat = np.reshape(d_grid, -1)
grid2_flat = np.reshape(b_grid_shift, -1)
mask_nan = ~np.isnan(grid1_flat) & ~np.isnan(grid2_flat)
grid1_flat_clean = grid1_flat[mask_nan]
grid2_flat_clean = grid2_flat[mask_nan]
ratio = grid2_flat_clean / grid1_flat_clean
diff = grid2_flat_clean - grid1_flat_clean

fig4 = plt.figure()
plt.axvline(x=1, alpha=0.4, color='red', linestyle='--', linewidth=2)
plt.hist(ratio, bins=120, alpha=1, range=(-1, 5))
plt.title("Distribution of real/calc depth")
plt.show()

fig5 = plt.figure()
plt.axvline(x=0, alpha=0.4, color='red', linestyle='--', linewidth=2)
plt.hist(diff, bins=100, alpha=1)
plt.title("Distribution of Δ depth")
plt.text(-6000, 2*10**5, "Deeper than calculation")
plt.text(2000, 2*10**5, "Shallower than calculation")
plt.show()
'''
