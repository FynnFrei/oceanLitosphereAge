from operator import is_not

from scipy.io import netcdf_file
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
from scripts import functions
from scripts.functions import get_land_grid

path = "grids/"
accuracy = 6  # how many minutes° are one gridpoint
age_file = netcdf_file(path + f'age.2020.1.GTS2012.{accuracy}m.grd', 'r', mmap=False)
bathymetry_file = netcdf_file(path + f'GRIDONE_2D_2008_{accuracy}m.nc', mmap=False)
#coastline_file = netcdf_file(path + f'binned_GSHHS_f.nc', mmap=False)
#print(coastline_file.variables.keys())
#print(bathymetry_file.variables.keys())

# dataset a -> age_file
a_lon_var = age_file.variables['lon']
a_lat_var = age_file.variables['lat']
a_grid_var = age_file.variables['z']
a_lon = a_lon_var[:]  # reads longitude values into list
a_lat = a_lat_var[:]
a_grid = a_grid_var[:]  # age in Ma

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

b_lon_shift, b_grid_shift = functions.shift_longitude(b_lon, b_grid, 200)  # shift lon to new center
a_lon_shift, a_grid_shift = functions.shift_longitude(a_lon, a_grid, 200)  # shift lon to new center

# adjust grids
a_lat_trim, a_lon_trim_shift, a_grid_trim_shift = functions.trim_grid(a_lat, a_lon_shift, a_grid_shift,
                                                                      b_grid)  # trim to same shape
b_grid_shift_ocean = functions.exclude_continental_plate(a_grid_trim_shift, b_grid_shift)  # make continents NaN

# calculate difference:
avr_ridge_depth = functions.calc_ridge_depth(a_grid_trim_shift, b_grid_shift_ocean, 0.01)
d_grid = functions.age_to_depth(a_grid_trim_shift, avr_ridge_depth)
depth_anomaly = abs(d_grid) - abs(b_grid_shift_ocean)  # array shows how much calculated depth deviates from actual depth
depth_anomaly_abs = abs(depth_anomaly)

# create subplots:

cmap1 = plt.cm.viridis
cmap2 = plt.cm.turbo_r
cmap3 = plt.cm.seismic
cmap4 = plt.cm.Greens
continent_color = (0.5, 0.5, 0.5)
transparent_nan = (0, 0, 0, 0)
cmap1.set_bad(color=continent_color)
cmap2.set_bad(color=continent_color)
cmap3.set_bad(color=continent_color)
cmap4.set_bad(color=transparent_nan)
font_axes_scale = 1.8  # axes values
font_header_scale = 4  # factor that scales the header font to the standard font size
font_label_scale = 2.5  # labels, source(GEBCO)
font_text_scale = 2  # normal text, description

# all maps-------------------------------------------------------------------------------------------------------------
'''
fig1, ax = plt.subplots(2, 2, figsize=(22, 10))

mesh1 = ax[0][0].imshow(
    a_grid_shift,
    extent=[a_lon_shift.min(), a_lon_shift.max(), a_lat.min(), a_lat.max()],
    origin="lower",
    cmap=cmap2,
    vmax=200,
    aspect=1
)
ax[0][0].set_title("Age of oceanic lithosphere")
ax[0][0].text(-170, -85, "Ogg 2012")
ax[0][0].set_aspect(1)
fig1.colorbar(mesh1, ax=ax[0][0], label="Age in Ma")
ax[0][0].grid()

mesh2 = ax[0][1].imshow(
    b_grid_shift_ocean,
    extent=[b_lon_shift.min(), b_lon_shift.max(), b_lat.min(), b_lat.max()],
    origin="lower",
    cmap=cmap1, vmin=-7000, vmax=0, aspect=1
)
ax[0][1].set_title("Bathymetric map")
ax[0][1].text(-170, -85, "GEBCO_2020 Grid")
ax[0][1].set_aspect(1)
fig1.colorbar(mesh2, ax=ax[0][1], label="Depth in m")
ax[0][1].grid()

mesh3 = ax[1][0].imshow(
    d_grid,
    extent=[b_lon_shift.min(), b_lon_shift.max(), b_lat.min(), b_lat.max()],
    origin="lower",
    cmap=cmap1, vmin=-7000, vmax=-1000, aspect=1
)
ax[1][0].set_title("Calculated depth")
ax[1][0].text(-170, -85, "Ogg 2012")
ax[1][0].set_aspect(1)
fig1.colorbar(mesh3, ax=ax[1][0], label="Depth in m")
ax[1][0].grid()

mesh4 = ax[1][1].imshow(
    depth_anomaly,
    extent=[b_lon_shift.min(), b_lon_shift.max(), b_lat.min(), b_lat.max()],
    origin="lower",
    cmap=cmap3, vmin=-3000, vmax=3000, aspect=1
)
ax[1][1].set_title("Calculated depth - bathymetric depth")
ax[1][1].text(-170, -85, "GEBCO_2020 Grid, Ogg 2012")
ax[1][1].set_aspect(1)
fig1.colorbar(mesh4, ax=ax[1][1], label="Δ Depth in m")
ax[1][1].grid()

plt.show()
'''

# depth anomaly map-----------------------------------------------------------------------------------------------------

fig2, ax = plt.subplots(figsize=(18, 9))
ax.imshow(
    depth_anomaly,
    extent=[a_lon_shift.min(), a_lon_shift.max(), a_lat.min(), a_lat.max()],
    origin="lower",
    cmap=cmap3,
    vmin=-3000,
    vmax=3000,
    aspect=1
)
land_grid = get_land_grid(b_grid_shift)
masked_land_grid = np.ma.masked_where(np.isnan(land_grid), land_grid)
ax.imshow(
    land_grid,
    extent=[a_lon_shift.min(), a_lon_shift.max(), a_lat.min(), a_lat.max()],
    origin="lower",
    cmap=cmap4,
    vmin=0,
    vmax=1,
    alpha=0.3,
    aspect=1
)
relative_font_size = fig2.get_size_inches()[0]*0.5  # gets font size relative to figure size
plt.title("Calculated Depth - Bathymetric Depth", fontsize=relative_font_size*font_header_scale)
plt.text(-170, -85, "GEBCO_2020 Grid, Ogg 1012", fontsize=relative_font_size*font_label_scale)
'''ax.cbar = plt.colorbar()
ax.cbar.ax.tick_params(labelsize=relative_font_size*font_axes_scale)   # modify font size of cbar values
ax.cbar.set_label(label="Δ depth (in meters)", fontsize=relative_font_size*font_label_scale)'''
plt.tick_params(labelsize=relative_font_size*font_axes_scale)
plt.xlabel("Longitude", fontsize=relative_font_size*font_label_scale)
plt.ylabel("Latitude", fontsize=relative_font_size*font_label_scale)
plt.grid()
plt.show()


# functions.plot_scatter(d_grid, b_grid_shift_ocean)

# depth anomaly histogram-----------------------------------------------------------------------------------------------
"""
grid1_flat = np.reshape(d_grid, -1)
grid2_flat = np.reshape(b_grid_shift_ocean, -1)
mask_nan = ~np.isnan(grid1_flat) & ~np.isnan(grid2_flat)
grid1_flat_clean = grid1_flat[mask_nan]
grid2_flat_clean = grid2_flat[mask_nan]
# ratio = grid2_flat_clean / grid1_flat_clean
diff = grid2_flat_clean - grid1_flat_clean

fig5 = plt.figure()
plt.axvline(x=0, alpha=0.4, color='red', linestyle='--', linewidth=2)
plt.hist(diff, bins=100, alpha=1, range=(-3000, 7000))
plt.title("Distribution of Δ depth", fontsize=font_size_header)
plt.text(5800, 1.95*10**6, "GEBCO_2020 Grid, Ogg 1012")
plt.xlabel("anomaly in m", fontsize=font_size_label)
plt.ylabel("amount of data points", fontsize=font_size_label)
plt.text(-3000, 8*10**5, "Deeper than calculation", color="blue", fontsize=font_size_text)
plt.text(2000, 8*10**5, "Shallower than calculation", color="red", fontsize=font_size_text)
plt.show()
"""
'''fig4 = plt.figure()
plt.axvline(x=1, alpha=0.4, color='red', linestyle='--', linewidth=2)
plt.hist(ratio, bins=120, alpha=1, range=(-1, 5))
plt.title("Distribution of real/calc depth")
plt.show()'''
'''
land_grid = get_land_grid(b_grid)
masked_land_grid = np.ma.masked_where(np.isnan(land_grid), land_grid)
fig2 = plt.figure(figsize=(18, 9))
plt.imshow(
    masked_land_grid,
    extent=[a_lon_shift.min(), a_lon_shift.max(), a_lat.min(), a_lat.max()],
    origin="lower",
    cmap=cmap3,
    vmin=0,
    vmax=1,
    aspect=1
)
plt.grid()
plt.show()
'''