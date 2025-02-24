import numpy as np
from scripts import data_prep as data
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors

# TODO hist and all maps typo and scales
# TODO plot actual depth over expected depth for all grid points in absolute
# TODO heat map for scatter
# variables

cmap1 = plt.cm.viridis
cmap2 = plt.cm.turbo_r
cmap3 = plt.cm.seismic
cmap_reds = plt.cm.Reds
cmap_black = mcolors.ListedColormap([0, 0, 0, 0.9])
continental_color = (0.5, 0.5, 0.5)
transparent = (0, 0, 0, 0)
cmap1.set_bad(color=continental_color)  # gives NaN pixels specific color
cmap2.set_bad(color=continental_color)
cmap3.set_bad(color=continental_color)
cmap_reds.set_bad(color=continental_color)
cmap_black.set_bad(color=transparent)
font_axes_scale = 1.8  # axes values
font_header_scale = 3.8  # factor that scales the header font to the standard font size
font_label_scale = 2.5  # labels, source(GEBCO)
font_text_scale = 2  # normal text, description


# all maps-------------------------------------------------------------------------------------------------------------
def all_maps():
    fig1, ax = plt.subplots(2, 2, figsize=(22, 10))

    mesh1 = ax[0][0].imshow(
        data.a_grid_shift,
        extent=[data.a_lon_shift.min(), data.a_lon_shift.max(), data.a_lat.min(), data.a_lat.max()],
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
        data.b_grid_shift_ocean,
        extent=[data.b_lon_shift.min(), data.b_lon_shift.max(), data.b_lat.min(), data.b_lat.max()],
        origin="lower",
        cmap=cmap1, vmin=-7000, vmax=0, aspect=1
    )
    ax[0][1].set_title("Bathymetric map")
    ax[0][1].text(-170, -85, "GEBCO_2020 Grid")
    ax[0][1].set_aspect(1)
    fig1.colorbar(mesh2, ax=ax[0][1], label="Depth in m")
    ax[0][1].grid()

    mesh3 = ax[1][0].imshow(
        data.d_grid,
        extent=[data.b_lon_shift.min(), data.b_lon_shift.max(), data.b_lat.min(), data.b_lat.max()],
        origin="lower",
        cmap=cmap1, vmin=-7000, vmax=-1000, aspect=1
    )
    ax[1][0].set_title("Calculated depth")
    ax[1][0].text(-170, -85, "Ogg 2012")
    ax[1][0].set_aspect(1)
    fig1.colorbar(mesh3, ax=ax[1][0], label="Depth in m")
    ax[1][0].grid()

    mesh4 = ax[1][1].imshow(
        data.depth_anomaly,
        extent=[data.b_lon_shift.min(), data.b_lon_shift.max(), data.b_lat.min(), data.b_lat.max()],
        origin="lower",
        cmap=cmap3, vmin=-3000, vmax=3000, aspect=1
    )
    ax[1][1].set_title("Calculated depth - bathymetric depth")
    ax[1][1].text(-170, -85, "GEBCO_2020 Grid, Ogg 2012")
    ax[1][1].set_aspect(1)
    fig1.colorbar(mesh4, ax=ax[1][1], label="Δ Depth in m")
    ax[1][1].grid()

    plt.show()


# depth anomaly map-----------------------------------------------------------------------------------------------------
def depth_anomaly_map(size):
    ratio = size * 16/9
    fig2, ax = plt.subplots(figsize=(ratio, size))
    depth_img = ax.imshow(
        data.depth_anomaly,
        extent=[data.a_lon_shift.min(), data.a_lon_shift.max(), data.a_lat.min(), data.a_lat.max()],
        origin="lower",
        cmap=cmap3,
        vmin=-4000,
        vmax=4000,
        aspect=1
    )
    ax.imshow(
        data.land_grid,
        extent=[data.a_lon_shift.min(), data.a_lon_shift.max(), data.a_lat.min(), data.a_lat.max()],
        origin="lower",
        cmap=cmap_black,
        vmin=0,
        vmax=1,
        alpha=0.3,
        aspect=1
    )
    relative_font_size = fig2.get_size_inches()[0] * 0.5  # gets font size relative to figure size
    plt.title("Difference between \nCalculated and Bathymetric Depth", fontsize=relative_font_size * font_header_scale,
              pad=10)
    plt.tick_params(labelsize=relative_font_size * font_axes_scale)
    plt.xlabel("Longitude", fontsize=relative_font_size * font_label_scale)
    plt.ylabel("Latitude", fontsize=relative_font_size * font_label_scale)
    cbar = plt.colorbar(depth_img, ax=ax, fraction=0.0232)
    cbar_ticks = np.linspace(-4000, 4000, num=5)  # Choose appropriate tick positions
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels([f"{tick * 0.001}km" for tick in cbar_ticks])  # Convert to km and format
    cbar.ax.tick_params(labelsize=relative_font_size * font_axes_scale)  # modify font size of cbar values
    # cbar.set_label(label="Δ Depth in m", fontsize=relative_font_size*font_label_scale)
    cbar.ax.text(-relative_font_size * 0.1, 0, "deeper than calculated",
                 fontsize=relative_font_size * font_text_scale,
                 rotation=90, ha="center", va="bottom", transform=cbar.ax.transAxes)
    cbar.ax.text(-relative_font_size * 0.1, 1, "shallower than calculated",
                 fontsize=relative_font_size * font_text_scale,
                 rotation=90, ha="center", va="top", transform=cbar.ax.transAxes)
    plt.grid()
    plt.text(-170, -85, "GEBCO_2020 Grid, Ogg 1012", fontsize=relative_font_size * font_label_scale)
    plt.show()


# depth absolute anomaly map--------------------------------------------------------------------------------------------
def depth_abs_anomaly_map(size):
    ratio = size * 16/9
    fig2, ax = plt.subplots(figsize=(ratio, size))
    depth_img = ax.imshow(
        data.depth_anomaly_abs,
        extent=[data.a_lon_shift.min(), data.a_lon_shift.max(), data.a_lat.min(), data.a_lat.max()],
        origin="lower",
        cmap=cmap_reds,
        vmin=0,
        vmax=4000,
        aspect=1
    )
    ax.imshow(
        data.land_grid,
        extent=[data.a_lon_shift.min(), data.a_lon_shift.max(), data.a_lat.min(), data.a_lat.max()],
        origin="lower",
        cmap=cmap_black,
        alpha=0.3,
        aspect=1
    )
    relative_font_size = fig2.get_size_inches()[0] * 0.5  # gets font size relative to figure size
    plt.title("Absolute Difference between\n Calculated and Bathymetric Depth",
              fontsize=relative_font_size * font_header_scale, pad=10)
    plt.tick_params(labelsize=relative_font_size * font_axes_scale)
    plt.xlabel("Longitude", fontsize=relative_font_size * font_label_scale)
    plt.ylabel("Latitude", fontsize=relative_font_size * font_label_scale)
    cbar = plt.colorbar(depth_img, ax=ax, fraction=0.0232)
    cbar_ticks = np.linspace(0, 4000, num=5)  # Choose appropriate tick positions
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels([f"{tick * 0.001}km" for tick in cbar_ticks])  # Convert to km and format
    cbar.ax.tick_params(labelsize=relative_font_size * font_axes_scale)  # modify font size of cbar values
    # cbar.set_label(label="Δ Depth in m", fontsize=relative_font_size*font_label_scale)
    cbar.ax.text(-relative_font_size * 0.1, 0.5, "Δ Depth in km",
                 fontsize=relative_font_size * font_label_scale,
                 rotation=90, ha="center", va="center", transform=cbar.ax.transAxes)
    plt.grid()
    plt.text(-170, -85, "GEBCO_2020 Grid, Ogg 1012", fontsize=relative_font_size * font_label_scale)
    plt.show()


# functions.plot_scatter(d_grid, b_grid_shift_ocean)

# depth anomaly histogram-----------------------------------------------------------------------------------------------
def depth_anomaly_hist():
    grid1_flat = np.reshape(data.d_grid, -1)
    grid2_flat = np.reshape(data.b_grid_shift_ocean, -1)
    mask_nan = ~np.isnan(grid1_flat) & ~np.isnan(grid2_flat)
    grid1_flat_clean = grid1_flat[mask_nan]
    grid2_flat_clean = grid2_flat[mask_nan]
    diff = grid2_flat_clean - grid1_flat_clean

    plt.figure()
    plt.axvline(x=0, alpha=0.4, color='red', linestyle='--', linewidth=2)
    plt.hist(diff, bins=100, alpha=1, range=(-3000, 7000))
    plt.title("Distribution of Δ depth", fontsize=font_header_scale)
    plt.text(5800, 1.95 * 10 ** 6, "GEBCO_2020 Grid, Ogg 1012")
    plt.xlabel("anomaly in m", fontsize=font_label_scale)
    plt.ylabel("amount of data points", fontsize=font_label_scale)
    plt.text(-3000, 8 * 10 ** 5, "Deeper than calculation", color="blue", fontsize=font_text_scale)
    plt.text(2000, 8 * 10 ** 5, "Shallower than calculation", color="red", fontsize=font_text_scale)
    plt.show()

