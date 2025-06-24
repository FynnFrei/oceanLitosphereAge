import numpy as np
from scripts import data_prep as data
from scripts import functions
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.stats import skew, kurtosis

# TODO fix the lon grid shift problem (it always shows 0° in center) right now o do it manually which is lazy code
# TODO Tiefseegräben identifizieren (Heatmap)
# TODO map ausschnitt von Tiefseegraben
# variables

font_path_regular = "C:/Users/fynnf/Documents/Weiteres/Fonts/static/Montserrat-Medium.ttf"
fm.fontManager.addfont(font_path_regular)  # Register font
montserrat_regular = fm.FontProperties(fname=font_path_regular).get_name()  # Get font name

# Set Montserrat as the default font
plt.rcParams["font.family"] = montserrat_regular


# colormaps-------------------------------------------------------------------------------------------------------------
#custom colormap:
cdict1 = {
    'red': (
        (0.0, 0.0, 0.0),
        (0.125, 0.0, 0.0),
        (0.48, 1.0, 1.0),
        (0.52, 1.0, 1.0),
        (0.875, 1.0, 1.0),
        (1.0, 0.5, 0.5),
    ),
    'green': (
        (0.0, 0.0, 0.0),
        (0.125, 0.0, 0.0),
        (0.48, 1.0, 1.0),
        (0.52, 1.0, 1.0),
        (0.875, 0.0, 0.0),
        (1.0, 0.0, 0.0),
    ),
    'blue': (
        (0.0, 0.5, 0.5),
        (0.125, 1.0, 1.0),
        (0.48, 1.0, 1.0),
        (0.52, 1.0, 1.0),
        (0.875, 0.0, 0.0),
        (1.0, 0.0, 0.0),
    )
}
reds_segmentdata = {
    'red': (
        (0.0,   1.0, 1.0),    # white
        (0.125, 1.0, 1.0),
        (0.25,  0.996, 0.996),
        (0.375, 0.988, 0.988),
        (0.5,   0.988, 0.988),
        (0.625, 0.984, 0.984),
        (0.75,  0.937, 0.937),
        (0.875, 0.843, 0.843),
        (1.0,   0.403, 0.403),
    ),
    'green': (
        (0.0,   1.0, 1.0),
        (0.125, 0.902, 0.902),
        (0.25,  0.816, 0.816),
        (0.375, 0.733, 0.733),
        (0.5,   0.604, 0.604),
        (0.625, 0.471, 0.471),
        (0.75,  0.318, 0.318),
        (0.875, 0.188, 0.188),
        (1.0,   0.0,   0.0),
    ),
    'blue': (
        (0.0,   1.0, 1.0),
        (0.125, 0.800, 0.800),
        (0.25,  0.643, 0.643),
        (0.375, 0.518, 0.518),
        (0.5,   0.404, 0.404),
        (0.625, 0.282, 0.282),
        (0.75,  0.220, 0.220),
        (0.875, 0.153, 0.153),
        (1.0,   0.051, 0.051),
    )
}
custom_Reds = LinearSegmentedColormap('CustomReds', reds_segmentdata)
new_seismic_colors = (
    plt.cm.Blues_r(np.linspace(-0.15, 1.05, 128)),
    np.array([1, 1, 1, 1]),     # one white value in the center
    custom_Reds(np.linspace(0,1,128))
)
colors_combined = np.vstack(new_seismic_colors)
cmap0 = LinearSegmentedColormap('BlueRed1', cdict1)
cmap1 = plt.cm.viridis
cmap2 = plt.cm.turbo_r
cmap3 = plt.cm.seismic
cmap_reds = plt.cm.Reds
cmap_new_seismic = LinearSegmentedColormap.from_list('BlueWhiteRed', colors_combined)
cmap_black = mcolors.ListedColormap([0, 0, 0, 0.9])
continental_color = (0.5, 0.5, 0.5)
transparent = (0, 0, 0, 0)
cmap0.set_bad(color=continental_color)  # gives NaN pixels specific color
cmap1.set_bad(color=continental_color)  # gives NaN pixels specific color
cmap2.set_bad(color=continental_color)
cmap3.set_bad(color=continental_color)
cmap_reds.set_bad(color=continental_color)
cmap_new_seismic.set_bad(color=continental_color)
cmap_black.set_bad(color=transparent)
font_axes_scale = 1.4  # axes values
font_header_scale = 3.4  # factor that scales the header font to the standard font size
font_label_scale = 1.9  # labels
font_text_scale = 1.7  # normal text, description, source(GEBCO)


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


# example map-----------------------------------------------------------------------------------------------------
def example_map(size):
    ratio = size * 16 / 9
    fig2, ax = plt.subplots(figsize=(ratio, size))
    min_lon = data.a_lon_shift.min()
    max_lon = data.a_lon_shift.max()
    min_lat = data.a_lat.min()
    max_lat = data.a_lat.max()
    depth_img = ax.imshow(
        data.a_grid_shift,
        extent=[min_lon, max_lon, min_lat, max_lat],
        origin="lower",
        cmap=cmap2,
        vmin=0,
        vmax=200,
        aspect=1
    )
    ax.imshow(
        data.land_grid,
        extent=[min_lon, max_lon, min_lat, max_lat],
        origin="lower",
        cmap=cmap_black,
        vmin=0,
        vmax=1,
        alpha=0.3,
        aspect=1
    )
    relative_font_size = fig2.get_size_inches()[0] * 0.5  # gets font size relative to figure size
    plt.title("Age of Oceanic Lithosphere", fontsize=relative_font_size * font_header_scale,
              pad=10)

    lon_ticks = [-140, -80, -20, 40, 100, 160]
    lon_labels = ['60°E', '120°E', '180°', '120°W', '60°W', '0°']
    ax.set_xticks(lon_ticks)
    ax.set_xticklabels(lon_labels)
    lat_ticks = [-75, -50, -25, 0, 25, 50, 75]
    lat_labels = ['75°S', '50°S', '25°S', '0°', '25°N', '50°N', '75°N']
    ax.set_yticks(lat_ticks)
    ax.set_yticklabels(lat_labels)

    plt.tick_params(labelsize=relative_font_size * font_axes_scale)
    plt.xlabel("Longitude", fontsize=relative_font_size * font_label_scale)
    plt.ylabel("Latitude", fontsize=relative_font_size * font_label_scale)
    cbar = plt.colorbar(depth_img, ax=ax, fraction=0.0232)
    cbar_ticks = np.linspace(0, 200, num=5)  # Choose appropriate tick positions
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels([f"{int(tick)}Ma" for tick in cbar_ticks])  # Convert to km and format
    cbar.ax.tick_params(labelsize=relative_font_size * font_axes_scale)  # modify font size of cbar values
    cbar.ax.text(-relative_font_size * 0.1, 0.5, "Age",
                 fontsize=relative_font_size * font_label_scale,
                 rotation=90, ha="center", va="center", transform=cbar.ax.transAxes)
    '''cbar.ax.text(-relative_font_size * 0.1, 0, "deeper than",
                 fontsize=relative_font_size * font_text_scale,
                 rotation=90, ha="center", va="bottom", transform=cbar.ax.transAxes)
    cbar.ax.text(-relative_font_size * 0.1, 0, "calculated",
                 fontsize=relative_font_size * font_text_scale,
                 rotation=90, ha="center", va="bottom", transform=cbar.ax.transAxes)
    cbar.ax.text(-relative_font_size * 0.1, 1, "shallower than",
                 fontsize=relative_font_size * font_text_scale,
                 rotation=90, ha="center", va="top", transform=cbar.ax.transAxes)
    cbar.ax.text(-relative_font_size * 0.1, 1, "calculated",
                 fontsize=relative_font_size * font_text_scale,
                 rotation=90, ha="center", va="top", transform=cbar.ax.transAxes)'''
    plt.grid()
    plt.text(-170, -85, "Ogg 2012", fontsize=relative_font_size * font_label_scale)
    plt.show()


# depth anomaly map-----------------------------------------------------------------------------------------------------
def depth_anomaly_map(size):
    ratio = size * 16 / 9

    anomaly_grid_window = functions.trim_grid(data.depth_anomaly, 0, 200, 180, data.accuracy)
    land_grid_window = functions.trim_grid(data.land_grid, 0, 200, 180, data.accuracy)

    fig2, ax = plt.subplots(figsize=(ratio, size))
    min_lon = data.a_lon_shift.min()
    max_lon = data.a_lon_shift.max()
    min_lat = data.a_lat.min()
    max_lat = data.a_lat.max()
    depth_img = ax.imshow(
        anomaly_grid_window,
        extent=[min_lon, max_lon, min_lat, max_lat],
        origin="lower",
        cmap=cmap_new_seismic,
        #norm=TwoSlopeNorm(vmin=-2000, vcenter=0, vmax=2000),
        vmin=-4000,
        vmax=4000,
        aspect=1
    )
    ax.imshow(
        land_grid_window,
        extent=[min_lon, max_lon, min_lat, max_lat],
        origin="lower",
        cmap=cmap_black,
        vmin=0,
        vmax=1,
        alpha=0.3,
        aspect=1
    )
    relative_font_size = fig2.get_size_inches()[0] * 0.45  # gets font size relative to figure size
    '''plt.title("Sea Floor Depth Difference between \nCalculation and Bathymetry", fontsize=relative_font_size * font_header_scale,
              pad=10)'''

    lon_ticks = [-140, -80, -20, 40, 100, 160]
    lon_labels = ['60°E', '120°E', '180°', '120°W', '60°W', '0°']
    ax.set_xticks(lon_ticks)
    ax.set_xticklabels(lon_labels)
    lat_ticks = [-75, -50, -25, 0, 25, 50, 75]
    lat_labels = ['75°S', '50°S', '25°S', '0°', '25°N', '50°N', '75°N']
    ax.set_yticks(lat_ticks)
    ax.set_yticklabels(lat_labels)

    plt.tick_params(labelsize=relative_font_size * font_axes_scale)
    plt.xlabel("Longitude", fontsize=relative_font_size * font_label_scale)
    plt.ylabel("Latitude", fontsize=relative_font_size * font_label_scale)
    cbar = plt.colorbar(depth_img, ax=ax, fraction=0.0232)
    cbar_ticks = np.linspace(-4000, 4000, num=5)  # Choose appropriate tick positions
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels([f"{tick * 0.001}km" for tick in cbar_ticks])  # Convert to km and format
    cbar.ax.tick_params(labelsize=relative_font_size * font_axes_scale)  # modify font size of cbar values
    # cbar.set_label(label="Δ Depth in m", fontsize=relative_font_size*font_label_scale)
    '''cbar.ax.text(-relative_font_size * 0.14, 0, "deeper than calculated",
                 fontsize=relative_font_size * font_text_scale,
                 rotation=90, ha="center", va="bottom", transform=cbar.ax.transAxes)'''
    cbar.ax.text(-relative_font_size * 0.05, 0, "deeper than calculated",
                 fontsize=relative_font_size * font_text_scale,
                 rotation=90, ha="center", va="bottom", transform=cbar.ax.transAxes)
    '''cbar.ax.text(-relative_font_size * 0.14, 1, "shallower than calculated",
                 fontsize=relative_font_size * font_text_scale,
                 rotation=90, ha="center", va="top", transform=cbar.ax.transAxes)'''
    cbar.ax.text(-relative_font_size * 0.05, 1, "shallower than calculated",
                 fontsize=relative_font_size * font_text_scale,
                 rotation=90, ha="center", va="top", transform=cbar.ax.transAxes)
    # plt.grid()
    plt.text(-170, -85, "GEBCO_2020 Grid, Ogg 2012", fontsize=relative_font_size * font_text_scale*0.8)
    plt.savefig("plots/depth_anomaly_map.png", dpi=600, bbox_inches='tight')
    plt.show()


# depth absolute anomaly map--------------------------------------------------------------------------------------------
def depth_abs_anomaly_map(size):
    ratio = size * 16 / 9
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
    relative_font_size = fig2.get_size_inches()[0] * 0.45  # gets font size relative to figure size
    '''plt.title("Absolute Sea Floor Depth Difference \nbetween Calculation and Bathymetry",
              fontsize=relative_font_size * font_header_scale, pad=10)'''

    lon_ticks = [-140, -80, -20, 40, 100, 160]
    lon_labels = ['60°E', '120°E', '180°', '120°W', '60°W', '0°']
    ax.set_xticks(lon_ticks)
    ax.set_xticklabels(lon_labels)
    lat_ticks = [-75, -50, -25, 0, 25, 50, 75]
    lat_labels = ['75°S', '50°S', '25°S', '0°', '25°N', '50°N', '75°N']
    ax.set_yticks(lat_ticks)
    ax.set_yticklabels(lat_labels)

    plt.tick_params(labelsize=relative_font_size * font_axes_scale)
    plt.xlabel("Longitude", fontsize=relative_font_size * font_label_scale)
    plt.ylabel("Latitude", fontsize=relative_font_size * font_label_scale)
    cbar = plt.colorbar(depth_img, ax=ax, fraction=0.0232)
    cbar_ticks = np.linspace(0, 4000, num=5)  # Choose appropriate tick positions
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels([f"{tick * 0.001}" for tick in cbar_ticks])  # Convert to km and format
    cbar.ax.tick_params(labelsize=relative_font_size * font_axes_scale)  # modify font size of cbar values
    # cbar.set_label(label="Δ Depth in m", fontsize=relative_font_size*font_label_scale)
    cbar.ax.text(-relative_font_size * 0.15, 0.5, "Absolute depth anomaly in km",
                 fontsize=relative_font_size * font_label_scale,
                 rotation=90, ha="center", va="center", transform=cbar.ax.transAxes)
    plt.grid()
    plt.text(-170, -85, "GEBCO_2020 Grid, Ogg 2012", fontsize=relative_font_size * font_label_scale)
    # plt.savefig("plots\depth_anomaly_hires.png", dpi=300, bbox_inches='tight')
    plt.show()


# functions.plot_scatter(d_grid, b_grid_shift_ocean)
def depth_anomaly_scatter(size):
    #mask_positive = data.grid2 > 0  # where bathymetry is > 0 -> land
    #mask_negative = data.grid2 < 0  # where bathymetry is subsurface
    grid_point_diffs, grid1, grid2 = functions.compare_grid_points(data.d_grid, data.b_grid_shift_ocean)

    optimum_line = np.linspace(-10000, 0, 2)
    x_depth_max = -10000
    x_depth_min = 0
    y_depth_max = 0     # range of real depth display
    y_depth_min = -10000
    x_depth_max_tick = -10000
    x_depth_min_tick = 0
    num_of_x_ticks = 6
    ax1_ticks = np.linspace(x_depth_min_tick, x_depth_max_tick, num_of_x_ticks)
    age_values = [None, None]  # stores age values for x_age-axis
    for i in range(2, num_of_x_ticks):
        age = functions.depth_to_age(ax1_ticks[i], data.avr_ridge_depth)
        age_values.append(age)
    '''# age axis manually:
    ax2_ticks = np.array([data.avr_ridge_depth, 4000, 6000, 8000])
    for i in range(len(ax2_ticks)):
        age = functions.depth_to_age(ax2_ticks[i], data.avr_ridge_depth)
        age_values.append(age)'''

    ratio = size * 12 / 9
    fig = plt.figure(figsize=(ratio, size))
    relative_font_size = fig.get_size_inches()[0] * 0.6  # gets font size relative to figure size
    ax1 = fig.add_subplot()

    hb = ax1.hexbin(grid1, grid2, gridsize=500, cmap=cmap1, bins='log', mincnt=1)
    plt.axhline(y=0, alpha=1, color='orange', linestyle='--', linewidth=1)
    plt.plot(optimum_line, optimum_line, color="red", linewidth=1.5)
    plt.axvline(x=data.avr_ridge_depth, alpha=0.7, color='orange', linestyle='--', linewidth=2)

    ax1.set_xlim(x_depth_min, x_depth_max)
    ax1.set_ylim(y_depth_min, y_depth_max)
    ax1.set_xticks(ax1_ticks)  # sets x_depth-axis labels
    '''plt.title("Bathymetric over calculated depth", fontsize=relative_font_size*font_header_scale,
              pad=int(relative_font_size))'''
    '''plt.text(0.98, 0.02, "GEBCO_2020 Grid, Ogg 2012",
             fontsize=relative_font_size*font_text_scale,
             ha="right", va="bottom", transform=ax1.transAxes)'''
    ax1.set_xlabel("Calculated depth in m", fontsize=relative_font_size*font_label_scale)
    ax1.set_ylabel("Bathymetric depth in m", fontsize=relative_font_size*font_label_scale)

    # Age axis
    ax2 = ax1.twiny()
    ax2.set_xlim(ax1.get_xlim())  # so that ax1 and ax2 ticks match positions
    ax2.set_xticks(ax1_ticks)
    ax2.set_xticklabels(["", "", f"{age_values[2]:.1f}", f"{age_values[3]:.1f}", f"{age_values[4]:.1f}", f"{age_values[5]:.1f}"])
    plt.tick_params(labelsize=relative_font_size * font_axes_scale)
    ax1.tick_params(labelsize=relative_font_size * font_axes_scale)
    ax2.set_xlabel("Age in Ma", fontsize=relative_font_size*font_label_scale)
    '''
    # Inset colorbar
    cax = inset_axes(ax1, width="30%", height="4%", loc="lower right", borderpad=2)
    cb = plt.colorbar(hb, cax=cax, orientation='horizontal')
    cb.set_label("log10(Density)", fontsize=10)
    cb.ax.xaxis.set_label_position('top')
    cb.ax.tick_params(labelsize=8)
    '''

    cb = fig.colorbar(hb, ax=ax1)
    cb.set_label("Hit count per hexagon", fontsize=relative_font_size*font_label_scale)
    cb.ax.tick_params(labelsize=relative_font_size*font_axes_scale)

    # legend
    legend_elements = [
        Line2D([0], [0], color='orange', linestyle='--', linewidth=2, label="lithosphere age = 0 Ma"),
        Line2D([0], [0], color='red', linewidth=2, label="calc. depth = bath. depth")
    ]
    ax1.legend(
        handles=legend_elements,
        loc="lower left",
        fontsize=relative_font_size*font_text_scale,
        frameon=True
    )
    ax1.grid()
    plt.savefig("plots/real_over_calc_depth_submarine_hires.png", dpi=300, bbox_inches='tight')
    #plt.show()


# depth anomaly histogram-----------------------------------------------------------------------------------------------
def depth_anomaly_hist(size, x_min, x_max, bins, cumulative=False):
    if data.accuracy == 6:  # how many data points do we have
        y_factor = 10
    else:
        y_factor = 1
    y_max = int(2500000000 / bins)
    ratio = size * 3/3  # figure dimensions

    # single grid point comparison (for histogram)
    grid_point_diffs, grid1, grid2 = functions.compare_grid_points(data.d_grid, data.b_grid_shift_ocean)

    # grid data preparation to area
    lon_num = len(data.depth_anomaly[0])  # 3600 for accuracy=6
    lat_num = len(data.depth_anomaly[:])  # 1800
    U_equator = 40075
    R_equator = 6378.137  # in km
    pixel_height = U_equator / lon_num  # height is independent of lat and equal for every pixel
    lat_idx = np.repeat(np.arange(lat_num), lon_num)
    lat_val = (- lat_idx + 0.5 * len(data.depth_anomaly)) / (lat_num / 180)
    pixel_width = 2 * np.pi * R_equator * np.cos(np.pi * lat_val / 180) / lon_num  # width of each pixel
    pixel_area = pixel_height * pixel_width
    grid_1D = np.reshape(data.depth_anomaly, -1)  # 1D bathymetric grid
    grid_1D_area = np.stack((grid_1D, pixel_area), axis=1)  # now every grid point contains an area value
    grid_1D_clean = grid_1D_area[~np.isnan(grid_1D_area[:, 0])]  # delete NaNs
    print(f'shape: {np.shape(grid_1D_clean)}')
    # statistics
    mean_val = np.average(grid_1D_clean[:, 0], weights=grid_1D_clean[:, 1])
    std_dev = np.sqrt(np.average((grid_1D_clean[:, 0] - mean_val) ** 2, weights=grid_1D_clean[:, 1]))
    skewness = np.average((grid_1D_clean[:, 0] - mean_val) ** 3, weights=grid_1D_clean[:, 1]) / (std_dev ** 3)
    kurt_excess = np.average((grid_1D_clean[:, 0] - mean_val) ** 4, weights=grid_1D_clean[:, 1]) / (
                std_dev ** 4) - 3  # difference to gauss curve
    stats_text = (
        f"Mean = {mean_val * 0.001:.3f} km\n"
        f"Std Dev = {std_dev * 0.001:.3f} km\n"
        f"Skewness = {skewness:.3f}\n"
        f"Excess Kurtosis = {kurt_excess:.3f}"
    )

    # Compute histogram data without plotting - this section is for colormap in histogram
    counts, bin_edges = np.histogram(grid_1D_clean[:,0], weights=grid_1D_clean[:,1], bins=bins, range=(x_min, x_max))
    # Normalize bin centers for colormap
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    norm = plt.Normalize(vmin=-4000, vmax=4000)
    colors = cmap_new_seismic(norm(bin_centers))

    fig3, ax1 = plt.subplots(figsize=(ratio, size))
    relative_font_size = fig3.get_size_inches()[0] * 0.8  # gets font size relative to figure size
    plt.axvline(x=0, alpha=0.4, color='black', linestyle='--', linewidth=2)
    #plt.hist(grid_point_diffs, bins=bins, alpha=1, range=(x_min, x_max), color=colors, edgecolor='black')
    plt.bar(bin_edges[:-1], counts, width=np.diff(bin_edges), color=colors, edgecolor='black', align='edge')
    plt.tick_params(labelsize=relative_font_size * font_axes_scale)
    x_ticks = np.linspace(x_min, x_max, num=5)  # Choose appropriate tick positions
    ax1.set_xlim(x_min, x_max+1)
    ax1.set_xticks(x_ticks)
    ax1.set_xticklabels([f"{int(tick * 0.001)}" for tick in x_ticks])  # Convert to km and format
    y_ticks = np.linspace(y_max / 5, y_max, num=5)  # Choose appropriate tick positions
    ax1.set_yticks(y_ticks)
    ax1.set_yticklabels([f"{int(tick * 0.000001)}" for tick in y_ticks])  # Convert to km and format
    # plt.title("Quantitative distribution \nof depth anomaly", fontsize=relative_font_size * font_header_scale, pad=10)
    '''plt.text(0.98, 0.98, "GEBCO_2020 Grid, Ogg 2012",
             fontsize=relative_font_size * font_text_scale,
             ha="right", va="top", transform=ax1.transAxes)'''
    plt.xlabel("Depth difference between Calculation and Bathymetry", fontsize=relative_font_size * font_label_scale)
    ax1.text(0.045, 0.55, "Deeper than\ncalculated",
            fontsize=relative_font_size * font_text_scale, color=cmap_new_seismic(0.1),
            ha="left", va="bottom", transform=ax1.transAxes)
    ax1.text(0.55, 0.55, "Shallower than\ncalculated",
            fontsize=relative_font_size * font_text_scale, color=cmap_new_seismic(0.9),
            ha="left", va="bottom", transform=ax1.transAxes)
    plt.ylabel(r'Area in million km$^2$', color='black', fontsize=relative_font_size * font_label_scale)
    ax1.grid(axis='x')

    # second y-axis for cumulative graph
    if cumulative:
        step_size = (abs(x_min) + abs(x_max)) / bins    # step size for cumulative x-Axis
        x_cumulative, y_cumulative = functions.cumulative_function(grid_point_diffs, x_min, x_max, step_size)
        ax2 = ax1.twinx()
        plt.tick_params(labelsize=relative_font_size * font_axes_scale)
        ax2.plot(x_cumulative, y_cumulative, color="purple", linewidth=3)  # cumulative plot (from 0-100%)
        ax2.fill_between(x_cumulative, y_cumulative, 0, color='purple', alpha=0.1)
        ax2.set_ylabel("Cumulative Distribution in %", color='purple', fontsize=relative_font_size * font_label_scale)
        ax2.set_ylim(0, 1.1)
        y_tick_vals = np.arange(0, 1.01, 0.2)
        ax2.set_yticks(y_tick_vals)
        ax2.set_yticklabels([f"{v * 100:.0f}" for v in y_tick_vals])
        ax2.grid(True)
    # plot for legend
    plt.plot([], [], ' ', label=stats_text)
    plt.legend(loc='upper right', fontsize=relative_font_size * font_text_scale)
    #plt.savefig("plots/histogram_depthdiff.png", dpi=300, bbox_inches='tight')
    #plt.show()


# bathymetric depth histogram-----------------------------------------------------------------------------------------------
def b_depth_hist(size, x_min, x_max, bins, cumulative=False):
    if data.accuracy == 2:  # how many data points do we have
        y_factor = 10
    else:
        y_factor = 1
    y_max = int(1000000000 / bins) * y_factor
    ratio = size * 3/3

    # grid data preparation to area
    lon_num = len(data.b_grid_shift_ocean[0])   # 3600 for accuracy=6
    lat_num = len(data.b_grid_shift_ocean[:])   # 1800
    U_equator = 40075
    R_equator = 6378.137    # in km
    pixel_height = U_equator / lon_num  # height is independent of lat and equal for every pixel
    lat_idx = np.repeat(np.arange(lat_num), lon_num)
    lat_val = (- lat_idx + 0.5 * len(data.b_grid_shift_ocean)) / (lat_num/180)
    pixel_width = 2*np.pi*R_equator*np.cos(np.pi*lat_val/180) / lon_num     # width of each pixel
    pixel_area = pixel_height * pixel_width
    b_1D = np.reshape(data.b_grid_shift_ocean, -1)  # 1D bathymetric grid
    b_1D_area = np.stack((b_1D, pixel_area), axis=1)   # now every grid point contains an area value
    b_1D_clean = b_1D_area[~np.isnan(b_1D_area[:, 0]) & (b_1D_area[:, 0] < 0)]  # delete NaNs and positive values
    print(np.shape(b_1D_clean))
    print(b_1D_clean[1500000,1])
    # statistics
    mean_val = np.average(b_1D_clean[:,0], weights=b_1D_clean[:,1])
    std_dev = np.sqrt(np.average((b_1D_clean[:,0] - mean_val)**2, weights=b_1D_clean[:,1]))
    skewness = np.average((b_1D_clean[:,0] - mean_val)**3, weights=b_1D_clean[:,1]) / (std_dev**3)
    kurt_excess = np.average((b_1D_clean[:,0] - mean_val)**4, weights=b_1D_clean[:,1]) / (std_dev**4) - 3  # difference to gauss curve
    kurt_normal = kurtosis(b_1D, fisher=False)
    stats_text = (
        f"Mean = {mean_val * 0.001:.3f} km\n"
        f"Std Dev = {std_dev * 0.001:.3f} km\n"
        f"Skewness = {skewness:.3f}\n"
        f"Excess Kurtosis = {kurt_excess:.3f}"
    )

    # plot
    fig3, ax1 = plt.subplots(figsize=(ratio, size))
    relative_font_size = fig3.get_size_inches()[0] * 0.8  # gets font size relative to figure size
    plt.hist(b_1D_clean[:,0], weights=b_1D_clean[:,1], bins=bins, alpha=1, range=(x_min, x_max), color=(0.05, 0.15, 0.5), edgecolor='black')
    plt.tick_params(labelsize=relative_font_size * font_axes_scale)
    x_ticks = np.linspace(x_min, x_max, num=int(abs(x_max-x_min)/1000 + 1))  # Choose appropriate tick positions
    ax1.set_xlim(x_min, x_max+1)
    ax1.set_xticks(x_ticks)
    ax1.set_xticklabels([f"{int(tick * 0.001)}.0" for tick in x_ticks])  # Convert to km and format
    y_ticks = np.linspace(y_max / 4, y_max, num=4)  # Choose appropriate tick positions
    ax1.set_yticks(y_ticks)
    ax1.set_yticklabels([f"{int(tick * 0.000001)}*10^6" for tick in y_ticks])  # Convert to km and format
    plt.text(0.98, 0.98, "GEBCO_2020 Grid",
             fontsize=relative_font_size * font_text_scale,
             ha="right", va="top", transform=ax1.transAxes)
    plt.xlabel("Bathymetric depth in km", fontsize=relative_font_size * font_label_scale)
    plt.ylabel("Area in km^2", color=(0.05, 0.15, 0.5), fontsize=relative_font_size * font_label_scale)
    ax1.grid(axis='x')

    # second y-axis for cumulative graph
    if cumulative:
        step_size = (abs(x_min) + abs(x_max)) / bins  # step size for cumulative x-Axis
        cumulative_x_axis_b, cumulative_y_axis_b = functions.cumulative_function(b_1D_clean, x_min, x_max, step_size)
        ax2 = ax1.twinx()
        plt.tick_params(labelsize=relative_font_size * font_axes_scale)
        ax2.plot(cumulative_x_axis_b, cumulative_y_axis_b, color="red", linewidth=3)  # cumulative plot (from 0-100%)
        ax2.fill_between(cumulative_x_axis_b, cumulative_y_axis_b, 0, color='red', alpha=0.1)
        ax2.set_ylabel("Cumulative Distribution in %", color="red", fontsize=relative_font_size * font_label_scale)
        ax2.set_ylim(0, 1.1)
        y_tick_vals = np.arange(0, 1.01, 0.2)
        ax2.set_yticks(y_tick_vals)
        ax2.set_yticklabels([f"{v * 100:.0f}" for v in y_tick_vals])
        ax2.grid(True)
    # Dummy plot for legend
    plt.plot([], [], ' ', label=stats_text)
    plt.legend(loc='center right', fontsize=relative_font_size * font_text_scale)

    plt.show()


def a_grid_hist(size, bins, cumulative=False):
    # grid data preparation
    a_1D = np.reshape(data.a_grid, -1)  # 1D age grid
    a_1D_clean = a_1D[~np.isnan(a_1D)]  # delete NaNs
    # statistics
    mean_val = np.mean(a_1D_clean)
    std_dev = np.std(a_1D_clean)
    skewness = skew(a_1D_clean)
    kurt_excess = kurtosis(a_1D_clean)  # difference to gauss curve
    kurt_normal = kurtosis(a_1D_clean, fisher=False)
    stats_text = (
        f"Mean = {mean_val * 0.001:.3f} km\n"
        f"Std Dev = {std_dev * 0.001:.3f} km\n"
        f"Skewness = {skewness:.3f}\n"
        f"Excess Kurtosis = {kurt_excess:.3f}"
    )

    x_min, x_max = 0, 200
    dx = x_max / bins
    y_max = int(10000000 / bins)
    ratio = size

    # plot
    fig, ax1 = plt.subplots(figsize=(ratio, size))
    relative_font_size = fig.get_size_inches()[0] * 0.8  # gets font size relative to figure size
    plt.hist(a_1D_clean, bins=bins, alpha=1, range=(x_min, x_max), color=(0.1, 0.2, 0.5), edgecolor='black')
    plt.tick_params(labelsize=relative_font_size * font_axes_scale)
    x_ticks = np.linspace(x_min, x_max, num=int(abs(x_max - x_min) / 40 + 1))  # Choose appropriate tick positions
    ax1.set_xlim(x_min, x_max + 1)
    ax1.set_ylim(0, y_max +1)
    ax1.set_xticks(x_ticks)
    y_ticks = np.linspace(y_max/4, y_max, num=4)  # Choose appropriate tick positions
    ax1.set_yticks(y_ticks)
    ax1.set_yticklabels([f"{int(tick * 0.001)}k" for tick in y_ticks])  # Convert to k and format
    plt.text(0.98, 0.98, "Ogg 2012",
             fontsize=relative_font_size * font_text_scale,
             ha="right", va="top", transform=ax1.transAxes)
    plt.xlabel("Age in Ma", fontsize=relative_font_size * font_label_scale)
    plt.ylabel("Amount of data points in #", color=(0.1, 0.2, 0.5), fontsize=relative_font_size * font_label_scale)
    ax1.grid(axis='x')

    # second y-axis for cumulative graph
    if cumulative:
        x, y2 = functions.cumulative_function(a_1D_clean, x_min, x_max, dx)
        ax2 = ax1.twinx()
        plt.tick_params(labelsize=relative_font_size * font_axes_scale)
        ax2.plot(x, y2, color="red", linewidth=3)  # cumulative plot (from 0-100%)
        ax2.fill_between(x, y2, 0, color='red', alpha=0.1)
        ax2.set_ylabel("Cumulative Distribution in %", color="red", fontsize=relative_font_size * font_label_scale)
        ax2.set_ylim(0, 1.1)
        y_tick_vals = np.arange(0, 1.01, 0.2)
        ax2.set_yticks(y_tick_vals)
        ax2.set_yticklabels([f"{v * 100:.0f}" for v in y_tick_vals])
        ax2.grid(True)
    # Dummy plot for legend
    plt.plot([], [], ' ', label=stats_text)
    plt.legend(loc='center right', fontsize=relative_font_size * font_text_scale)
    plt.show()


'''
# functions.plot_scatter(d_grid, b_grid_shift_ocean)
def old_depth_anomaly_scatter():
    mask_positive = data.grid2 > 0  # where bathymetry is > 0 -> land

    optimum_line = np.linspace(-8000, -2000, 2)
    x_depth_max = -8500
    x_depth_min = data.avr_ridge_depth
    x_depth_max_tick = -8000
    x_depth_min_tick = -3000
    num_of_x_ticks = 6
    ax1_ticks = np.linspace(x_depth_min_tick, x_depth_max_tick, num_of_x_ticks)
    age_values = []  # stores age values for x_age-axis
    for i in range(num_of_x_ticks):
        if ax1_ticks[i] > data.avr_ridge_depth:
            age = 0  # avoid negative age values
        else:
            age = functions.depth_to_age(ax1_ticks[i], data.avr_ridge_depth)
        age_values.append(age)
    print(age_values)

    fig = plt.figure(figsize=(15, 8))
    ax1 = fig.add_subplot()
    plt.axvline(x=data.avr_ridge_depth, alpha=0.4, color='orange', linestyle='--', linewidth=2)
    plt.scatter(data.grid1[~mask_positive], data.grid2[~mask_positive], s=0.01, color="blue")
    plt.scatter(data.grid1[mask_positive], data.grid2[mask_positive], s=0.01, color="green")  # above sea level
    plt.plot(optimum_line, optimum_line, color="red")
    ax1.set_xlim(x_depth_min, x_depth_max)
    ax1.set_xticks(ax1_ticks)  # sets x_depth-axis labels
    #plt.gca().invert_xaxis()
    plt.title("Real Depth Over Calculated Depth", fontsize=28)
    plt.text(-2000, -10000, "GEBCO_2020 Grid, Ogg 2012")
    ax1.set_xlabel("calculated depth", fontsize=16)
    ax1.set_ylabel("real depth", fontsize=16)

    # Age axis
    ax2 = ax1.twiny()
    ax2.set_xlim(ax1.get_xlim())  # so that ax1 and ax2 ticks match positions
    ax2.set_xticks(ax1_ticks)
    ax2.set_xticklabels([f"{age:.1f} Ma" for age in age_values])
    # Custom legend markers for larger dots
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=10, label="Below Sea Level"),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='green', markersize=10, label="Above Sea Level"),
        Line2D([0], [0], color='red', linewidth=2, label="calc. depth = real depth")
    ]
    plt.legend(handles=legend_elements, loc="upper right", fontsize=14)
    plt.grid()
    plt.show()
'''


def sqrt_t():
    # Define the x range, avoiding x=0 to prevent division by zero
    x = np.linspace(0.1, 5, 400)
    y = -np.sqrt(x)

    # Plotting
    fig3, ax1 = plt.subplots(figsize=(10, 6))
    relative_font_size = fig3.get_size_inches()[0]  # gets font size relative to figure size
    plt.plot(x, y, color=cmap_new_seismic(0.9))
    plt.xlabel('Age', fontsize=relative_font_size * font_label_scale)
    plt.ylabel('Sea floor depth', fontsize=relative_font_size * font_label_scale)
    plt.text(0.65, 0.6, r'$depth \propto \sqrt{t}$',
             fontsize=relative_font_size * font_text_scale * 1.5,
             ha="right", va="top", transform=ax1.transAxes)
    plt.xticks([])
    plt.yticks([])
    plt.grid(False)
    plt.show()
