import numpy as np
from scripts import data_prep as data
from scripts import functions
from matplotlib import pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# TODO all_maps typo and scales
# TODO fix the lon grid shift problem (it always shows 0° in center) right now o do it manually which is bad code
# variables

font_path_regular = "C:/Users/fynnf/Documents/Weiteres/Fonts/static/Montserrat-Medium.ttf"
fm.fontManager.addfont(font_path_regular)  # Register font
montserrat_regular = fm.FontProperties(fname=font_path_regular).get_name()  # Get font name

# Set Montserrat as the default font
plt.rcParams["font.family"] = montserrat_regular

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
font_header_scale = 3.4  # factor that scales the header font to the standard font size
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
    fig2, ax = plt.subplots(figsize=(ratio, size))
    min_lon = data.a_lon_shift.min()
    max_lon = data.a_lon_shift.max()
    min_lat = data.a_lat.min()
    max_lat = data.a_lat.max()
    depth_img = ax.imshow(
        data.depth_anomaly,
        extent=[min_lon, max_lon, min_lat, max_lat],
        origin="lower",
        cmap=cmap3,
        vmin=-4000,
        vmax=4000,
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
    relative_font_size = fig2.get_size_inches()[0] * 0.78  # gets font size relative to figure size
    plt.title("Sea Floor Depth Difference between \nCalculation and Bathymetry", fontsize=relative_font_size * font_header_scale,
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
    cbar_ticks = np.linspace(-4000, 4000, num=5)  # Choose appropriate tick positions
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels([f"{tick * 0.001}km" for tick in cbar_ticks])  # Convert to km and format
    cbar.ax.tick_params(labelsize=relative_font_size * font_axes_scale)  # modify font size of cbar values
    # cbar.set_label(label="Δ Depth in m", fontsize=relative_font_size*font_label_scale)
    cbar.ax.text(-relative_font_size * 0.3, 0, "deeper than",
                 fontsize=relative_font_size * font_text_scale,
                 rotation=90, ha="center", va="bottom", transform=cbar.ax.transAxes)
    cbar.ax.text(-relative_font_size * 0.1, 0, "calculated",
                 fontsize=relative_font_size * font_text_scale,
                 rotation=90, ha="center", va="bottom", transform=cbar.ax.transAxes)
    cbar.ax.text(-relative_font_size * 0.3, 1, "shallower than",
                 fontsize=relative_font_size * font_text_scale,
                 rotation=90, ha="center", va="top", transform=cbar.ax.transAxes)
    cbar.ax.text(-relative_font_size * 0.1, 1, "calculated",
                 fontsize=relative_font_size * font_text_scale,
                 rotation=90, ha="center", va="top", transform=cbar.ax.transAxes)
    plt.grid()
    plt.text(-170, -85, "GEBCO_2020 Grid, Ogg 2012", fontsize=relative_font_size * font_label_scale)
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
    relative_font_size = fig2.get_size_inches()[0] * 0.7  # gets font size relative to figure size
    plt.title("Absolute Sea Floor Depth Difference \nbetween Calculation and Bathymetry",
              fontsize=relative_font_size * font_header_scale, pad=10)

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
    cbar.set_ticklabels([f"{tick * 0.001}km" for tick in cbar_ticks])  # Convert to km and format
    cbar.ax.tick_params(labelsize=relative_font_size * font_axes_scale)  # modify font size of cbar values
    # cbar.set_label(label="Δ Depth in m", fontsize=relative_font_size*font_label_scale)
    cbar.ax.text(-relative_font_size * 0.2, 0.5, "Δ Depth",
                 fontsize=relative_font_size * font_label_scale,
                 rotation=90, ha="center", va="center", transform=cbar.ax.transAxes)
    plt.grid()
    plt.text(-170, -85, "GEBCO_2020 Grid, Ogg 2012", fontsize=relative_font_size * font_label_scale)
    plt.show()


# functions.plot_scatter(d_grid, b_grid_shift_ocean)
def depth_anomaly_scatter(size):
    mask_positive = data.grid2 > 0  # where bathymetry is > 0 -> land

    optimum_line = np.linspace(-8500, -2500, 2)
    x_depth_max = -8500
    x_depth_min = data.avr_ridge_depth
    x_depth_max_tick = -8000
    x_depth_min_tick = -3000
    num_of_x_ticks = 6
    ax1_ticks = np.linspace(x_depth_min_tick, x_depth_max_tick, num_of_x_ticks)
    age_values = []  # stores age values for x_age-axis
    for i in range(num_of_x_ticks):
        age = functions.depth_to_age(ax1_ticks[i], data.avr_ridge_depth)
        age_values.append(age)

    ratio = size * 12 / 9
    fig = plt.figure(figsize=(ratio, size))
    relative_font_size = fig.get_size_inches()[0] * 0.6  # gets font size relative to figure size
    ax1 = fig.add_subplot()

    hb = ax1.hexbin(data.grid1, data.grid2, gridsize=1000, cmap=cmap1, bins='log', mincnt=1)
    plt.axhline(y=0, alpha=1, color='orange', linestyle='--', linewidth=1)
    plt.plot(optimum_line, optimum_line, color="red", linewidth=1.5)

    ax1.set_xlim(x_depth_min, x_depth_max)
    ax1.set_xticks(ax1_ticks)  # sets x_depth-axis labels
    plt.title("Bathymetric over calculated depth", fontsize=relative_font_size*font_header_scale,
              pad=int(relative_font_size))
    plt.text(0.98, 0.98, "GEBCO_2020 Grid, Ogg 2012",
             fontsize=relative_font_size*font_text_scale,
             ha="right", va="top", transform=ax1.transAxes)
    ax1.set_xlabel("Calculated depth (in km)", fontsize=relative_font_size*font_label_scale)
    ax1.set_ylabel("Bathymetric depth (in km)", fontsize=relative_font_size*font_label_scale)

    # Age axis
    ax2 = ax1.twiny()
    ax2.set_xlim(ax1.get_xlim())  # so that ax1 and ax2 ticks match positions
    ax2.set_xticks(ax1_ticks)
    ax2.set_xticklabels([f"{age:.1f}" for age in age_values])
    plt.tick_params(labelsize=relative_font_size * font_axes_scale)
    ax1.tick_params(labelsize=relative_font_size * font_axes_scale)
    ax2.set_xlabel("Age (in million years)", fontsize=relative_font_size*font_label_scale)
    '''
    # Inset colorbar
    cax = inset_axes(ax1, width="30%", height="4%", loc="lower right", borderpad=2)
    cb = plt.colorbar(hb, cax=cax, orientation='horizontal')
    cb.set_label("log10(Density)", fontsize=10)
    cb.ax.xaxis.set_label_position('top')
    cb.ax.tick_params(labelsize=8)
    '''

    cb = fig.colorbar(hb, ax=ax1)
    cb.set_label("Density of data points", fontsize=relative_font_size*font_label_scale)
    cb.ax.tick_params(labelsize=relative_font_size*font_axes_scale)

    # legend
    legend_elements = [
        Line2D([0], [0], color='orange', linestyle='--', linewidth=1, label="sea surface"),
        Line2D([0], [0], color='red', linewidth=1.5, label="calc. depth = bath. depth")
    ]
    ax1.legend(
        handles=legend_elements,
        loc="lower left",
        fontsize=relative_font_size*font_text_scale,
        frameon=True
    )
    plt.grid()
    plt.show()


# depth anomaly histogram-----------------------------------------------------------------------------------------------
def depth_anomaly_hist(size, bins):
    if data.accuracy == 2:  # how many data points do we have
        y_factor = 10
    else:
        y_factor = 1
    y_max = int(20000000 / bins) * y_factor
    ratio = size * 3/3
    fig3, ax = plt.subplots(figsize=(ratio, size))
    relative_font_size = fig3.get_size_inches()[0] * 0.9  # gets font size relative to figure size
    plt.axvline(x=10, alpha=0.4, color='red', linestyle='--', linewidth=2)
    plt.axvline(x=-10, alpha=0.4, color='blue', linestyle='--', linewidth=2)
    plt.hist(data.grid_point_diffs, bins=bins, alpha=1, range=(-2000, 6000), color=(0.9, 0.5, 0.1))
    plt.tick_params(labelsize=relative_font_size * font_axes_scale)
    x_ticks = np.linspace(-2000, 6000, num=5)  # Choose appropriate tick positions
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([f"{int(tick * 0.001)} km" for tick in x_ticks])  # Convert to km and format
    y_ticks = np.linspace(y_max / 4, y_max, num=4)  # Choose appropriate tick positions
    ax.set_yticks(y_ticks)
    ax.set_yticklabels([f"{int(tick * 0.001)}k" for tick in y_ticks])  # Convert to km and format
    plt.title("Quantitative distribution \nof depth anomaly", fontsize=relative_font_size * font_header_scale, pad=10)
    plt.text(0.98, 0.98, "GEBCO_2020 Grid, Ogg 2012",
             fontsize=relative_font_size * font_text_scale,
             ha="right", va="top", transform=ax.transAxes)
    plt.xlabel("Calculated depth - Bathymetric depth", fontsize=relative_font_size * font_label_scale)
    ax.text(0.02, 0.6, "Deeper than\ncalculated",
            fontsize=relative_font_size * font_text_scale, color='b',
            ha="left", va="bottom", transform=ax.transAxes)
    ax.text(0.52, 0.6, "Shallower than\ncalculated",
            fontsize=relative_font_size * font_text_scale, color='r',
            ha="left", va="bottom", transform=ax.transAxes)
    plt.ylabel("Amount of data points", fontsize=relative_font_size * font_label_scale)
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