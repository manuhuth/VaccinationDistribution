import pickle
import pandas as pd
import numpy as np

from models.vaccination.optimization_functions_estimagic import get_sum_of_states
from visualization.plot_gridspec import plot_gridspec

import matplotlib.pyplot as plt

plt.style.use("default")
# plt.style.use('ggplot')
plt.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica"],
        "grid.color": "white",
    }
)
plt.rcParams["lines.linewidth"] = 1.5
plt.rc("text", usetex=True)
plt.rcParams["text.latex.preamble"] = [r"\usepackage{amsmath}"]
plt.rcParams["axes.facecolor"] = "#E6E6E6"


def sum_states(df, names):
    for name in names:
        columns = [col for col in df.columns if all(x in col for x in names)]
    df_reduced = df[columns]
    result = df_reduced @ np.repeat(1, df_reduced.shape[1])

    return result


def get_input_dicts(names, results_type, deterministic):
    dict_optimal = {}
    for index in range(len(results_type)):
        df = results_type[index]
        df_sum = sum_states(df, names) / 10 ** 6

        key = f"key{index}"

        dict_optimal[key] = {
            "y": df_sum,
            "label": "Stochastic",
            "color": "steelblue",
            "linestyle": "-",
            "linewidth": 0.2,
        }

    if deterministic is not None:
        dict_optimal["z_deterministic"] = {
            "y": deterministic,
            "label": "Deterministic",
            "color": "firebrick",
            "linestyle": "--",
            "linewidth": 1.5,
        }
    return dict_optimal


function = "splines"
# model_type = "current"

plot_path = (
    "/home/manuel/Documents/VaccinationDistribution/paper/images/"
    + function
    + "_stochastic"
)

time = np.linspace(0, 140, 6000) / 7

types = ["current", "optimal", "unrestricted"]

results = {}

for index in types:
    path = (
        f"/home/manuel/Documents/VaccinationDistribution/code/objects/output_stochastic_{function}_{index}_rawSim"
        + ".pkl"
    )
    with open(
        path,
        "rb",
    ) as input:
        results[index] = pickle.load(input)

load = (
    "/home/manuel/Documents/VaccinationDistribution/code/objects/output_"
    + function
    + ".pkl"
)

with open(
    load,
    "rb",
) as input:
    dict_out = pickle.load(input)
df_optimal_results = dict_out["df_optimal_results"]
df_unrestricted_results = dict_out["df_unrestricted_results"]
model_optimal = dict_out["model_optimal"]
model_unrestricted = dict_out["model_seperated"]
model_current = dict_out["model_current"]
time_points = dict_out["time_points"]


names_A = ["infectious", "countryA"]
names_B = ["infectious", "countryB"]

until = 500
# Pareto Optimal
results_optimal = results["optimal"][0:until]
deterministic_optimal_A = sum_states(model_optimal["states"], names_A) / 10 ** 6
optimal_countryA = get_input_dicts(names_A, results_optimal, deterministic_optimal_A)

deterministic_optimal_B = sum_states(model_optimal["states"], names_B) / 10 ** 6
optimal_countryB = get_input_dicts(names_B, results_optimal, deterministic_optimal_B)

# unrestricted
results_unrestricted = results["unrestricted"][0:until]
deterministic_unrestricted_A = (
    sum_states(model_unrestricted["states"], names_A) / 10 ** 6
)
unrestricted_countryA = get_input_dicts(
    names_A, results_unrestricted, deterministic_unrestricted_A
)

deterministic_unrestricted_B = (
    sum_states(model_unrestricted["states"], names_B) / 10 ** 6
)
unrestricted_countryB = get_input_dicts(
    names_B, results_unrestricted, deterministic_unrestricted_B
)

# current
results_current = results["current"][0:until]
deterministic_current_A = sum_states(model_current["states"], names_A) / 10 ** 6
current_countryA = get_input_dicts(names_A, results_current, deterministic_current_A)

deterministic_current_B = sum_states(model_current["states"], names_B) / 10 ** 6
current_countryB = get_input_dicts(names_B, results_current, deterministic_current_B)


y = [
    optimal_countryA,
    unrestricted_countryA,
    optimal_countryB,
    unrestricted_countryB,
    current_countryA,
    current_countryA,
]

path_vaccine = plot_path + "_infectious"

plot_gridspec(
    y=y,
    time=time,
    title="Number of infectious individuals in million",
    legend_next_to_plot=True,
    legend_location="lower center",
    ylim=[0, 35],
    xlim=None,
    bbox=(0.515, -0.05),
    line_thickness=0.6,
    ylabel=["Infectious", "Infectious"],
    xlabel="Weeks",
    row_labels=["Country A", "Country B"],
    title_col1="Pareto strategy",
    title_col2="Unrestricted strategy",
    title_col3="Current strategy",
    xticks=np.linspace(0, 20, 6),
    path=path_vaccine,
)


# -------------------------------------------------------------------------------------------------------------------------------
def sum_over_all(results2, names, end_sum=True, list_input=True):

    output = []
    if list_input is True:
        iterate = results2
    else:
        iterate = results2.keys()
    for key in iterate:
        if list_input is True:
            df = key
        else:
            df = results2[key]
        values = sum_states(df, names)

        if end_sum is True:
            output.append(list(values)[-1])
        else:
            output.append(values)
    if end_sum is False:
        colnames = ["run" + str(i) for i in (range(len(output)))]
        output = pd.DataFrame(np.array(output).transpose(), columns=colnames)
    else:
        output = np.array(output)
    return output


names_A = ["dead", "countryA"]
names_B = ["dead", "countryB"]
mod_type = "optimal"
dead_countryA_optimal = (
    sum_over_all(results[mod_type], names_A, end_sum=True, list_input=True) / 10 ** 6
)
dead_countryB_optimal = (
    sum_over_all(results[mod_type], names_B, end_sum=True, list_input=True) / 10 ** 6
)

mod_type = "unrestricted"
dead_countryA_unrestricted = (
    sum_over_all(results[mod_type], names_A, end_sum=True, list_input=True) / 10 ** 6
)
dead_countryB_unrestricted = (
    sum_over_all(results[mod_type], names_B, end_sum=True, list_input=True) / 10 ** 6
)

mod_type = "current"
dead_countryA_current = (
    sum_over_all(results[mod_type], names_A, end_sum=True, list_input=True) / 10 ** 6
)
dead_countryB_current = (
    sum_over_all(results[mod_type], names_B, end_sum=True, list_input=True) / 10 ** 6
)


# plt.hist2d(dead_countryA_unrestricted, dead_countryB_unrestricted, bins=(70, 70))
# plt.colorbar()
# plt.show()
from matplotlib.ticker import NullFormatter, MaxNLocator
from numpy import linspace

plt.ion()


def ellipse(ra, rb, ang, x0, y0, Nb=100):
    xpos, ypos = x0, y0
    radm, radn = ra, rb
    an = ang
    co, si = np.cos(an), np.sin(an)
    the = linspace(0, 2 * np.pi, Nb)
    X = radm * np.cos(the) * co - si * radn * np.sin(the) + xpos
    Y = radm * np.cos(the) * si + co * radn * np.sin(the) + ypos
    return X, Y


def plot_2d_historgam(
    x,
    y,
    path,
    title="Number of deceased individuals in million",
    xlabel="Country B",
    ylabel="Country A",
    fontsize=20,
    ellipse=False,
):
    # Set up default x and y limits
    xlims = [min(x), max(x)]
    ylims = [min(y), max(y)]

    # Set up your x and y labels

    # Define the locations for the axes
    left, width = 0.12, 0.55
    bottom, height = 0.12, 0.55
    bottom_h = left_h = left + width + 0.02

    # Set up the geometry of the three plots
    rect_temperature = [left, bottom, width, height]  # dimensions of temp plot
    rect_histx = [left, bottom_h, width, 0.25]  # dimensions of x-histogram
    rect_histy = [left_h, bottom, 0.25, height]  # dimensions of y-histogram

    # Set up the size of the figure
    fig = plt.figure(1, figsize=(9.5, 9))
    fig.suptitle(title, fontsize=fontsize)

    # Make the three plots
    axTemperature = plt.axes(rect_temperature)  # temperature plot
    axHistx = plt.axes(rect_histx)  # x histogram
    axHisty = plt.axes(rect_histy)  # y histogram

    # Remove the inner axes numbers of the histograms
    nullfmt = NullFormatter()
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # Find the min/max of the data
    xmin = min(xlims)
    xmax = max(xlims)
    ymin = min(ylims)
    ymax = max(y)

    # Make the 'main' temperature plot
    # Define the number of bins
    nxbins = 50
    nybins = 50
    nbins = 100

    xbins = linspace(start=xmin, stop=xmax, num=nxbins)
    ybins = linspace(start=ymin, stop=ymax, num=nybins)
    xcenter = (xbins[0:-1] + xbins[1:]) / 2.0
    ycenter = (ybins[0:-1] + ybins[1:]) / 2.0
    aspectratio = 1.0 * (xmax - 0) / (1.0 * ymax - 0)

    H, xedges, yedges = np.histogram2d(y, x, bins=(ybins, xbins))
    X = xcenter
    Y = ycenter
    Z = H

    # Plot the temperature data
    cax = axTemperature.imshow(
        H,
        extent=[xmin, xmax, ymin, ymax],
        interpolation="none",
        origin="lower",
        aspect="auto",
    )

    # Plot the temperature plot contours
    contourcolor = "white"
    xcenter = np.mean(x)
    ycenter = np.mean(y)
    ra = np.std(x)
    rb = np.std(y)
    ang = 0

    if ellipse is True:
        X, Y = ellipse(ra, rb, ang, xcenter, ycenter)
        axTemperature.plot(X, Y, "k:", ms=1, linewidth=2.0, color=contourcolor)
        axTemperature.annotate(
            "$1\\sigma$",
            xy=(X[15], Y[15]),
            xycoords="data",
            xytext=(10, 10),
            textcoords="offset points",
            horizontalalignment="right",
            verticalalignment="bottom",
            fontsize=25,
            color=contourcolor,
        )

        X, Y = ellipse(2 * ra, 2 * rb, ang, xcenter, ycenter)
        axTemperature.plot(X, Y, "k:", color=contourcolor, ms=1, linewidth=2.0)
        axTemperature.annotate(
            "$2\\sigma$",
            xy=(X[15], Y[15]),
            xycoords="data",
            xytext=(10, 10),
            textcoords="offset points",
            horizontalalignment="right",
            verticalalignment="bottom",
            fontsize=25,
            color=contourcolor,
        )

    # X,Y=ellipse(3*ra,3*rb,ang,xcenter,ycenter)
    # axTemperature.plot(X,Y,"k:",color = contourcolor, ms=1,linewidth=2.0)
    # axTemperature.annotate('$3\\sigma$', xy=(X[15], Y[15]), xycoords='data',xytext=(10, 10),
    #                       textcoords='offset points',horizontalalignment='right',
    #                       verticalalignment='bottom',fontsize=25, color = contourcolor)

    # Plot the axes labels
    axTemperature.set_xlabel(xlabel, fontsize=25)
    axTemperature.set_ylabel(ylabel, fontsize=25)

    # Make the tickmarks pretty
    ticklabels = axTemperature.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(18)
        label.set_family("serif")

    ticklabels = axTemperature.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(18)
        label.set_family("serif")

    # Set up the plot limits
    axTemperature.set_xlim(xlims)
    axTemperature.set_ylim(ylims)

    # Set up the histogram bins
    xbins = np.arange(xmin, xmax, (xmax - xmin) / nbins)
    ybins = np.arange(ymin, ymax, (ymax - ymin) / nbins)

    # Plot the histograms
    axHistx.hist(x, bins=xbins, color="steelblue")
    axHisty.hist(y, bins=ybins, orientation="horizontal", color="seagreen")

    # Set up the histogram limits
    axHistx.set_xlim(min(x), max(x))
    axHisty.set_ylim(min(y), max(y))

    axHistx.spines["top"].set_alpha(0)
    # axHistx.spines["bottom"].set_alpha(0)
    axHistx.spines["right"].set_alpha(0)
    # axHistx.spines["left"].set_alpha(0)
    axHistx.grid(alpha=0.5)

    axHisty.spines["top"].set_alpha(0)
    # axHistx.spines["bottom"].set_alpha(0)
    axHisty.spines["right"].set_alpha(0)
    # axHistx.spines["left"].set_alpha(0)
    axHisty.grid(alpha=0.5)

    # Make the tickmarks pretty
    ticklabels = axHistx.get_yticklabels()
    for label in ticklabels:
        label.set_fontsize(12)
        label.set_family("serif")

    # Make the tickmarks pretty
    ticklabels = axHisty.get_xticklabels()
    for label in ticklabels:
        label.set_fontsize(12)
        label.set_family("serif")

    # Cool trick that changes the number of tickmarks for the histogram axes
    axHisty.xaxis.set_major_locator(MaxNLocator(4))
    axHistx.yaxis.set_major_locator(MaxNLocator(4))

    # Show the plot
    plt.draw()

    plt.savefig(path)

    return axTemperature, axHistx, axHisty


# unrestricted
y = dead_countryA_unrestricted
x = dead_countryB_unrestricted
total_unrestricted = y + x
path_unrestricted = plot_path + "_histogram_deceased_unrestricted"
axT_unrestricted, axx_unrestricted, axy_unrestricted = plot_2d_historgam(
    x, y, path_unrestricted
)

# optimal
y = dead_countryA_optimal
x = dead_countryB_optimal
total_optimal = y + x
path_optimal = plot_path + "_histogram_deceased_optimal"
faxT_optimal, axx_optimal, axy_optimal = plot_2d_historgam(x, y, path_optimal)

# current
y = dead_countryA_current
x = dead_countryB_current
total_current = y + x
path_current = plot_path + "_histogram_deceased_current"
faxT_current, axx_current, axy_current = plot_2d_historgam(x, y, path_current)


# ----------------------Histogram---------------------------------------------
names = ["dead"]


# Pareto Optimal
deterministic_optimal = sum_states(model_optimal["states"], names) / 10 ** 6
deterministic_unrestricted = sum_states(model_unrestricted["states"], names) / 10 ** 6
deterministic_current = sum_states(model_current["states"], names) / 10 ** 6

np.mean(total_optimal) - np.array(deterministic_optimal)[-1]
np.mean(total_unrestricted) - np.array(deterministic_unrestricted)[-1]
np.mean(total_current) - np.array(deterministic_current)[-1]


fig, ax = plt.subplots()
# bins = numpy.linspace(-10, 10, 100)
bbox = (0.515, -0.05)
legend_location = "lower center"
ax.hist(total_unrestricted, alpha=0.5, label="Unrestricted strategy")
ax.hist(total_optimal, alpha=0.5, label="Pareto optimal strategy")
ax.hist(total_current, alpha=0.5, label="Current strategy")
handles, labels = ax.get_legend_handles_labels()
legend = fig.legend(
    handles,
    labels,
    bbox_to_anchor=bbox,
    loc=legend_location,
    ncol=3,
    facecolor="white",
    framealpha=1,
)
frame = legend.get_frame()
frame.set_linewidth(0)
ax.set_xlabel("Total number of deceased individuals in million")
ax.set_ylabel("Absolute frequency")
ax.spines["top"].set_alpha(0)
# ax.spines["bottom"].set_alpha(0)
ax.spines["right"].set_alpha(0)
# ax.spines["left"].set_alpha(0)
ax.vlines(np.mean(total_unrestricted), ymin=0, ymax=250, color="C0", linestyle="dotted")
ax.vlines(np.mean(total_optimal), ymin=0, ymax=250, color="C1", linestyle="dotted")
ax.vlines(np.mean(total_current), ymin=0, ymax=250, color="C2", linestyle="dotted")
ax.grid(alpha=0.5)
# fig.suptitle("Deceased individuals (Splines)")
plt.tight_layout()
path_hist = plot_path + "_histogram.png"
plt.savefig(path_hist, bbox_inches="tight")
