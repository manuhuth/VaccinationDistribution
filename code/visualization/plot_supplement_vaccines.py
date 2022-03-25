import pickle
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

from functions.plot_tools import get_spline
from functions.plot_tools import plot_best_strategy
from functions.plot_tools import stacked_bar, get_spline, plot_best_splines
from functions.plot_tools import plot_bars_vac, compute_incidences_countries, plot_incidences_country, plot_incidences, compute_incidences

from numpy import trapz

from functions.plot_tools import plot_horizontal_bars_annotated


font = {"family": "normal", "weight": "normal", "size": 10}

matplotlib.rc("font", **font)
import matplotlib.patches as mpatches

from functions.plot_tools import plot_pareto_front
from functions.plot_tools import plot_bars_deaths

# load outputs 2 countries, ordinary set-up
path = "/home/manuel/Documents/VaccinationDistribution/code/objects/"
names = [
    "initalEqual_vacEqual",
    "initalUnequal_vacEqual",
    "initalEqual_vacUnequal",
    "initalUnequal_vacUnequal_nvacc_60000",
]
dicts = {}

for i in names:
    path1 = path + i + ".pkl"
    with open(
        path1,
        "rb",
    ) as input:
        loaded_object = pickle.load(input)
        dicts[i] = loaded_object


y_position = 1.35
y_size = 28
x_position = -0.3
size = 19
font = {"family": "normal", "weight": "normal", "size": size}
matplotlib.rc("font", **font)
matplotlib.rcParams['xtick.labelsize'] = size - 2 
matplotlib.rcParams['ytick.labelsize'] = size - 2 

count_plot = 65
color_vaccines = "steelblue"

fig = plt.figure(constrained_layout=True, figsize = (16,10))
gs = GridSpec(2, 2, figure=fig)
gs.update(wspace=0.04, hspace=0.035)



data_all = dicts["initalUnequal_vacUnequal_nvacc_60000"]["all_strategies"]
data_pareto = dicts["initalUnequal_vacUnequal_nvacc_60000"]["pareto_improvements"]

periods=10
length=14
grid_points=6000
n=20
doses = 30000
alpha = 0.6
color = "steelblue"
ylabel = "Doses allocated to Country 1"
xlabel = "Weeks"
add_periods = [0.5, 0.5]
title = "Optimal - Vaccine 1"






ax = fig.add_subplot(gs[0, 0])

plot_best_splines(ax, data_all, "vac1", periods, length, grid_points, n, doses,
                      alpha, color, ylabel, "", add_periods, "Optimal strategy - Vaccine 1")
ax.text(
        x_position,
        y_position,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=y_size,
)
count_plot += 1


ax = fig.add_subplot(gs[1, 0])
plot_best_splines(ax, data_all, "vac2", periods, length, grid_points, n, doses,
                      alpha, color, ylabel, xlabel, add_periods, "Optimal strategy - Vaccine 2")



ax = fig.add_subplot(gs[0, 1])

plot_best_splines(ax, data_pareto, "vac1", periods, length, grid_points, n, doses,
                      alpha, color, "", "", add_periods, "Pareto strategy - Vaccine 1")
ax.text(
        x_position,
        y_position,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=y_size,
)
count_plot += 1


ax = fig.add_subplot(gs[1, 1])
plot_best_splines(ax, data_pareto, "vac2", periods, length, grid_points, n, doses,
                      alpha, color, "", xlabel, add_periods, "Pareto strategy - Vaccine 2")



fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/plot_supplements_best_vaccines",
    bbox_inches="tight",
)



#waterfall plot 50 best results, pareto and optimal startegy, horizontal line fpr pop based value
n = 20
data_all = dicts["initalUnequal_vacUnequal_nvacc_60000"]["all_strategies"]
data_pareto = dicts["initalUnequal_vacUnequal_nvacc_60000"]["pareto_improvements"]


smallest_all = (data_all["fval"].nsmallest(n))
smallest_pareto = (data_pareto["fval"].nsmallest(n))

fig = plt.figure(constrained_layout=True, figsize = (16,10))
gs = GridSpec(1, 2, figure=fig)

ax = fig.add_subplot(gs[0, 0])
ax.scatter(np.linspace(0, len(smallest_all)-1, len(smallest_all)), smallest_all )

