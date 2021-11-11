import pickle
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

from functions.plot_tools import get_spline
from functions.plot_tools import plot_best_strategy

font = {"family": "normal", "weight": "normal", "size": 30}

matplotlib.rc("font", **font)
import matplotlib.patches as mpatches

from functions.plot_tools import plot_pareto_front
from functions.plot_tools import plot_bars_deaths

# load outputs 2 countries, ordinary set-up
path = "/home/manuel/Documents/VaccinationDistribution/code/objects/"

names = [
    "R_value_8_countryA_0.6",
    "R_value_8_countryA_0.4",
    "R_value_8_countryB_0.6",
    "R_value_8_countryB_0.4",
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


Titles = ["Country A; 0.6", "Country A; 0.4", "Country B; 0.6", "Country B; 0.4"]


# start plot
fig, axs = plt.subplots(4, 4, figsize=(65, 60))
count_plot = 97
for i in range(len(dicts)):
    ax = plot_bars_deaths(
        dict_output=dicts,
        ax=axs[0][i],
        case=names[i],
        total=False,
        unit=10 ** 6,
        title=Titles[i],
        xlabel="Areas",
        ylabel="Difference in %",
        label_optimal="Optimal strategy",
        label_Pareto="Pareto optimal strategy",
        label_population="Population",
        ylim=[-75, 410],
    )
    ax.text(
        -0.05,
        1.03,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=35,
    )
    if i == 0:
        ax.text(
            -0.25,
            0.5,
            "Deaths compared to \n population based strategy",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            weight="bold",
            size=35,
            rotation="vertical",
        )
    count_plot += 1

handles1, labels1 = ax.get_legend_handles_labels()
fig.legend(handles1, labels1, loc="lower center", ncol=2, bbox_to_anchor=[0.5, 0.695])


for i in range(len(dicts)):
    ax = plot_pareto_front(
        dict_output=dicts,
        fig=fig,
        ax=axs[1][i],
        case=names[i],
        size_points=8.7,
        size_optimal=22,
        title=Titles[i],
        alpha=0.3,
        color_pareto="C0",
        linewidth_pareto=4.8,
        xlabel="Deaths Country A",
        ylabel="Deaths Country B",
    )
    ax.text(
        -0.05,
        1.03,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=35,
    )
    if i == 0:
        ax.text(
            -0.25,
            0.5,
            "Pareto front and \n Pareto set",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            weight="bold",
            size=35,
            rotation="vertical",
        )
    count_plot += 1

red_patch = mpatches.Patch(color="firebrick", label="Optimal point", linestyle="dashed")
green_patch = mpatches.Patch(color="seagreen", label="Pareto Optimal point")
blue_patch = mpatches.Patch(color="C0", label="Pareto set")
handles, labels = ax.get_legend_handles_labels()
fig.legend(
    handles=[red_patch, green_patch, blue_patch],
    loc="lower center",
    ncol=3,
    bbox_to_anchor=[0.5, 0.49],
)

for i in range(len(dicts)):
    ax = plot_best_strategy(
        dict_output=dicts,
        fig=fig,
        var_interest="unconstrained",
        ax=axs[2][i],
        case=names[i],
        title=Titles[i],
        periods=10,
        xlabel="Days",
        ylabel="% of vaccine one",
        col_vac1="C4",
        col_vac2="C5",
        linewidth=4,
        linewidth_vac=3,
    )
    ax.text(
        -0.05,
        1.03,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=35,
    )
    if i == 0:
        ax.text(
            -0.25,
            0.5,
            "% of vaccines assigned to \n country A & % of population being \n immune (optimal)",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            weight="bold",
            size=35,
            rotation="vertical",
        )
    count_plot += 1

for i in range(len(dicts)):
    ax = plot_best_strategy(
        dict_output=dicts,
        fig=fig,
        var_interest="pareto",
        ax=axs[3][i],
        case=names[i],
        title=Titles[i],
        periods=10,
        length=14,
        total_length=140,
        grid_points=6000,
        xlabel="Days",
        ylabel="% of vaccine two",
        col_vac1="C4",
        col_vac2="C5",
        linewidth=4,
        linewidth_vac=3,
    )
    ax.text(
        -0.05,
        1.03,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=35,
    )
    if i == 0:
        ax.text(
            -0.25,
            0.5,
            "% of vaccines assigned to \n country A & % of population being \n immune (Pareto optimal) ",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            weight="bold",
            size=35,
            rotation="vertical",
        )
    count_plot += 1
handles4, labels4 = ax.get_legend_handles_labels()
fig.legend(
    handles=handles4,
    labels=labels4,
    loc="lower center",
    ncol=4,
    bbox_to_anchor=[0.5, 0.29],
)

fig.subplots_adjust(hspace=0.4)
fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/scenario4",
    bbox_inches="tight",
)
