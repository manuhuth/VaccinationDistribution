import pickle
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

from functions.plot_tools import get_spline
from functions.plot_tools import plot_best_strategy

font = {"family": "normal", "weight": "normal", "size": 18}

matplotlib.rc("font", **font)
import matplotlib.patches as mpatches

from functions.plot_tools import plot_pareto_front
from functions.plot_tools import plot_bars_deaths

# load outputs 2 countries, ordinary set-up
path = "/home/manuel/Documents/VaccinationDistribution/code/objects/"

names = [
    "reaction_4",
    "reaction_8",
    "reaction_12",
    "reaction_16",
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


Titles = [
    "Reaction in week 4",
    "Reaction in week 8",
    "Reaction in week 12",
    "Reaction in week 16",
]


# start plot
fig, axs = plt.subplots(1, 4, figsize=(40, 10))
count_plot = 97
for i in range(len(dicts)):
    ax = plot_bars_deaths(
        dict_output=dicts,
        ax=axs[i],
        case=names[i],
        total=False,
        unit=10 ** 6,
        title=Titles[i],
        xlabel="Areas",
        ylabel="Difference in %",
        label_optimal="Optimal strategy",
        label_Pareto="Pareto optimal strategy",
        label_population="Population",
        ylim=[-25, 25],
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
            "",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            weight="bold",
            size=35,
            rotation="vertical",
        )
    count_plot += 1

handles1, labels1 = ax.get_legend_handles_labels()
fig.legend(handles1, labels1, loc="lower center", ncol=2, bbox_to_anchor=[0.5, -0.06])
fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/presentation_bar",
    bbox_inches="tight",
)

fig, axs = plt.subplots()
count_plot = 97
for i in range(0, 1):
    ax = plot_pareto_front(
        dict_output=dicts,
        fig=fig,
        ax=axs,
        case=names[i],
        size_points=8.7,
        size_optimal=22,
        title="",
        alpha=0.3,
        color_pareto="C0",
        linewidth_pareto=1,
        xlabel="Deaths Country A",
        ylabel="Deaths Country B",
    )
    # ax.text(-0.05, 1.03,chr(count_plot), horizontalalignment='center',
    #        verticalalignment='center',
    #        transform = ax.transAxes, weight="bold", size= 35)
    if i == 0:
        ax.text(
            -0.25,
            0.5,
            "",
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
    bbox_to_anchor=[0.5, -0.23],
)
fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/presentation_pareto",
    bbox_inches="tight",
)

fig, axs = plt.subplots(1, 4, figsize=(40, 10))
count_plot = 97
for i in range(len(dicts)):
    ax = plot_best_strategy(
        dict_output=dicts,
        fig=fig,
        var_interest="unconstrained",
        ax=axs[i],
        case=names[i],
        title=Titles[i],
        periods=10,
        xlabel="Days",
        ylabel="% of vaccine allocated to country A",
        col_vac1="C4",
        col_vac2="C5",
        linewidth=4,
        linewidth_vac=0,
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
            "",
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
    handles=handles4[0:2],
    labels=labels4[0:2],
    loc="lower center",
    ncol=4,
    bbox_to_anchor=[0.5, -0.06],
)
fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/presentation_vaccine",
    bbox_inches="tight",
)

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
    "/home/manuel/Documents/VaccinationDistribution/paper/images/presentation",
    bbox_inches="tight",
)
