import pickle
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

from functions.plot_tools import get_spline
from functions.plot_tools import plot_best_strategy
from functions.plot_tools import stacked_bar
from functions.plot_tools import plot_bars_vac

font = {"family": "normal", "weight": "normal", "size": 30}

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
    "initalUnequal_vacUnequal",
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
    "Equal initial states \nEqual vaccines",
    "Unequal initial states \nEqual vaccines",
    "Equal initial states \nUnequal vaccines",
    "Unequal initial states \nUnequal vaccines",
]

y_position = 1.08
y_size = 45


# start plot
fig, axs = plt.subplots(6, 4, figsize=(65, 55), gridspec_kw={'height_ratios': [0.2,0.3, 0.8, 0.8, 0.5, 0.5]})
category_names = ["Virus wild type", "Virus mutant"]
for i in range(len(dicts)):
    ax = axs[0][i]
    if i == 0:
        ylab = True
        legend=True
    else:
        ylab = False
        legend=True
    if i in [0,2]:
        results = {
            'Country A': [5, 5],
            'Country B': [5, 5],
        }
    else:
        results = {
            'Country A': [10, 0],
            'Country B': [0, 10],
        }
    stacked_bar(results, category_names, ax=ax, ylabel=ylab, legend=legend, map_col = "Set1")
    
    if i == 0:
        ax.text(
            -0.25,
            0.5,
            "Initially \ninfected",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            weight="bold",
            size=35,
            rotation="vertical",
        )
    ax.text(
            0.5,
            1.5,
            Titles[i],
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            weight="bold",
            size=35,
        )

for i in range(len(dicts)):
    ax = axs[1][i]
   
    if i in [1, 3]:
        vaccine1 = [
                0.95,
                0.6,
                0.9,
                0.9,
            ]
        vaccine2 = [
                0.6,
                0.95,
                0.9,
                0.9,
            ]
    else:
        vaccine1 = [
                0.95,
                0.6,
                0.9,
                0.9,
            ]
        vaccine2 = [
                0.95,
                0.6,
                0.9,
                0.9,
            ]
    if i ==0:
        title = "Reduction in %"
    else:
        title = ""
    plot_bars_vac(ax, vaccine1, vaccine2, barWidth=0.25, title = title)
    if i == 0:
        ax.text(
            -0.25,
            0.5,
            "Vaccine \nproperties",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            weight="bold",
            size=35,
            rotation="vertical",
        )

count_plot = 97
for i in range(len(dicts)):
    if i == 0:
        ylab = "Difference in %"
    else:
        ylab = ""
    ax = plot_bars_deaths(
        dict_output=dicts,
        ax=axs[2][i],
        case=names[i],
        total=False,
        unit=10 ** 6,
        title="",
        xlabel="Areas",
        ylabel=ylab,
        label_optimal="Optimal",
        label_Pareto="Pareto",
        label_population="Population",
        ylim=None,
    )
    ax.text(
        -0.05,
        y_position,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=y_size,
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
        ax.legend()
    count_plot += 1

# handles1, labels1 = ax.get_legend_handles_labels()
# fig.legend(handles1, labels1, loc="lower center", ncol=2, bbox_to_anchor=[0.5, 0.695])


for i in range(len(dicts)):
    if i == 0:
        ylab = "Deaths Country B"
    else:
        ylab = ""
    ax = plot_pareto_front(
        dict_output=dicts,
        fig=fig,
        ax=axs[3][i],
        case=names[i],
        size_points=90,
        size_optimal=22,
        title="",
        alpha=0.3,
        color_pareto="grey",
        color_min="C0",
        color_pareto_improvement="C1",
        linewidth_pareto=4.8,
        xlabel="Deaths Country A",
        ylabel=ylab,
    )
    ax.text(
        -0.05,
        y_position,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=y_size,
    )

    if i ==0:
        ax.text(
            -0.25,
            0.5,
            "Pareto front and \n Pareto set",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            weight="bold",
            size=35,
            rotation="vertical")
           
        red_patch = mpatches.Patch(color="C0", label="Minimum", linestyle="dashed")
        green_patch = mpatches.Patch(color="C1", label="Pareto minimum")
        blue_patch = mpatches.Patch(color="grey", label="Set of Pareto \nimprovements")
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles=[red_patch, green_patch, blue_patch])
    count_plot += 1


for i in range(len(dicts)):
    if i == 0:
        ylab = "% of vaccine in Country A"
    else:
        ylab = ""
    ax = plot_best_strategy(
        dict_output=dicts,
        fig=fig,
        ax=axs[4][i],
        case=names[i],
        x_scatter=np.linspace(0, 140, 9) / 7,
        vac_interest="vac1",
        periods=8,
        length=14,
        total_length=140,
        grid_points=6000,
        col_unconstrained="C0",
        label_unconstrained="Optimal",
        col_pareto="C1",
        label_pareto="Pareto optimal",
        xlabel="Weeks",
        ylabel=ylab,
        title="",
        linewidth=4,
        s_scatter=0,
        label_scatter="",
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
            "Allocation of \nVaccine one",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            weight="bold",
            size=35,
            rotation="vertical",
        )
        ax.text(
            -0.25,
            0.5,
            "Allocation of \nVaccine one",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            weight="bold",
            size=35,
            rotation="vertical",
        )
        ax.legend()
    count_plot += 1

for i in range(2,len(dicts)):
    if i == 2:
        ylab = "% of vaccine in Country A"
    else:
        ylab = ""
    ax = plot_best_strategy(
        dict_output=dicts,
        fig=fig,
        ax=axs[5][i],
        case=names[i],
        x_scatter=np.linspace(0, 140, 9) / 7,
        vac_interest="vac2",
        periods=8,
        length=14,
        total_length=140,
        grid_points=6000,
        col_unconstrained="C0",
        label_unconstrained="Optimal",
        col_pareto="C1",
        label_pareto="Pareto optimal",
        xlabel="weeks",
        ylabel=ylab,
        title="",
        linewidth=4,
        s_scatter=0,
        label_scatter="",
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
    if i == 2:
        ax.text(
            -0.25,
            0.5,
            "Allocation of \nVaccine two",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            weight="bold",
            size=35,
            rotation="vertical",
        )
        ax.legend()
    count_plot += 1


fig.subplots_adjust(hspace=0.4)

for i in range(2):
    axs[5][i].get_yaxis().set_visible(False)
    axs[5][i].get_xaxis().set_visible(False)
    axs[5][i].spines['left'].set_visible(False)
    axs[5][i].spines['bottom'].set_visible(False)
fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/scenario1",
    bbox_inches="tight",
)
