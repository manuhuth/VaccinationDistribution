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
    "reaction_2",
    "reaction_6",
    "reaction_8",
    "reaction_12",
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
    "Reaction from week 0",
    "Reaction from week 4",
    "Reaction from week 6",
    "Reaction from week 10",
]


# start plot
fig, axs = plt.subplots(5, 4, figsize=(65, 57), gridspec_kw={'height_ratios': [0.2, 0.8, 0.8, 0.5, 0.5]})

count_plot = 97
starts = [0,4,6,10]
length = 200

color = "seagreen"

for i in range(len(dicts)):
    texts = {"virus wild type" : [0, 200/7], 
         "Virus mutant" : [4, 200/7],
         "Optimize vaccinations" : [starts[i], 200/7]}

    ax = axs[0][i]
    for j in range(len(texts.keys())):
        key = list(texts.keys())[j]
        y_up = np.repeat(1 - 1 / len(texts.keys()) * j, 2)
        y_down = np.repeat(1 - 1 / len(texts.keys()) * (j+1),2)
        alpha = 1 - 1 / len(texts.keys()) * j
        ax.fill_between(
                texts[key], y_up, y_down, step="pre", alpha=alpha, color=color
            )
        ax.set_ylim([0,1])
        ax.set_xlabel("Weeks")
        ax.text(np.sum(texts[key])/2, (y_up+y_down)[0]/2, key, ha="center", va="center")
        ax.get_yaxis().set_visible(False)
        ax.spines['left'].set_visible(False)
    ax.text(
            0.5,
            1.3,
            Titles[i],
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            weight="bold",
            size=35,
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
    count_plot += 1



for i in range(len(dicts)):
    if i == 0:
        ylab = "Difference in %"
    else:
        ylab = ""
    ax = plot_bars_deaths(
        dict_output=dicts,
        ax=axs[1][i],
        case=names[i],
        total=False,
        unit=10 ** 6,
        title="",
        xlabel="Areas",
        ylabel=ylab,
        label_optimal="Optimal",
        label_Pareto="Pareto optimal",
        label_population="Population",
        ylim=[-65, 30],
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
        ax.legend()
    count_plot += 1


for i in range(len(dicts)):
    if i == 0:
        ylab = "Deaths Country B"
    else:
        ylab = ""
    ax = plot_pareto_front(
        dict_output=dicts,
        fig=fig,
        ax=axs[2][i],
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
        ax.annotate('Pareto front', xy=(0.75*10**5, 0.3*10**5), xytext=(0.8*10**5, 0.6*10**5),
                    arrowprops=dict(facecolor='black', shrink=0.1))
        
        
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
        ax=axs[3][i],
        case=names[i],
        x_scatter=np.linspace(0, 200, 15) / 7,
        vac_interest="vac1",
        periods=14,
        length=14,
        total_length=200,
        grid_points=1000,
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
        ax.legend()
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
        x_scatter=np.linspace(0, 200, 15) / 7,
        vac_interest="vac2",
        periods=14,
        length=14,
        total_length=200,
        grid_points=1000,
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
    if i == 0:
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


fig.subplots_adjust(hspace=0.3)
fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/scenario3",
    bbox_inches="tight",
)
