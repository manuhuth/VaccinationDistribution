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

y_position = 1.08
y_size = 45

label1 = 0.8
label2 = 0.6

names = [
    f"R_value_2_countryA_{label1}",
    f"R_value_2_countryA_{label2}",
    f"R_value_2_countryB_{label1}",
    f"R_value_2_countryB_{label2}",
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


xlims = [1*10**4, 0.5*10**5, 0.7*10**5, 0.75*10**5]
ylims = [1.1*10**5, 1.2*10**5, 1.7*10**4, 1*10**5]
xlims = [   0.5*10**5, 1*10**4,  0.75*10**5,  0.7*10**5]                          
ylims = [   1.2*10**5, 1.1*10**5, 1*10**5,      1.7*10**4]                        

# start plot
fig, axs = plt.subplots(5, 4, figsize=(65, 46), gridspec_kw={'height_ratios': [0.3, 0.8, 0.8, 0.5, 0.5]})
count_plot = 97


i = 0
alphas = [0.8, 0.8, 0.4, 0.4]
color_npi = "C2"
xlabel = "Weeks"
ylabel = "Degree of NPIs"
labels = [label1, label2]
countries = ["country A", "country A", "country B", "country B"]

for i in range(len(dicts)):
    ax = axs[0][i]
    if i == 0:
        ylab = ylabel
    else:
        ylab = ""
    
    j = i % 2
    ax.plot([0, 200/7], np.repeat(1 - labels[j], 2), color = color_npi)
    ax.fill_between([0, 200/7], np.repeat(1 - labels[j], 2), color = color_npi, alpha=alphas[i])
    ax.set_ylim([0, 0.6])
    ax.set_xlabel(xlabel)
    ax.text(200/7 / 2,  (1 - labels[j]) / 2, f"Degree of NPIs in {countries[i]}",
            horizontalalignment="center",
            verticalalignment="center",)
    ax.text(
        -0.05,
        1.17,
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
            "Degree of NPIs",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            weight="bold",
            size=35,
            rotation="vertical",
        )
    if i % 2 == 0:
        ax.text(
                1.2,
                1.4,
                f"NPIs only in {countries[i]}",
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
    ax.set_xlim([0, xlims[i]])
    ax.set_ylim([0, ylims[i]])
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
        ax.annotate('Pareto front', xy=(3*10**4, 0.32*10**5), xytext=(3.4*10**4, 0.42*10**5),
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
        x_scatter=np.linspace(0, 140, 9) / 7,
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

for i in range(0,len(dicts)):
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


fig.subplots_adjust(hspace=0.35)
fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/scenario4",
    bbox_inches="tight",
)
