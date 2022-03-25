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
from functions.plot_tools import stacked_bar
from functions.plot_tools import plot_bars_vac, compute_incidences_countries, plot_incidences_country, plot_incidences, compute_incidences

from numpy import trapz

from functions.plot_tools import plot_horizontal_bars_annotated


size = 19
font = {"family": "normal", "weight": "normal", "size": size}
matplotlib.rc("font", **font)
matplotlib.rcParams['xtick.labelsize'] = size - 2 
matplotlib.rcParams['ytick.labelsize'] = size - 2 

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
count_plot = 65
color_vaccines = "seagreen"
x_position = -0.3

fig = plt.figure(constrained_layout=True, figsize = (20,14))
gs = GridSpec(4, 3, figure=fig, height_ratios=[0.5, 1, 1, 1], width_ratios = [1,1,1.2])
gs.update(wspace=0.02, hspace=0.05)

#plot one 
ax = fig.add_subplot(gs[1, 0:2])
plot_horizontal_bars_annotated(ax, dict_use = dicts["initalUnequal_vacUnequal_nvacc_60000"],
                               scale = 1,
                               color_vline = "black",
                               linestyle_vline = "dashed",
                               title = "Number of deceased individuals", font_size = size - 3)


ax.text(
        x_position + 0.13,
        y_position,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=y_size,
)
count_plot += 1


#plot two

ax = fig.add_subplot(gs[0:2, 2:3])

plot_pareto_front(
        dict_output=dicts,
        fig=fig,
        ax=ax,
        case="initalUnequal_vacUnequal_nvacc_60000",
        size_points=20,
        size_optimal=12,
        title="",
        alpha=0.3,
        color_pareto="grey",
        color_min="steelblue",
        color_pareto_improvement="seagreen",
        linewidth_pareto=3,
        xlabel="Deceased in Country 1",
        ylabel="Deceased in Country 2",
    )
ax.annotate('Pareto front', xy=(1.5*10**5, 0.8*10**5), xytext=(1.6*10**5, 1.1*10**5),
                            arrowprops=dict(facecolor='black', shrink=0.1))
ax.get_yaxis().set_major_formatter(
matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
ax.get_xaxis().set_major_formatter(
matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
ax.text(
        x_position,
        y_position*0.85,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=y_size,
)
red_patch = mpatches.Patch(color="steelblue", label="Minimum", linestyle="dashed")
green_patch = mpatches.Patch(color="seagreen", label="Pareto minimum")
blue_patch = mpatches.Patch(color="grey", label="Set of Pareto \nimprovements")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=[red_patch, green_patch, blue_patch], loc = "upper right")
count_plot += 1



#plot three 
add_additional = {"integers" : [9,10],
                 "number" : 0.5}
ax = fig.add_subplot(gs[2, 0])
plot_best_strategy(
        dict_output=dicts,
        fig=fig,
        ax=ax,
        case="initalUnequal_vacUnequal_nvacc_60000",
        x_scatter=np.linspace(0, 140, 11) / 7,
        vac_interest="vac1",
        periods=10,
        length=14,
        total_length=140,
        grid_points=6000,
        col_unconstrained=color_vaccines,
        label_unconstrained="",
        col_pareto=color_vaccines,
        label_pareto="Pareto optimal",
        xlabel="",
        ylabel="Vaccine 1 allocated\nto Country 1",
        title="Optimal strategy",
        linewidth=4,
        s_scatter=0,
        label_scatter="",
        plot="optimal",
        n_vacc = 30000,
        x_total=17,
        y_total=25500,
        scale_total=10*3,
        add_additional=add_additional,
        font_size = size-4,
)
ax.yaxis.grid(alpha=0.6)
ax.set_ylim([0,30500])
ax.get_yaxis().set_major_formatter(
matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
ax.get_xaxis().set_major_formatter(
matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
ax.legend(loc="lower right")
ax.text(
        x_position-0.1,
        y_position,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=y_size,
)

#plot four
ax = fig.add_subplot(gs[2, 1])
plot_best_strategy(
        dict_output=dicts,
        fig=fig,
        ax=ax,
        case="initalUnequal_vacUnequal_nvacc_60000",
        x_scatter=np.linspace(0, 140, 11) / 7,
        vac_interest="vac1",
        periods=10,
        length=14,
        total_length=140,
        grid_points=6000,
        col_unconstrained=color_vaccines,
        label_unconstrained="Optimal",
        col_pareto=color_vaccines,
        label_pareto="Pareto optimal",
        xlabel="",
        ylabel="",
        title="Pareto strategy",
        linewidth=4,
        s_scatter=0,
        label_scatter="",
        plot="pareto",
        n_vacc = 30000,
        x_total=17,
        y_total=25500,
        scale_total=10*3,
        add_additional=add_additional,
        font_size = size-4,
)
ax.set_ylim([0,30500])
ax.get_yaxis().set_major_formatter(
matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
ax.get_xaxis().set_major_formatter(
matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
count_plot += 1
ax.yaxis.grid(alpha=0.6)

#plot six
ax = fig.add_subplot(gs[3, 0])
plot_best_strategy(
        dict_output=dicts,
        fig=fig,
        ax=ax,
        case="initalUnequal_vacUnequal_nvacc_60000",
        x_scatter=np.linspace(0, 112, 11) / 7,
        vac_interest="vac2",
        periods=10,
        length=14,
        total_length=140,
        grid_points=6000,
        col_unconstrained=color_vaccines,
        label_unconstrained="",
        col_pareto=color_vaccines,
        label_pareto="Pareto optimal",
        xlabel="Weeks",
        ylabel="Vaccine 2 allocated\nto Country 1",
        title="",
        linewidth=4,
        s_scatter=0,
        label_scatter="",
        plot="optimal",
        n_vacc = 30000,
        x_total=6.5,
        y_total=7000,
        scale_total=10*3,
        add_additional=add_additional,
        font_size = size-4,
)
ax.set_ylim([0,30500])
#ax.legend()
ax.yaxis.grid(alpha=0.6)
ax.get_yaxis().set_major_formatter(
matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
ax.get_xaxis().set_major_formatter(
matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))


#plot seven
ax = fig.add_subplot(gs[3, 1])
plot_best_strategy(
        dict_output=dicts,
        fig=fig,
        ax=ax,
        case="initalUnequal_vacUnequal_nvacc_60000",
        x_scatter=np.linspace(0, 140, 11) / 7,
        vac_interest="vac2",
        periods=10,
        length=14,
        total_length=140,
        grid_points=6000,
        col_unconstrained=color_vaccines,
        label_unconstrained="",
        col_pareto=color_vaccines,
        label_pareto="Pareto optimal",
        xlabel="Weeks",
        ylabel="",
        title="",
        linewidth=4,
        s_scatter=0,
        label_scatter="",
        plot="pareto",
        n_vacc = 30000,
        x_total=6.5,
        y_total=7000,
        scale_total=1,
        add_additional=add_additional,
        font_size = size-4,
)
ax.set_ylim([0,30500])
ax.yaxis.grid(alpha=0.6)
ax.get_yaxis().set_major_formatter(
matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
ax.get_xaxis().set_major_formatter(
matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))



#plot incidences
time =  np.linspace(0, 140, 6000) 
trajectories = ["trajectories_best",
                "trajectories_pareto",
                "trajectories_pop",
                ]

incidences = compute_incidences_countries(dicts, time,  trajectories,
                                          name = "initalUnequal_vacUnequal_nvacc_60000",
                                 viruses=["virus1", "virus2"],
                                 countries = ["countryA", "countryB"],
                                 lambda1 = 0.1, 
                                 habitant_scale = 0.01)

#plot six first due to ylims
ax = fig.add_subplot(gs[3, 2])
plot_incidences_country(ax=ax,trajectories=trajectories, time=time,
                        incidences=incidences, viruses = ["virus1", "virus2"],
                            index_country = "countryB",
                            label_type = ["Optimal", "Paretol", "Pop.-based"],
                            colors = ["steelblue", "seagreen", "black"])

ax.set_title("")
ax.set_ylabel("Country 2")
ax.set_xlabel("Weeks")
#ax.legend()
ax.yaxis.grid(alpha=0.6)

ylim = ax.get_ylim()
ax.get_yaxis().set_major_formatter(
matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
ax.get_xaxis().set_major_formatter(
matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))



ax = fig.add_subplot(gs[2, 2])
plot_incidences_country(ax=ax,trajectories=trajectories, time=time,
                        incidences=incidences, viruses = ["virus1", "virus2"],
                            index_country = "countryA",
                            label_type = ["Optimal", "Pareto", "Pop.-based"],
                            colors = ["steelblue", "seagreen", "black"])
ax.yaxis.grid(alpha=0.6)
ax.set_title("7-day incidences")
ax.set_ylabel("Country 1")
ax.set_xlabel("")
ax.set_ylim(ylim)
ax.legend()
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
ax.get_yaxis().set_major_formatter(
matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
ax.get_xaxis().set_major_formatter(
matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))




fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/plot_one", bbox_inches='tight'
)

