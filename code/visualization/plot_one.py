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


y_position = 1.08
y_size = 12
count_plot = 97
color_vaccines = "seagreen"

fig = plt.figure(constrained_layout=True, figsize = (16,10))
gs = GridSpec(3, 3, figure=fig)


#plot one 
ax = fig.add_subplot(gs[0, 0:2])
plot_horizontal_bars_annotated(ax, dict_use = dicts["initalUnequal_vacUnequal_nvacc_60000"],
                               scale = 10**5,
                               color_vline = "black",
                               linestyle_vline = "dashed",
                               title = "Number of deaths per 100,000 inhabitants")
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
count_plot += 1


#plot two

ax = fig.add_subplot(gs[0, 2:3])

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
        color_min="C0",
        color_pareto_improvement="C1",
        linewidth_pareto=3,
        xlabel="Deaths Country A",
        ylabel="Deaths Country B",
    )
ax.annotate('Pareto front', xy=(1.5*10**5, 0.8*10**5), xytext=(1.6*10**5, 1.1*10**5),
                            arrowprops=dict(facecolor='black', shrink=0.1))
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
red_patch = mpatches.Patch(color="C0", label="Minimum", linestyle="dashed")
green_patch = mpatches.Patch(color="C1", label="Pareto minimum")
blue_patch = mpatches.Patch(color="grey", label="Set of Pareto \nimprovements")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=[red_patch, green_patch, blue_patch])
count_plot += 1



#plot three 
add_additional = {"integers" : [9,10],
                 "number" : 0.5}
ax = fig.add_subplot(gs[1, 0])
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
        xlabel="Weeks",
        ylabel="Doses of vaccine one allocated \n to country A in 1,000",
        title="Optimal strategy",
        linewidth=4,
        s_scatter=0,
        label_scatter="",
        plot="optimal",
        n_vacc = 30,
        x_total=3.6,
        y_total=7,
        scale_total=10*3,
        add_additional=add_additional,
)
ax.set_ylim([0,30.5])
ax.legend()
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

#plot four
ax = fig.add_subplot(gs[1, 1])
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
        xlabel="Weeks",
        ylabel="",
        title="Pareto optimal strategy",
        linewidth=4,
        s_scatter=0,
        label_scatter="",
        plot="pareto",
        n_vacc = 30,
        x_total=17,
        y_total=23,
        scale_total=10*3,
        add_additional=add_additional,
)
ax.set_ylim([0,30.5])
count_plot += 1


#plot six
ax = fig.add_subplot(gs[2, 0])
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
        ylabel="Doses of vaccine two allocated \n to country A in 10,000",
        title="Optimal strategy",
        linewidth=4,
        s_scatter=0,
        label_scatter="",
        plot="optimal",
        n_vacc = 30,
        x_total=5.5,
        y_total=7,
        scale_total=10*3,
        add_additional=add_additional,
)
ax.set_ylim([0,30.5])
#ax.legend()


#plot seven
ax = fig.add_subplot(gs[2, 1])
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
        title="Pareto optimal strategy",
        linewidth=4,
        s_scatter=0,
        label_scatter="",
        plot="pareto",
        n_vacc = 30,
        x_total=6.5,
        y_total=6,
        scale_total=10*3,
        add_additional=add_additional,
)
ax.set_ylim([0,30.5])





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
ax = fig.add_subplot(gs[2, 2])
plot_incidences_country(ax=ax,trajectories=trajectories, time=time,
                        incidences=incidences, viruses = ["virus1", "virus2"],
                            index_country = "countryB",
                            label_type = ["Optimal", "Pareto optimal", "Population\nbased"],
                            colors = ["C0", "C1", "black"])

ax.set_title("Country B")
ax.set_ylabel("7-day incidence per\n100,000 habtiants")
ax.set_xlabel("Weeks")
ax.legend()


ylim = ax.get_ylim()


ax = fig.add_subplot(gs[1, 2])
plot_incidences_country(ax=ax,trajectories=trajectories, time=time,
                        incidences=incidences, viruses = ["virus1", "virus2"],
                            index_country = "countryA",
                            label_type = ["Optimal", "Pareto optimal", "Population\nbased"],
                            colors = ["C0", "C1", "black"])

ax.set_title("Country A")
ax.set_ylabel("7-day incidence per\n100,000 habtiants")
ax.set_xlabel("Weeks")
ax.set_ylim(ylim)
ax.legend()
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
count_plot += 1




fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/plot_one",
    bbox_inches="tight",
)

