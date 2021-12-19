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
from functions.plot_tools import plot_bars_vac

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


y_position = 1.08
y_size = 12
count_plot = 97
color_vaccines = "C2"

fig = plt.figure(constrained_layout=True, figsize = (16,10))
gs = GridSpec(3, 3, figure=fig)


#plot one 
ax = fig.add_subplot(gs[0, 0:2])
plot_horizontal_bars_annotated(ax, dict_use = dicts["initalUnequal_vacUnequal"],
                               scale = 10**5,
                               color_vline = "black",
                               linestyle_vline = "dashed",
                               title = "Number of deaths per 100.000 inhabitants")

fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/plot_one_bars",
    bbox_inches="tight",
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
count_plot += 1


#plot two
fig = plt.figure(constrained_layout=True, figsize = (16,10))
gs = GridSpec(3, 3, figure=fig)
ax = fig.add_subplot(gs[0, 0:1])

plot_pareto_front(
        dict_output=dicts,
        fig=fig,
        ax=ax,
        case="initalUnequal_vacUnequal",
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


red_patch = mpatches.Patch(color="C0", label="Minimum", linestyle="dashed")
green_patch = mpatches.Patch(color="C1", label="Pareto minimum")
blue_patch = mpatches.Patch(color="grey", label="Set of Pareto \nimprovements")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles=[red_patch, green_patch, blue_patch])
count_plot += 1

fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/plot_one_pareto_front",
    bbox_inches="tight",
)



#plot three 
fig = plt.figure(constrained_layout=True, figsize = (16,10))
gs = GridSpec(3, 3, figure=fig)


add_additional = {"integers" : [9,10],
                 "number" : 0.5}
ax = fig.add_subplot(gs[0, 0])
plot_best_strategy(
        dict_output=dicts,
        fig=fig,
        ax=ax,
        case="initalUnequal_vacUnequal",
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
        ylabel="Doses of vaccine one allocated \n to country A in 10,000",
        title="Optimal strategy",
        linewidth=4,
        s_scatter=0,
        label_scatter="",
        plot="optimal",
        n_vacc = 0.6,
        x_total=3.7,
        y_total=0.15,
        scale_total=10*5,
        add_additional=add_additional,
)
ax.set_ylim([0,0.65])
ax.legend()

#plot four
ax = fig.add_subplot(gs[0, 1])
plot_best_strategy(
        dict_output=dicts,
        fig=fig,
        ax=ax,
        case="initalUnequal_vacUnequal",
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
        n_vacc = 0.6,
        x_total=17,
        y_total=0.5,
        scale_total=10*5,
        add_additional=add_additional,
)
ax.set_ylim([0,0.65])
count_plot += 1




#plot six
ax = fig.add_subplot(gs[1, 0])
plot_best_strategy(
        dict_output=dicts,
        fig=fig,
        ax=ax,
        case="initalUnequal_vacUnequal",
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
        n_vacc = 0.6,
        x_total=5.5,
        y_total=0.15,
        scale_total=10*5,
        add_additional=add_additional,
)
ax.set_ylim([0,0.65])
#ax.legend()



#plot seven
ax = fig.add_subplot(gs[1, 1])
plot_best_strategy(
        dict_output=dicts,
        fig=fig,
        ax=ax,
        case="initalUnequal_vacUnequal",
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
        n_vacc = 0.6,
        x_total=7.5,
        y_total=0.12,
        scale_total=10*5,
        add_additional=add_additional,
)
ax.set_ylim([0,0.65])


fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/plot_allocation",
    bbox_inches="tight",
)




#plot five
fig = plt.figure(constrained_layout=True, figsize = (16,10))
gs = GridSpec(3, 3, figure=fig)
ax = fig.add_subplot(gs[1, 2])
time = np.linspace(0, 140, 6000) / 7
ax.plot(time, dicts["initalUnequal_vacUnequal"]["optimal_frac_vaccinated"]["total"] * 100)
ax.set_xlabel("Weeks")
ax.set_ylabel("% of vaccinated or recovered")

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