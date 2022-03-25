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
from functions.plot_tools import plot_bars_vac, compute_deceased, compute_splines_from_results

from numpy import trapz



n_vaccs = list(np.linspace(0.00001, 0.01, 100)) + list(np.linspace(0.015, 0.167, 48)) + list(np.linspace(0.2, 1, 20))
n_vaccs = np.array(n_vaccs)
results = []
for index in n_vaccs:
    path = (
            f"/home/manuel/Documents/VaccinationDistribution/code/objects/distance_{index}.pkl"
            )
    
    with open(
          path,
          "rb",
    ) as input:
          loaded_object = pickle.load(input)
    results.append(loaded_object)


type_opti = "pareto_improvements"
vac = "vac2"
periods = 10
length = 14
total_length = 140
grid_points = 6000
y_position = 1.35
y_size = 28
x_position = -0.3
size = 19
font = {"family": "normal", "weight": "normal", "size": size}
matplotlib.rc("font", **font)
matplotlib.rcParams['xtick.labelsize'] = size - 2 
matplotlib.rcParams['ytick.labelsize'] = size - 2 
count_plot = 65


add_additional = {"integers" : [9,10],
                 "number" : 0.5}


fig = plt.figure(constrained_layout=True, figsize = (23,15))
gs = GridSpec(3, 2, figure=fig)




ax = fig.add_subplot(gs[0, 0])    
dec_A = compute_deceased(results, "countryA")
dec_B = compute_deceased(results, "countryB")
colors = ["C1", "steelblue", "black"]
labels = ["Pareto optimal", "Optimal", "Population-size\nbased"]
linestyles = ["solid", "solid", "dashed"]

for index in range(len(dec_A.keys())):
    keys = list(dec_A.keys())
    ax.plot(1-n_vaccs, np.array(dec_A[keys[index]]), label = labels[index], color=colors[index],
            linestyle=linestyles[index])
ax.set_title("Deaths in Country 1")
ax.set_xlabel("Distance")
ax.set_ylabel("Deaths")
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

ax = fig.add_subplot(gs[0, 1])  
for index in range(len(dec_B.keys())):
    keys = list(dec_B.keys())
    ax.plot(1-n_vaccs, np.array(dec_B[keys[index]]), label = labels[index], color=colors[index],
            linestyle=linestyles[index])
ax.set_title("Deaths in Country 2")
ax.set_xlabel("Distance")
ax.set_ylabel("Deaths")



strategy_identifier = ["all_strategies", "pareto_improvements"]
name_strategy = ["Optimal strategy", "Pareto optimal strategy"]
colors = ["steelblue", "C1"]
paper_value = 1-  1/1000
ylim = (0.2, 0.95)

index = 1
ax = fig.add_subplot(gs[1, 0])
fractions = compute_splines_from_results(type_opti = "all_strategies",
                                      vac = "vac1",results=results, n_vaccs=n_vaccs)
ax.scatter(1-n_vaccs, fractions, label = "Optimal", color = "steelblue")
#ax.fill_between(n_vaccs, np.repeat(0, len(fractions)), fractions, color = colors[index], alpha=0.3)
ax.set_ylabel("Share assigned to Country 1")
ax.plot(1-n_vaccs, np.repeat(0.5, len(fractions)), color= "black", linestyle = "dashed", label = "Population-size\nbased")
ax.set_xlabel("Vaccines per week")
ax.axvline(paper_value, 0, 1, color = "firebrick", linestyle = "dotted", label = "From main part")
ax.set_title("Optimal strategy - Vaccine 1")
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


ax = fig.add_subplot(gs[2, 0])
fractions = compute_splines_from_results(type_opti = "all_strategies",
                                      vac = "vac2",results=results, n_vaccs=n_vaccs)
ax.scatter(1-n_vaccs, fractions, label = "Optimal", color = "steelblue")
#ax.fill_between(n_vaccs, np.repeat(0, len(fractions)), fractions, color = colors[index], alpha=0.3)
ax.set_ylabel("Share assigned to Country 1")
ax.plot(1-n_vaccs, np.repeat(0.5, len(fractions)), color= "black", linestyle = "dashed", label = "Population-size\nbased")
ax.set_xlabel("Distance")
ax.axvline(paper_value, 0, 1, color = "firebrick", linestyle = "dotted", label = "Main analysis")
ax.set_title("Optimal strategy - Vaccine 2")
ax.set_ylim(ylim)

                

index = 1
ax = fig.add_subplot(gs[1, 1])
fractions = compute_splines_from_results(type_opti = "pareto_improvements",
                                      vac = "vac1",results=results, n_vaccs=n_vaccs)
ax.scatter(1-n_vaccs, fractions, label = "Pareto optimal", color = "steelblue")
#ax.fill_between(n_vaccs, np.repeat(0, len(fractions)), fractions, color = colors[index], alpha=0.3)
ax.set_ylabel("Share assigned to Country 1")
ax.plot(1-n_vaccs, np.repeat(0.5, len(fractions)), color= "black", linestyle = "dashed", label = "Population-size\nbased")
ax.set_xlabel("Distance")
ax.axvline(paper_value, 0, 1, color = "firebrick", linestyle = "dotted", label = "Main analysis")
ax.set_title("Pareto strategy - Vaccine 1")
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
         
ax = fig.add_subplot(gs[2, 1])
fractions = compute_splines_from_results(type_opti = "pareto_improvements",
                                      vac = "vac2",results=results, n_vaccs=n_vaccs)
ax.scatter(1-n_vaccs, fractions, label = "Pareto optimal", color = "steelblue")
#ax.fill_between(n_vaccs, np.repeat(0, len(fractions)), fractions, color = colors[index], alpha=0.3)
ax.set_ylabel("Share assigned to Country 1")
ax.plot(1-n_vaccs, np.repeat(0.5, len(fractions)), color= "black", linestyle = "dashed", label = "Population-size\nbased")
ax.set_xlabel("Distance")
ax.axvline(paper_value, 0, 1, color = "firebrick", linestyle = "dotted", label = "Main analysis")
ax.set_title("Pareto strategy - Vaccine 2")
ax.set_ylim(ylim)
plt.tight_layout()

fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/plot_supplements_distances",
    bbox_inches="tight",
)