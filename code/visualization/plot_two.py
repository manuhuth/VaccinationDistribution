import pickle
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

from matplotlib.gridspec import GridSpec
import matplotlib.cm as cm
from functions.plot_tools import get_spline
from functions.plot_tools import plot_best_strategy
from functions.plot_tools import stacked_bar
from functions.plot_tools import plot_bars_vac, compute_deceased, compute_splines_from_results
from functions.plot_tools import compute_splines_from_results_initials

from numpy import trapz

font = {"family": "normal", "weight": "normal", "size": 18}

matplotlib.rc("font", **font)

initials = [0,1,2,3,4,5]#125, 150, 175, 200, 225, 250
vac = "Unequal" #only 
initial_states = "Unequal"  
specification = f"inital{initial_states}_vac{vac}"

results = []
for index in initials:
    path = (
        f"/home/manuel/Documents/VaccinationDistribution/code/objects/{specification}_initial_{index}.pkl"
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
y_position = 1.08
y_size = 29
count_plot = 65
lwidth = 2.5

add_additional = {"integers" : [9,10],
                 "number" : 0.5}


fig = plt.figure(constrained_layout=True, figsize = (30, 18))
gs = GridSpec(3, 3, figure=fig, height_ratios=[1.5,2,2])

#data = [[10,9,8,7,6,5,4,3,2,1,0],
#        [0,1,2,3,4,5,6,7,8,9,10],
#        [0,1,2,3,4,5,6,7,8,9,10],
#        [10,9,8,7,6,5,4,3,2,1,0],
#        ]

data = [[10,9,8,7,6,5],
        [0,1,2,3,4,5],
        [0,1,2,3,4,5],
        [10,9,8,7,6,5],
        ]

ax = fig.add_subplot(gs[0, 0:2]) 
cmap = cm.get_cmap("Greens", 14)
hex_col = []
for i in range(cmap.N):
    rgba = cmap(i)
    # rgb2hex accepts rgb or rgba
    hex_col.append(matplotlib.colors.rgb2hex(rgba))

color_array = pd.DataFrame(columns=range(6), index = range(4))
for index in range(np.array(data).shape[0]):
    for index_col in range(np.array(data).shape[1]):
        color_array.loc[index, index_col] = hex_col[np.array(data)[index, index_col]]

row_labels=["Country A\nwild-type", "Country A\nMutant",
            "Country B\nwild-type", "Country B\nMutant",]
col_labels = [f"Case {i}" for i in range(11)]
x_table = ax.table(cellText=data, rowLabels=row_labels, colLabels= col_labels, loc="center", bbox=[0, 0.13, 1, 0.8],
                   cellLoc = "center", cellColours=np.array(color_array),)
ax.set_title("Initial virus distributions")
ax.axis("off")
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
ax.set_ylim([0.2, 0.9])
count_plot += 1



   
dec_A = compute_deceased(results, "countryA")
dec_B = compute_deceased(results, "countryB")
colors = ["C1", "C0", "C0"]
labels = ["Pareto optimal", "Optimal", "Population\nbased"]
linestyles = ["solid", "solid", "dashed"]
keys = list(dec_B.keys())

ax = fig.add_subplot(gs[1, 0])  
for index in range(len(dec_B.keys())):
    if keys[index] in ["optimal", "pop"]:
        if keys[index] == "pop":
            marker = "o"
        else:
            marker = "X"
        ax.scatter(initials, (np.array(dec_A[keys[index]]) + np.array(dec_B[keys[index]]))/1000,
                label = labels[index], color=colors[index],
                linestyle=linestyles[index], s=190, marker=marker)

pop = (np.array(dec_A["pop"]) + np.array(dec_B["pop"]))/1000
optimal = (np.array(dec_A["optimal"]) + np.array(dec_B["optimal"]))/1000
pareto = (np.array(dec_A["pareto"]) + np.array(dec_B["pareto"]))/1000
for index in range(len(initials)):
    ax.plot(np.repeat(initials[index],2), [pop[index], optimal[index]], color = "C0",
            linewidth = lwidth, linestyle = "solid" )
    diff = np.round(- pop[index] + optimal[index],2)
    ax.annotate(f"{diff}", (initials[index]+0.1, optimal[index] -diff/2 ) )
    
    
ax.set_title("Optimal and population-size based strategy")

ax.set_ylabel("Deaths in 1,000")
ylim = ax.get_ylim()
ax.set_xlim(-0.5, 5.5)
col_label_sym = [f"Case {i}" for i in range(5)] + ["Case 5"]
ax.set_xticklabels([""] + col_label_sym)
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




ax = fig.add_subplot(gs[2, 0])
for index in range(len(dec_B.keys())):
    if keys[index] in ["pareto", "pop"]:
        if keys[index] == "pop":
            marker = "o"
        else:
            marker = "X"
        ax.scatter(initials, (np.array(dec_A[keys[index]]) + np.array(dec_B[keys[index]]))/1000,
                label = labels[index], color="C1",
                linestyle=linestyles[index], s=190, marker=marker)
for index in range(len(initials)):
    ax.plot(np.repeat(initials[index],2), [pop[index], pareto[index]], color = "C1", linewidth = lwidth )
    diff = np.round(- pop[index] + pareto[index],2)
    ax.annotate(f"{diff}", (initials[index]+0.1, pareto[index] -diff/2 ) )
    
    
ax.set_title("Improvement through Pareto optimal strategy")
ax.set_ylabel("")
ax.set_xlim(-0.5, 5.5)
ax.set_ylabel("Total deaths in 1,000")
ax.set_xticklabels([""] + col_label_sym )
ax.set_ylim(ylim)

ax.legend()



strategy_identifier = ["all_strategies", "pareto_improvements"]
name_strategy = ["Optimal strategy", "Pareto optimal strategy"]
colors = ["C0", "C1"]
paper_value = 30000
ylim = (0.2, 0.95)

index = 1
ax = fig.add_subplot(gs[1, 1])
fractions = compute_splines_from_results_initials(initials=initials, type_opti = "all_strategies",
                                      vac = "vac1", results=results)

ax.scatter(initials, fractions,
           label = "Optimal", color="seagreen",
           s=190, marker="X")
ax.scatter(initials, np.repeat(0.5, len(initials)),
           label = "Population\nbased", color="seagreen",
           s=190, marker="o")
for index in range(len(initials)):
    ax.plot(np.repeat(initials[index],2), [0.5, fractions[index]], color = "seagreen", linewidth = lwidth )
    diff = np.round(-0.5 + fractions[index],2)
    ax.annotate(f"{diff}", (initials[index]+0.1, 0.5  + diff/2 ) )
ax.set_title("Optimal - Difference in\nvaccine one allocation")
ax.set_ylim([0.27, 0.93])
ax.set_ylabel("Fraction allocated to country A")
ax.legend()
cols_labels_small = [f"Case {i}" for i in range(len(initials))]
ax.set_xticklabels([""] + cols_labels_small )
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


ax = fig.add_subplot(gs[1, 2])
fractions = compute_splines_from_results_initials(initials=initials, type_opti = "pareto_improvements",
                                      vac = "vac1", results=results)

ax.scatter(initials, fractions,
           label = "Pareto optimal", color="seagreen",
           s=190, marker="X")
ax.scatter(initials, np.repeat(0.5, len(initials)),
           label = "Population\nbased", color="seagreen",
           s=190, marker="o")
for index in range(len(initials)):
    ax.plot(np.repeat(initials[index],2), [0.5, fractions[index]], color = "seagreen", linewidth = lwidth )
    diff = np.round(-0.5 + fractions[index],2)
    ax.annotate(f"{diff}", (initials[index]+0.1, 0.5  + diff/2 ) )
ax.set_title("Pareto optimal - Difference in\nvaccine one allocation")
ax.set_ylim([0.27, 0.93])
ax.legend()
ax.set_xticklabels([""] + cols_labels_small )
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

ax = fig.add_subplot(gs[2, 1])
fractions = compute_splines_from_results_initials(initials=initials, type_opti = "all_strategies",
                                      vac = "vac2", results=results)

ax.scatter(initials, fractions,
           label = "Optimal", color="seagreen",
           s=190, marker="X")
ax.scatter(initials, np.repeat(0.5, len(initials)),
           label = "Population\nbased", color="seagreen",
           s=190, marker="o")
for index in range(len(initials)):
    ax.plot(np.repeat(initials[index],2), [0.5, fractions[index]], color = "seagreen", linewidth = lwidth )
    diff = np.round(-0.5 + fractions[index],2)
    ax.annotate(f"{diff}", (initials[index]+0.1, 0.5  + diff/2 ) )
ax.set_title("Optimal - Difference in\nvaccine two allocation")
ax.set_ylabel("Fraction allocated to country A")
ax.set_ylim([0.27, 0.93])
ax.set_xticklabels([""] + cols_labels_small )
ax.legend()


ax = fig.add_subplot(gs[2, 2])
fractions = compute_splines_from_results_initials(initials=initials, type_opti = "pareto_improvements",
                                      vac = "vac1", results=results)

ax.scatter(initials, fractions,
           label = "Pareto optimal", color="seagreen",
           s=190, marker="X")
ax.scatter(initials, np.repeat(0.5, len(initials)),
           label = "Population\nbased", color="seagreen",
           s=190, marker="o")
for index in range(len(initials)):
    ax.plot(np.repeat(initials[index],2), [0.5, fractions[index]], color = "seagreen", linewidth = lwidth )
    diff = np.round(-0.5 + fractions[index],2)
    ax.annotate(f"{diff}", (initials[index]+0.1, 0.5  + diff/2 ) )
ax.set_title("Pareto optimal - Difference in\nvaccine two allocation")
ax.set_ylim([0.27, 0.93])
ax.legend()
ax.set_xticklabels([""] + cols_labels_small )

fig.subplots_adjust(wspace=1.3)
fig.tight_layout(pad=2)


fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/plot_two",
    bbox_inches="tight",
)