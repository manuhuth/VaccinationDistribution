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
y_size = 16
count_plot = 65


add_additional = {"integers" : [9,10],
                 "number" : 0.5}


fig = plt.figure(constrained_layout=True, figsize = (23,16))
gs = GridSpec(2, 3, figure=fig)




   
dec_A = compute_deceased(results, "countryA")
dec_B = compute_deceased(results, "countryB")
colors = ["C1", "C0", "black"]
labels = ["Pareto optimal", "Optimal", "Population\nbased"]
linestyles = ["solid", "solid", "dashed"]
keys = list(dec_B.keys())

ax = fig.add_subplot(gs[0, 0])  
for index in range(len(dec_B.keys())):
    if keys[index] in ["optimal", "pop"]:
        ax.scatter(initials, (np.array(dec_A[keys[index]]) + np.array(dec_B[keys[index]]))/100000,
                label = labels[index], color=colors[index],
                linestyle=linestyles[index], s=40)

pop = (np.array(dec_A["pop"]) + np.array(dec_B["pop"]))/100000
optimal = (np.array(dec_A["optimal"]) + np.array(dec_B["optimal"]))/100000 
pareto = (np.array(dec_A["pareto"]) + np.array(dec_B["pareto"]))/100000 
for index in range(len(initials)):
    ax.plot(np.repeat(initials[index],2), [pop[index], optimal[index]], color = "firebrick" )
    diff = np.round(- pop[index] + optimal[index],2)
    ax.annotate(f"{diff}", (initials[index]+0.1, optimal[index] -diff/2 ) )
    
    
ax.set_title("Total deaths")
ax.set_xlabel("Initial mutant (wild-type) infected\nat $t=0$ in country A (B)")
ax.set_ylabel("Deaths in 100,000")
ylim = ax.get_ylim()
ax.legend()


ax = fig.add_subplot(gs[0, 1])
for index in range(len(dec_B.keys())):
    if keys[index] in ["pareto", "pop"]:
        ax.scatter(initials, (np.array(dec_A[keys[index]]) + np.array(dec_B[keys[index]]))/100000,
                label = labels[index], color=colors[index],
                linestyle=linestyles[index], s=40)
for index in range(len(initials)):
    ax.plot(np.repeat(initials[index],2), [pop[index], pareto[index]], color = "firebrick" )
    diff = np.round(- pop[index] + pareto[index],2)
    ax.annotate(f"{diff}", (initials[index]+0.1, pareto[index] -diff/2 ) )
    
    
ax.set_title("Total deaths")
ax.set_xlabel("Initial mutant (wild-type) infected\nat $t=0$ in country A (B)")
ax.set_ylabel("Deaths in 100,000")
ax.set_ylim(ylim)
ax.legend()






ax = fig.add_subplot(gs[0, 1])  
for index in range(len(dec_A.keys())):
    keys = list(dec_A.keys())
    ax.scatter(initials, np.array(dec_A[keys[index]])/100000, label = labels[index], color=colors[index],
            linestyle=linestyles[index])
ax.set_title("Deaths in country A")
ax.set_xlabel("Initial mutant (wild-type) infected\nat $t=0$ in country A (B)")
ax.set_ylabel("Deaths in 100,000")
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

ax = fig.add_subplot(gs[0, 2])  
for index in range(len(dec_B.keys())):
    keys = list(dec_B.keys())
    ax.scatter(initials, np.array(dec_B[keys[index]])/100000, label = labels[index], color=colors[index],
            linestyle=linestyles[index])
ax.set_title("Deaths in country B")
ax.set_xlabel("Initial mutant (wild-type) infected\nat $t=0$ in country A (B)")
ax.set_ylabel("Deaths in 100,000")





strategy_identifier = ["all_strategies", "pareto_improvements"]
name_strategy = ["Optimal strategy", "Pareto optimal strategy"]
colors = ["C0", "C1"]
paper_value = 30000
ylim = (0.2, 0.95)

index = 1
ax = fig.add_subplot(gs[1, 0])
fractions = compute_splines_from_results(type_opti = "all_strategies",
                                      vac = "vac1",)
ax.plot(initials, fractions, label = "Optimal", color = "C0")
#ax.fill_between(n_vaccs/2, np.repeat(0, len(fractions)), fractions, color = colors[index], alpha=0.3)
ax.set_ylabel("Fraction assigned\nto country A")
ax.plot(initials, np.repeat(0.5, len(fractions)), color= "black", linestyle = "dashed", label = "Population\nbased")
ax.set_xlabel("Vaccines per week")
ax.axvline(paper_value, 0, 1, color = "firebrick", linestyle = "dotted", label = "Used in the\nmain part")
ax.set_title("Optimal strategy (Vaccine one)")
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