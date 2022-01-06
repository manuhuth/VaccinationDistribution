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

def autolabel(rects, rounding = 0):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        if rounding == 0:
            number = '{}'.format(int(np.round(height)))
        else:
            number = '{}'.format(np.round(height,rounding))
            if np.round(height,rounding) == int(np.round(height,rounding)):
                number = '{}'.format(int(np.round(height,rounding)))
        ax.annotate(number,
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

appearence = list(range(1, 12))
vac = "Unequal" #only 
initial_states = "Unequal"  
specification = f"inital{initial_states}_vac{vac}"

results = []
for index in appearence:
    specification = f"appears_week_{index}"
    path = (
        f"/home/manuel/Documents/VaccinationDistribution/code/objects/{specification}.pkl"
    )
    
    with open(
          path,
          "rb",
    ) as input:
          loaded_object = pickle.load(input)
    results.append(loaded_object)

appearence = [0] + appearence
    
deaths_countryA = pd.DataFrame(data=None, index=range(len(appearence)), columns=["optimal", "pareto", "population"])
deaths_countryB = deaths_countryA.copy()

deaths_countryA.iloc[0,:] = np.array([0.24, 0.94, 0.98]) * 10**5
deaths_countryB.iloc[0,:] = np.array([1.68, 1.14, 1.19]) * 10**5
for index in range(1, len(appearence)):
    r = results[index-1]
    index_optimal = np.argmin(r["all_strategies"]["fval"])
    index_pareto = np.argmin(r["pareto_improvements"]["fval"])
    
    
    A = [r["all_strategies"].iloc[index_optimal]["countryA"],
         r["pareto_improvements"].iloc[index_pareto]["countryA"],
         r["population_based"][0]]

    B = [r["all_strategies"].iloc[index_optimal]["countryB"],
         r["pareto_improvements"].iloc[index_pareto]["countryB"],
         r["population_based"][1]]
    deaths_countryA.iloc[index] = A
    deaths_countryB.iloc[index] = B

total_deaths = deaths_countryA + deaths_countryB

y_position = 1.08
y_size = 30
count_plot = 65


fig = plt.figure(constrained_layout=True, figsize = (75,17))
gs = GridSpec(3, 3, figure=fig)



xlabel = "Week of mutant appearence"
#------------------------------------------------------------------------------
ax = fig.add_subplot(gs[0, 2])
x = np.array(appearence)  # the label locations
width = 0.37  # the width of the bars
rects1 = ax.bar(x, total_deaths["population"]/10**3, width, label = "Population based", color = "C2", alpha = 0.7)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Deaths in 1,000')
ax.set_title('Deaths by population-size based vaccine allocation')
ax.set_xticks(x)
ax.set_xticklabels(appearence)
ax.legend()
ax.set_xlabel(xlabel)
autolabel(rects1, rounding=1)

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
#-------------------------------------------------------------------------------
ax = fig.add_subplot(gs[1, 2])

rects1 = ax.bar(x - width/2,(total_deaths["population"]- total_deaths["optimal"])/10**3, width, label='Optimal', color = "steelblue", alpha = 0.7)
rects2 = ax.bar(x + width/2, (total_deaths["population"]-total_deaths["pareto"])/10**3, width, label='Pareto optimal', color = "darkorange", alpha=0.7)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Saved lifes in 1,000')
ax.set_title('Number of saved lifes compared to population-size based vaccine allocation')
ax.set_xticks(x)
ax.set_xlabel(xlabel)

ax.set_xticklabels(appearence)
ax.legend()
autolabel(rects1, rounding=1)
autolabel(rects2, rounding=1)
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


#------------------------------------------------------------------------------
ax = fig.add_subplot(gs[2, 2])
optimal_percentage = -(total_deaths["optimal"] / total_deaths["population"] - 1)*100

pareto_percentage = -(total_deaths["pareto"] / total_deaths["population"] - 1)*100

rects1 = ax.bar(x - width/2, optimal_percentage, width, label='Optimal', color = "steelblue", alpha = 0.7)
rects2 = ax.bar(x + width/2,pareto_percentage, width, label='Pareto optimal', color = "darkorange", alpha=0.7)

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Relative number of saved lifes in %')
ax.set_title('Saved lifes relative to population-size based vaccine allocation')
ax.set_xticks(x)
ax.set_xticklabels(appearence)
ax.legend()
ax.set_ylim(0, 100)
ax.set_xlabel(xlabel)

autolabel(rects1)
autolabel(rects2)
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



fig.tight_layout(pad=3.5)

fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/plot_three",
    bbox_inches="tight",
)




