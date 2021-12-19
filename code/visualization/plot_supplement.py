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



n_vaccs = np.linspace(10, 100, 19)*10**3 #125, 150, 175, 200, 225, 250

results = []
for index in n_vaccs:
    path = (
            f"/home/manuel/Documents/VaccinationDistribution/code/objects/initalUnequal_vacUnequal_nvacc_{int(index)}.pkl"
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
gs = GridSpec(3, 2, figure=fig)



def compute_deceased(results, country):
    pareto_deceased = []
    optimal_deceased = []
    pop_deceased = []
    for index in range(len(results)):
        result = results[index]
        pareto_result = result["trajectories_pareto"]
        optimal_result = result["trajectories_best"]
        population_result = result["trajectories_pop"]
        cols = [x for x in pareto_result.columns if "dead" in x and country in x]
        
        deceased_pareto = list(pareto_result[cols].sum(axis=1))
        deceased_optimal = list(optimal_result[cols].sum(axis=1))
        deceased_pop = list(population_result[cols].sum(axis=1))
        
        pareto_deceased.append(deceased_pareto[-1])
        optimal_deceased.append(deceased_optimal[-1])
        pop_deceased.append(deceased_pop[-1])
        
    out = {"pareto" : pareto_deceased,
           "optimal": optimal_deceased,
           "pop" : pop_deceased,}
    return out
def compute_splines_from_results(type_opti = "pareto_improvements",
                                 vac = "vac2",
                                 periods = 10,
                                 length = 14,
                                 total_length = 140,
                                 grid_points = 6000,
                                 add_additional = {"integers" : [9,10],
                                                     "number" : 0.5},):
    fractions = []
    for index in range(len(n_vaccs)):
        result = results[index]
        n_vacc = n_vaccs[index]
        df_optimal = result[type_opti]
        if not(add_additional is None):
            for index in add_additional["integers"]:
                for index2 in ["countryA"]:
                    for index3 in ["vac1", "vac2"]:
                        name = f"yy_{index2}_{index3}_{index}"
                        df_optimal[name] = add_additional["number"]
        
        cols = [x for x in df_optimal.columns if vac in x]
        argmin = np.argmin(df_optimal["fval"])
        
        yy_points = pd.Series(list(df_optimal.iloc[argmin][cols]))
        time = np.linspace(0, total_length, grid_points) / 7
        spline = get_spline(yy_points,
                            periods=periods,
                            length=length,
                            total_length=total_length,
                            grid_points=grid_points,
                        )
        area = trapz(spline*n_vacc/2, dx=(time[1] - time[0]))
        total_area = n_vacc*10
        fraction_country_A = area/total_area
        fractions.append(fraction_country_A)
        
    return fractions

ax = fig.add_subplot(gs[0, 0])    
dec_A = compute_deceased(results, "countryA")
dec_B = compute_deceased(results, "countryB")
colors = ["C1", "C0", "black"]
labels = ["Pareto optimal", "Optimal", "Population\nbased"]
linestyles = ["solid", "solid", "dashed"]

for index in range(len(dec_A.keys())):
    keys = list(dec_A.keys())
    ax.plot(n_vaccs/2, np.array(dec_A[keys[index]])/100000, label = labels[index], color=colors[index],
            linestyle=linestyles[index])
ax.set_title("Deaths in country A")
ax.set_xlabel("Vaccine doses per week")
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

ax = fig.add_subplot(gs[0, 1])  
for index in range(len(dec_B.keys())):
    keys = list(dec_B.keys())
    ax.plot(n_vaccs/2, np.array(dec_B[keys[index]])/100000, label = labels[index], color=colors[index],
            linestyle=linestyles[index])
ax.set_title("Deaths in country B")
ax.set_xlabel("Vaccine doses per week")
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
ax.plot(n_vaccs/2, fractions, label = "Optimal", color = "C0")
#ax.fill_between(n_vaccs/2, np.repeat(0, len(fractions)), fractions, color = colors[index], alpha=0.3)
ax.set_ylabel("Fraction assigned\nto country A")
ax.plot(n_vaccs/2, np.repeat(0.5, len(fractions)), color= "black", linestyle = "dashed", label = "Population\nbased")
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


ax = fig.add_subplot(gs[2, 0])
fractions = compute_splines_from_results(type_opti = "all_strategies",
                                      vac = "vac2",)
ax.plot(n_vaccs/2, fractions, label = "Optimal", color = "C0")
#ax.fill_between(n_vaccs/2, np.repeat(0, len(fractions)), fractions, color = colors[index], alpha=0.3)
ax.set_ylabel("Fraction assigned\nto country A")
ax.plot(n_vaccs/2, np.repeat(0.5, len(fractions)), color= "black", linestyle = "dashed", label = "Population\nbased")
ax.set_xlabel("Vaccines per week")
ax.axvline(paper_value, 0, 1, color = "firebrick", linestyle = "dotted", label = "Used in the\nmain part")
ax.set_title("Optimal strategy (Vaccine two)")
ax.set_ylim(ylim)

                

index = 1
ax = fig.add_subplot(gs[1, 1])
fractions = compute_splines_from_results(type_opti = "pareto_improvements",
                                      vac = "vac1",)
ax.plot(n_vaccs/2, fractions, label = "Pareto optimal", color = "C0")
#ax.fill_between(n_vaccs/2, np.repeat(0, len(fractions)), fractions, color = colors[index], alpha=0.3)
ax.set_ylabel("Fraction assigned\nto country A")
ax.plot(n_vaccs/2, np.repeat(0.5, len(fractions)), color= "black", linestyle = "dashed", label = "Population\nbased")
ax.set_xlabel("Vaccines per week")
ax.axvline(paper_value, 0, 1, color = "firebrick", linestyle = "dotted", label = "Used in the\nmain part")
ax.set_title("Pareto optimal strategy (Vaccine one)")
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
               
ax = fig.add_subplot(gs[2, 1])
fractions = compute_splines_from_results(type_opti = "pareto_improvements",
                                      vac = "vac2",)
ax.plot(n_vaccs/2, fractions, label = "Pareto optimal", color = "C0")
#ax.fill_between(n_vaccs/2, np.repeat(0, len(fractions)), fractions, color = colors[index], alpha=0.3)
ax.set_ylabel("Fraction assigned\nto country A")
ax.plot(n_vaccs/2, np.repeat(0.5, len(fractions)), color= "black", linestyle = "dashed", label = "Population\nbased")
ax.set_xlabel("Vaccines per week")
ax.axvline(paper_value, 0, 1, color = "firebrick", linestyle = "dotted", label = "Used in the\nmain part")
ax.set_title("Pareto optimal strategy (Vaccine two)")
ax.set_ylim(ylim)
plt.tight_layout()
fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/plot_supplements_nvaccines",
    bbox_inches="tight",
)