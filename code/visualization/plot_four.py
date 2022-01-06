import pickle
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns

from functions.plot_tools import get_spline
from functions.plot_tools import plot_best_strategy
from functions.plot_tools import stacked_bar
from functions.plot_tools import plot_bars_vac, compute_incidences_countries, plot_incidences_country, plot_incidences, compute_incidences
from functions.plot_tools import compute_splines_from_Rval_results

from numpy import trapz

from functions.plot_tools import plot_horizontal_bars_annotated


font = {"family": "normal", "weight": "normal", "size": 20}

matplotlib.rc("font", **font)
import matplotlib.patches as mpatches

from functions.plot_tools import plot_pareto_front
from functions.plot_tools import plot_bars_deaths


y_position = 1.08
y_size = 29
count_plot = 65

dict_R_cases = {}
case = 0
for index1 in range(3, 11):
    for index2 in range(3,11):
        dict_R_cases[f"case{case}"] = np.array([index1, index2])/10
        case += 1

# load outputs 2 countries
path = "/home/manuel/Documents/VaccinationDistribution/code/objects/"

dicts = {}

for i in range(len(dict_R_cases.keys())):
    path1 = path + f"degree_R_{i}" + ".pkl"
    with open(
        path1,
        "rb",
    ) as input:
        loaded_object = pickle.load(input)
        dicts[i] = loaded_object

vac_pareto_vac2 = compute_splines_from_Rval_results(dicts, type_opti = "pareto_improvements",
                                 vac = "vac2",
                                 periods = 10,
                                 length = 14,
                                 total_length = 140,
                                 grid_points = 6000,
                                 add_additional = {"integers" : [9,10],
                                                     "number" : 0.5},)
vac_pareto_vac1 = compute_splines_from_Rval_results(dicts, type_opti = "pareto_improvements",
                                 vac = "vac1",
                                 periods = 10,
                                 length = 14,
                                 total_length = 140,
                                 grid_points = 6000,
                                 add_additional = {"integers" : [9,10],
                                                     "number" : 0.5},)


vac_optimal_vac2 = compute_splines_from_Rval_results(dicts, type_opti = "all_strategies",
                                 vac = "vac2",
                                 periods = 10,
                                 length = 14,
                                 total_length = 140,
                                 grid_points = 6000,
                                 add_additional = {"integers" : [9,10],
                                                     "number" : 0.5},)
vac_optimal_vac1 = compute_splines_from_Rval_results(dicts, type_opti = "all_strategies",
                                 vac = "vac1",
                                 periods = 10,
                                 length = 14,
                                 total_length = 140,
                                 grid_points = 6000,
                                 add_additional = {"integers" : [9,10],
                                                     "number" : 0.5},)

df_R_pop = pd.DataFrame(np.nan, columns = np.round(1-np.linspace(0.3,1, 8),2), index =np.round(1-np.linspace(0.3,1, 8),2))
df_R_optimal = df_R_pop.copy()
df_R_pareto = df_R_pop.copy()

df_vac1_optimal = df_R_pop.copy()
df_vac1_pareto = df_R_pop.copy()

df_vac2_optimal = df_R_pop.copy()
df_vac2_pareto = df_R_pop.copy()

for index in range(len(dict_R_cases.keys())):
    key = f"case{index}"
    R_vals = dict_R_cases[key]
    
    optimal_df = dicts[index]["all_strategies"]
    index_optimal = np.argmin(optimal_df["fval"])
    deaths_optimal = optimal_df.loc[index_optimal,:]["fval"]
    
    pareto_df = dicts[index]["pareto_improvements"].reset_index(drop=True)
    index_pareto = np.argmin(pareto_df["fval"])
    deaths_pareto = pareto_df.loc[index_pareto,:]["fval"]
    
    deaths_pop = dicts[index]["population_based"]
    #death_pop_A = deaths_pop[0]
    #death_pop_B = deaths_pop[1]
    
    df_R_pop.loc[np.round(1-R_vals[0],2), np.round(1-R_vals[1],2)] = np.sum(deaths_pop)
    df_R_optimal.loc[np.round(1-R_vals[0],2), np.round(1-R_vals[1],2)] = deaths_optimal
    df_R_pareto.loc[np.round(1-R_vals[0],2), np.round(1-R_vals[1],2)] = deaths_pareto
    
    df_vac1_optimal.loc[np.round(1-R_vals[0],2), np.round(1-R_vals[1],2)] = vac_optimal_vac1[index]
    df_vac1_pareto.loc[np.round(1-R_vals[0],2), np.round(1-R_vals[1],2)] = vac_pareto_vac1[index]
    df_vac2_optimal.loc[np.round(1-R_vals[0],2), np.round(1-R_vals[1],2)] = vac_optimal_vac2[index]
    df_vac2_pareto.loc[np.round(1-R_vals[0],2), np.round(1-R_vals[1],2)] = vac_pareto_vac2[index]

d_optimal = - df_R_optimal + df_R_pop 
d_pareto = - df_R_pareto + df_R_pop 

df_vac_pareto = (df_vac1_pareto + df_vac2_pareto)/2
df_vac_optimal = (df_vac1_optimal + df_vac2_optimal)/2

maximum = np.max(d_optimal.to_numpy().flatten())
minimum = np.min(d_pareto.to_numpy().flatten())
scale=1000
fig = plt.figure(figsize = (39, 28))
gs = GridSpec(3, 3, figure=fig)



ax =  fig.add_subplot(gs[0, 0]) 
sns.heatmap(np.round(d_optimal.iloc[:,::-1] / scale ,0), cmap= "Blues", annot=True, fmt='g',
            ax=ax, vmin=minimum/scale, vmax=maximum/scale, cbar = False)
ax.set_xlabel("Degree of NPIs in country B")
ax.set_ylabel("Degree of NPIs in country A")
ax.set_title("Optimal strategy - Saved lives in 1,0000")
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


ax =  fig.add_subplot(gs[0, 1]) 
sns.heatmap(np.round(d_pareto.iloc[:,::-1]/scale,0), cmap= "Blues", annot=True, fmt='g', ax=ax, vmin=minimum/scale, vmax=maximum/scale, cbar=True)
ax.set_xlabel("Degree of NPIs in country B")
ax.set_ylabel("Degree of NPIs in country A")
ax.set_title("Pareto optimal  strategy - Saved lives in 1,0000")
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


ax =  fig.add_subplot(gs[1, 0]) 
sns.heatmap(np.round(df_vac_optimal.iloc[:,::-1]*100 ,0), cmap= "Greens", annot=True, fmt='g',
            ax=ax, vmin=20, vmax=100, cbar=False)
ax.set_xlabel("Degree of NPIs in country B")
ax.set_ylabel("Degree of NPIs in country A")
ax.set_title("Optimal strategy - % of vaccine assigend to country A")
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


ax =  fig.add_subplot(gs[1, 1]) 
sns.heatmap(np.round(df_vac_pareto.iloc[:,::-1]*100 ,0), cmap= "Greens", annot=True, fmt='g', ax=ax, vmin=20, vmax=100, cbar=True)
ax.set_xlabel("Degree of NPIs in country B")
ax.set_ylabel("Degree of NPIs in country A")
ax.set_title("Pareto strategy - % of vaccine assigend to country A")
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


fig.tight_layout(pad=3)

fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/plot_four",
    bbox_inches="tight",
)