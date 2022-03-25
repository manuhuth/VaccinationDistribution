import pickle
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.cbook import get_sample_data

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


def imscatter(x, y, image, ax=None, zoom=1):
    if ax is None:
        ax = plt.gca()
    try:
        image = plt.imread(image)
    except TypeError:
        # Likely already an array...
        pass
    im = OffsetImage(image, zoom=zoom)
    x, y = np.atleast_1d(x, y)
    artists = []
    for x0, y0 in zip(x, y):
        ab = AnnotationBbox(im, (x0, y0), xycoords='data', frameon=False)
        artists.append(ax.add_artist(ab))
    ax.update_datalim(np.column_stack([x, y]))
    ax.autoscale()
    return artists

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

count_plot = 65
lwidth = 2.5
y_position = 1.2
y_size = 42
x_position = -0.15
font_size_annot = 21

size = 32
font = {"family": "normal", "weight": "normal", "size": size}
matplotlib.rc("font", **font)
matplotlib.rcParams['xtick.labelsize'] = size - 4
matplotlib.rcParams['ytick.labelsize'] = size - 4

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

d_optimal = (- df_R_optimal + df_R_pop) / df_R_pop *100
d_pareto = (- df_R_pareto + df_R_pop) / df_R_pop *100

df_vac_pareto = (df_vac1_pareto + df_vac2_pareto)/2
df_vac_optimal = (df_vac1_optimal + df_vac2_optimal)/2

maximum = np.max(d_optimal.to_numpy().flatten())
minimum = np.min(d_pareto.to_numpy().flatten())
scale=1


#----------------------------------------------------------------------------------
fig = plt.figure(figsize = (39, 36))
gs = GridSpec(3, 3, figure=fig, width_ratios=[1, 1.2, 1])
gs.update(wspace=0.05, hspace=0.5)

NPIs = np.linspace(0, 0.7, 8)
wild = 0.24
mutant = wild*2

R0_wild = (1-NPIs) * wild / 0.1
R0_mutant = (1-NPIs)*mutant / 0.1
ax =  fig.add_subplot(gs[0, 0]) 
ax.yaxis.grid(alpha=0.6)
ax.plot(NPIs, R0_wild, color = "seagreen", label = "Wild-type")
ax.plot(NPIs, R0_mutant, color = "steelblue", label = "Variant")
for index in range(len(NPIs)):
    imscatter([NPIs[index]], [R0_mutant[index]],
              "/home/manuel/Documents/VaccinationDistribution/paper/images/virus_blue.png",
              ax=ax, zoom = R0_mutant[index]/1.8*0.08)
for index in range(len(NPIs)):
    imscatter([NPIs[index]], [R0_wild[index]],
              "/home/manuel/Documents/VaccinationDistribution/paper/images/virus_green.png",
              ax=ax, zoom = R0_wild[index]/1.8*0.08)

#ax.scatter(NPIs, R0_wild, color = "C0",  s= 140)
#ax.scatter(NPIs, R0_mutant, color = "C1",  s=140)
ax.legend()
ax.set_title("Effective reproduction rate\nof unvaccinated individuals")
ax.set_xlabel("Stringency Index")
ax.set_ylabel("Effective reproduction rate")
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
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

ax =  fig.add_subplot(gs[0, 1]) 
sns.heatmap(np.round(df_vac_optimal.iloc[:,::-1]*100 ,0), cmap= "Greens", annot=True, fmt='g',
            ax=ax, vmin=20, vmax=100, cbar=False, square=True, annot_kws={"fontsize":font_size_annot})
for t in ax.texts: t.set_text(t.get_text() + "%")
ax.set_xlabel("Stringency Index in Country 2")
ax.set_ylabel("Stringency Index in Country 1")
ax.set_title("Vaccine assigend to Country 1\n- Optimal strategy")
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


ax =  fig.add_subplot(gs[0, 2]) 
sns.heatmap(np.round(df_vac_pareto.iloc[:,::-1]*100 ,0), cmap= "Greens", annot=True,
            fmt='g', ax=ax, vmin=20, vmax=100, cbar=True, square=True, annot_kws={"fontsize":font_size_annot})
for t in ax.texts: t.set_text(t.get_text() + "%")
ax.set_xlabel("Stringency Index in Country 2")
ax.set_ylabel("")
ax.set_title("Vaccine assigend to Country 2\n- Pareto strategy")



count_plot += 1


ax =  fig.add_subplot(gs[1, 0]) 
sns.heatmap(np.round(df_R_pop .iloc[:,::-1] / scale ,0), cmap= "Blues", annot=False, fmt='f',
            ax=ax, square=True)#, vmin=minimum/scale, vmax=maximum/scale, cbar = False)

ax.set_xlabel("Stringency Index in Country 2")
ax.set_ylabel("Stringency Index in Country 1")
ax.set_title("Deceasd individuals\n- Pop.-based strategy")
ax.text(
        x_position-0.07,
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
sns.heatmap(np.round(d_optimal.iloc[:,::-1] ,0), cmap= "Blues", annot=True, fmt='g',
            ax=ax, vmin=minimum, vmax=maximum, cbar = False, square=True, annot_kws={"fontsize":font_size_annot})
for t in ax.texts: t.set_text(t.get_text() + "%")
ax.set_xlabel("Stringency Index in Country 2")
ax.set_ylabel("Stringency Index in Country 1")
ax.set_title("Reduction in deceased individuals\n- Optimal strategy")
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


ax =  fig.add_subplot(gs[1, 2]) 
sns.heatmap(np.round(d_pareto.iloc[:,::-1],0), cmap= "Blues", annot=True,
            fmt='g', ax=ax, vmin=minimum, vmax=maximum, cbar=True, square=True, annot_kws={"fontsize":font_size_annot})
for t in ax.texts: t.set_text(t.get_text() + "%")
ax.set_xlabel("Stringency Index in Country 2")
ax.set_ylabel("")
ax.set_title("Reduction in deceased individuals\n- Pareto strategy")


fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/plot_four",
    bbox_inches="tight",
)