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

from numpy import trapz

from functions.plot_tools import plot_horizontal_bars_annotated


font = {"family": "normal", "weight": "normal", "size": 10}

matplotlib.rc("font", **font)
import matplotlib.patches as mpatches

from functions.plot_tools import plot_pareto_front
from functions.plot_tools import plot_bars_deaths


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

df_R_pop = pd.DataFrame(np.nan, columns = np.round(1-np.linspace(0.3,1, 8),2), index =np.round(1-np.linspace(0.3,1, 8),2))

for index in range(len(dict_R_cases.keys())):
    key = f"case{index}"
    R_vals = dict_R_cases[key]
    
    
    deaths_pop = dicts[index]["population_based"]
    #death_pop_A = deaths_pop[0]
    #death_pop_B = deaths_pop[1]
    
    df_R_pop.loc[np.round(1-R_vals[0],2), np.round(1-R_vals[1],2)] = np.sum(deaths_pop)



sns.heatmap(df_R_pop/1000)

