import matplotlib.pyplot as plt
import pickle
import numpy as np
from pymoo.visualization.radviz import Radviz
from models.vaccination.spline_example import get_spline
from matplotlib.gridspec import GridSpec
import matplotlib as mpl


from functions.plot_tools import plot_horizontal_bars_annotated


cmap = mpl.cm.get_cmap("tab10")
font = {"family": "normal", "weight": "normal", "size": 10}

mpl.rc("font", **font)
path = "/home/manuel/Documents/VaccinationDistribution/code/objects/results_trust_constr_4countries.pkl"

with open(
    path,
    "rb",
) as input:
    dict_out = pickle.load(input)

spline_xx = {
    "xx0": 0.0,
    "xx1": 95.33333333333333,
    "xx2": 190.66666666666666,
    "xx3": 286.0,
    "xx4": 381.3333333333333,
    "xx5": 420.0,
}
scale = 10 ** 6

pop_deaths = np.array(dict_out["pop_trajectories"]["deaths"])
unconstr_deaths = np.array(dict_out["best_trajectories"]["deaths"])
constr_deaths = np.array(dict_out["pareto_trajectories"]["deaths"])

results = {"pop" : pop_deaths, "optimal" : unconstr_deaths, "pareto" : constr_deaths}
cat_names = ["Optimal\nStrategy", "Population\nStrategy", "Pareto optimal\nStrategy"]

letter_y = 1.09
letter_size = 16
fig = plt.figure(constrained_layout=True, figsize=(18, 10))
gs = GridSpec(4, 4, figure=fig)
count_plot = 65


# Fig1 - Deaths relative to pop
ax = fig.add_subplot(gs[0, 0:1])


plot_horizontal_bars_annotated(results=results,
                               category_names=cat_names,
                               ax=ax,)
