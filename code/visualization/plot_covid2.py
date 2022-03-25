import matplotlib.pyplot as plt
import pickle
import pandas as pd
import numpy as np
from models.vaccination.spline_example import get_spline
from matplotlib.gridspec import GridSpec
import matplotlib as mpl


from functions.plot_tools import plot_horizontal_bars_annotated_many
from functions.plot_tools import compute_incidences

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

sizes = [11455.5 * 1000, 67012.9 * 1000, 83019.2 * 1000, 66647.1 * 1000]
countries = ["countryA", "countryB", "countryC", "countryD"]

tra = ["best_trajectories", "pareto_trajectories", "pop_trajectories"]

incidences = {}
for index2 in tra:
    incidence = {}
    for index in range(len(countries)):
        incidence[countries[index]] = compute_incidences(trajectories=dict_out[index2]["trajectories"].reset_index(drop =True),
                               viruses = ["virus1", "virus2"],
                               countries = [countries[index]] ,
                               time = np.linspace(0, 420, 6000),
                               lambda1 = 0.1,
                               habitant_scale = 100000 / sizes[index],)
    incidences[index2] = incidence



country = ["Belgium", "France", "Germany", "United Kingdom"]





y_position = 1.08
y_size = 16
count_plot = 65
fig = plt.figure(constrained_layout=True, figsize=(12, 9))
gs = GridSpec(3, 2, figure=fig)

# Fig1 - Deaths relative to pop
ax = fig.add_subplot(gs[0, 0:2])


plot_horizontal_bars_annotated_many(ax,
                               global_optimum=unconstr_deaths[1:5], population_based=pop_deaths[1:5], pareto_optimum=constr_deaths[1:5],
                               scale = 10**5,
                               color_vline = "black",
                               linestyle_vline = "dashed", title="",
                               category_names=country,
                               category_colors=["steelblue", "C1", "seagreen", "firebrick"], bbox_to_anchor=(0.75,0))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
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
ax.set_title("Deaths in millions")



time = pd.date_range("2021-01-01", "2022-02-25", freq="6050s")#np.linspace(0,140, 5999)
labels=["Optimal", "Pareto optimal", "Population\nbased"]
colors = ["steelblue", "C1", "black"]
linestyles=["solid", "solid", "dashed"]
alphas = [0.7,0.7,0.8]
y_lim = [0, 1800]
y_label= "7-day incidences per\n100,000 inhabtiants"

for index2 in range(len(countries)):
    g = [gs[1, 0], gs[1, 1], gs[2,0], gs[2,1]]

    ax = fig.add_subplot(g[index2])
    
    for index in range(len(tra)):
        incidence= incidences[tra[index]]
        country_tra = incidence[countries[index2]]
    
        country_sum = country_tra["virus1"] + country_tra["virus2"]
        ax.plot(time,country_sum, label = labels[index], color=colors[index], linestyle=linestyles[index], alpha=alphas[index])
        if index == 10:
            incidence_plus = incidences[tra[index+1]]
            country_tra_plus = incidence_plus[countries[index2]]
            country_sum_plus = country_tra_plus["virus1"] + country_tra_plus["virus2"]
            ax.fill_between(time, country_sum.to_numpy().flatten(), country_sum_plus.to_numpy().flatten(),
                            where = country_sum.to_numpy().flatten() < country_sum_plus.to_numpy().flatten(), facecolor=colors[index+1],
                            alpha=alphas[index])
            ax.fill_between(time, country_sum.to_numpy().flatten(), country_sum_plus.to_numpy().flatten(),
                            where = country_sum.to_numpy().flatten() > country_sum_plus.to_numpy().flatten(), facecolor=colors[index],
                            alpha=alphas[index])
            
            ax.fill_between(time, np.repeat(0,len(country_sum_plus.to_numpy().flatten())), country_sum.to_numpy().flatten(),
                            where = country_sum.to_numpy().flatten() < country_sum_plus.to_numpy().flatten(), facecolor=colors[index],
                            alpha=alphas[index])
            ax.fill_between(time, np.repeat(0,len(country_sum_plus.to_numpy().flatten())), country_sum_plus.to_numpy().flatten(),
                            where = country_sum.to_numpy().flatten() > country_sum_plus.to_numpy().flatten(), facecolor=colors[index+1],
                            alpha=alphas[index])
        if index in [10]:
            ax.fill_between(time, np.repeat(0,len(country_sum.to_numpy().flatten())), country_sum.to_numpy().flatten(),
                            facecolor=colors[index],
                            alpha=alphas[index])
        ax.yaxis.grid(alpha=0.6)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
    if index2==0:
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

    ax.set_title(country[index2])
    ax.set_ylim(y_lim)
    if index2 in [0,2]:
        ax.set_ylabel(y_label)
    ax.xaxis.set_tick_params(rotation=45)

    
fig.tight_layout(pad=1.5)

fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/covid_plot_two",
    bbox_inches="tight",
)