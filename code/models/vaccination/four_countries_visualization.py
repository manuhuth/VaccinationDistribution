import matplotlib.pyplot as plt
import pickle
import numpy as np
from pymoo.visualization.radviz import Radviz
from models.vaccination.spline_example import get_spline
from matplotlib.gridspec import GridSpec
import matplotlib as mpl

from functions.plot_tools import plot_bars_multiple
from functions.plot_tools import plot_trajectories_aggregated
from functions.plot_tools import plot_trajectories_aggregated_vac
from functions.plot_tools import plot_vac_allocated

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
letter_y = 1.09
letter_size = 16
fig = plt.figure(constrained_layout=True, figsize=(18, 10))
gs = GridSpec(4, 4, figure=fig)
count_plot = 97


# Fig1 - Deaths relative to pop
ax = fig.add_subplot(gs[0, 0:1])
pop_deaths = np.array(dict_out["pop_trajectories"]["deaths"])
unconstr_deaths = np.array(dict_out["best_trajectories"]["deaths"])
constr_deaths = np.array(dict_out["pareto_trajectories"]["deaths"])

plot_bars_multiple(
    ax,
    unconstr_deaths,
    pop_deaths,
    constr_deaths,
    label_optimal="Optimal",
    label_Pareto="Pareto",
    X=["Total"] + ["Belgium", "France", "Germany", "United \nKingdom"],
    xlabel="Countries",
    ylabel="Difference in %",
    title="Number of deaths per strategy and \ncountry compared to population strategy",
    color_good="seagreen",
    color_bad="darkred",
    alpha=0.4,
    xlim=[-0.42, 5.7],
)
ax.text(
    -0.05,
    letter_y,
    chr(count_plot),
    horizontalalignment="center",
    verticalalignment="center",
    transform=ax.transAxes,
    weight="bold",
    size=letter_size,
)

count_plot += 1

# Fig 2,3,4 - Global Deaths, Infections,
targets = ["dead", "infectious", "susceptible"]
label_names = ["deaths", "active infectious cases", "susceptible"]
for index in range(len(targets)):
    x_pos = 1 + index
    ax = fig.add_subplot(gs[0, x_pos])
    plot_trajectories_aggregated(
        ax,
        length=420,
        pop_trajectory=dict_out["pop_trajectories"]["trajectories"],
        unconstr_trajectory=dict_out["best_trajectories"]["trajectories"],
        constr_trajectory=dict_out["pareto_trajectories"]["trajectories"],
        labels=["Optimal", "Pareto", "Population"],
        colors=[cmap(0), cmap(1), cmap(2)],
        alphas=[0.2, 0.3, 0.3],
        xlabel="Weeks",
        ylabel=f"{label_names[index].capitalize()} \nin millions",
        title=f"Total number of {label_names[index]}",
        target=[targets[index]],
        scale=scale,
    )
    ax.text(
        -0.05,
        letter_y,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=letter_size,
    )
    count_plot += 1


# Fig 5 - 12
areas = ["countryA", "countryB", "countryC", "countryD"]
countries = ["Belgium", "France", "Germany", "the United Kingdom"]
for index_y in range(2):
    y_pos = index_y + 1
    for index in range(len(countries)):
        x_pos = index
        ax = fig.add_subplot(gs[y_pos, x_pos])
        plot_legend = False
        if index == 0:
            plot_legend = True
        title = f"{label_names[index_y].capitalize()} in {countries[index]}"
        if index_y == 1:
            title = f"{label_names[index_y].capitalize()} \nin {countries[index]}"
        plot_trajectories_aggregated(
            ax,
            length=420,
            pop_trajectory=dict_out["pop_trajectories"]["trajectories"],
            unconstr_trajectory=dict_out["best_trajectories"]["trajectories"],
            constr_trajectory=dict_out["pareto_trajectories"]["trajectories"],
            labels=["Optimal", "Pareto", "Population"],
            colors=[cmap(0), cmap(1), cmap(2)],
            alphas=[0.1, 0.1, 0.1],
            xlabel="Weeks",
            ylabel=f"{label_names[index_y].capitalize()} \nin millions",
            title=title,
            target=[targets[index_y], areas[index]],
            scale=scale,
            plot_legend=plot_legend,
            fill_between=False,
        )
        ax.text(
            -0.05,
            letter_y,
            chr(count_plot),
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            weight="bold",
            size=letter_size,
        )
        count_plot += 1


# Fig 13-16
# for index in range(len(countries)):
#        x_pos = index
#        ax = fig.add_subplot(gs[3, x_pos])
#        plot_legend = False
#        if index == 0:
#            plot_legend = True
#        plot_trajectories_aggregated_vac(ax, length = 420,
#                                         pop_trajectory = dict_out["pop_trajectories"]["trajectories"],
#                                         unconstr_trajectory = dict_out["best_trajectories"]["trajectories"],
#                                         constr_trajectory = dict_out["pareto_trajectories"]["trajectories"],
#                                         labels = ["Population", "Optimal", "Pareto"],
#                                         colors= ["C0", "C1", "C2"],
#                                         alphas= [0.6, 0.4, 0.2],
#                                         xlabel = "Weeks",
#                                         ylabel = "% immune",
#                                        title = f"Proportion of immune population in {countries[index]}",
#                                         scale=10**6,
#                                         fill_between = False, plot_legend = plot_legend)
#        ax.text(-0.05, letter_y, chr(count_plot), horizontalalignment='center',
#                verticalalignment='center',
#                transform = ax.transAxes, weight="bold", size= letter_size)
#        count_plot += 1
fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/4countries_deaths",
    bbox_inches="tight",
)


# Figures Vaccines
fig = plt.figure(constrained_layout=True, figsize=(18, 10))
gs = GridSpec(2, 4, figure=fig)
count_plot = 97
vac = ["vac1", "vac2"]
vac_type = ["mRNA", "Vector"]

# Fig - 1 - 16
for index_vac in range(len(vac)):
    for index_areas in range(len(areas)):

        ax = fig.add_subplot(gs[index_vac, index_areas])

        plot_vac_allocated(
            ax,
            colors=[cmap(0), cmap(1), cmap(2)],
            time=dict_out["vaccine"]["t"],
            dict_out=dict_out,
            index_vac=index_vac,
            index_areas=index_areas,
            areas=areas,
            scale=scale,
            countries=countries,
            col_vac1="C7",
            col_vac2="C8",
            label_vac1="mRNA",
            label_vac2="Vector",
            vac=["vac1", "vac2"],
            types=["unconstrained", "constrained", "pop"],
            ylabel="% of total \navailable vaccine",
            xlabel="Weeks",
            labels=["Optimal", "Pareto", "Population"],
            alphas=[0.1, 0.1, 0.1],
            axvline_x=40,
            ylim=[-0.05, 0.9],
            title="Percentage of available "
            + vac_type[index_vac]
            + " \nvaccine allocated to ",
            total=False,
            spline_xx=spline_xx,
            numb_xx=4,
            s=30,
        )

        ax.text(
            -0.05,
            letter_y,
            chr(count_plot),
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            weight="bold",
            size=letter_size,
        )
        count_plot += 1


fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/4countries_vaccines",
    bbox_inches="tight",
)


# add one plot with sensitivities for all vaccines (say they are all within an improvement of maybe colour by improvement)


fig = plt.figure(constrained_layout=True, figsize=(18, 10))
gs = GridSpec(2, 4, figure=fig)
count_plot = 97
vac = ["vac1", "vac2"]
vac_type = ["mRNA", "Vector"]
linewidth = 0.5
color_best = "firebrick"
linewidth_best = 1
color_rest = "C0"
linewidth_rest = 0.5
ylabel = "% of total \nvaccine allocated"
xlabel = "Weeks"
title = "Vaccine allocation for best \n10 runs "

df = dict_out["vaccine_all_paths_best"]
for index_vac in range(len(vac)):
    for index_areas in range(len(areas)):

        ax = fig.add_subplot(gs[index_vac, index_areas])
        for index_position in range(0, 10):

            cols = [
                x
                for x in df.columns
                if areas[index_areas] in x
                and vac[index_vac] in x
                and f"a{index_position}_" in x
            ]
            for index_cols in cols:
                if index_position == 0:
                    color = color_best
                    linewidth = linewidth_best
                else:
                    color = color_rest
                    linewidth = linewidth_rest
                ax.plot(df["t"] / 7, df[index_cols], linewidth=linewidth, color=color)
                ax.set_title(title + f"in {countries[index_areas].capitalize()}")
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        ax.text(
            -0.05,
            letter_y,
            chr(count_plot),
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            weight="bold",
            size=letter_size,
        )
        count_plot += 1


path = "/home/manuel/Documents/VaccinationDistribution/code/objects/manyoo_dict.pkl"

with open(
    path,
    "rb",
) as input:
    di = pickle.load(input)

f = di["f"]
x = di["x"]
p = di["p"]

a = np.concatenate(
    [
        f,
        [
            p[
                0,
            ]
        ],
    ]
)
best = a[np.argmin(np.sum(a, axis=1)), :]

# plot Pareto front -> 2 dimensional scatter plots
best_point = best
names = ["Belgium", "France", "Germany", "UK"]
scale = 10 ** 6
name_pop = "Population based point"
color_improvements = "C1"
color_pareto_improvements = "C2"
color_pop = "C4"
s_improvements = 22
s_pareto_improvements = 22
s_pop = 22
label_pareto_improvements = "Pareto improvements"
label_improvements = "Non-Pareto improvements"
label_best = "Optimal point"
s_best = 22
color_best = "C3"
color_all = "C0"
label_all = "No improvement"
s_all = 14


n_countries = len(names)
index_improvement = np.sum(f, axis=1) <= np.sum(p, axis=1)
index_pareto_improvement = np.sum((f - p) < 0, axis=1) == f.shape[1]
pareto_improvements = f.loc[index_pareto_improvement]
improvements = f.loc[index_improvement]
fig, axs = plt.subplots(n_countries, n_countries, figsize=(10, 10))
for i in range(len(names)):
    for j in range(len(names)):
        if j <= i:
            ax = axs[i][j]
            if i == j:
                ax.text(
                    0.5, 0.5, names[i], ha="center", va="center", zorder=10, fontsize=20
                )
                ax.get_xaxis().set_visible(False)
                ax.get_yaxis().set_visible(False)
            else:

                ax.scatter(
                    f[i] / scale,
                    f[j] / scale,
                    s=s_all,
                    color=color_all,
                    label=label_all,
                    alpha=0.5,
                )
                ax.scatter(
                    improvements[i] / scale,
                    improvements[j] / scale,
                    s=s_improvements,
                    color=color_improvements,
                    label=label_improvements,
                )
                ax.scatter(
                    pareto_improvements[i] / scale,
                    pareto_improvements[j] / scale,
                    s=s_pareto_improvements,
                    color=color_pareto_improvements,
                    label=label_pareto_improvements,
                )
                ax.scatter(
                    best_point[i] / scale,
                    best_point[j] / scale,
                    s=s_best,
                    color=color_best,
                    label=label_best,
                )

                ax.set_xlabel(names[i])
                ax.set_ylabel(names[j])
                ax.scatter(
                    p[0, :][i] / scale,
                    p[0, :][j] / scale,
                    s=s_pop,
                    label=name_pop,
                    color=color_pop,
                )

# fig.subplots_adjust(hspace=0.4)
# fig.subplots_adjust(vspace=0.4)
fig.tight_layout()

handles4, labels4 = axs[1][0].get_legend_handles_labels()
fig.legend(
    handles=handles4,
    labels=labels4,
    loc="lower center",
    ncol=4,
    bbox_to_anchor=[0.5, -0.05],
)
