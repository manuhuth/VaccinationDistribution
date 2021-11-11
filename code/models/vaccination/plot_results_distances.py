import pickle

import matplotlib
import matplotlib.pyplot as plt

from functions.plot_tools import plot_distance_curves

font = {"family": "normal", "weight": "normal", "size": 30}

matplotlib.rc("font", **font)


# load outputs 2 countries, ordinary set-up
path = "/home/manuel/Documents/VaccinationDistribution/code/objects/"
names = [
    "initalEqual_vacEqual",
    "initalUnequal_vacEqual",
    "initalEqual_vacUnequal",
    "initalUnequal_vacUnequal",
]
dicts = {}

for i in names:
    path1 = path + i + "_distance.pkl"
    with open(
        path1,
        "rb",
    ) as input:
        loaded_object = pickle.load(input)
        dicts[i] = loaded_object


Titles = [
    "Eq. initial states; Eq. vac.",
    "Uneq. initial states; Eq. vac.",
    "Eq. initial states; Uneq. vac.",
    "Uneq. initial states; Uneq. vac.",
]

max_i = None
fig, axs = plt.subplots(3, 4, figsize=(65, 60))

count_plot = 97
for i in range(len(dicts)):
    ax = plot_distance_curves(
        dict_use=dicts[names[i]],
        max_index=max_i,
        linewidth=5,
        color_optimal="C0",
        color_pareto="C1",
        color_pop="C2",
        label_optimal="Optimal",
        label_pareto="pareto",
        label_pop="pop",
        var="fval",
        relative=True,
        ax=axs[0][i],
        x_label="Distance parameter",
        y_label="% deaths",
        title=Titles[i],
        vline=None,
        vline_color="C3",
        vline_width=4,
        vline_label="Previous parameter",
        ylim=[-26, 0.1],
    )

    ax.text(
        -0.05,
        1.03,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=35,
    )
    if i == 0:
        ax.text(
            -0.25,
            0.5,
            "Whole parameter space",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            weight="bold",
            size=35,
            rotation="vertical",
        )
    count_plot += 1

for i in range(len(dicts)):
    ax = plot_distance_curves(
        dict_use=dicts[names[i]],
        max_index=14,
        linewidth=5,
        color_optimal="C0",
        color_pareto="C1",
        color_pop="C2",
        label_optimal="Optimal",
        label_pareto="pareto",
        label_pop="pop",
        var="fval",
        relative=True,
        ax=axs[1][i],
        x_label="Distance parameter",
        y_label="% deaths",
        title=Titles[i],
        vline=1 - 1 / 5000,
        vline_color="C3",
        vline_width=4,
        vline_label="Parameter used before",
        ylim=[-26, 0.1],
        v_ymin=-25,
        v_ymax=0,
    )
    ax.text(
        -0.05,
        1.03,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=35,
    )
    if i == 0:
        ax.text(
            -0.25,
            0.5,
            "Relevant interval",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            weight="bold",
            size=35,
            rotation="vertical",
        )
    count_plot += 1
handles, labels = ax.get_legend_handles_labels()

for i in range(len(dicts)):
    ax = plot_distance_curves(
        dict_use=dicts[names[i]],
        max_index=max_i,
        linewidth=5,
        color_optimal="C0",
        color_pareto="C1",
        color_pop="C2",
        label_optimal="Optimal strategy",
        label_pareto="Pareto optimal strategy",
        label_pop="Population based strategy",
        var="fval",
        relative=False,
        ax=axs[2][i],
        x_label="Distance parameter",
        y_label="Deaths",
        title=Titles[i],
        vline=None,
        vline_color="C3",
        vline_width=4,
        vline_label="Previous parameter",
        ylim=[200000 - 10, 360000 + 10],
    )

    ax.text(
        -0.05,
        1.03,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=35,
    )
    if i == 0:
        ax.text(
            -0.25,
            0.5,
            "Total numbers",
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            weight="bold",
            size=35,
            rotation="vertical",
        )
    count_plot += 1


handles4, labels4 = ax.get_legend_handles_labels()
fig.legend(
    handles=handles4 + [handles[2]],
    labels=labels4 + [labels[2]],
    loc="lower center",
    ncol=4,
    bbox_to_anchor=[0.5, 0.08],
)

fig.subplots_adjust(hspace=0.3)
fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/scenario2",
    bbox_inches="tight",
)
