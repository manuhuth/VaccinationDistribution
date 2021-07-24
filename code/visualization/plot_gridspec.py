import pickle
from models.vaccination.optimization_functions_estimagic import get_sum_of_states

import numpy as np
import matplotlib.pyplot as plt

plt.style.use("default")
# plt.style.use('ggplot')
plt.rcParams.update(
    {
        "text.usetex": True,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica"],
        "grid.color": "white",
    }
)
plt.rcParams["lines.linewidth"] = 1.5
plt.rc("text", usetex=True)
plt.rcParams["text.latex.preamble"] = [r"\usepackage{cmbright}"]
plt.rcParams["axes.facecolor"] = "#E6E6E6"


from functions.run_sbml import get_model_and_solver_from_sbml

type_optim = "piecewise"
if type_optim == "piecewise":
    plot_name = "Piecewise constant"
else:
    plot_name = "Spline"
# ----------Load optimization data---------------------------------------------
load = (
    "/home/manuel/Documents/VaccinationDistribution/code/objects/output_"
    + type_optim
    + ".pkl"
)
with open(
    load,
    "rb",
) as input:
    dict_out = pickle.load(input)
df_optimal_results = dict_out["df_optimal_results"]
df_unrestricted_results = dict_out["df_unrestricted_results"]
model_optimal = dict_out["model_optimal"]
model_unrestricted = dict_out["model_seperated"]
model_current = dict_out["model_current"]
model_directory = dict_out["model_directory"]
path_sbml = dict_out["path_sbml"]
model_name = dict_out["model_name"]
compare_total = dict_out["compare_total"]
compare_A = dict_out["compare_A"]
compare_B = dict_out["compare_B"]
time_points = dict_out["time_points"]
type_method = dict_out["type"]
print(compare_total)
print(compare_A)
print(compare_B)

# ------Set plot path----------------------------------------------------------
plot_path = (
    "/home/manuel/Documents/VaccinationDistribution/paper/images/" + type_method + "_"
)

model_solver = model_and_solver = get_model_and_solver_from_sbml(
    path_sbml=path_sbml,
    model_name=model_name,
    model_directory=model_directory,
    only_import=True,
)
model = model_solver["model"]
observables = model.getObservableNames()


def sum_states(df, names):
    for name in names:
        columns = [col for col in df.columns if all(x in col for x in names)]
    df_reduced = df[columns]
    result = df_reduced @ np.repeat(1, df_reduced.shape[1])

    return result


def plot_gridspec(
    y,
    time,
    title="Amount of vaccine doses (in mil.)",
    legend_next_to_plot=True,
    legend_location="lower center",
    ylim=[0, 1.4],
    xlim=None,
    bbox=(0.515, -0.05),
    line_thickness=0.6,
    ylabel=["Vaccine doses", "Vaccine_doses"],
    xlabel="Weeks",
    row_labels=["Country A", "Country B"],
    title_col1="Pareto strategy",
    title_col2="Unrestricted strategy",
    title_col3="Current strategy",
    xticks=np.linspace(0, 20, 6),
    path=None,
    twin_axes=None,
):

    cols = [title_col1, title_col2, title_col3]

    fig = plt.figure()
    gs = fig.add_gridspec(2, 3, hspace=0.1, wspace=0.1)
    (ax1, ax2, ax5), (ax3, ax4, ax6) = gs.subplots(sharex="col", sharey="row")
    fig.suptitle(title)

    axes = [ax1, ax2, ax3, ax4, ax5, ax6]

    for index in range(len(axes)):
        ax = axes[index]
        y_i = y[index]
        count = 0
        for index2 in list(y_i.keys()):
            y_plot = y_i[index2]
            if not ("linewidth" in y_plot):
                y_plot["linewidth"] = 1.5
            if count == 0 or count == 1  or index2 == "z_deterministic":
                ax.plot(
                    time,
                    y_plot["y"],
                    label=y_plot["label"],
                    color=y_plot["color"],
                    linestyle=y_plot["linestyle"],
                    linewidth=y_plot["linewidth"],
                )
            else:
                ax.plot(
                    time,
                    y_plot["y"],
                    color=y_plot["color"],
                    linestyle=y_plot["linestyle"],
                    linewidth=y_plot["linewidth"],
                )
            count += 1

    ax1.spines["top"].set_alpha(0)
    ax1.spines["bottom"].set_alpha(0)
    ax1.spines["right"].set_alpha(0)
    ax1.spines["left"].set_alpha(0)
    ax1.grid(alpha=line_thickness)
    ax1.set_ylabel(ylabel[0])
    # ax1.set_title(title_col1)

    ax2.spines["top"].set_alpha(0)
    ax2.spines["bottom"].set_alpha(0)
    ax2.spines["right"].set_alpha(0)
    ax2.spines["left"].set_alpha(0)
    ax2.grid(alpha=line_thickness)

    # ax2.set_title(title_col2)

    ax3.spines["top"].set_alpha(0)
    ax3.spines["bottom"].set_alpha(0)
    ax3.spines["right"].set_alpha(0)
    ax3.spines["left"].set_alpha(0)
    ax3.grid(alpha=line_thickness)
    ax3.set_xticks(xticks)
    ax3.set_ylabel(ylabel[1])
    ax3.set_xlabel(xlabel)

    ax4.spines["top"].set_alpha(0)
    ax4.spines["bottom"].set_alpha(0)
    ax4.spines["right"].set_alpha(0)
    ax4.spines["left"].set_alpha(0)
    ax4.grid(alpha=line_thickness)
    ax4.set_xticks(xticks)
    ax4.set_xlabel(xlabel)

    ax5.spines["top"].set_alpha(0)
    ax5.spines["bottom"].set_alpha(0)
    ax5.spines["right"].set_alpha(0)
    ax5.spines["left"].set_alpha(0)
    ax5.grid(alpha=line_thickness)

    ax6.spines["top"].set_alpha(0)
    ax6.spines["bottom"].set_alpha(0)
    ax6.spines["right"].set_alpha(0)
    ax6.spines["left"].set_alpha(0)
    ax6.grid(alpha=line_thickness)
    ax6.set_xticks(xticks)
    ax6.set_xlabel(xlabel)

    if ylim is not None:
        for index in axes:
            index.set_ylim(ylim)

    if xlim is not None:
        for index in axes:
            index.set_xlim(xlim)

    if twin_axes is not None:
        ax1_twin = ax1.twinx()
        ax2_twin = ax2.twinx()
        ax3_twin = ax3.twinx()
        ax4_twin = ax4.twinx()
        ax5_twin = ax5.twinx()
        ax6_twin = ax6.twinx()
        twins = [ax1_twin, ax2_twin, ax3_twin, ax4_twin, ax5_twin, ax6_twin]

        for index in range(len(twins)):
            ax_twin = twins[index]
            twin_y_i = twin_axes[index]
            for index2 in twin_y_i.keys():
                twin_y_plot = twin_y_i[index2]
                ax_twin.plot(
                    time,
                    twin_y_plot["y"],
                    label=twin_y_plot["label"],
                    color=twin_y_plot["color"],
                    linestyle=twin_y_plot["linestyle"],
                )
                ax_twin.spines["top"].set_alpha(0)
                ax_twin.spines["bottom"].set_alpha(0)
                ax_twin.spines["right"].set_alpha(0)
                ax_twin.spines["left"].set_alpha(0)
                ax_twin.xaxis.set_ticks_position("none")
                ax_twin.set_ylabel(twin_y_plot["ylabel"])

                ax_twin.set_ylim(twin_y_plot["ylim"])
                if index in [0, 1, 2, 3]:
                    ax_twin.set_yticklabels([])
                    ax_twin.yaxis.set_ticks_position("none")

    ax6.yaxis.set_ticks_position("none")
    ax5.yaxis.set_ticks_position("none")
    ax5.xaxis.set_ticks_position("none")
    ax4.yaxis.set_ticks_position("none")
    # ax3.xaxis.set_ticks_position("none")
    # ax3.yaxis.set_ticks_position("none")
    ax2.xaxis.set_ticks_position("none")
    ax2.yaxis.set_ticks_position("none")
    ax1.xaxis.set_ticks_position("none")
    # ax1.yaxis.set_ticks_position("none")
    pad = 2
    ax_cols = [ax1, ax2, ax5]
    for ax, col in zip(ax_cols, cols):
        ax.annotate(
            col,
            xy=(0.5, 1.05),
            xytext=(0, pad),
            xycoords="axes fraction",
            textcoords="offset points",
            size="large",
            ha="center",
            va="baseline",
        )
    ax_rows = [ax1, ax3]
    for ax, row in zip(ax_rows, row_labels):
        ax.annotate(
            row,
            xy=(-0.05, 0.5),
            xytext=(-ax.yaxis.labelpad - pad, 0),
            xycoords=ax.yaxis.label,
            textcoords="offset points",
            size="large",
            ha="right",
            va="center",
        )

    if legend_next_to_plot is True:
        if twin_axes is not None:
            handles1, labels1 = ax4.get_legend_handles_labels()
            handles_twin, labels_twin = ax4_twin.get_legend_handles_labels()
            handles = handles1 + handles_twin
            labels = labels1 + labels_twin
        else:
            handles, labels = ax4.get_legend_handles_labels()
        legend = fig.legend(
            handles,
            labels,
            bbox_to_anchor=bbox,
            loc=legend_location,
            ncol=3,
            facecolor="white",
            framealpha=1,
        )
        frame = legend.get_frame()
        frame.set_linewidth(0)

    plt.tight_layout()

    if path is not None:
        fig.savefig(path, bbox_inches="tight")

    return fig, axes


# Plot absolute vaccine doses
y11 = {
    "y": model_optimal["observables"]["observable_quantity_countryA_vac1"] / 10 ** 6,
    "label": "mRNA",
    "color": "steelblue",
    "linestyle": "-",
}
y12 = {
    "y": model_optimal["observables"]["observable_quantity_countryA_vac2"] / 10 ** 6,
    "label": "Vector",
    "color": "coral",
    "linestyle": "--",
}
y1 = {"first": y11, "second": y12}

y21 = {
    "y": model_unrestricted["observables"]["observable_quantity_countryA_vac1"]
    / 10 ** 6,
    "label": "mRNA",
    "color": "steelblue",
    "linestyle": "-",
}
y22 = {
    "y": model_unrestricted["observables"]["observable_quantity_countryA_vac2"]
    / 10 ** 6,
    "label": "Vector",
    "color": "coral",
    "linestyle": "--",
}
y2 = {"first": y21, "second": y22}

y31 = {
    "y": model_optimal["observables"]["observable_quantity_countryB_vac1"] / 10 ** 6,
    "label": "mRNA",
    "color": "steelblue",
    "linestyle": "-",
}
y32 = {
    "y": model_optimal["observables"]["observable_quantity_countryB_vac2"] / 10 ** 6,
    "label": "Vector",
    "color": "coral",
    "linestyle": "--",
}
y3 = {"first": y31, "second": y32}

y41 = {
    "y": model_unrestricted["observables"]["observable_quantity_countryB_vac1"]
    / 10 ** 6,
    "label": "mRNA",
    "color": "steelblue",
    "linestyle": "-",
}
y42 = {
    "y": model_unrestricted["observables"]["observable_quantity_countryB_vac2"]
    / 10 ** 6,
    "label": "Vector",
    "color": "coral",
    "linestyle": "--",
}
y4 = {"first": y41, "second": y42}

y51 = {
    "y": model_current["observables"]["observable_quantity_countryA_vac1"] / 10 ** 6,
    "label": "mRNA",
    "color": "steelblue",
    "linestyle": "-",
}
y52 = {
    "y": model_current["observables"]["observable_quantity_countryA_vac2"] / 10 ** 6,
    "label": "Vector",
    "color": "coral",
    "linestyle": "--",
}
y5 = {"first": y51, "second": y52}

y61 = {
    "y": model_current["observables"]["observable_quantity_countryB_vac1"] / 10 ** 6,
    "label": "mRNA",
    "color": "steelblue",
    "linestyle": "-",
}
y62 = {
    "y": model_current["observables"]["observable_quantity_countryB_vac2"] / 10 ** 6,
    "label": "Vector",
    "color": "coral",
    "linestyle": "--",
}
y6 = {"first": y61, "second": y62}

y = [y1, y2, y3, y4, y5, y6]


path_vaccine = plot_path + "vaccine_total_quantity"
plot_gridspec(
    y=y,
    time=np.array(time_points) / 7.0,
    title=f"Amount of vaccine doses in millions",
    legend_next_to_plot=True,
    legend_location="lower center",
    ylim=[0, 1.4],
    xlim=None,
    bbox=(0.515, -0.05),
    line_thickness=0.6,
    ylabel=["Vaccine doses", "Vaccine doses"],
    xlabel="Weeks",
    row_labels=["Country A", "Country B"],
    title_col1="Pareto strategy",
    title_col2="Unrestricted strategy",
    title_col3="Current strategy",
    xticks=np.linspace(0, 20, 6),
    path=path_vaccine,
)

# ----------------------------------------------
# Plot relative vaccines
y11 = {
    "y": model_optimal["observables"]["proportion_countryA_vac1"],
    "label": "mRNA",
    "color": "steelblue",
    "linestyle": "-",
}
y12 = {
    "y": model_optimal["observables"]["proportion_countryA_vac2"],
    "label": "Vector",
    "color": "coral",
    "linestyle": "--",
}
y1 = {"first": y11, "second": y12}

y21 = {
    "y": model_unrestricted["observables"]["proportion_countryA_vac1"],
    "label": "mRNA",
    "color": "steelblue",
    "linestyle": "-",
}
y22 = {
    "y": model_unrestricted["observables"]["proportion_countryA_vac2"],
    "label": "Vector",
    "color": "coral",
    "linestyle": "--",
}
y2 = {"first": y21, "second": y22}

y31 = {
    "y": model_optimal["observables"]["proportion_countryB_vac1"],
    "label": "mRNA",
    "color": "steelblue",
    "linestyle": "-",
}
y32 = {
    "y": model_optimal["observables"]["proportion_countryB_vac2"],
    "label": "Vector",
    "color": "coral",
    "linestyle": "--",
}
y3 = {"first": y31, "second": y32}

y41 = {
    "y": model_unrestricted["observables"]["proportion_countryB_vac1"],
    "label": "mRNA",
    "color": "steelblue",
    "linestyle": "-",
}
y42 = {
    "y": model_unrestricted["observables"]["proportion_countryB_vac2"],
    "label": "Vector",
    "color": "coral",
    "linestyle": "--",
}
y4 = {"first": y41, "second": y42}


y51 = {
    "y": model_current["observables"]["proportion_countryA_vac1"],
    "label": "mRNA",
    "color": "steelblue",
    "linestyle": "-",
}
y52 = {
    "y": model_current["observables"]["proportion_countryA_vac2"],
    "label": "Vector",
    "color": "coral",
    "linestyle": "--",
}
y5 = {"first": y51, "second": y52}

y61 = {
    "y": model_current["observables"]["proportion_countryB_vac1"],
    "label": "mRNA",
    "color": "steelblue",
    "linestyle": "-",
}
y62 = {
    "y": model_current["observables"]["proportion_countryB_vac2"],
    "label": "Vector",
    "color": "coral",
    "linestyle": "--",
}
y6 = {"first": y61, "second": y62}

y = [y1, y2, y3, y4, y5, y6]


path_vaccine = plot_path + "vaccine_fractions"
plot_gridspec(
    y=y,
    time=np.array(time_points) / 7.0,
    title=f"Fraction of vaccine doses",
    legend_next_to_plot=True,
    legend_location="lower center",
    ylim=[0, 1.1],
    xlim=None,
    bbox=(0.515, -0.05),
    line_thickness=0.6,
    ylabel=["Fraction", "Fraction"],
    xlabel="Weeks",
    row_labels=["Country A", "Country B"],
    title_col1="Pareto strategy",
    title_col2="Unrestricted strategy",
    title_col3="Current strategy",
    xticks=np.linspace(0, 20, 6),
    path=path_vaccine,
)


# ----------------------------------------------
# Plot absolute number of infectious individuals
col1 = "seagreen"
col2 = "firebrick"
y11 = {
    "y": sum_states(model_optimal["states"], ["infectious", "countryA", "virW"])
    / 10 ** 6,
    "label": "Infectious wild type",
    "color": col1,
    "linestyle": "-",
}
y12 = {
    "y": sum_states(model_optimal["states"], ["infectious", "countryA", "virM"])
    / 10 ** 6,
    "label": "Infectious  mutant",
    "color": col2,
    "linestyle": "--",
}
y1 = {"first": y11, "second": y12}

y21 = {
    "y": sum_states(model_unrestricted["states"], ["infectious", "countryA", "virW"])
    / 10 ** 6,
    "label": "Infectious  wild type",
    "color": col1,
    "linestyle": "-",
}
y22 = {
    "y": sum_states(model_unrestricted["states"], ["infectious", "countryA", "virM"])
    / 10 ** 6,
    "label": "Infectious  mutant",
    "color": col2,
    "linestyle": "--",
}
y2 = {"first": y21, "second": y22}

y31 = {
    "y": sum_states(model_optimal["states"], ["infectious", "countryB", "virW"])
    / 10 ** 6,
    "label": "Infectious wild type",
    "color": col1,
    "linestyle": "-",
}
y32 = {
    "y": sum_states(model_optimal["states"], ["infectious", "countryB", "virM"])
    / 10 ** 6,
    "label": "Infectious mutant",
    "color": col2,
    "linestyle": "--",
}
y3 = {"first": y31, "second": y32}

y41 = {
    "y": sum_states(model_unrestricted["states"], ["infectious", "countryB", "virW"])
    / 10 ** 6,
    "label": "Infectious  wild type",
    "color": col1,
    "linestyle": "-",
}
y42 = {
    "y": sum_states(model_unrestricted["states"], ["infectious", "countryB", "virM"])
    / 10 ** 6,
    "label": "Infectious  mutant",
    "color": col2,
    "linestyle": "--",
}
y4 = {"first": y41, "second": y42}

y51 = {
    "y": sum_states(model_current["states"], ["infectious", "countryA", "virW"])
    / 10 ** 6,
    "label": "mRNA",
    "color": col1,
    "linestyle": "-",
}
y52 = {
    "y": sum_states(model_current["states"], ["infectious", "countryA", "virM"])
    / 10 ** 6,
    "label": "Vector",
    "color": col2,
    "linestyle": "--",
}
y5 = {"first": y51, "second": y52}

y61 = {
    "y": sum_states(model_current["states"], ["infectious", "countryB", "virW"])
    / 10 ** 6,
    "label": "mRNA",
    "color": col1,
    "linestyle": "-",
}
y62 = {
    "y": sum_states(model_current["states"], ["infectious", "countryB", "virM"])
    / 10 ** 6,
    "label": "Vector",
    "color": col2,
    "linestyle": "--",
}
y6 = {"first": y61, "second": y62}

y = [y1, y2, y3, y4, y5, y6]


path_vaccine = plot_path + "infectious"
plot_gridspec(
    y=y,
    time=np.array(time_points) / 7.0,
    title=f"Number of infectious individuals in millions",
    legend_next_to_plot=True,
    legend_location="lower center",
    ylim=[0, 27],
    xlim=None,
    bbox=(0.515, -0.05),
    line_thickness=0.6,
    ylabel=["Infectious", "Infectious"],
    xlabel="Weeks",
    row_labels=["Country A", "Country B"],
    title_col1="Pareto strategy",
    title_col2="Unrestricted strategy",
    title_col3="Current strategy",
    xticks=np.linspace(0, 20, 6),
    path=path_vaccine,
)


# ---------------------------------------------
# Plot absolute number of deceased individuals
color = "steelblue"
tw11 = {
    "y": (
        model_optimal["states"]["dead_countryA_vac2_virW"]
        + model_optimal["states"]["dead_countryA_vac2_virM"]
        + model_optimal["states"]["dead_countryA_vac1_virW"]
        + model_optimal["states"]["dead_countryA_vac1_virM"]
        + model_optimal["states"]["dead_countryA_vac0_virW"]
        + model_optimal["states"]["dead_countryA_vac0_virM"]
    )
    / 10 ** 6,
    "label": "Deceased",
    "color": color,
    "linestyle": "dotted",
    "ylim": [0, 2.0],
    "ylabel": "",
}
tw1 = {"first": tw11}

tw21 = {
    "y": (
        model_unrestricted["states"]["dead_countryA_vac2_virW"]
        + model_unrestricted["states"]["dead_countryA_vac2_virM"]
        + model_unrestricted["states"]["dead_countryA_vac1_virW"]
        + model_unrestricted["states"]["dead_countryA_vac1_virM"]
        + model_unrestricted["states"]["dead_countryA_vac0_virW"]
        + model_unrestricted["states"]["dead_countryA_vac0_virM"]
    )
    / 10 ** 6,
    "label": "Deceased",
    "color": color,
    "linestyle": "dotted",
    "ylim": [0, 2.0],
    "ylabel": "",
}
tw2 = {"first": tw21}

tw31 = {
    "y": (
        model_optimal["states"]["dead_countryB_vac2_virW"]
        + model_optimal["states"]["dead_countryB_vac2_virM"]
        + model_optimal["states"]["dead_countryB_vac1_virW"]
        + model_optimal["states"]["dead_countryB_vac1_virM"]
        + model_optimal["states"]["dead_countryB_vac0_virW"]
        + model_optimal["states"]["dead_countryB_vac0_virM"]
    )
    / 10 ** 6,
    "label": "Deceased",
    "color": color,
    "linestyle": "dotted",
    "ylim": [0, 2.0],
    "ylabel": "",
}
tw3 = {"first": tw31}

tw41 = {
    "y": (
        model_unrestricted["states"]["dead_countryB_vac2_virW"]
        + model_unrestricted["states"]["dead_countryB_vac2_virM"]
        + model_unrestricted["states"]["dead_countryB_vac1_virW"]
        + model_unrestricted["states"]["dead_countryB_vac1_virM"]
        * model_unrestricted["states"]["dead_countryB_vac0_virW"]
        + model_unrestricted["states"]["dead_countryB_vac0_virM"]
    )
    / 10 ** 6,
    "label": "Deceased",
    "color": color,
    "linestyle": "dotted",
    "ylim": [0, 2.0],
    "ylabel": "",
}
tw4 = {"first": tw41}

tw51 = {
    "y": (
        model_current["states"]["dead_countryA_vac2_virW"]
        + model_current["states"]["dead_countryA_vac2_virM"]
        + model_current["states"]["dead_countryA_vac1_virW"]
        + model_current["states"]["dead_countryA_vac1_virM"]
        + model_current["states"]["dead_countryA_vac0_virW"]
        + model_current["states"]["dead_countryA_vac0_virM"]
    )
    / 10 ** 6,
    "label": "Deceased",
    "color": color,
    "linestyle": "dotted",
    "ylim": [0, 2.0],
    "ylabel": "Deceased",
}
tw5 = {"first": tw51}

tw61 = {
    "y": (
        model_current["states"]["dead_countryB_vac2_virW"]
        + model_current["states"]["dead_countryB_vac2_virM"]
        + model_current["states"]["dead_countryB_vac1_virW"]
        + model_current["states"]["dead_countryB_vac1_virM"]
        + model_current["states"]["dead_countryB_vac0_virW"]
        + model_current["states"]["dead_countryB_vac0_virM"]
    )
    / 10 ** 6,
    "label": "Deceased",
    "color": color,
    "linestyle": "dotted",
    "ylim": [0, 2.0],
    "ylabel": "Deceased",
}
tw6 = {"first": tw61}

twin_axes = [tw1, tw2, tw3, tw4, tw5, tw6]


path_vaccine = plot_path + "infectious_dead"
plot_gridspec(
    y=y,
    time=np.array(time_points) / 7.0,
    title=f"Number of infectious and deceased individuals in millions",
    legend_next_to_plot=True,
    legend_location="lower center",
    ylim=[0, 27],
    xlim=None,
    bbox=(0.515, -0.05),
    line_thickness=0.6,
    ylabel=["Infectious", "Infectious"],
    xlabel="Weeks",
    row_labels=["Country A", "Country B"],
    title_col1="Pareto strategy",
    title_col2="Unrestricted strategy",
    title_col3="Current strategy",
    xticks=np.linspace(0, 20, 6),
    path=path_vaccine,
    twin_axes=twin_axes,
)


# create waterfall plot
fvalues1 = df_optimal_results["fval"][0:20]
fvalues2 = df_unrestricted_results["fval"][0:20]

fval1_rel = fvalues1 / fvalues1[0]
fval2_rel = fvalues2 / fvalues2[0]

fig = plt.figure()
gs = fig.add_gridspec(1, 2, hspace=0.1, wspace=0.1)
(ax1, ax2) = gs.subplots(sharex="col", sharey="row")

relative = [fval1_rel, fval2_rel]
axes = [ax1, ax2]
color = ["steelblue", "coral"]
label = ["Pareto optimal", "Unrestricted optimal"]

for index in range(len(axes)):
    ax = axes[index]
    ax.spines["top"].set_alpha(0)
    ax.spines["bottom"].set_alpha(0)
    ax.spines["right"].set_alpha(0)
    ax.spines["left"].set_alpha(0)
    ax.grid(alpha=0.6)
    ax.set_xticks(np.linspace(0, 20, 5))
    ax.set_xlim([-2, 20])
    if index == 1:
        ax.yaxis.set_ticks_position("none")
    if index == 0:
        ax.set_ylabel("Relative function value")
    length = len(relative[index])
    ax.scatter(
        np.linspace(0, length - 1, length),
        relative[index],
        color=color[index],
        s=25,
        label=label[index],
    )
    ax.set_xlabel("Optimizer run")

fig.suptitle("Ordered multi-start function values relative to best value")

bbox = (0.515, -0.07)
legend_location = "lower center"
handles, labels = ax.get_legend_handles_labels()
# legend = ax.legend(handles, labels, bbox_to_anchor = bbox, loc=legend_location, ncol = 2, facecolor='white', framealpha=1)
fig.legend(loc=legend_location, bbox_to_anchor=bbox, ncol=2, frameon=False)
path_wf = plot_path + "waterfall"
fig.savefig(path_wf, bbox_inches="tight")


# -----------------------------------------------------------------------------
var_interest = "dead"
# Deceased plots
unrestricted_deceased_A = get_sum_of_states(
    model,
    model_unrestricted,
    state_type=[var_interest, "countryA"],
    final_amount=True,
)
unrestricted_deceased_B = get_sum_of_states(
    model,
    model_unrestricted,
    state_type=[var_interest, "countryB"],
    final_amount=True,
)
unrestricted_deceased = unrestricted_deceased_A + unrestricted_deceased_B


current_deceased_A = get_sum_of_states(
    model, model_current, state_type=[var_interest, "countryA"], final_amount=True
)
current_deceased_B = get_sum_of_states(
    model, model_current, state_type=[var_interest, "countryB"], final_amount=True
)
current_deceased = current_deceased_A + current_deceased_B

dead_A_current = get_sum_of_states(
    model, model_current, state_type=[var_interest, "countryA"], final_amount=True
)
dead_B_current = get_sum_of_states(
    model, model_current, state_type=[var_interest, "countryB"], final_amount=True
)

dead_A_optimal = get_sum_of_states(
    model, model_optimal, state_type=[var_interest, "countryA"], final_amount=True
)
dead_B_optimal = get_sum_of_states(
    model, model_optimal, state_type=[var_interest, "countryB"], final_amount=True
)

optimal_deceased = get_sum_of_states(
    model, model_optimal, state_type=[var_interest], final_amount=True
)


# deviations from current policy
dev_pareto_A = (dead_A_optimal / current_deceased_A - 1) * 100
dev_pareto_B = (dead_B_optimal / current_deceased_B - 1) * 100
dev_pareto_AB = (optimal_deceased / current_deceased - 1) * 100

dev_unr_A = (unrestricted_deceased_A / current_deceased_A - 1) * 100
dev_unr_B = (unrestricted_deceased_B / current_deceased_B - 1) * 100
dev_unr_AB = (unrestricted_deceased / current_deceased - 1) * 100

c1 = [dev_pareto_AB, dev_pareto_A, dev_pareto_B]
c2 = [dev_unr_AB, dev_unr_A, dev_unr_B]

data = [c1, c2]

barWidth = 0.2
X = np.array([0, 0.5, 1])
# [left, bottom, width, height]
fig = plt.figure()
ax = fig.add_axes([0, 0, 1, 1])
ax.bar(
    X + 0.00,
    data[0],
    color="steelblue",
    width=barWidth,
    label="Pareto optimal strategy",
)
ax.bar(X + 0.2, data[1], color="coral", width=barWidth, label="Unrestricted strategy")
ax.spines["top"].set_alpha(0)
ax.spines["bottom"].set_alpha(0)
ax.spines["right"].set_alpha(0)
ax.spines["left"].set_alpha(0)
ax.grid(alpha=0.3)
# ax.set_xticklabels([])
ax.xaxis.set_ticks([0.1, 0.6, 1.1])
ax.xaxis.set_ticklabels(["Total", "Country A", "Country B"])
ax.set_ylabel("Deviation in \%")
ax.set_title("Deceased individuals relative deviation to the curent strategy")
ax.hlines(0, -0.20, 1.4, color="black", linewidth=0.5)
ax.set_xlim(-0.1, 1.3)
bbox = (0.515, -0.1)
legend_location = "lower center"
handles, labels = ax.get_legend_handles_labels()
# legend = ax.legend(handles, labels, bbox_to_anchor = bbox, loc=legend_location, ncol = 2, facecolor='white', framealpha=1)
fig.legend(loc=legend_location, bbox_to_anchor=bbox, ncol=2, frameon=False)
path_wf = plot_path + "percentage_deviation"
fig.savefig(path_wf, bbox_inches="tight")

if type_optim == "piecewise":

    def survey(results, category_names, title):
        """
        Parameters
        ----------
        results : dict
            A mapping from question labels to a list of answers per category.
            It is assumed all lists contain the same number of entries and that
            it matches the length of *category_names*.
        category_names : list of str
            The category labels.
        """
        labels = list(results.keys())
        data = np.array(list(results.values()))
        data_cum = data.cumsum(axis=1)

        category_colors = [[0.184, 0.49, 0.333, 1], [0.70, 0.13, 0.13, 1]]

        fig, ax = plt.subplots(figsize=(9.2, 3))
        ax.invert_yaxis()
        ax.xaxis.set_visible(True)
        ax.set_xlim(0, np.sum(data, axis=1).max() + 0.01)
        ax.set_ylim(-0.5, 2.5)

        for i, (colname, color) in enumerate(zip(category_names, category_colors)):

            widths = data[:, i]
            starts = data_cum[:, i] - widths
            ax.barh(labels, widths, left=starts, height=0.5, label=colname, color=color)
            xcenters = starts + widths / 2

            r, g, b, _ = color
            text_color = "white" if r * g * b < 0.5 else "darkgrey"
            for y, (x, c) in enumerate(zip(xcenters, widths)):
                ax.text(
                    x,
                    y,
                    str(np.round(c, 2)),
                    ha="center",
                    va="center",
                    color=text_color,
                )
        ax.legend(
            ncol=len(category_names),
            bbox_to_anchor=(0.22, -0.13),
            loc="lower left",
            fontsize="small",
        )
        ax.spines["top"].set_alpha(0)
        ax.spines["bottom"].set_alpha(0)
        ax.spines["right"].set_alpha(0)
        ax.spines["left"].set_alpha(0)
        ax.grid(alpha=0)
        ax.set_title(title)
        end_values = np.round(
            [
                np.round(unrestricted_deceased_A / 10 ** 6, 2)
                + np.round(unrestricted_deceased_B / 10 ** 6, 2),
                np.round(dead_A_optimal / 10 ** 6, 2)
                + np.round(dead_B_optimal / 10 ** 6, 2),
                np.round(current_deceased_A / 10 ** 6, 2)
                + np.round(current_deceased_B / 10 ** 6, 2),
            ],
            2,
        )
        ax.xaxis.set_ticks(end_values)
        ax.vlines(
            end_values[0], ymin=-0.5, ymax=0.25, color="dimgrey", linestyle="dotted"
        )
        ax.vlines(
            end_values[1], ymin=-0.5, ymax=1.25, color="dimgrey", linestyle="dotted"
        )
        ax.vlines(
            end_values[2], ymin=-0.5, ymax=2.26, color="dimgrey", linestyle="dotted"
        )
        # ax.annotate("", xy=(0.5, 0.5), xytext=(0, 0), arrowprops={"arrowstyle":"->", "relpos" :(1, 0.5)})
        # ax.arrow(1.86, -0.4, 2.33, 1.1)
        # ax.annotate('', xy=(1.86/2.23, -0.2), xycoords='axes fraction', xytext=(1, -0.20),
        #        arrowprops=dict(arrowstyle="->", color='black'))
        # ax.text(2.01, -1.03, f"{int(np.round(1.86/2.23 - 1,2)*100)}\%")
        ax.annotate(
            "",
            xy=(0.826, 0.16),
            xycoords="axes fraction",
            xytext=(1, 0.16),
            arrowprops=dict(arrowstyle="->", color="black"),
        )
        ax.text(1.95, 0.055, f"{int(np.round(1.86/2.23 - 1,2)*100)}\%")
        # ax.annotate('', xy=(2.12/2.23, -0.12), xycoords='axes fraction', xytext=(1, -0.12),
        #        arrowprops=dict(arrowstyle="->", color='black'))
        # ax.text(2.15, -0.8, f"{int(np.round(2.12/2.23 - 1,2)*100)}\%")
        ax.annotate(
            "",
            xy=(0.94, 0.49),
            xycoords="axes fraction",
            xytext=(1, 0.49),
            arrowprops=dict(arrowstyle="->", color="black"),
        )
        ax.text(2.138, 1.043, f"{int(np.round(2.12/2.23 - 1,2)*100)}\%")

        return fig, ax

    category_names = ["Deceased country A", "Deceased country B"]
    results = {
        "Unrestricted strategy": np.round(
            [unrestricted_deceased_A / 10 ** 6, unrestricted_deceased_B / 10 ** 6], 2
        ),
        "Pareto strategy": np.round(
            [dead_A_optimal / 10 ** 6, dead_B_optimal / 10 ** 6], 2
        ),
        "Current stretagy": np.round(
            [current_deceased_A / 10 ** 6, current_deceased_B / 10 ** 6], 2
        ),
    }

    fig, ax = survey(
        results,
        category_names,
        f"Number of deceased individuals in millions",
    )
    path_wf = plot_path + "percentage_deviation"
    fig.savefig(path_wf, bbox_inches="tight")
else:

    def survey(results, category_names, title):
        """
        Parameters
        ----------
        results : dict
            A mapping from question labels to a list of answers per category.
            It is assumed all lists contain the same number of entries and that
            it matches the length of *category_names*.
        category_names : list of str
            The category labels.
        """
        labels = list(results.keys())
        data = np.array(list(results.values()))
        data_cum = data.cumsum(axis=1)
        category_colors = [[0.184, 0.49, 0.333, 1], [0.70, 0.13, 0.13, 1]]

        fig, ax = plt.subplots(figsize=(9.2, 3))
        ax.invert_yaxis()
        ax.xaxis.set_visible(True)
        ax.set_xlim(0, np.sum(data, axis=1).max() + 0.01)
        ax.set_ylim(-0.5, 2.5)

        for i, (colname, color) in enumerate(zip(category_names, category_colors)):
            widths = data[:, i]
            starts = data_cum[:, i] - widths
            ax.barh(labels, widths, left=starts, height=0.5, label=colname, color=color)
            xcenters = starts + widths / 2

            r, g, b, _ = color
            text_color = "white" if r * g * b < 0.5 else "darkgrey"
            for y, (x, c) in enumerate(zip(xcenters, widths)):
                ax.text(
                    x,
                    y,
                    str(np.round(c, 2)),
                    ha="center",
                    va="center",
                    color=text_color,
                )
        ax.legend(
            ncol=len(category_names),
            bbox_to_anchor=(0.22, -0.13),
            loc="lower left",
            fontsize="small",
        )
        ax.spines["top"].set_alpha(0)
        ax.spines["bottom"].set_alpha(0)
        ax.spines["right"].set_alpha(0)
        ax.spines["left"].set_alpha(0)
        ax.grid(alpha=0)
        ax.set_title(title)
        end_values = np.round(
            [
                np.round(unrestricted_deceased_A / 10 ** 6, 2)
                + np.round(unrestricted_deceased_B / 10 ** 6, 2),
                np.round(dead_A_optimal / 10 ** 6, 2)
                + np.round(dead_B_optimal / 10 ** 6, 2),
                np.round(current_deceased_A / 10 ** 6, 2)
                + np.round(current_deceased_B / 10 ** 6, 2),
            ],
            2,
        )
        ax.xaxis.set_ticks(end_values)
        ax.vlines(
            end_values[0], ymin=-0.5, ymax=0.25, color="dimgrey", linestyle="dotted"
        )
        ax.vlines(
            end_values[1], ymin=-0.5, ymax=1.25, color="dimgrey", linestyle="dotted"
        )
        ax.vlines(
            end_values[2], ymin=-0.5, ymax=2.26, color="dimgrey", linestyle="dotted"
        )
        # ax.annotate("", xy=(0.5, 0.5), xytext=(0, 0), arrowprops={"arrowstyle":"->", "relpos" :(1, 0.5)})
        # ax.arrow(1.86, -0.4, 2.33, 1.1)
        # ax.annotate('', xy=(1.86/2.23, -0.2), xycoords='axes fraction', xytext=(1, -0.20),
        #        arrowprops=dict(arrowstyle="->", color='black'))
        # ax.text(2.01, -1.03, f"{int(np.round(1.86/2.23 - 1,2)*100)}\%")
        ax.annotate(
            "",
            xy=(0.826, 0.16),
            xycoords="axes fraction",
            xytext=(1, 0.16),
            arrowprops=dict(arrowstyle="->", color="black"),
        )
        ax.text(1.95, 0.055, f"{int(np.round(1.86/2.23 - 1,2)*100)}\%")
        # ax.annotate('', xy=(2.12/2.23, -0.12), xycoords='axes fraction', xytext=(1, -0.12),
        #        arrowprops=dict(arrowstyle="->", color='black'))
        # ax.text(2.15, -0.8, f"{int(np.round(2.12/2.23 - 1,2)*100)}\%")
        ax.annotate(
            "",
            xy=(0.94, 0.49),
            xycoords="axes fraction",
            xytext=(1, 0.49),
            arrowprops=dict(arrowstyle="->", color="black"),
        )
        ax.text(2.138, 1.043, f"{int(np.round(2.12/2.23 - 1,2)*100)}\%")

        return fig, ax

    category_names = ["Deceased country A", "Deceased country B"]
    results = {
        "Unrestricted strategy": np.round(
            [unrestricted_deceased_A / 10 ** 6, unrestricted_deceased_B / 10 ** 6], 2
        ),
        "Pareto strategy": np.round(
            [dead_A_optimal / 10 ** 6, dead_B_optimal / 10 ** 6], 2
        ),
        "Current stretagy": np.round(
            [current_deceased_A / 10 ** 6, current_deceased_B / 10 ** 6], 2
        ),
    }

    fig, ax = survey(
        results,
        category_names,
        f"Number of deceased individuals in millions",
    )
    path_wf = plot_path + "percentage_deviation"
    fig.savefig(path_wf, bbox_inches="tight")
