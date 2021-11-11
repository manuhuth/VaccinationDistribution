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
plt.rcParams["lines.linewidth"] = 0.5
plt.rc("text", usetex=True)
plt.rcParams["text.latex.preamble"] = [r"\usepackage{cmbright}"]
plt.rcParams["axes.facecolor"] = "#E6E6E6"


load = "/home/manuel/Documents/VaccinationDistribution/code/objects/output_splines_sensitivity_distance_toy.pkl"
with open(
    load,
    "rb",
) as input:
    df = pickle.load(input)

df_sum = df

for index in ["pareto", "unrestricted", "current"]:
    A = f"A_{index}"
    B = f"B_{index}"
    df_sum[index] = df[A] + df[B]


df = df / 10 ** 6
df["distance"] = df["distance"] * 10 ** 6
df = df.sort_values("distance")
colnames = df.columns
x = 1 / (1 + np.array(df["distance"]))
line_thickness = 0.5
bbox = (0.515, -0.05)
xlabel = r"$b$"
ylabels = ["Deaths", "Deaths", "Deaths"]
titles = ["Country A", "Country B", "Total"]

legend_location = "lower center"
strs = ["A_", "B_"]

fig, ax = plt.subplots(3, 1)


for index in range(len(ax) - 1):

    ax_it = ax[index]
    ax_it.spines["top"].set_alpha(0)
    ax_it.spines["bottom"].set_alpha(0)
    ax_it.spines["right"].set_alpha(0)
    ax_it.spines["left"].set_alpha(0)
    ax_it.grid(alpha=line_thickness)
    ax_it.set_xlabel(xlabel)
    ax_it.set_ylabel(ylabels[index])

    ax_it.set_title(titles[index])

    for index_2 in df.columns:
        if strs[index] in index_2:
            if "current" in index_2:
                label = "Current strategy"
                color = "C0"
            elif "pareto" in index_2:
                label = "Pareto optimal strategy"
                color = "C1"
            else:
                label = "Unrestricted strategy"
                color = "C2"
            ax_it.plot(x, df[index_2], label=label, color=color)

    ax_it = ax[2]
    ax_it.plot(x, df["current"], color="C0")
    ax_it.plot(x, df["pareto"], color="C1")
    ax_it.plot(x, df["unrestricted"], color="C2")
    ax_it.set_title("Total")
    ax_it.set_xlabel(r"$b$")
    ax_it.set_ylabel("Deaths")
    ax_it.spines["top"].set_alpha(0)
    ax_it.spines["bottom"].set_alpha(0)
    ax_it.spines["right"].set_alpha(0)
    ax_it.spines["left"].set_alpha(0)
    ax_it.grid(alpha=line_thickness)

    handles, labels = ax[0].get_legend_handles_labels()
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
fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/sensitivity_distance_toy",
    bbox_inches="tight",
)
