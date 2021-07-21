import numpy as np
from scipy.interpolate import CubicHermiteSpline

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
plt.rcParams["text.latex.preamble"] = [r"\usepackage{amsmath}"]
plt.rcParams["axes.facecolor"] = "#E6E6E6"


# piecewise linear
pl = np.random.uniform(0, 1, 10)
pl_grid = np.repeat(pl, 6000 / 10)

# spline
np.random.seed(1234)
x = np.linspace(0, 20, 21)
y = np.random.uniform(-8, 8, 21)
y_df = np.gradient(y)
spline = CubicHermiteSpline(x, y, y_df)
x_grid = np.linspace(0, 20, 6000)
spline_values = spline(x_grid)
fractions = 1.0 / (1.0 + np.exp(-np.array(spline_values)))

y_plot = [spline_values, fractions, pl_grid]
label = ["Spline value", "$f_l$", "$f_l$"]
title = [
    "Cubic Hermite spline",
    "Logistically transformed cubic Hermite spline",
    "Stepwise function",
]
color = ["steelblue", "steelblue", "coral"]


fig, ax_g = plt.subplots(3)
ax1 = ax_g[0]
ax2 = ax_g[1]
ax3 = ax_g[2]
axes = [ax1, ax2, ax3]
for index in range(len(axes)):
    ax = axes[index]
    ax.plot(x_grid, y_plot[index], color=color[index])
    ax.set_ylabel(label[index])
    ax.set_xlabel("Weeks")
    ax.spines["top"].set_alpha(0)
    ax.spines["bottom"].set_alpha(0)
    ax.spines["right"].set_alpha(0)
    ax.spines["left"].set_alpha(0)
    ax.grid(alpha=0.6)
    ax.set_xticks(np.linspace(0, 20, 6))
    ax.set_title(title[index])
    if index == 1:
        # ax.set_yticklabels([])
        ax.yaxis.set_ticks_position("none")
    if index in [1, 2]:
        ax.set_ylim(0, 1)

fig.tight_layout()
path_wf = (
    "/home/manuel/Documents/VaccinationDistribution/paper/images/" + "example_methods"
)
fig.savefig(path_wf, bbox_inches="tight")
