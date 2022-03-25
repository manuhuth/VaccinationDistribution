import numpy as np
from scipy.interpolate import CubicHermiteSpline

import matplotlib.pyplot as plt

plt.style.use("default")
# plt.style.use('ggplot')



# piecewise linear
np.random.seed(1234)
pl = np.random.uniform(0, 1, 10)
pl_grid = np.repeat(pl, 6000 / 10)

# spline

x = np.linspace(0, 20, 21)
y = np.random.uniform(-5, 5.2, 21)

y_odds = 1 / (1+ np.exp(-y))

y_df = np.gradient(y)
spline = CubicHermiteSpline(x, y, y_df)
x_grid = np.linspace(0, 20, 6000)
spline_values = spline(x_grid)
fractions = 1.0 / (1.0 + np.exp(-np.array(spline_values)))

y_plot = [pl_grid, spline_values, fractions]
label = ["$f_l$", "Spline value", "$f_l$"]
title = [
    "Stepwise function",
    "Cubic Hermite spline",
    "Logistically transformed cubic Hermite spline",
]
color = ["coral", "steelblue", "steelblue"]


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
    #ax.spines["bottom"].set_alpha(0)
    ax.spines["right"].set_alpha(0)
    #ax.spines["left"].set_alpha(0)
    #ax.grid(alpha=0.6)
    ax.set_xticks(np.linspace(0, 20, 6))
    ax.set_title(title[index])
    if index == 1:
        # ax.set_yticklabels([])
        ax.yaxis.set_ticks_position("none")
    if index in [0, 2]:
        ax.set_ylim(0, 1)

fig.tight_layout()
path_wf = (
    "/home/manuel/Documents/VaccinationDistribution/paper/images/" + "example_methods"
)
fig.savefig(path_wf, bbox_inches="tight")


import numpy as np
from scipy.interpolate import CubicHermiteSpline

import matplotlib.pyplot as plt
import matplotlib


size = 19
font = {"family": "normal", "weight": "normal", "size": size}
matplotlib.rc("font", **font)
matplotlib.rcParams['xtick.labelsize'] = size - 2 
matplotlib.rcParams['ytick.labelsize'] = size - 2 


fig, ax_g = plt.subplots(1, figsize = (8,3))
ax1 = ax_g
axes = [ax1, ax1, ax1]
for index in [2]:
    ax = axes[index]
    ax.plot(x_grid[0:3000], y_plot[index][0:3000], color=color[index])
    
    ax.scatter(x[0:11], y_odds[0:11] , color=color[index])
    ax.set_ylabel("Share of vaccine")
    ax.set_xlabel("Time")
    ax.spines["top"].set_alpha(0)
    #ax.spines["bottom"].set_alpha(0)
    ax.spines["right"].set_alpha(0)
    #ax.spines["left"].set_alpha(0)
    #ax.grid(alpha=0.6)
    ax.set_xticks(np.linspace(0, 10, 11))
    ax.set_title("")
    ax.xaxis.set_ticklabels([])
    if index == 1:
        # ax.set_yticklabels([])
        ax.xaxis.set_ticks_position("none")
    if index in [0, 2]:
        ax.set_ylim(0, 1)

for index in range(len(x[0:11])):
    ax.plot([x[index], x[index]], [0, y_odds[index]], color = "darkgrey", linestyle = "dotted" )

fig.tight_layout()

path_wf = (
    "/home/manuel/Documents/VaccinationDistribution/paper/images/" + "vacc_pathways"
)
fig.savefig(path_wf, bbox_inches="tight")