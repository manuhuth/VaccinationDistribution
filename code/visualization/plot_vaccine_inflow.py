import numpy as np
from functions.create_vaccine_doses_inflow import create_inflow_from_data
import matplotlib.pyplot as plt
plt.style.use('default')
#plt.style.use('ggplot')
plt.rcParams.update(
    {"text.usetex": True, "font.family": "sans-serif", "font.sans-serif": ["Helvetica"],
     "grid.color": "white"}
)
plt.rcParams['lines.linewidth'] = 1.5
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rcParams['axes.facecolor'] = "#E6E6E6"

inflow = create_inflow_from_data(number_decision_periods = 10, format_df = "wide")


fig, ax = plt.subplots()
vaccines = ["vac1", "vac2"]
label = ["mRNA", "Vector"]
color = ["steelblue", "coral"]
linestyle = ["solid", "dashed"]

for index in range(len(vaccines)):
    ax.plot(np.linspace(0, 20, 6000), np.repeat(inflow[vaccines[index]], 600)/10**6,
            color = color[index], label = label[index], linestyle=linestyle[index],)

ax.spines["top"].set_alpha(0)
ax.spines["bottom"].set_alpha(0)
ax.spines["right"].set_alpha(0)
ax.spines["left"].set_alpha(0)
ax.grid(alpha=0.6)
ax.set_xticks(np.linspace(0, 20, 6))
ax.set_xlabel("Weeks")
ax.set_ylabel("Available vaccine doses")
ax.set_title("Total available vaccine doses (in mil.)")

bbox = (0.515, -0.22)
legend_location = 'lower center'
handles, labels = ax.get_legend_handles_labels()
#legend = ax.legend(handles, labels, bbox_to_anchor = bbox, loc=legend_location, ncol = 2, facecolor='white', framealpha=1)
ax.legend(loc=legend_location, bbox_to_anchor = bbox, ncol=2, frameon=False)

path="/home/manuel/Documents/VaccinationDistribution/paper/images/available_vaccine"
fig.savefig(path, bbox_inches="tight")
          
          