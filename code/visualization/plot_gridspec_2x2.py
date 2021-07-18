import pickle

import numpy as np
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


from functions.run_sbml import get_model_and_solver_from_sbml

type_optim = "splines"

# ----------Load optimization data---------------------------------------------
load = "/home/manuel/Documents/VaccinationDistribution/code/objects/output_" + type_optim + ".pkl"
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
plot_path = "/home/manuel/Documents/VaccinationDistribution/paper/images/" + type_method + "_"

model_solver = model_and_solver = get_model_and_solver_from_sbml(
    path_sbml=path_sbml,
    model_name=model_name,
    model_directory=model_directory,
    only_import=True,
)
model = model_solver["model"]
observables = model.getObservableNames()


def plot_gridspec(y, time, title = "Amount of vaccine doses (in mil.)",
                    legend_next_to_plot = True,
                    legend_location = 'lower center',
                    ylim = [0, 1.4],
                    xlim = None,
                    bbox = (0.515, -0.05),
                    line_thickness = 0.6,
                    ylabel = ["Vaccine doses", "Vaccine_doses"],
                    xlabel = "Weeks",
                    row_labels = ["Country A", "Country B"], 
                    title_col1 = "Pareto strategy",
                    title_col2 = "Unrestricted strategy",
                    xticks = np.linspace(0, 20, 6),
                    path = None,
                    twin_axes = None):
    
    cols = [title_col1, title_col2]
    
    fig = plt.figure()
    gs = fig.add_gridspec(2, 2, hspace=0.1, wspace=0.1)
    (ax1, ax2), (ax3, ax4) = gs.subplots(sharex='col', sharey='row')
    fig.suptitle(title)
    
    axes = [ax1, ax2, ax3, ax4]
    
    for index in range(len(axes)):
        ax = axes[index]
        y_i = y[index] 
        for index2 in y_i.keys():
            y_plot = y_i[index2]
            ax.plot(time,
                    y_plot["y"],
                    label=y_plot["label"],
                    color=y_plot["color"],
                    linestyle = y_plot["linestyle"]
                    )
    
    
    ax1.spines["top"].set_alpha(0)
    ax1.spines["bottom"].set_alpha(0)
    ax1.spines["right"].set_alpha(0)
    ax1.spines["left"].set_alpha(0)
    ax1.grid(alpha=line_thickness)
    ax1.set_ylabel(ylabel[0])
    #ax1.set_title(title_col1)

    
    
    ax2.spines["top"].set_alpha(0)
    ax2.spines["bottom"].set_alpha(0)
    ax2.spines["right"].set_alpha(0)
    ax2.spines["left"].set_alpha(0)
    ax2.grid(alpha=line_thickness)

    #ax2.set_title(title_col2)
    
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
    
    if ylim is not None:
        for index in axes:
            index.set_ylim(ylim)
    
    if xlim is not None:
        for index in axes:
            index.set_xlim(xlim)
            
    if twin_axes is not None:
        ax1_twin=ax1.twinx()
        ax2_twin=ax2.twinx()
        ax3_twin=ax3.twinx()
        ax4_twin=ax4.twinx()
        twins = [ax1_twin, ax2_twin, ax3_twin, ax4_twin]
        
        for index in range(len(twins)):
            ax_twin = twins[index]
            twin_y_i = twin_axes[index] 
            for index2 in twin_y_i.keys():
                twin_y_plot = twin_y_i[index2]
                ax_twin.plot(time,
                        twin_y_plot["y"],
                        label=twin_y_plot["label"],
                        color=twin_y_plot["color"],
                        linestyle = twin_y_plot["linestyle"]
                        )
                ax_twin.spines["top"].set_alpha(0)
                ax_twin.spines["bottom"].set_alpha(0)
                ax_twin.spines["right"].set_alpha(0)
                ax_twin.spines["left"].set_alpha(0)
                ax_twin.xaxis.set_ticks_position("none")
                ax_twin.set_ylabel(twin_y_plot["ylabel"])
                
                ax_twin.set_ylim(twin_y_plot["ylim"])
                if index in [0,2]:
                    ax_twin.set_yticklabels([])
                    ax_twin.yaxis.set_ticks_position("none")

        
    #ax4.xaxis.set_ticks_position("none")
    ax4.yaxis.set_ticks_position("none")
    #ax3.xaxis.set_ticks_position("none")
    #ax3.yaxis.set_ticks_position("none")
    ax2.xaxis.set_ticks_position("none")
    ax2.yaxis.set_ticks_position("none")
    ax1.xaxis.set_ticks_position("none")
    #ax1.yaxis.set_ticks_position("none")
    pad = 2
    ax_cols = [ax1, ax2]
    for ax, col in zip(ax_cols, cols):
        ax.annotate(col, xy=(0.5, 1.05), xytext=(0, pad),
                    xycoords='axes fraction', textcoords='offset points',
                    size='large', ha='center', va='baseline')
    ax_rows = [ax1, ax3]    
    for ax, row in zip(ax_rows, row_labels):
        ax.annotate(row, xy=(-0.05, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                    xycoords=ax.yaxis.label, textcoords='offset points',
                    size='large', ha='right', va='center')

    
    if legend_next_to_plot is True:
        if twin_axes is not None:
            handles1, labels1 = ax4.get_legend_handles_labels()
            handles_twin, labels_twin = ax4_twin.get_legend_handles_labels()
            handles = handles1 + handles_twin
            labels = labels1 + labels_twin
        else:
            handles, labels = ax4.get_legend_handles_labels()
        legend = fig.legend(handles, labels, bbox_to_anchor = bbox, loc=legend_location, ncol = 3, facecolor='white', framealpha=1)
        frame = legend.get_frame()
        frame.set_linewidth(0)
    
    plt.tight_layout()
    
    if path is not None:
        fig.savefig(path, bbox_inches="tight")
    
    return fig, axes




#Plot absolute vaccine doses
y11 = {"y" : model_optimal["observables"]["observable_quantity_countryA_vac1"]/10**6, "label" : "mRNA", "color" : "steelblue", "linestyle" : "-"}
y12 = {"y" : model_optimal["observables"]["observable_quantity_countryA_vac2"]/10**6, "label" : "Vector", "color" : "coral", "linestyle" : "--"}
y1 = {"first" : y11, "second" : y12}

y21 = {"y" : model_unrestricted["observables"]["observable_quantity_countryA_vac1"]/10**6, "label" : "mRNA", "color" : "steelblue", "linestyle" : "-"}
y22 = {"y" : model_unrestricted["observables"]["observable_quantity_countryA_vac2"]/10**6, "label" : "Vector", "color" :"coral", "linestyle" : "--"}
y2 = {"first" : y21, "second" : y22}

y31 = {"y" : model_optimal["observables"]["observable_quantity_countryB_vac1"]/10**6, "label" : "mRNA", "color" : "steelblue", "linestyle" : "-"}
y32 = {"y" : model_optimal["observables"]["observable_quantity_countryB_vac2"]/10**6, "label" : "Vector", "color" : "coral", "linestyle" : "--"}
y3 = {"first" : y31, "second" : y32}

y41 = {"y" : model_unrestricted["observables"]["observable_quantity_countryB_vac1"]/10**6, "label" : "mRNA", "color" : "steelblue", "linestyle" : "-"}
y42 = {"y" : model_unrestricted["observables"]["observable_quantity_countryB_vac2"]/10**6, "label" : "Vector", "color" :"coral", "linestyle" : "--"}
y4 = {"first" : y41, "second" : y42}

y = [y1, y2, y3, y4]



path_vaccine = plot_path +  "vaccine_total_quantity"
plot_gridspec(y=y, time=np.array(time_points)/7.0, title = "Amount of vaccine doses (in mil.)",
                    legend_next_to_plot = True,
                    legend_location = 'lower center',
                    ylim = [0, 1.4],
                    xlim = None,
                    bbox = (0.515, -0.05),
                    line_thickness = 0.6,
                    ylabel = ["Vaccine doses", "Vaccine doses"],
                    xlabel = "Weeks",
                    row_labels = ["Country A", "Country B"], 
                    title_col1 = "Pareto strategy",
                    title_col2 = "Unrestricted strategy",
                    xticks = np.linspace(0, 20, 6),
                    path = path_vaccine)

#----------------------------------------------
#Plot relative vaccines
y11 = {"y" : model_optimal["observables"]["proportion_countryA_vac1"], "label" : "mRNA", "color" : "steelblue", "linestyle" : "-"}
y12 = {"y" : model_optimal["observables"]["proportion_countryA_vac2"], "label" : "Vector", "color" : "coral", "linestyle" : "--"}
y1 = {"first" : y11, "second" : y12}

y21 = {"y" : model_unrestricted["observables"]["proportion_countryA_vac1"], "label" : "mRNA", "color" : "steelblue", "linestyle" : "-"}
y22 = {"y" : model_unrestricted["observables"]["proportion_countryA_vac2"], "label" : "Vector", "color" :"coral", "linestyle" : "--"}
y2 = {"first" : y21, "second" : y22}

y31 = {"y" : model_optimal["observables"]["proportion_countryB_vac1"], "label" : "mRNA", "color" : "steelblue", "linestyle" : "-"}
y32 = {"y" : model_optimal["observables"]["proportion_countryB_vac2"], "label" : "Vector", "color" : "coral", "linestyle" : "--"}
y3 = {"first" : y31, "second" : y32}

y41 = {"y" : model_unrestricted["observables"]["proportion_countryB_vac1"], "label" : "mRNA", "color" : "steelblue", "linestyle" : "-"}
y42 = {"y" : model_unrestricted["observables"]["proportion_countryB_vac2"], "label" : "Vector", "color" :"coral", "linestyle" : "--"}
y4 = {"first" : y41, "second" : y42}

y = [y1, y2, y3, y4]



path_vaccine = plot_path +  "vaccine_fractions"
plot_gridspec(y=y, time=np.array(time_points)/7.0, title = "Fraction of vaccine doses",
                    legend_next_to_plot = True,
                    legend_location = 'lower center',
                    ylim = [0, 1.1],
                    xlim = None,
                    bbox = (0.515, -0.05),
                    line_thickness = 0.6,
                    ylabel = ["Fraction", "Fraction"],
                    xlabel = "Weeks",
                    row_labels = ["Country A", "Country B"], 
                    title_col1 = "Pareto strategy",
                    title_col2 = "Unrestricted strategy",
                    xticks = np.linspace(0, 20, 6),
                    path = path_vaccine)


#----------------------------------------------
#Plot absolute number of infectious individuals
col1 ="seagreen"
col2 = "firebrick"
y11 = {"y" : model_optimal["states"]["infectious_countryA_vac0_virW"]/10**6, "label" : "Infectious wild type", "color" : col1, "linestyle" : "-"}
y12 = {"y" : model_optimal["states"]["infectious_countryA_vac0_virM"]/10**6, "label" : "Infectious  mutant", "color" : col2, "linestyle" : "--"}
y1 = {"first" : y11, "second" : y12}

y21 = {"y" : model_unrestricted["states"]["infectious_countryA_vac0_virW"]/10**6, "label" : "Infectious  wild type", "color" :col1, "linestyle" : "-"}
y22 = {"y" : model_unrestricted["states"]["infectious_countryA_vac0_virM"]/10**6, "label" : "Infectious  mutant", "color" :col2, "linestyle" : "--"}
y2 = {"first" : y21, "second" : y22}

y31 = {"y" : model_optimal["states"]["infectious_countryB_vac0_virW"]/10**6, "label" : "Infectious wild type", "color" : col1, "linestyle" : "-"}
y32 = {"y" : model_optimal["states"]["infectious_countryB_vac0_virM"]/10**6, "label" : "Infectious mutant", "color" : col2, "linestyle" : "--"}
y3 = {"first" : y31, "second" : y32}

y41 = {"y" : model_unrestricted["states"]["infectious_countryB_vac0_virW"]/10**6, "label" : "Infectious  wild type", "color" : col1, "linestyle" : "-"}
y42 = {"y" : model_unrestricted["states"]["infectious_countryB_vac0_virM"]/10**6, "label" : "Infectious  mutant", "color" :col2, "linestyle" : "--"}
y4 = {"first" : y41, "second" : y42}

y = [y1, y2, y3, y4]


path_vaccine = plot_path +  "infectious"
plot_gridspec(y=y, time=np.array(time_points)/7.0, title = "Number of infectious individuals (in mil.)",
                    legend_next_to_plot = True,
                    legend_location = 'lower center',
                    ylim = [0,25],
                    xlim = None,
                    bbox = (0.515, -0.05),
                    line_thickness = 0.6,
                    ylabel = ["Infectious", "Infectious"],
                    xlabel = "Weeks",
                    row_labels = ["Country A", "Country B"], 
                    title_col1 = "Pareto strategy",
                    title_col2 = "Unrestricted strategy",
                    xticks = np.linspace(0, 20, 6),
                    path = path_vaccine)



#---------------------------------------------
#Plot absolute number of deceased individuals 
color ="steelblue"
tw11 = {"y" : (model_optimal["states"]["dead_countryA_vac2_virW"] + model_optimal["states"]["dead_countryA_vac2_virM"]+model_optimal["states"]["dead_countryA_vac1_virW"] + model_optimal["states"]["dead_countryA_vac1_virM"]+model_optimal["states"]["dead_countryA_vac0_virW"] + model_optimal["states"]["dead_countryA_vac0_virM"])/10**6,
        "label" : "Deceased", "color" : color, "linestyle" : "dotted", "ylim" : [0, 2.0], "ylabel" : ""}
tw1 = {"first" : tw11}

tw21 = {"y" : (model_unrestricted["states"]["dead_countryA_vac2_virW"] + model_unrestricted["states"]["dead_countryA_vac2_virM"]+model_unrestricted["states"]["dead_countryA_vac1_virW"] + model_unrestricted["states"]["dead_countryA_vac1_virM"]+model_unrestricted["states"]["dead_countryA_vac0_virW"] + model_unrestricted["states"]["dead_countryA_vac0_virM"])/10**6,
        "label" : "Deceased", "color" : color, "linestyle" : "dotted", "ylim" : [0, 2.0], "ylabel" : "Deceased"}
tw2 = {"first" : tw21}

tw31 = {"y" : (model_optimal["states"]["dead_countryB_vac2_virW"] + model_optimal["states"]["dead_countryB_vac2_virM"]+model_optimal["states"]["dead_countryB_vac1_virW"] + model_optimal["states"]["dead_countryB_vac1_virM"]+model_optimal["states"]["dead_countryB_vac0_virW"] + model_optimal["states"]["dead_countryB_vac0_virM"])/10**6, 
        "label" : "Deceased", "color" : color, "linestyle" : "dotted", "ylim" : [0, 2.0], "ylabel" : ""}
tw3 = {"first" : tw31}

tw41 = {"y" : (model_unrestricted["states"]["dead_countryB_vac2_virW"] + model_unrestricted["states"]["dead_countryB_vac2_virM"]+model_unrestricted["states"]["dead_countryB_vac1_virW"] + model_unrestricted["states"]["dead_countryB_vac1_virM"]*model_unrestricted["states"]["dead_countryB_vac0_virW"] + model_unrestricted["states"]["dead_countryB_vac0_virM"])/10**6,
        "label" : "Deceased", "color" : color, "linestyle" : "dotted", "ylim" : [0, 2.0], "ylabel" : "Deceased"}
tw4 = {"first" : tw41}

twin_axes = [tw1, tw2, tw3, tw4]


path_vaccine = plot_path +  "infectious_dead"
plot_gridspec(y=y, time=np.array(time_points)/7.0, title = "Number of infectious and deceased individuals (in mil.)",
                    legend_next_to_plot = True,
                    legend_location = 'lower center',
                    ylim = [0,25],
                    xlim = None,
                    bbox = (0.515, -0.05),
                    line_thickness = 0.6,
                    ylabel = ["Infectious", "Infectious"],
                    xlabel = "Weeks",
                    row_labels = ["Country A", "Country B"], 
                    title_col1 = "Pareto strategy",
                    title_col2 = "Unrestricted strategy",
                    xticks = np.linspace(0, 20, 6),
                    path = path_vaccine,
                    twin_axes=twin_axes)





#create waterfall plot
fvalues1 = df_optimal_results["fval"][0:20]
fvalues2 = df_unrestricted_results["fval"][0:20]

fval1_rel = fvalues1 / fvalues1[0]
fval2_rel = fvalues2 / fvalues2[0]

fig = plt.figure()
gs = fig.add_gridspec(1, 2, hspace=0.1, wspace=0.1)
(ax1, ax2) = gs.subplots(sharex='col', sharey='row')

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
    ax.set_xlim([-2,20])
    if index == 1:
        ax.yaxis.set_ticks_position("none")
    if index == 0:
        ax.set_ylabel("Relative function value")
    length = len(relative[index])
    ax.scatter(np.linspace(0, length-1, length), relative[index], color = color[index], s=25, label=label[index] )
    ax.set_xlabel("Optimizer run")
    
fig.suptitle("Ordered multi-start function values relative to best value")

bbox = (0.515, -0.07)
legend_location = 'lower center'
handles, labels = ax.get_legend_handles_labels()
#legend = ax.legend(handles, labels, bbox_to_anchor = bbox, loc=legend_location, ncol = 2, facecolor='white', framealpha=1)
fig.legend(loc=legend_location, bbox_to_anchor = bbox, ncol=2, frameon=False)
path_wf = plot_path  + "waterfall"
fig.savefig(path_wf, bbox_inches="tight")




