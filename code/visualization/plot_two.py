import pickle
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import matplotlib.cm as cm
from functions.plot_tools import get_spline
from functions.plot_tools import plot_best_strategy
from functions.plot_tools import stacked_bar
from functions.plot_tools import plot_bars_vac, compute_deceased, compute_splines_from_results
from functions.plot_tools import compute_splines_from_results_initials

from numpy import trapz

def imscatter(x, y, image, ax=None, zoom=1):
    if ax is None:
        ax = plt.gca()
    try:
        image = plt.imread(image)
    except TypeError:
        # Likely already an array...
        pass
    im = OffsetImage(image, zoom=zoom)
    x, y = np.atleast_1d(x, y)
    artists = []
    for x0, y0 in zip(x, y):
        ab = AnnotationBbox(im, (x0, y0), xycoords='data', frameon=False)
        artists.append(ax.add_artist(ab))
    ax.update_datalim(np.column_stack([x, y]))
    ax.autoscale()
    return artists

size = 19
font = {"family": "normal", "weight": "normal", "size": size}
matplotlib.rc("font", **font)
matplotlib.rcParams['xtick.labelsize'] = size - 2 
matplotlib.rcParams['ytick.labelsize'] = size - 2 

initials = [0,1,2,3,4,5]#125, 150, 175, 200, 225, 250
vac = "Unequal" #only 
initial_states = "Unequal"  
specification = f"inital{initial_states}_vac{vac}"

results = []
for index in initials:
    path = (
        f"/home/manuel/Documents/VaccinationDistribution/code/objects/{specification}_initial_{index}.pkl"
    )
    
    with open(
          path,
          "rb",
    ) as input:
          loaded_object = pickle.load(input)
    results.append(loaded_object)
    
type_opti = "pareto_improvements"
vac = "vac2"
periods = 10
length = 14
total_length = 140
grid_points = 6000
count_plot = 65
lwidth = 2.5
y_position = 1.35
y_size = 28
x_position = -0.3

add_additional = {"integers" : [9,10],
                 "number" : 0.5}


fig = plt.figure(constrained_layout=True, figsize = (20,14))
gs = GridSpec(4, 3, figure=fig, height_ratios=[1.5, 0.4 ,2,2])
gs.update(wspace=0.04, hspace=0.035)

#data = [[10,9,8,7,6,5,4,3,2,1,0],
#        [0,1,2,3,4,5,6,7,8,9,10],
#        [0,1,2,3,4,5,6,7,8,9,10],
#        [10,9,8,7,6,5,4,3,2,1,0],
#        ]

data = [[10,9,8,7,6,5],
        [0,1,2,3,4,5],
        [0,1,2,3,4,5],
        [10,9,8,7,6,5],
        ]

ax = fig.add_subplot(gs[0, 0:3]) 
#cmap = cm.get_cmap("Greens", 14)
#hex_col = []
#for i in range(cmap.N):
#    rgba = cmap(i)
    # rgb2hex accepts rgb or rgba
#    hex_col.append(matplotlib.colors.rgb2hex(rgba))

#color_array = pd.DataFrame(columns=range(6), index = range(4))
#for index in range(np.array(data).shape[0]):
#    for index_col in range(np.array(data).shape[1]):
#        color_array.loc[index, index_col] = hex_col[np.array(data)[index, index_col]]

#row_labels=["Country A\nwild-type", "Country A\nMutant",
#            "Country B\nwild-type", "Country B\nMutant",]
#col_labels = [f"Case {i}" for i in range(11)]
#x_table = ax.table(cellText=data, rowLabels=row_labels, colLabels= col_labels, loc="center", bbox=[0, 0.13, 1, 0.8],
#                   cellLoc = "center", cellColours=np.array(color_array),)
#ax.set_title("Initial virus distributions")
#ax.axis("off")
#ax.text(
#        -0.05,
#        y_position,
#        chr(count_plot),
#        horizontalalignment="center",
#        verticalalignment="center",
#        transform=ax.transAxes,
#        weight="bold",
#        size=y_size,
#)
#ax.set_ylim([0.2, 0.9])
#count_plot += 1

ax.set_xlim(-1, 15)
ax.set_ylim(0,12)



width = 1.5
height = 3
c_y = 5
c_x = 3

for j in range(2):
    for i in range(6):
        left = 0 + c_x*i
        bottom = 0.5 + c_y*j
        left_bottom = (left, bottom)
        if j < 1:
            left_bottom_big = (left - 0.35, bottom - 0.5)
            width_big = width + 2*0.35
            height_big = height + 1 + c_y
            ax.add_patch(Rectangle(left_bottom_big, width_big,height_big, fc = "lightgray" , color="lightgray", alpha=0.6))
            ax.annotate(f"Case {i}", (left_bottom_big[0] + width_big/2, left_bottom_big[1] - 1 ),
                        ha='center', va= "center")


for j in range(2):
    for i in range(6):
        left = 0 + c_x*i
        bottom = 0.5 + c_y*j
        left_bottom = (left, bottom)

        ax.add_patch(Rectangle(left_bottom, width, height,fc="white", color="black"))
        
        if bottom > 1:
            country = "Country 1"
        else:
            country = "Country 2"
        
        y_coord = bottom + height/2    
        
        if i < 1: 
            ax.annotate(country, (left - 1.1, y_coord), ha='center', va= "center")#, rotation=90)
        



#case 1
left = 0
bottom = 0.5
delta = 0.2
xs = np.linspace(left + delta, left + width - delta, 5 )


for i in [1,2]:
    imscatter(xs, np.repeat(bottom + c_y + height * i/3,5) ,
              "/home/manuel/Documents/VaccinationDistribution/paper/images/virus_green.png",
               ax=ax, zoom = 0.07)
    imscatter(xs, np.repeat(bottom + height * i/3,5) ,
          "/home/manuel/Documents/VaccinationDistribution/paper/images/virus_blue.png",
           ax=ax, zoom = 0.07)

#case 2-4

for n in range(1, 6):
    left = 3*n
    bottom = 0.5
    delta = 0.2
    
    xs = np.linspace(left + delta, left + width - delta, 5 )
    
    
    for i in [1,2]:
        imscatter(xs, np.repeat(bottom + c_y + height * i/3,5) ,
                  "/home/manuel/Documents/VaccinationDistribution/paper/images/virus_green.png",
                   ax=ax, zoom = 0.07)
        imscatter(xs, np.repeat(bottom + height * i/3,5) ,
              "/home/manuel/Documents/VaccinationDistribution/paper/images/virus_blue.png",
               ax=ax, zoom = 0.07)
        
    
    imscatter(xs[0:n], np.repeat(bottom + height * 2/3,5)[0:n] ,
                  "/home/manuel/Documents/VaccinationDistribution/paper/images/virus_green.png",
                   ax=ax, zoom = 0.07)
    imscatter(xs[0:n], np.repeat(bottom + c_y + height * 2/3,5)[0:n] ,
              "/home/manuel/Documents/VaccinationDistribution/paper/images/virus_blue.png",
               ax=ax, zoom = 0.07)
ax.set_ylim(-1,10)
ax.yaxis.set_visible(False)
ax.xaxis.set_visible(False)
ax.spines['left'].set_visible(False)   
ax.spines['bottom'].set_visible(False)   
ax.set_title("Initial virus distributions")

ax.text(
        -0.075,
        y_position,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=y_size,
)
count_plot += 1

ax = fig.add_subplot(gs[1, 0:4])
col = "seagreen"
ax.plot([-0.35, 14.35], [0,0], color = col, linewidth= 5)
ax.plot([-0.35,-0.35], [0, 1], color=col, linewidth=3)
ax.plot([-0.35, 14.35], [1, 0], color=col, linewidth=3)
ax.fill_between([-0.35, 14.35], [0,0], [1,0], color = col, alpha = 0.4)
ax.set_xlim(-1, 15)
ax.set_ylim([0, 1])
ax.spines['left'].set_visible(False)   
ax.spines['bottom'].set_visible(False) 
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.tick_params(bottom = False, left=False) 
ax.set_xlabel("Degree of unequal initial distribution")

#--------------------------------------------------------------------------------------------




   
dec_A = compute_deceased(results, "countryA")
dec_B = compute_deceased(results, "countryB")
colors = ["steelblue", "steelblue", "steelblue"]
labels = ["Pareto optimal", "Optimal", "Pop.-based"]
linestyles = ["solid", "solid", "dashed"]
keys = list(dec_B.keys())

ax = fig.add_subplot(gs[2, 0])  
#for index in range(len(dec_B.keys())):
#    if keys[index] in ["optimal", "pop"]:
#        if keys[index] == "pop":
#            marker = "o"
#        else:
#            marker = "X"
#        ax.scatter(initials, (np.array(dec_A[keys[index]]) + np.array(dec_B[keys[index]])),
#                label = labels[index], color=colors[index],
#                linestyle=linestyles[index], s=190, marker=marker)

pop = (np.array(dec_A["pop"]) + np.array(dec_B["pop"]))
optimal = (np.array(dec_A["optimal"]) + np.array(dec_B["optimal"]))
pareto = (np.array(dec_A["pareto"]) + np.array(dec_B["pareto"]))
diffs = ((-optimal + pop) / pop) 
for index in range(len(initials)):
    ax.plot(np.repeat(initials[index],2), [0, diffs[index]], color = "steelblue",
            linewidth = lwidth, linestyle = "dotted" )
    #diff = np.round((pop[index] - optimal[index]) / pop[index] * 100,1)
    #diff_plot = (pop[index] - optimal[index])
    #ax.annotate(f"{-diff}%", (initials[index], optimal[index] - 5000), ha='center' )

ax.scatter(initials, diffs,
                color="steelblue",
                s=200)
   
ax.yaxis.grid(alpha=0.6)    
ax.set_title("Reduction in deceased individuals")

ax.set_ylabel("Optimal strategy")
#ax.set_ylim([180000, 235000])
ylim = ax.get_ylim()
ax.set_xlim(-0.5, 5.5)
col_label_sym = [f"{i}" for i in range(5)] + ["5"]
ax.set_xticklabels([])
vals = ax.get_yticks()
ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
ax.legend()


ax.text(
        x_position,
        y_position,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=y_size,
)
count_plot += 1



#----------------------------------------------------------------------------------------------------
ax = fig.add_subplot(gs[3, 0])
#for index in range(len(dec_B.keys())):
#    if keys[index] in ["pareto", "pop"]:
#        if keys[index] == "pop":
#            marker = "o"
#        else:
#            marker = "X"
#        ax.scatter(initials, (np.array(dec_A[keys[index]]) + np.array(dec_B[keys[index]])),
#                label = labels[index], color="C0",
#                linestyle=linestyles[index], s=190, marker=marker)
diffs = ((-pareto + pop) / pop) 
for index in range(len(initials)):
    ax.plot(np.repeat(initials[index],2), [0, diffs[index]], color = "steelblue",
            linewidth = lwidth, linestyle = "dotted" )
#    diff = np.round((pop[index] - pareto[index]) / pop[index] * 100,1)
#    diff_plot = (pop[index] - pareto[index])
#    ax.annotate(f"{-diff}%", (initials[index]+0.1, pareto[index] + diff_plot / 3 ) )
ax.scatter(initials, diffs,
                color="steelblue",
                s=200)
ax.yaxis.grid(alpha=0.6)    
ax.set_title("")
ax.set_xlim(-0.5, 5.5)
ax.set_ylabel("Pareto strategy")
ax.set_xlabel("Case")
ax.set_xticklabels([""] + col_label_sym )
ax.set_ylim(ylim)
vals = ax.get_yticks()
ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])

ax.legend()


#----------------------------------------------------------------------------
strategy_identifier = ["all_strategies", "pareto_improvements"]
name_strategy = ["Optimal strategy", "Pareto optimal strategy"]
colors = ["steelblue", "steelblue"]
paper_value = 30000
ylim = [-0.4, 0.79]

index = 1
ax = fig.add_subplot(gs[2, 1])
fractions = compute_splines_from_results_initials(initials=initials, type_opti = "all_strategies",
                                      vac = "vac1", results=results)


#ax.scatter(initials, fractions,
#           label = "Optimal", color="C1",
#           s=190, marker="X")
#ax.scatter(initials, np.repeat(0.5, len(initials)),
#           label = "Population\nbased", color="C1",
#           s=190, marker="o")
#for index in range(len(initials)):
#    ax.plot(np.repeat(initials[index],2), [0.5, fractions[index]], color = "C1", linewidth = lwidth )
 #   diff = np.round(-0.5 + fractions[index],2)
 #   diff_plot = (pop[index] - optimal[index])
 #   ax.annotate(f"{diff}", (initials[index]+0.1, 0.5  + diff_plot/2 ) )
diffs = np.array(fractions) / 0.5 - 1
for index in range(len(initials)):
   ax.plot(np.repeat(initials[index],2), [0, diffs[index]], color = "steelblue",
           linewidth = lwidth, linestyle = "dotted" )
#    diff = np.round((pop[index] - pareto[index]) / pop[index] * 100,1)
#    diff_plot = (pop[index] - pareto[index])
#    ax.annotate(f"{-diff}%", (initials[index]+0.1, pareto[index] + diff_plot / 3 ) )
ax.scatter(initials, diffs,
                color="steelblue",
                s=200)
ax.yaxis.grid(alpha=0.6)
ax.set_title("Change of Vaccine 1 doses\nallocated to Country 1")
ax.set_ylim(ylim)
ax.set_ylabel("Opimal strategy")
ax.legend()
cols_labels_small = [f"{i}" for i in range(len(initials))]
#ax.set_xticklabels([""] + cols_labels_small )
ax.set_xticklabels([""] + col_label_sym)
#ax.set_xlabel("Case")
vals = ax.get_yticks()
ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
ax.text(
        x_position,
        y_position,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=y_size,
)
count_plot += 1

#--------------------------------------------------------------------------------------------------------
ax = fig.add_subplot(gs[3, 1])
fractions = compute_splines_from_results_initials(initials=initials, type_opti = "pareto_improvements",
                                      vac = "vac1", results=results)

#ax.scatter(initials, fractions,
#           label = "Pareto optimal", color="C1",
#           s=190, marker="X")
#ax.scatter(initials, np.repeat(0.5, len(initials)),
#           label = "Population\nbased", color="C1",
#           s=190, marker="o")
#ax.yaxis.grid(alpha=0.6)
#for index in range(len(initials)):
#    ax.plot(np.repeat(initials[index],2), [0.5, fractions[index]], color = "C1", linewidth = lwidth )
#    diff = np.round(-0.5 + fractions[index],2)
#    ax.annotate(f"{diff}", (initials[index]+0.1, 0.5  + diff/2 ) )
diffs = np.array(fractions) / 0.5 - 1
for index in range(len(initials)):
   ax.plot(np.repeat(initials[index],2), [0, diffs[index]], color = "steelblue",
           linewidth = lwidth, linestyle = "dotted" )
#    diff = np.round((pop[index] - pareto[index]) / pop[index] * 100,1)
#    diff_plot = (pop[index] - pareto[index])
#    ax.annotate(f"{-diff}%", (initials[index]+0.1, pareto[index] + diff_plot / 3 ) )
ax.scatter(initials, diffs,
                color="steelblue",
                s=200)
ax.set_ylim(ylim)
ax.yaxis.grid(alpha=0.6)
vals = ax.get_yticks()
ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
ax.set_title("")
ax.set_ylabel("Pareto strategy")
#ax.set_ylim([0.27, 0.93])
ax.legend()
ax.set_xticklabels([""] + col_label_sym )
ax.set_xlabel("Case")
ax.text(
        x_position,
        y_position*0.82,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=y_size,
)

#--------------------------------------------------------------------------------------------------------
ax = fig.add_subplot(gs[2, 2])
fractions = compute_splines_from_results_initials(initials=initials, type_opti = "all_strategies",
                                      vac = "vac2", results=results)

#ax.scatter(initials, fractions,
#           label = "Optimal", color="C1",
#           s=190, marker="X")
#ax.scatter(initials, np.repeat(0.5, len(initials)),
#           label = "Population\nbased", color="C1",
#           s=190, marker="o")
#ax.yaxis.grid(alpha=0.6)
#for index in range(len(initials)):
#    ax.plot(np.repeat(initials[index],2), [0.5, fractions[index]], color = "C1", linewidth = lwidth )
#    diff = np.round(-0.5 + fractions[index],2)
#    ax.annotate(f"{diff}", (initials[index]+0.1, 0.5  + diff/2 ) )
diffs = np.array(fractions) / 0.5 - 1
for index in range(len(initials)):
   ax.plot(np.repeat(initials[index],2), [0, diffs[index]], color = "steelblue",
           linewidth = lwidth, linestyle = "dotted" )
#    diff = np.round((pop[index] - pareto[index]) / pop[index] * 100,1)
#    diff_plot = (pop[index] - pareto[index])
#    ax.annotate(f"{-diff}%", (initials[index]+0.1, pareto[index] + diff_plot / 3 ) )
ax.scatter(initials, diffs,
                color="steelblue",
                s=200)
ax.set_ylim(ylim)
ax.yaxis.grid(alpha=0.6)
ax.set_title("Change of Vaccine 2 doses\nallocated to Country 1")
ax.set_ylabel("")
#ax.set_xlabel("Case")
vals = ax.get_yticks()
ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
#ax.set_ylim([0.27, 0.93])
ax.set_xticklabels( [""] + cols_labels_small)
ax.legend()

#--------------------------------------------------------------------------------------------------------
ax = fig.add_subplot(gs[3, 2])
fractions = compute_splines_from_results_initials(initials=initials, type_opti = "pareto_improvements",
                                      vac = "vac2", results=results)

#ax.scatter(initials, fractions,
#           label = "Pareto optimal", color="C1",
#           s=190, marker="X")
#ax.scatter(initials, np.repeat(0.5, len(initials)),
#           label = "Population\nbased", color="C1",
#           s=190, marker="o")
diffs = np.array(fractions) / 0.5 - 1
for index in range(len(initials)):
   ax.plot(np.repeat(initials[index],2), [0, diffs[index]], color = "steelblue",
           linewidth = lwidth, linestyle = "dotted" )
#    diff = np.round((pop[index] - pareto[index]) / pop[index] * 100,1)
#    diff_plot = (pop[index] - pareto[index])
#    ax.annotate(f"{-diff}%", (initials[index]+0.1, pareto[index] + diff_plot / 3 ) )
ax.scatter(initials, diffs,
                color="steelblue",
                s=200)
ax.set_ylim(ylim)
ax.yaxis.grid(alpha=0.6)


#for index in range(len(initials)):
#    ax.plot(np.repeat(initials[index],2), [0.5, fractions[index]], color = "C1", linewidth = lwidth )
#    diff = np.round(-0.5 + fractions[index],2)
#    ax.annotate(f"{diff}", (initials[index]+0.1, 0.5  + diff/2 ) )
ax.set_title("")
vals = ax.get_yticks()
ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
#ax.set_ylim([0.27, 0.93])
#ax.legend()
ax.set_xticklabels([""] + cols_labels_small )
ax.set_xlabel("Case")



fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/plot_two",
    bbox_inches="tight",
)


#---------------------------------subplots------------------------------------
type_opti = "pareto_improvements"
vac = "vac2"
periods = 10
length = 14
total_length = 140
grid_points = 6000
count_plot = 65
lwidth = 2.5
y_position = 1.35
y_size = 28
x_position = -0.3

add_additional = {"integers" : [9,10],
                 "number" : 0.5}


fig = plt.figure(constrained_layout=True, figsize = (20,14))
gs = GridSpec(4, 3, figure=fig)
gs.update(wspace=0.04, hspace=0.035)

dict_out = results

types = ["vac1", "vac2"]
index = 0
for index_1 in range(len(types)):
    for index_2 in range(len(results)):
        x_plot = int(index % 3)
        y_plot = int(np.floor(index / 3))

        ax = fig.add_subplot(gs[y_plot, x_plot])
        
        color_vaccines = "steelblue"
        ylab = ""
        if (index % 3) == 0:
            ylab = f"Vaccine {index_1 + 1} allocated\nto Country 1"
        
        print(index in [0,6])
            
        plot_best_strategy(
                dict_output=dict_out,
                fig=fig,
                ax=ax,
                case=index_2,
                x_scatter=np.linspace(0, 140, 11) / 7,
                vac_interest=types[index_1],
                periods=10,
                length=14,
                total_length=140,
                grid_points=6000,
                col_unconstrained=color_vaccines,
                label_unconstrained="",
                col_pareto=color_vaccines,
                label_pareto="",
                xlabel="",
                ylabel=ylab,
                title=f"Case {index_2}",
                linewidth=4,
                s_scatter=0,
                label_scatter="",
                plot="optimal",
                n_vacc = 30000,
                x_total=17,
                y_total=25500,
                scale_total=10*3,
                add_additional=add_additional,
        )
        ax.yaxis.grid(alpha=0.6)
        ax.set_ylim([0,30500])
       
        if index in [0, 6]:
            ax.text(
                    x_position,
                    y_position*0.82,
                    chr(count_plot),
                    horizontalalignment="center",
                    verticalalignment="center",
                    transform=ax.transAxes,
                    weight="bold",
                    size=y_size,
                    )
            count_plot += 1
        index = index + 1

fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/plot_two_timings_optimal",
    bbox_inches="tight",
)

fig = plt.figure(constrained_layout=True, figsize = (20,14))
gs = GridSpec(4, 3, figure=fig)
gs.update(wspace=0.04, hspace=0.035)
count_plot = 65
dict_out = results


index = 0
for index_1 in range(len(types)):
    for index_2 in range(len(results)):
        x_plot = int(index % 3)
        y_plot = int(np.floor(index / 3))

        ax = fig.add_subplot(gs[y_plot, x_plot])
        
        color_vaccines = "steelblue"
        ylab = ""
        if (index % 3) == 0:
            ylab = f"Vaccine {index_1 + 1} allocated\nto Country 1"
            
            
        plot_best_strategy(
                dict_output=dict_out,
                fig=fig,
                ax=ax,
                case=index_2,
                x_scatter=np.linspace(0, 140, 11) / 7,
                vac_interest=types[index_1],
                periods=10,
                length=14,
                total_length=140,
                grid_points=6000,
                col_unconstrained=color_vaccines,
                label_unconstrained="",
                col_pareto=color_vaccines,
                label_pareto="",
                xlabel="",
                ylabel=ylab,
                title=f"Case {index_2}",
                linewidth=4,
                s_scatter=0,
                label_scatter="",
                plot="pareto",
                n_vacc = 30000,
                x_total=17,
                y_total=25500,
                scale_total=10*3,
                add_additional=add_additional,
        )
        ax.yaxis.grid(alpha=0.6)
        ax.set_ylim([0,30500])
        if index in [0, 6]:
            ax.text(
                    x_position,
                    y_position*0.82,
                    chr(count_plot),
                    horizontalalignment="center",
                    verticalalignment="center",
                    transform=ax.transAxes,
                    weight="bold",
                    size=y_size,
                    )
            count_plot += 1
        index = index + 1


fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/plot_two_timings_pareto",
    bbox_inches="tight",
)