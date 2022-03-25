import pickle
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from matplotlib.cbook import get_sample_data
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

def autolabel(rects, rounding = 0):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        if rounding == 0:
            number = '{}'.format(int(np.round(height)))
        else:
            number = '{}'.format(np.round(height,rounding))
            if np.round(height,rounding) == int(np.round(height,rounding)):
                number = '{}'.format(int(np.round(height,rounding)))
        ax.annotate(number,
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

appearence = [0.05] + list(range(1, 13))
vac = "Unequal" #only 
initial_states = "Unequal"  
specification = f"inital{initial_states}_vac{vac}"

results = []
for index in appearence:
    specification = f"appears_week_{index}"
    path = (
        f"/home/manuel/Documents/VaccinationDistribution/code/objects/knowledge_optimal_{specification}.pkl"
    )
    
    with open(
          path,
          "rb",
    ) as input:
          loaded_object = pickle.load(input)
    results.append(loaded_object)

appearence = np.floor(appearence)
    
deaths_countryA = pd.DataFrame(data=None, index=range(len(appearence)), columns=["optimal", "pareto", "population"])
deaths_countryB = deaths_countryA.copy()

deaths_countryA.iloc[0,:] = np.array([0.24, 0.94, 0.98]) * 10**5
deaths_countryB.iloc[0,:] = np.array([1.68, 1.14, 1.19]) * 10**5
for index in range(0, len(appearence)):
    r = results[index]
    index_optimal = np.argmin(r["all_strategies"]["fval"])
    index_pareto = np.argmin(r["pareto_improvements"]["fval"])
    
    
    A = [r["all_strategies"].iloc[index_optimal]["countryA"],
         r["pareto_improvements"].iloc[index_pareto]["countryA"],
         r["population_based"][0]]

    B = [r["all_strategies"].iloc[index_optimal]["countryB"],
         r["pareto_improvements"].iloc[index_pareto]["countryB"],
         r["population_based"][1]]
    deaths_countryA.iloc[index] = A
    deaths_countryB.iloc[index] = B

total_deaths = deaths_countryA + deaths_countryB
optimal = total_deaths["optimal"]
pareto = total_deaths["pareto"]
pop = total_deaths["population"]

#-----------------------------------------------------------------------------------
y_position = 1.35
x_position = -0.3
y_size = 28
count_plot = 65
lwidth = 2.5
initials = list(range(13))

size = 19
font = {"family": "normal", "weight": "normal", "size": size}
matplotlib.rc("font", **font)
matplotlib.rcParams['xtick.labelsize'] = size - 2 
matplotlib.rcParams['ytick.labelsize'] = size - 2 

x = np.array(appearence)  # the label locations
width = 0.37  # the width of the bars


fig = plt.figure(constrained_layout=True, figsize = (20,14))
gs = GridSpec(3, 3, figure=fig, height_ratios=[0.8,1,1])
gs.update(wspace=0.04, hspace=0.04)

xlabel = "Week of variant appearence"


#--------------------------------------------------------------------------------------
ax = fig.add_subplot(gs[0, 0:3])

top_stripe = 0.52
bottom_stripe = 0.48
color_stripe = "gray"
alpha_color = 0.6
distance_stripe = 0.2
distance_bar = 0.15
color_bar = "seagreen"
color_bar2 = "steelblue"
length_x = 1
length_y = 1
alpha_bar = 0.55
epsilon= 0
text_im = 0.055
plus_y_text = 0.07


#ax.fill_between([0,length_x], [bottom_stripe, bottom_stripe], [top_stripe, top_stripe],
#                color = color_stripe, alpha = alpha_color)
#ax.plot([0,length_x], [bottom_stripe, bottom_stripe], color = color_stripe)
#ax.plot([0,length_x], [top_stripe, top_stripe], color = color_stripe)
#ax.plot([0,0], [bottom_stripe, top_stripe], color=color_stripe)
#ax.plot([length_x, length_x], [bottom_stripe, top_stripe], color=color_stripe)

imscatter([0 + text_im], top_stripe + distance_stripe + plus_y_text,
          "/home/manuel/Documents/VaccinationDistribution/paper/images/virus_green.png",
          ax=ax, zoom = 0.14)

ax.annotate("Wild-type\nis present",
            xy=(0, 0), xycoords='data',
            xytext=(0, top_stripe + distance_stripe), textcoords='data', ha='center', 
            arrowprops=dict(arrowstyle="-",
                            connectionstyle="arc3"),
            )


pos_blue = length_x / 3
imscatter([pos_blue + text_im + 0.025], top_stripe + distance_stripe + plus_y_text ,
          "/home/manuel/Documents/VaccinationDistribution/paper/images/virus_blue.png",
          ax=ax, zoom = 0.14)

ax.annotate("Variant appears\nin Country 2",
            xy=(pos_blue,0), xycoords='data',
            xytext=(pos_blue, top_stripe + distance_stripe), textcoords='data', ha='center',
            arrowprops=dict(arrowstyle="-",
                            connectionstyle="arc3"),
            )

pos_top_bar = bottom_stripe - distance_stripe - distance_bar
ax.fill_between([0,pos_blue], [0, 0], [pos_top_bar, pos_top_bar],
                color = color_bar, alpha = alpha_bar)
#ax.plot([0,pos_blue], [0, 0], color = color_bar)
#ax.plot([0,pos_blue], [pos_top_bar, pos_top_bar], color = color_bar)
#ax.plot([0,0], [0, pos_top_bar], color=color_bar)
#ax.plot([pos_blue-epsilon, pos_blue-epsilon], [0, pos_top_bar], color=color_bar)

ax.text(pos_blue/2, 0*pos_top_bar+0.2  , "Optimization according to wild-type", ha="center", va="center")

ax.fill_between([pos_blue, length_x], [0, 0], [pos_top_bar, pos_top_bar],
                color = color_bar, alpha = alpha_bar)

ax.fill_between([pos_blue, length_x], [pos_top_bar, pos_top_bar], [2*pos_top_bar, 2*pos_top_bar],
                color = color_bar2, alpha = alpha_bar)

ax.text(pos_blue + (length_x - pos_blue)/2, pos_top_bar + 0.2, "Optimization according to wild-type & variant", ha="center", va="center")



#ax.set_xlabel("Week")
ax.set_ylim([0,length_y])
ax.set_xlim([0,length_x])
ax.spines['left'].set_visible(False) 
ax.yaxis.set_visible(False)

ax.set_xticks([0, pos_blue, length_x])
ax.set_xticklabels(["Start of\nsimulation", "Week of variant\nappearence", "End of\n simulation"])

#rects1 = ax.bar(x, total_deaths["population"]/10**3, width, color = "white", alpha = 0)
#delta = 0.15
#x_m = np.array(range(0, 13)) - delta
#x_p = x_m + 2*delta
#imscatter(x_p, range(0, 13),
#          "/home/manuel/Documents/VaccinationDistribution/paper/images/virus_blue.png",
#          ax=ax, zoom = 0.08)
#imscatter(x_m, np.repeat(0, 13),
#          "/home/manuel/Documents/VaccinationDistribution/paper/images/virus_green.png",
#          ax=ax, zoom = 0.08)
#ax.set_ylim(-2, 42)
#ax.set_ylabel("")
#ax.get_yaxis().set_visible(False)
##ax.set_xlim(0,11)
#ax.set_xlabel(xlabel)
#ax.xaxis.grid(alpha=1)
#ax.spines['left'].set_visible(False)   
#ax.set_xticks(x)

#width = 0.37  # the width of the bars
#ax.yaxis.set_ticks_position('none') 

#ax.annotate('Wild-type', xy=(1-delta, 2), xytext=(0, 14),
#                            arrowprops=dict(facecolor='black', shrink=0.1))

#ax.annotate('Variant', xy=(5 + delta - 0.1, 6), xytext=(4.3, 15),
#                            arrowprops=dict(facecolor='black', shrink=0.1))

#ax.set_title("Scenarios of variant appearence")
#ax.yaxis.set_visible(False)
#ax.spines['left'].set_visible(False)
ax.text(
       x_position + 0.220,
        y_position,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=y_size,
)
count_plot += 1
#-----------------------------------------------------------------------------
ax = fig.add_subplot(gs[1, 0])
diffs = ((-optimal + pop) / pop) 
for index in range(len(initials)):
    ax.plot(np.repeat(initials[index],2), [0, diffs[index]], color = "C0",
            linewidth = lwidth, linestyle = "dotted" )
    #diff = np.round((pop[index] - optimal[index]) / pop[index] * 100,1)
    #diff_plot = (pop[index] - optimal[index])
    #ax.annotate(f"{-diff}%", (initials[index], optimal[index] - 5000), ha='center' )

ax.scatter(initials, diffs,
                color="C0",
                s=200)
   
ax.yaxis.grid(alpha=0.6)    
ax.set_title("Reduction in deceased individuals")

ax.set_ylabel("Optimal strategy")
#ax.set_ylim([180000, 235000])
ylim = ax.get_ylim()
ax.set_xlim(-0.5, 12.5)
col_label_sym = [f"{i}" for i in range(5)] + ["5"]
ax.set_xticklabels([])
vals_1 = ax.get_yticks()
ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals_1])

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

#-----------------------------------------------------------------------------
ax = fig.add_subplot(gs[2, 0])
diffs = ((-pareto + pop) / pop) 
for index in range(len(initials)):
    ax.plot(np.repeat(initials[index],2), [0, diffs[index]], color = "C0",
            linewidth = lwidth, linestyle = "dotted" )
    #diff = np.round((pop[index] - optimal[index]) / pop[index] * 100,1)
    #diff_plot = (pop[index] - optimal[index])
    #ax.annotate(f"{-diff}%", (initials[index], optimal[index] - 5000), ha='center' )

ax.scatter(initials, diffs,
                color="C0",
                s=200)
   
ax.yaxis.grid(alpha=0.6)    
ax.set_title("")

ax.set_ylabel("Pareto strategy")
#ax.set_ylim([180000, 235000])
ylim = ax.get_ylim()
ax.set_xlim(-0.5, 13)
col_label_sym = [f"{i}" for i in range(5)] + ["5"]
#ax.set_xticklabels([])
ax.set_xlabel("Week of variant appearence")

ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals_1])
ax.set_ylim([-0.05, 0.65])
x_range = 2*np.array(range(0, 7))
ax.set_xticks(x_range)
ax.legend()


#------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------
ylim = [-0.49, 0.49]
ax = fig.add_subplot(gs[1, 1])
fractions = compute_splines_from_results_initials(initials=initials, type_opti = "all_strategies",
                                      vac = "vac1", results=results, periods = 14,
                                      length = 21,
                                      total_length = 294, add_additional = None,
                                      grid_points = 1000)


diffs = np.array(fractions) / 0.5 - 1
for index in range(len(initials)):
   ax.plot(np.repeat(initials[index],2), [0, diffs[index]], color = "steelblue",
           linewidth = lwidth, linestyle = "dotted" )

ax.scatter(initials, diffs,
                color="steelblue",
                s=200)
ax.yaxis.grid(alpha=0.6)
ax.set_title("Change of Vaccine 1 doses\nallocated to Country 1")
#ax.set_ylim(ylim)
ax.set_ylabel("Opimal strategy")
ax.legend()
ax.set_xlim(-0.5, 13)
#cols_labels_small = [f"{i}" for i in range(len(initials))]
#ax.set_xticklabels([""] + cols_labels_small )
ax.set_xticklabels([])
#ax.set_xlabel("Case")
ax.set_ylim(ylim)
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
ax = fig.add_subplot(gs[2, 1])
fractions = compute_splines_from_results_initials(initials=initials, type_opti = "pareto_improvements",
                                      vac = "vac1", results=results, periods = 14,
                                      length = 21,
                                      total_length = 294, add_additional = None,
                                      grid_points = 1000)


diffs = np.array(fractions) / 0.5 - 1
for index in range(len(initials)):
   ax.plot(np.repeat(initials[index],2), [0, diffs[index]], color = "steelblue",
           linewidth = lwidth, linestyle = "dotted" )

ax.scatter(initials, diffs,
                color="steelblue",
                s=200)
ax.yaxis.grid(alpha=0.6)
ax.set_title("")
ax.set_ylim(ylim)
#ax.set_ylim(ylim)
ax.set_ylabel("Pareto strategy")
ax.legend()
ax.set_xlim(-0.5, 13)
#cols_labels_small = [f"{i}" for i in range(len(initials))]
#ax.set_xticklabels([""] + cols_labels_small )
#ax.set_xticklabels([""] + col_label_sym)
#ax.set_xlabel("Case")
vals = ax.get_yticks()
ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
ax.set_xticks(x_range)
ax.set_xlabel("Week of variant appearence")

#--------------------------------------------------------------------------------------------------------
ax = fig.add_subplot(gs[1, 2])
fractions = compute_splines_from_results_initials(initials=initials, type_opti = "all_strategies",
                                      vac = "vac2", results=results, periods = 14,
                                      length = 21,
                                      total_length = 294, add_additional = None,
                                      grid_points = 1000)


diffs = np.array(fractions) / 0.5 - 1
for index in range(len(initials)):
   ax.plot(np.repeat(initials[index],2), [0, diffs[index]], color = "steelblue",
           linewidth = lwidth, linestyle = "dotted" )

ax.scatter(initials, diffs,
                color="steelblue",
                s=200)
ax.yaxis.grid(alpha=0.6)
ax.set_title("Change of Vaccine 2 doses\nallocated to Country 1")
#ax.set_ylim(ylim)
ax.set_ylabel("")
ax.legend()
ax.set_ylim(ylim)
ax.set_xlim(-0.5, 13)
#cols_labels_small = [f"{i}" for i in range(len(initials))]
#ax.set_xticklabels([""] + cols_labels_small )
ax.set_xticklabels([])
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
ax = fig.add_subplot(gs[2, 2])
fractions = compute_splines_from_results_initials(initials=initials, type_opti = "pareto_improvements",
                                      vac = "vac2", results=results, periods = 14,
                                      length = 21,
                                      total_length = 294,
                                      grid_points = 1000, add_additional = None)


diffs = np.array(fractions) / 0.5 - 1
for index in range(len(initials)):
   ax.plot(np.repeat(initials[index],2), [0, diffs[index]], color = "steelblue",
           linewidth = lwidth, linestyle = "dotted" )

ax.scatter(initials, diffs,
                color="steelblue",
                s=200)
ax.yaxis.grid(alpha=0.6)
ax.set_title("")
#ax.set_ylim(ylim)
ax.set_ylabel("")

ax.legend()
ax.set_ylim(ylim)
ax.set_xlim(-0.5, 13)
#cols_labels_small = [f"{i}" for i in range(len(initials))]
#ax.set_xticklabels([""] + cols_labels_small )
#ax.set_xticklabels([""] + col_label_sym)
#ax.set_xlabel("Case")
vals = ax.get_yticks()
ax.set_yticklabels(['{:,.0%}'.format(x) for x in vals])
ax.set_xticks(x_range)
ax.set_xlabel("Week of variant appearence")











fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/plot_three_reaction",
    bbox_inches="tight",
)



#------------------------------------------------------------------------------
#ax = fig.add_subplot(gs[1, 2])

#rects1 = ax.bar(x, total_deaths["population"]/10**3, width, label = "Population based", color = "steelblue", alpha = 0.7)

# Add some text for labels, title and custom x-axis tick labels, etc.
#ax.set_ylabel('Deaths in 1,000')
#ax.set_title('Deaths by population-size based vaccine allocation')
#ax.set_xticks(x)
#ax.set_xticklabels(appearence)
#ax.legend()
#ax.set_xlabel(xlabel)
#autolabel(rects1, rounding=1)
#ax.yaxis.grid(alpha=0.6)
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
#count_plot += 1
#-------------------------------------------------------------------------------
#ax = fig.add_subplot(gs[2, 2])

#rects1 = ax.bar(x - width/2,(total_deaths["population"]- total_deaths["optimal"])/10**3, width, label='Optimal', color = "steelblue", alpha = 0.7)
#rects2 = ax.bar(x + width/2, (total_deaths["population"]-total_deaths["pareto"])/10**3, width, label='Pareto optimal', color = "darkorange", alpha=0.7)

## Add some text for labels, title and custom x-axis tick labels, etc.
#ax.set_ylabel('Saved lifes in 1,000')
#ax.set_title('Number of saved lifes compared to population-size based vaccine allocation')
#ax.set_xticks(x)
#ax.set_xlabel(xlabel)#

#ax.set_xticklabels(appearence)
#ax.legend()
#ax.yaxis.grid(alpha=0.6)
#autolabel(rects1, rounding=1)
#autolabel(rects2, rounding=1)
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
#count_plot += 1


#------------------------------------------------------------------------------
#ax = fig.add_subplot(gs[2, 2])
#optimal_percentage = -(total_deaths["optimal"] / total_deaths["population"] - 1)*100

#pareto_percentage = -(total_deaths["pareto"] / total_deaths["population"] - 1)*100

#rects1 = ax.bar(x - width/2, optimal_percentage, width, label='Optimal', color = "steelblue", alpha = 0.7)
#rects2 = ax.bar(x + width/2,pareto_percentage, width, label='Pareto optimal', color = "darkorange", alpha=0.7)

# Add some text for labels, title and custom x-axis tick labels, etc.
#ax.set_ylabel('Relative number of saved lifes in %')
#ax.set_title('Saved lifes relative to population-size based vaccine allocation')
#ax.set_xticks(x)
#ax.set_xticklabels(appearence)
#ax.legend()
#ax.set_ylim(0, 100)
#ax.set_xlabel(xlabel)
#ax.yaxis.grid(alpha=0.6)
#autolabel(rects1)
#autolabel(rects2)
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
#count_plot += 1



#fig.tight_layout(pad=3.5)



#fig.savefig(
#    "/home/manuel/Documents/VaccinationDistribution/paper/images/plot_three_reaction",
#    bbox_inches="tight",
#)




