import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.gridspec import GridSpec
import matplotlib
import geopandas as gpd
from descartes import PolygonPatch
import datetime as dt
from functions.create_vaccine_doses_inflow import create_inflow_from_data_sum


def plotCountryPatch( axes, country_name, fcolor, alpha ):
    # plot a country on the provided axes
    nami = world[world.name == country_name]
    namigm = nami.__geo_interface__['features']  # geopandas's geo_interface
    namig0 = {'type': namigm[0]['geometry']['type'], \
              'coordinates': namigm[0]['geometry']['coordinates']}
    axes.add_patch(PolygonPatch( namig0, fc=fcolor, ec="black", alpha=alpha, zorder=2, label = country_name ))


world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
# Remove French Guiana from France.
shape = world[world['name'] == 'France']['geometry'].all()

# Multipolygon ValueError Workaround.
fr_df = pd.Series(['France', 'France'], name='country')
fr_df = gpd.GeoDataFrame(fr_df, geometry=[shape[1], shape[2]])
fr_df = fr_df.dissolve(by='country')
world.at[world['name'] == 'France', 'geometry'] = fr_df['geometry'].values 
#world.plot()

europe = world[world.continent=="Europe"]
#europe.plot()



bfgu = europe[(europe.name == "Belgium") |  (europe.name == "France") | (europe.name == "Germany") | (europe.name == "United Kingdom") | (europe.name == "Netherlands") 
              | (europe.name == "Luxembourg") | (europe.name == "Switzerland") | (europe.name == "Austria") | (europe.name == "Italy") 
              | (europe.name == "Poland") | (europe.name == "Czechia") |  (europe.name == "Spain") | (europe.name == "Denmark") ] 


countries = ["Belgium", "France", "Germany", "United Kingdom"]
color_countries_map = "blues" #"seagreen"
alpha_countries_map = 0.5


#----------------------------read stringency index----------------------------
stringency_path = "/home/manuel/Documents/VaccinationDistribution/stringency.csv"

stringencies = pd.read_csv(stringency_path)
coun = ["Belgium", "France", "Germany", "Uk"]

#---------Define for time course plot-----------------------------------------
length = 420
end_data = 280
number_yy = 6
interventionA = {
    "t": 46,  # first detected at day 53
    "parameter": {
        "infectious_countryD_vac0_virus2_t0": 20,
        "infectious_countryB_vac0_virus2_t0": 1,
        "infectious_countryC_vac0_virus2_t0": 20,
    },
}
interventionB = {
    "t": 76,  # first detected at day 53
    "parameter": {
        "infectious_countryA_vac0_virus2_t0": 20,
    },
}
interventions = [interventionA, interventionB]
spline_xx = {}
for index in range(8):
    
    spline_xx[f"xx{index}"] = index * (47 + 2/3)
spline_xx["xx8"] = length
#---------------------TRue vaccine doses--------------------------------------
vaccine_inflow = create_inflow_from_data_sum()
europe_vaccine = (vaccine_inflow["europe"].iloc[1:(-1)]).reset_index(drop=True)
uk_vaccine = vaccine_inflow["uk"]

#-------------------Start plot------------------------------------------------
y_position = 1.35
y_size = 28
count_plot = 65
x_position = -0.35

size = 19
font = {"family": "normal", "weight": "normal", "size": size}
matplotlib.rc("font", **font)
matplotlib.rcParams['xtick.labelsize'] = size - 2
matplotlib.rcParams['ytick.labelsize'] = size - 2



color =["Steelblue", "C1", "Seagreen", "Firebrick"]
fig = plt.figure(figsize = (20, 20))
gs = GridSpec(4, 4, figure=fig)
gs.update(wspace=0.4, hspace=0.9)



time = [dt.datetime.strptime(d,'%Y-%m-%d').date() for d in stringencies["Time"]] 
#map
ax = fig.add_subplot(gs[0:2, 0:2]) 
#ax.set_facecolor('lightsteelblue')
bfgu.plot(ax=ax, edgecolor=u'gray', color='whitesmoke', legend=True)
#color = iter(cm.Pastel1([0,1,2,3]))
for index in range(len(countries)):
    c = color[index]
    plotCountryPatch(ax, countries[index], c, alpha = alpha_countries_map )
ax.legend(loc="lower right")
ax.get_xaxis().set_visible(False)
ax.get_yaxis().set_visible(False)
ax.set_title("Spatial allocation of countries")
ax.text(
        x_position+0.14,
        1.04,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=y_size,
)
count_plot += 1

#color = iter(cm.Pastel1([0,1,2,3]))


for index in range(len(countries)):
    x_plot = index % 2
    y_plot = int(np.floor(index/2))

    ax=fig.add_subplot(gs[y_plot, 2+x_plot]) 
    if index == 0:
        ax.text(
                x_position,
                y_position*0.86,
                chr(count_plot),
                horizontalalignment="center",
                verticalalignment="center",
                transform=ax.transAxes,
                weight="bold",
                size=y_size,
        )
        count_plot += 1

    c = color[index]
    ax.plot(time,  stringencies[coun[index]], alpha=alpha_countries_map, color=c)
    if x_plot == 0:
        ax.set_ylabel("Stringency Index")
    if y_plot== 0:
        ax.tick_params(labelbottom=False) 
    ax.set_title(countries[index])
    ax.fill_between(stringencies["Time"], np.repeat(0,stringencies.shape[0]), stringencies[coun[index]], alpha= alpha_countries_map, color=c)
    ax.set_ylim(0,100)
    ax.xaxis.set_tick_params(rotation=45)
    #ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%Y'))
    #ax.xaxis.set_major_locator(mdates.DayLocator())
    ax.yaxis.grid(alpha=0.6)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
#---------------------Vaccines------------------------------------
ax=fig.add_subplot(gs[2, 0])
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

alpha_vacc = 0.5
scale= 10**6
vacs = list(uk_vaccine.columns)
sdate = dt.date(2021,1,1)
edate = dt.date(2022,12,24)
days = pd.date_range(sdate,edate-dt.timedelta(days=1),freq='d')
weeks = days[::7][0:uk_vaccine.shape[0]]
labels = ["mRNA", "Vector"]
colors = ["steelblue", "C1"]

for index in range(len(vacs)):
    ax.plot(weeks, europe_vaccine[vacs[index]]/scale, label=labels[index], color=colors[index], alpha=0.6)
    if index == 0:
        bottom=europe_vaccine[vacs[index+1]]/scale
    else:
        bottom = np.repeat(0, len(weeks))
    ax.fill_between(weeks, bottom, europe_vaccine[vacs[index]]/scale, alpha = alpha_vacc, color=colors[index])
ax.set_ylabel("Vaccinated doses\nin millions")
ax.set_title("Europe")
ax.xaxis.set_tick_params(rotation=45)
ax.yaxis.grid(alpha=0.6)
#ax.legend()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

ax=fig.add_subplot(gs[2, 1])

for index in range(len(vacs)):
    ax.plot(weeks, uk_vaccine[vacs[index]]/scale, label=labels[index], color=colors[index], alpha=0.6)
    if index == 1:
        bottom=uk_vaccine[vacs[index-1]]/scale
    else:
        bottom = np.repeat(0, len(weeks))
    ax.fill_between(weeks, bottom, uk_vaccine[vacs[index]]/scale, alpha = alpha_vacc, color=colors[index])
#ax.set_ylabel("Vaccinated doses in millions")
ax.set_title("United Kingdom")
ax.xaxis.set_tick_params(rotation=45)
ax.yaxis.grid(alpha=0.6)
ax.legend()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

#-----------------------time course-------------------------------------
ax=fig.add_subplot(gs[2, 2:4])

ax.set_xlim([0, 60])
ax.text(
        x_position+0.2,
        y_position,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=y_size,
)
count_plot += 1

ax.set_ylim([0, 1])
ax.get_yaxis().set_visible(False)
    # ax.spines['left'].set_visible(False)
    # ax.spines['bottom'].set_position('center')

color_tl = "seagreen"
ax.fill_between(
        [0, length / 7], [1, 1], [0.8, 0.8], step="pre", alpha=0.95, color=color_tl
)
ax.text(length / 7 / 2, 0.9, "Alpha variant", ha="center", va="center", size = size-5)

ax.fill_between(
        [interventionA["t"] / 7, length / 7],
        [0.8, 0.8],
        [0.6, 0.6],
        step="pre",
        alpha=0.75,
        color=color_tl,
)
#ax.plot([interventionA["t"] / 7, interventionA["t"] / 7], [0.8, 0.6], color="lightgrey", linewidth=0.7, alpha=0.)
ax.text(
        ((length - interventionA["t"]) / 2 + interventionA["t"]) / 7,
        0.7,
        "Delta variant",
        ha="center",
        va="center", size = size - 5,
)

for i in range(8):
    key1 = f"xx{i}"
    key2 = f"xx{i+1}"
    ax.fill_between(
            [spline_xx[key1] / 7, spline_xx[key2] / 7],
            [0.6, 0.6],
            [0.4, 0.4],
            step="pre",
            alpha=0.55,
            color=color_tl,
    )
    ax.plot([spline_xx[key1] / 7, spline_xx[key1] / 7], [0.6, 0.4], color="lightgrey", linewidth=0.7, alpha=0.5)
    ax.text(
            ((spline_xx[key2] - spline_xx[key1]) / 2 + spline_xx[key1]) / 7,
            0.5,
            f"Spline {i+1}",
            ha="center",
            va="center", size = size-5
    )

ax.fill_between(
        [0, length / 7], [0.4, 0.4], [0.2, 0.2], step="pre", alpha=0.35, color=color_tl
)
ax.text(length / 7 / 2, 0.3, "Optimize vaccinations", ha="center", va="center", size = size -5)

#ax.fill_between(
#        [end_data / 7, length / 7],
#        [0.4, 0.4],
#        [0.2, 0.2],
#        step="pre",
#        alpha=0.35,
#        color=color_tl,
#)
#ax.plot([end_data / 7, end_data / 7], [0.4, 0.2], color="lightgrey", linewidth=0.7, alpha=0.8)
#ax.text(
#        ((length - end_data) / 2 + end_data) / 7,
#        0.3,
#        "Pop. based \nallocation",
#        ha="center",
#        va="center", size = size - 4,
#)

ax.fill_between(
        [0, length / 7], [0.2, 0.2], [0.0, 0.0], step="pre", alpha=0.15, color=color_tl
)
ax.text(length / 7 / 2, 0.1, "NPIs active", ha="center", va="center", size = size - 4)
ax.set_xticklabels("Weeks")
ax.set_title("Time course") 

for index in np.linspace(0,1,6):
    ax.plot([0, length/7], [index, index], color="lightgrey", linewidth=0.7, alpha=0.8)
ax.set_xticklabels(["2021-01", "2021-03", "2021-05", "2021-07", "2021-09", "2021-11", "2022-01"])
ax.xaxis.set_tick_params(rotation=45)


fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/covid_plot_one",
    bbox_inches="tight",
)