import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.interpolate import CubicHermiteSpline as cbs
from matplotlib.gridspec import GridSpec
from numpy import trapz
import matplotlib as mpl

mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.spines.right"] = False
font = {"size": 20}

mpl.rc("font", **font)


def plot_bars_deaths(
    dict_output,
    ax=None,
    case="initalEqual_vacEqual",
    total=False,
    unit=10 ** 6,
    title="Equal initial states; Equal vaccines",
    xlabel="Areas",
    ylabel="Number of deaths",
    label_optimal="Optimal",
    label_Pareto="Pareto",
    label_population="Population",
    ylim=None,
):

    # preprocess
    output = dict_output[case]
    appended_df = output["optimal_strategies"].append(output["pareto_frontier"])
    country_names = [
        x for x in appended_df.columns if "country" in x and not ("_" in x)
    ]
    add_row = dict(
        zip(
            country_names + ["fval"],
            output["population_based"] + [np.sum(output["population_based"])],
        )
    )

    appended_df = appended_df.append(add_row, ignore_index=True)

    unrestricted_min = appended_df.iloc[np.argmin(appended_df["fval"])][
        ["fval"] + country_names
    ]

    pareto_df = appended_df[appended_df["countryA"] <= output["population_based"][0]]
    for i in range(len(country_names)):
        pareto_df = pareto_df[
            pareto_df[country_names[i]] <= output["population_based"][i]
        ]

    pareto_optimal = pareto_df.iloc[np.argmin(pareto_df["fval"])][
        ["fval"] + country_names
    ]
    population_based = pd.Series(
        [np.sum(output["population_based"])] + output["population_based"],
        index=["fval"] + country_names,
    )

    if total is True:
        X = ["Total"] + country_names
        unrestricted = np.round(list(unrestricted_min / unit), 2)
        pareto = np.round(list(pareto_optimal / unit), 2)
        X_axis = np.arange(len(X))

        if ax is None:
            fig, ax = plt.subplots()
        rects1 = ax.bar(X_axis - 0.2, unrestricted, 0.2, label=label_optimal)
        rects2 = ax.bar(X_axis, pareto, 0.2, label=label_Pareto)
        rects3 = ax.bar(X_axis + 0.2, pareto, 0.2, label=label_population)
        ax.set_xticks(X_axis)
        ax.set_xticklabels(X)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(
            loc="upper center",
            bbox_to_anchor=(0.5, -0.14),
            fancybox=True,
            shadow=True,
            ncol=2,
        )
    elif total is False:
        X = ["Total"] + country_names
        unrestricted = np.round(
            list(((unrestricted_min - population_based) / population_based) * 100), 2
        )
        pareto = np.round(
            list(((pareto_optimal - population_based) / population_based) * 100), 2
        )

        X_axis = np.arange(len(X))

        if ax is None:
            fig, ax = plt.subplots()
        rects1 = ax.bar(
            X_axis - 0.2, unrestricted, 0.4, label=label_optimal, edgecolor="black"
        )
        rects2 = ax.bar(
            X_axis + 0.2, pareto, 0.4, label=label_Pareto, edgecolor="black"
        )
        ax.set_xticks(X_axis)
        ax.set_xticklabels(X)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)

        ax.bar_label(rects1, padding=5, fmt="%.1f%%")
        ax.bar_label(rects2, padding=5, fmt="%.1f%%")
        mini = np.min([unrestricted, pareto])
        maxi = np.max([unrestricted, pareto])

        ax.set_ylim([1.2 * mini, 1.2 * maxi])

        # ax.legend(
        #    loc="upper center",
        #    bbox_to_anchor=(0.5, -0.14),
        #    fancybox=True,
        #    shadow=True,
        #    ncol=2,
        # )
        if not (ylim is None):
            ax.set_ylim(ylim)
        return ax


def plot_bars_deaths2(
    dict_output,
    ax=None,
    case="initalEqual_vacEqual",
    total=False,
    unit=10 ** 6,
    title="Equal initial states; Equal vaccines",
    xlabel="Areas",
    ylabel="Number of deaths",
    label_optimal="Optimal",
    label_Pareto="Pareto",
    label_population="Population",
    ylim=None,
):

    # preprocess
    output = dict_output[case]
    appended_df = output["optimal_strategies"].append(output["pareto_frontier"])
    country_names = [
        x for x in appended_df.columns if "country" in x and not ("_" in x)
    ]

    appended_df = output["all_strategies"]

    unrestricted_min = appended_df.iloc[np.argmin(appended_df["fval"])][
        ["fval"] + country_names
    ]

    pareto_df = output["pareto_improvements"]

    pareto_optimal = pareto_df.iloc[np.argmin(pareto_df["fval"])][
        ["fval"] + country_names
    ]
    population_based = pd.Series(
        [np.sum(output["population_based"])] + output["population_based"],
        index=["fval"] + country_names,
    )

    if total is True:
        X = ["Total"] + country_names
        unrestricted = list(unrestricted_min / unit)
        pareto = list(pareto_optimal / unit)
        X_axis = np.arange(len(X))

        if ax is None:
            fig, ax = plt.subplots()
        ax.bar(X_axis - 0.2, unrestricted, 0.2, label=label_optimal)
        ax.bar(X_axis, pareto, 0.2, label=label_Pareto)
        ax.bar(X_axis + 0.2, pareto, 0.2, label=label_population)
        ax.set_xticks(X_axis)
        ax.set_xticklabels(X)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(
            loc="upper center",
            bbox_to_anchor=(0.5, -0.14),
            fancybox=True,
            shadow=True,
            ncol=2,
        )
    elif total is False:
        X = ["Total"] + country_names
        unrestricted = list(
            ((unrestricted_min - population_based) / population_based) * 100
        )
        pareto = list(((pareto_optimal - population_based) / population_based) * 100)

        X_axis = np.arange(len(X))

        if ax is None:
            fig, ax = plt.subplots()
        ax.bar(X_axis - 0.2, unrestricted, 0.4, label=label_optimal)
        ax.bar(X_axis + 0.2, pareto, 0.4, label=label_Pareto)
        ax.set_xticks(X_axis)
        ax.set_xticklabels(X)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        # for i in range(len(X)):
        # ax.annotate(np.round(unrestricted[0],4), (-0.25 + X_axis[0], unrestricted[0] - 0.1 ))
        # ax.legend(
        #    loc="upper center",
        #    bbox_to_anchor=(0.5, -0.14),
        #    fancybox=True,
        #    shadow=True,
        #    ncol=2,
        # )
        if not (ylim is None):
            ax.set_ylim(ylim)
        return ax


def plot_pareto_front(
    dict_output,
    case,
    fig=None,
    ax=None,
    size_optimal=10,
    size_points=3,
    title="",
    alpha=0.3,
    color_pareto="C2",
    linewidth_pareto=0.8,
    xlabel="Country A",
    ylabel="Country B",
    color_min="C0",
    color_pareto_improvement="C1",
):

    pareto_x = dict_output[case]["pareto_frontier"]["countryA"]
    pareto_y = dict_output[case]["pareto_frontier"]["countryB"]
    pareto = dict_output[case]["population_based"]
    output = dict_output[case]
    appended_df = output["optimal_strategies"].append(output["pareto_frontier"])
    country_names = [
        x for x in appended_df.columns if "country" in x and not ("_" in x)
    ]

    add_row = dict(
        zip(
            country_names + ["fval"],
            output["population_based"] + [np.sum(output["population_based"])],
        )
    )

    appended_df = appended_df.append(add_row, ignore_index=True)

    amin = np.argmin(appended_df["fval"])

    pareto_optimal_df = appended_df[
        (appended_df["countryA"] <= pareto[0]) & (appended_df["countryB"] <= pareto[1])
    ].reset_index()
    pareto_amin = np.argmin(pareto_optimal_df["fval"])

    if ax is None:
        fig, ax = plt.subplots()

    im = ax.scatter(
        pareto_x,
        pareto_y,
        c=pareto_x + pareto_y,
        s=size_points,
    )
    minA = np.min(pareto_x)
    minB = np.min(pareto_y)
    ax.fill_between(
        x=[minA, pareto[0]],
        y1=[pareto[1], pareto[1]],
        y2=[minB, minB],
        alpha=alpha,
        color=color_pareto,
    )
    ax.hlines(
        appended_df["countryB"][amin],
        xmin=minA,
        xmax=appended_df["countryA"][amin],
        linewidth=linewidth_pareto,
        linestyle="dashed",
        color=color_min,
    )
    ax.vlines(
        appended_df["countryA"][amin],
        ymin=minB,
        ymax=appended_df["countryB"][amin],
        linewidth=linewidth_pareto,
        linestyle="dashed",
        color=color_min,
    )

    ax.hlines(
        pareto_optimal_df["countryB"][pareto_amin],
        xmin=minA,
        xmax=pareto_optimal_df["countryA"][pareto_amin],
        linewidth=linewidth_pareto,
        linestyle="dashed",
        color=color_pareto_improvement,
    )
    ax.vlines(
        pareto_optimal_df["countryA"][pareto_amin],
        ymin=minB,
        ymax=pareto_optimal_df["countryB"][pareto_amin],
        linewidth=linewidth_pareto,
        linestyle="dashed",
        color=color_pareto_improvement,
    )
    cbar = fig.colorbar(im, ax=ax)
    cbar.formatter.set_powerlimits((0, 0))
    ax.scatter(
        appended_df["countryA"][amin],
        appended_df["countryB"][amin],
        color="firebrick",
        s=size_optimal,
    )
    ax.scatter(
        pareto_optimal_df["countryA"][pareto_amin],
        pareto_optimal_df["countryB"][pareto_amin],
        color="seagreen",
        s=size_optimal,
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.ticklabel_format(axis="both", style="sci", scilimits=(0, 0))
    ax.set_title(title)

    return ax


def plot_pareto_front2(
    dict_output,
    case,
    fig=None,
    ax=None,
    size_optimal=10,
    size_points=3,
    title="",
    alpha=0.3,
    color_pareto="C0",
    linewidth_pareto=0.8,
    xlabel="Country A",
    ylabel="Country B",
):

    pareto_x = dict_output[case]["pareto_frontier"]["countryA"]
    pareto_y = dict_output[case]["pareto_frontier"]["countryB"]
    pareto = dict_output[case]["population_based"]
    output = dict_output[case]

    appended_df = output["all_strategies"]

    amin = np.argmin(appended_df["fval"])

    pareto_optimal_df = output["pareto_improvements"]
    pareto_amin = np.argmin(pareto_optimal_df["fval"])

    if ax is None:
        fig, ax = plt.subplots()

    im = ax.scatter(
        pareto_x,
        pareto_y,
        c=pareto_x + pareto_y,
        s=size_points,
    )
    minA = np.min(pareto_x)
    minB = np.min(pareto_y)
    ax.fill_between(
        x=[minA, pareto[0]],
        y1=[pareto[1], pareto[1]],
        y2=[minB, minB],
        alpha=alpha,
        color=color_pareto,
    )
    ax.hlines(
        appended_df["countryB"][amin],
        xmin=minA,
        xmax=appended_df["countryA"][amin],
        linewidth=linewidth_pareto,
        linestyle="dashed",
        color="firebrick",
    )
    ax.vlines(
        appended_df["countryA"][amin],
        ymin=minB,
        ymax=appended_df["countryB"][amin],
        linewidth=linewidth_pareto,
        linestyle="dashed",
        color="firebrick",
    )

    ax.hlines(
        pareto_optimal_df["countryB"][pareto_amin],
        xmin=minA,
        xmax=pareto_optimal_df["countryA"][pareto_amin],
        linewidth=linewidth_pareto,
        linestyle="dashed",
        color="seagreen",
    )
    ax.vlines(
        pareto_optimal_df["countryA"][pareto_amin],
        ymin=minB,
        ymax=pareto_optimal_df["countryB"][pareto_amin],
        linewidth=linewidth_pareto,
        linestyle="dashed",
        color="seagreen",
    )
    cbar = fig.colorbar(im, ax=ax)
    cbar.formatter.set_powerlimits((0, 0))
    ax.scatter(
        appended_df["countryA"][amin],
        appended_df["countryB"][amin],
        color="firebrick",
        s=size_optimal,
    )
    ax.scatter(
        pareto_optimal_df["countryA"][pareto_amin],
        pareto_optimal_df["countryB"][pareto_amin],
        color="seagreen",
        s=size_optimal,
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.ticklabel_format(axis="both", style="sci", scilimits=(0, 0))
    ax.set_title(title)

    return ax


# From Lorenzos Code------------------------------------------------------------------------
def finite_differences(xx, yy):
    dd = []

    fd = onesidedFD(yy[0], yy[1], xx[1] - xx[0])
    dd.append(fd)

    for i in range(1, len(xx) - 1):
        dd.append(
            centeredFD(
                yy[i - 1], yy[i], yy[i + 1], xx[i] - xx[i - 1], xx[i + 1] - xx[i]
            )
        )

    fd = onesidedFD(yy[-2], yy[-1], xx[-1] - xx[-2])
    dd.append(fd)

    return np.asarray(dd)


def onesidedFD(y0, y1, h):

    return (y1 - y0) * 1 / h


def centeredFD(ym1, y0, yp1, hm, hp):
    if hm == hp:
        return (yp1 - ym1) / (2 * hm)
    else:
        return ((yp1 - y0) / hp + (y0 - ym1) / hm) / 2


# ----------------------------------------------------------------------------------------
def get_spline(array, periods, length, total_length, grid_points=6000, transform=True):

    y = np.array(np.log(array / (1 - array)))
    x = np.linspace(0, periods * length, periods + 1)

    fd = finite_differences(x, y)
    # raise ValueError(fd)
    spline = cbs(x, y, fd)
    grid = np.linspace(0, total_length, grid_points)

    spline_vals = spline(grid)
    if transform == True:
        logistic = 1 / (1 + np.exp(-spline_vals))

    return logistic


def plot_best_strategy(
    dict_output,
    vac_interest,
    case,
    x_scatter=None,
    fig=None,
    ax=None,
    periods=8,
    length=14,
    total_length=140,
    grid_points=6000,
    col_unconstrained="C0",
    label_unconstrained="Optimal",
    col_pareto="C1",
    label_pareto="Pareto optimal",
    xlabel="Time",
    ylabel="% of vaccine in Country A",
    title="",
    linewidth=1,
    s_scatter=4,
    label_scatter="",
    plot=None,
    n_vacc=1,
    x_total=16,
    y_total=0.2,
    scale_total=1,
    add_additional=None,
):

    pareto = dict_output[case]["population_based"]
    output = dict_output[case]
    appended_df = output["optimal_strategies"].append(output["pareto_frontier"])
    country_names = [
        x for x in appended_df.columns if "country" in x and not ("_" in x)
    ]
    if not(add_additional is None):
        for index in add_additional["integers"]:
            for index2 in ["country A"]:
                for index3 in ["vac1", "vac2"]:
                    name = f"yy_{index2}_{index3}_{index}"
                    appended_df[name] = add_additional["number"]
        
        
    pars = [x for x in appended_df.columns if "yy_" in x]
    pars1 = [x for x in pars if "vac1" in x]
    pars2 = [x for x in pars if "vac2" in x]

    add_row = dict(
        zip(
            country_names + ["fval"],
            output["population_based"] + [np.sum(output["population_based"])],
        )
    )

    appended_df = appended_df.append(add_row, ignore_index=True)

    amin = np.argmin(appended_df["fval"])

    pareto_optimal_df = appended_df[
        (appended_df["countryA"] <= pareto[0]) & (appended_df["countryB"] <= pareto[1])
    ].reset_index()
    pareto_amin = np.argmin(pareto_optimal_df["fval"])

    optimal_vacc_strategy1 = appended_df[pars1].iloc[amin]
    optimal_pareto_strategy1 = pareto_optimal_df[pars1].iloc[pareto_amin]

    if optimal_pareto_strategy1.isnull().values.any():
        optimal_pareto_strategy1 = pd.Series(
            np.repeat(0.5, len(optimal_pareto_strategy1)),
            index=pareto_optimal_df[pars1].iloc[pareto_amin].index,
        )

    if len(pars2) > 0:
        optimal_vacc_strategy2 = appended_df[pars2].iloc[amin]
        optimal_pareto_strategy2 = pareto_optimal_df[pars2].iloc[pareto_amin]

        if optimal_pareto_strategy2.isnull().values.any():
            optimal_pareto_strategy2 = pd.Series(
                np.repeat(0.5, len(optimal_pareto_strategy2)),
                index=pareto_optimal_df[pars2].iloc[pareto_amin].index,
            )

    if vac_interest == "vac1":
        vac_pareto = get_spline(
            optimal_pareto_strategy1,
            periods=periods,
            length=length,
            total_length=total_length,
            grid_points=grid_points,
        )
        vac_unconstrained = get_spline(
            optimal_vacc_strategy1,
            periods=periods,
            length=length,
            total_length=total_length,
            grid_points=grid_points,
        )
        scatter_vac_pareto = optimal_pareto_strategy1
        scatter_vac_unconstr = optimal_vacc_strategy1
    else:
        vac_pareto = get_spline(
            optimal_pareto_strategy2,
            periods=periods,
            length=length,
            total_length=total_length,
            grid_points=grid_points,
        )
        vac_unconstrained = get_spline(
            optimal_vacc_strategy2,
            periods=periods,
            length=length,
            total_length=total_length,
            grid_points=grid_points,
        )
        scatter_vac_pareto = optimal_pareto_strategy2
        scatter_vac_unconstr = optimal_vacc_strategy2
    if ax is None:
        fig, ax = plt.subplots()
    
    if plot == "pareto" or plot is None:

        ax.plot(
            np.linspace(0, total_length, grid_points) / 7,
            vac_pareto*n_vacc,
            color=col_pareto,
            linewidth=linewidth,
            label=label_pareto,
        )
        
        ax.fill_between(np.linspace(0, total_length, grid_points) / 7, vac_pareto*n_vacc, 
                        np.repeat(0, len(vac_pareto))*n_vacc, color=col_pareto, alpha = 0.3)
        ax.plot(np.linspace(0, total_length, grid_points) / 7,
                np.repeat(0.5, len(vac_pareto))*n_vacc,
                color="black", linestyle="dashed", label="Population \nallocation")
        time = np.linspace(0, total_length, grid_points) / 7
        area = trapz(vac_pareto*n_vacc, dx=(time[1] - time[0]))
        ax.text(x_total, y_total, f"Total doses: \n{np.round(area, 2)}",
                horizontalalignment="center",
                verticalalignment="center",
                bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))
        
    if plot == "optimal" or plot is None:
        ax.plot(
            np.linspace(0, total_length, grid_points) / 7,
            vac_unconstrained*n_vacc,
            color=col_unconstrained,
            linewidth=linewidth,
            label=label_unconstrained,
        )
        ax.plot(np.linspace(0, total_length, grid_points) / 7,
                np.repeat(0.5, len(vac_unconstrained))*n_vacc,
                color="black", linestyle="dashed", label="Population \nallocation")
        time = np.linspace(0, total_length, grid_points-1) / 7 
        area = trapz(vac_unconstrained*n_vacc, dx=(time[1] - time[0]))
        ax.fill_between(np.linspace(0, total_length, grid_points) / 7, vac_unconstrained*n_vacc, 
                        np.repeat(0, len(vac_unconstrained))*n_vacc, color=col_unconstrained,
                        alpha = 0.3)
        ax.text(x_total, y_total, f"Total doses: \n{np.round(area, 2)}",
                horizontalalignment="center",
                verticalalignment="center",
                bbox=dict(facecolor='none', edgecolor='black', boxstyle='round,pad=1'))

    #ax.scatter(x_scatter, scatter_vac_pareto, s=s_scatter, label=label_scatter,
    #           color=col_pareto)
    #ax.scatter(x_scatter, scatter_vac_unconstr, s=s_scatter, label=label_scatter,
    #           color=col_pareto)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_ylim([-0.05, 1.05])

    return ax


def plot_distance_curves(
    dict_use,
    max_index=None,
    linewidth=1,
    color_optimal="C0",
    color_pareto="C1",
    color_pop="C2",
    label_optimal="Optimal",
    label_pareto="pareto",
    label_pop="pop",
    var="fval",
    relative=True,
    ax=None,
    x_label="Distance parameter",
    y_label="% deaths compared to the Population based strategy",
    title="",
    vline=None,
    vline_color="C3",
    vline_width=4,
    vline_label="Previous parameter",
    v_ymin=0,
    v_ymax=-10,
    ylim=None,
):

    if ax is None:
        fig, ax = plt.subplots()

    optimal = (
        dict_use["optimal"]
        .sort_values(by=["distance"])[var][
            (dict_use["optimal"]["fval"] > 0)
            & (dict_use["optimal"]["countryA"] > 0)
            & (dict_use["optimal"]["countryB"] > 0)
            & (dict_use["pareto"]["fval"] > 0)
            & (dict_use["pareto"]["countryA"] > 0)
            & (dict_use["pareto"]["countryB"] > 0)
        ]
        .reset_index(drop=True)
    )

    pareto = (
        dict_use["pareto"]
        .sort_values(by=["distance"])[var][
            (dict_use["pareto"]["fval"] > 0)
            & (dict_use["pareto"]["countryA"] > 0)
            & (dict_use["pareto"]["countryB"] > 0)
            & (dict_use["optimal"]["fval"] > 0)
            & (dict_use["optimal"]["countryA"] > 0)
            & (dict_use["optimal"]["countryB"] > 0)
        ]
        .reset_index(drop=True)
    )

    pop = (
        dict_use["pop_based"]
        .sort_values(by=["distance"])[var][
            (dict_use["pop_based"]["fval"] > 0)
            & (dict_use["pop_based"]["countryA"] > 0)
            & (dict_use["pop_based"]["countryB"] > 0)
            & (dict_use["pareto"]["fval"] > 0)
            & (dict_use["pareto"]["countryA"] > 0)
            & (dict_use["pareto"]["countryB"] > 0)
            & (dict_use["optimal"]["fval"] > 0)
            & (dict_use["optimal"]["countryA"] > 0)
            & (dict_use["optimal"]["countryB"] > 0)
        ]
        .reset_index(drop=True)
    )

    distance = 1 - dict_use["optimal"].sort_values(by=["distance"])["distance"][
        (dict_use["pareto"]["fval"] > 0)
        & (dict_use["pareto"]["countryA"] > 0)
        & (dict_use["pareto"]["countryB"] > 0)
        & (dict_use["optimal"]["fval"] > 0)
        & (dict_use["optimal"]["countryA"] > 0)
        & (dict_use["optimal"]["countryB"] > 0)
    ].reset_index(drop=True)

    if max_index is None:
        max_index = len(pop)

    a = (optimal[0:max_index] - pop[0:max_index]) / pop[0:max_index]
    b = (pareto[0:max_index] - pop[0:max_index]) / pop[0:max_index]
    if relative is True:
        ax.plot(
            distance[0:max_index],
            a * 100,
            color=color_optimal,
            linewidth=linewidth,
            label=label_optimal,
        )
        ax.plot(
            distance[0:max_index],
            b * 100,
            color=color_pareto,
            linewidth=linewidth,
            label=label_pareto,
        )

    elif relative is False:
        ax.plot(
            distance[0:max_index],
            optimal[0:max_index],
            color=color_optimal,
            linewidth=linewidth,
            label=label_optimal,
        )
        ax.plot(
            distance[0:max_index],
            pareto[0:max_index],
            color=color_pareto,
            linewidth=linewidth,
            label=label_pareto,
        )
        ax.plot(
            distance[0:max_index],
            pop[0:max_index],
            color=color_pop,
            linewidth=linewidth,
            label=label_pop,
        )
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    if not (vline is None):
        ax.vlines(
            vline,
            ymin=v_ymin,
            ymax=v_ymax,
            color=vline_color,
            linewidth=vline_width,
            linestyles="dashed",
            label=vline_label,
        )
    if not (ylim is None):
        ax.set_ylim(ylim)
    return ax


# ---------------------------------------------------------------------------------------------------------------------
def plot_bars_multiple(
    ax,
    unconstr_deaths,
    pop_deaths,
    constr_deaths,
    label_optimal="Optimal strategy",
    label_Pareto="Pareto strategy",
    X=["Total"] + ["Belgium", "France", "Germany", "United \nKingdom"],
    xlabel="Countries",
    ylabel="Difference in %",
    title="Number of deaths per strategy and country compared to population strategy",
    color_good="seagreen",
    color_bad="firebrick",
    alpha=0.3,
    label_good="Improvement",
    label_bad="Deterioration",
    xlim=None,
):
    unrestricted = (unconstr_deaths / pop_deaths - 1) * 100
    pareto = (constr_deaths / pop_deaths - 1) * 100
    minimum = np.min([pareto, unrestricted]) * 1.05
    maximum = np.max([pareto, unrestricted]) * 1.05
    X_axis = np.arange(len(X))

    ax.bar(X_axis - 0.3, unrestricted, 0.3, label=label_optimal, edgecolor="grey")
    ax.bar(X_axis, pareto, 0.3, label=label_Pareto, edgecolor="grey")
    ax.set_xticks(X_axis - 0.1)
    ax.set_xticklabels(X)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    a = None
    if not (xlim is None):
        ax.set_xlim(xlim)
        a = xlim[1]
    if a is None:
        a = len(X)
    # ax.fill_between([-0.4, a  + 0.2], [0, 0], [minimum,minimum], step="pre", alpha=alpha, color = color_good,
    #                label=label_good)
    # ax.fill_between([-0.4, a  + 0.2], [maximum, maximum], [0,0], step="pre", alpha=alpha, color = color_bad,
    #                label=label_bad)

    if not (xlim is None):
        ax.set_xlim(xlim)
    ax.set_ylim([minimum, maximum])
    ax.legend()


def plot_trajectories_aggregated(
    ax,
    length,
    pop_trajectory,
    unconstr_trajectory,
    constr_trajectory,
    labels=["Population", "Optimal", "Pareto"],
    colors=["C0", "C1", "C2"],
    alphas=[0.6, 0.4, 0.2],
    xlabel="Weeks",
    ylabel="Infected individuals \nin millions",
    title="Total number of infected individuals",
    target="infectious",
    scale=10 ** 6,
    fill_between=False,
    plot_legend=True,
):
    index_axis = np.array(list(pop_trajectory.reset_index(drop=True).index))
    x_axis = index_axis / index_axis[-1] * length / 7
    trajectories = [unconstr_trajectory, constr_trajectory, pop_trajectory]
    sum_infectious = {}
    for index in range(len(trajectories)):
        df = trajectories[index]
        if len(target) == 1:
            states_infectious = [x for x in df.columns if target[0] in x]
        if len(target) == 2:
            states_infectious = [
                x for x in df.columns if (target[0] in x) and (target[1] in x)
            ]
        sum_infectious = df[states_infectious].sum(axis=1)
        if labels[index] == "Population":
            linestyle = "dashed"
            alpha = 0.8
        else:
            linestyle = "solid"
            alpha = 1
        ax.plot(
            x_axis,
            sum_infectious / scale,
            label=labels[index],
            color=colors[index],
            linestyle=linestyle,
            alpha=alpha,
        )
        if fill_between is True:
            ax.fill_between(
                x_axis,
                sum_infectious / scale,
                np.repeat(0, len(x_axis)),
                step="pre",
                alpha=alphas[index],
                color=colors[index],
            )
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        if plot_legend is True:
            ax.legend()


def plot_trajectories_aggregated_vac(
    ax,
    length,
    pop_trajectory,
    unconstr_trajectory,
    constr_trajectory,
    labels=["Population", "Optimal", "Pareto"],
    colors=["C0", "C1", "C2"],
    alphas=[0.6, 0.4, 0.2],
    xlabel="Weeks",
    ylabel="",
    title="",
    scale=10 ** 6,
    fill_between=False,
    plot_legend=True,
):
    index_axis = np.array(list(pop_trajectory.reset_index(drop=True).index))
    x_axis = index_axis / index_axis[-1] * length / 7
    trajectories = [unconstr_trajectory, constr_trajectory, pop_trajectory]
    sum_infectious = {}
    for index in range(len(trajectories)):
        df = trajectories[index]
        states_vaccinated = [
            x
            for x in df.columns
            if (("vac1" in x) or ("vac2" in x) or ("recoverede" in x))
            and not ("dead" in x)
        ]
        states_alive = [x for x in df.columns if not ("dead" in x)]
        sum_vaccinated = df[states_vaccinated].sum(axis=1)
        sum_alive = df[states_alive].sum(axis=1)
        prop_vac = sum_vaccinated / sum_alive

        ax.plot(
            x_axis,
            prop_vac,
            label=labels[index],
            color=colors[index],
        )
        if fill_between is True:
            ax.fill_between(
                x_axis,
                sum_infectious / scale,
                np.repeat(0, len(x_axis)),
                step="pre",
                alpha=alphas[index],
                color=colors[index],
            )
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        if plot_legend is True:
            ax.legend()


def plot_vac_allocated(
    ax,
    colors,
    time,
    dict_out,
    index_vac,
    index_areas,
    areas,
    scale,
    countries,
    col_vac1="C7",
    col_vac2="C8",
    label_vac1="mRNA",
    label_vac2="Vector",
    vac=["vac1", "vac2"],
    types=["unconstrained", "constrained", "pop"],
    ylabel="% received",
    xlabel="Weeks",
    labels=["Optimal", "Pareto", "Population"],
    alphas=[0.1, 0.1, 0.1],
    axvline_x=40,
    ylim=[-0.05, 0.9],
    title="Vaccine received in ",
    total=True,
    spline_xx=None,
    numb_xx=4,
    s=5,
):
    for index_type in range(len(types)):
        vac_available = dict_out["vaccine"][vac[index_vac]]

        name = f"{types[index_type]}_{areas[index_areas]}_{vac[index_vac]}"
        vac_prop = dict_out["allocated_best"][name]
        vac_allocated = vac_available * vac_prop

        if total is True:
            y = vac_allocated / scale
        else:
            y = vac_prop
        if types[index_type] == "pop":
            alpha = 0.8
            linestyle = "dashed"
        else:
            alpha = 1
            linestyle = "solid"
        ax.plot(
            time / 7,
            y,
            color=colors[index_type],
            label=labels[index_type],
            linestyle=linestyle,
            alpha=alpha,
        )
        if not (spline_xx is None):
            if types[index_type] != "pop":
                xx = np.array(list(spline_xx.values()))
                y_index = xx / xx[-1] * (len(y) - 1)
                ax.scatter(
                    xx[0:numb_xx] / 7,
                    y.loc[np.round(y_index[0:numb_xx])],
                    s=s,
                )
        # y_lim = [ax.get_yticks()[0], ax.get_yticks()[-1]]
        # ax.fill_between([mini, axvline_x], y_lim,
        #                    color="grey", step="pre", alpha=0.5)
        # ax.fill_between([mini, 60], y_lim,
        #                    color="grey", step="pre", alpha=0.2)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    ax.set_title(title + f"{countries[index_areas]}")

    # ax.axvline(axvline_x ,0,  1, color = "firebrick",
    #           linestyle = "dashed", label = "Last optimitaion \npoint", linewidth = 0.7)
    if [index_vac, index_areas] == [0, 0]:
        ax.legend()


# ----------------------------------------------------------------------------------------------------
def plot_four_country_overview(
    spline_xx,
    vaccine_inflow,
    number_yy,
    length,
    interventionA,
    interventionB,
    end_data,
    start_population,
    infectious_t0,
    recovered_t0,
    delta,
    omega,
    countries,
    areas,
    par_R,
    number_xx_R,
    total_grid,
    df_inf_true,
    df_infected,
    grid_data,
    grid_sim,
    scale,
    text_x,
    text_y,
    ylim,
    text_str,
    text_lockdown_x,
    text_lockdown_y=0.02,
    text_lockdown_str="Constant \nNPIs",
    color_prop=["C4", "C5", "C6", "C7"],
    label_vac1="mRNA",
    label_vac2="vector",
    color_vac1="C7",
    color_vac2="C8",
    title_vac="Available vaccines",
    title_setup="Set-up",
    position_start_vac=[0, 0.5],
    height_start_vac=0.3,
    letter_size=16,
    letter_y=1.06,
    size=(18, 16),
):
    linspace = []
    for j in range(len(spline_xx.values()) - 1):
        new = list(np.linspace(0, list(spline_xx.values())[j + 1], 1000))
        linspace += new

    vaccine_available = pd.DataFrame(
        {
            "vac1": np.repeat(
                np.array(list(vaccine_inflow.values())[0 : (number_yy - 1)]), 1000
            ),
            "vac2": np.repeat(
                np.array(
                    list(vaccine_inflow.values())[(number_yy - 1) : (2 * number_yy)]
                ),
                1000,
            ),
            "t": np.linspace(0, length, len(linspace)),
        }
    )

    fig = plt.figure(constrained_layout=True, figsize=size)
    gs = GridSpec(3, 4, figure=fig)
    count_plot = 97

    ax = fig.add_subplot(gs[0, :1])
    ax.set_xlim([0, 60])
    ax.set_ylim([0, 1])
    ax.get_yaxis().set_visible(False)
    # ax.spines['left'].set_visible(False)
    # ax.spines['bottom'].set_position('center')
    ax.text(
        -0.05,
        letter_y,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=letter_size,
    )

    color_tl = "seagreen"
    ax.fill_between(
        [0, length / 7], [1, 1], [0.8, 0.8], step="pre", alpha=0.6, color=color_tl
    )
    ax.text(length / 7 / 2, 0.9, "Alpha variant", ha="center", va="center")

    ax.fill_between(
        [interventionA["t"] / 7, length / 7],
        [0.8, 0.8],
        [0.6, 0.6],
        step="pre",
        alpha=0.5,
        color=color_tl,
    )
    ax.text(
        ((length - interventionA["t"]) / 2 + interventionA["t"]) / 7,
        0.7,
        "Delta variant",
        ha="center",
        va="center",
    )

    for i in range(3):
        key1 = f"xx{i}"
        key2 = f"xx{i+1}"
        ax.fill_between(
            [spline_xx[key1] / 7, spline_xx[key2] / 7],
            [0.6, 0.6],
            [0.4, 0.4],
            step="pre",
            alpha=0.4,
            color=color_tl,
        )
        ax.text(
            ((spline_xx[key2] - spline_xx[key1]) / 2 + spline_xx[key1]) / 7,
            0.5,
            f"Spline {i+1}",
            ha="center",
            va="center",
        )

    ax.fill_between(
        [0, end_data / 7], [0.4, 0.4], [0.2, 0.2], step="pre", alpha=0.3, color=color_tl
    )
    ax.text(end_data / 7 / 2, 0.3, "Optimize vaccinations", ha="center", va="center")

    ax.fill_between(
        [end_data / 7, length / 7],
        [0.4, 0.4],
        [0.2, 0.2],
        step="pre",
        alpha=0.3,
        color=color_tl,
    )
    ax.text(
        ((length - end_data) / 2 + end_data) / 7,
        0.3,
        "Pop. based \nallocation",
        ha="center",
        va="center",
    )

    ax.fill_between(
        [0, length / 7], [0.2, 0.2], [0.0, 0.0], step="pre", alpha=0.2, color=color_tl
    )
    ax.text(length / 7 / 2, 0.1, "NPIs active", ha="center", va="center")
    ax.set_xlabel("Weeks")
    ax.set_title("Time course")
    # ax = fig.add_subplot(gs[0, :2])
    # ax.set_title(title_setup)
    # ax.get_xaxis().set_visible(False)
    # ax.get_yaxis().set_visible(False)
    # ax.table(cellText=[[1,1],[2,2]], loc='upper center',
    #               rowLabels=['Alpha \nvariant','Delta \nvariant'],
    #               colLabels=['Vaccine mRNA','Vaccine \n vector'],
    #               colLoc="center", rowLoc = "center", colWidths=[0.2,0.2],
    #              )

    count_plot += 1
    ax = fig.add_subplot(gs[0, 1])

    susceptible = (
        start_population["susceptible_countryA_vac0_t0"]
        + start_population["susceptible_countryB_vac0_t0"]
        + start_population["susceptible_countryC_vac0_t0"]
        + start_population["susceptible_countryD_vac0_t0"]
    )
    infected = (
        infectious_t0["infectious_countryA_vac0_virus1_t0"]
        + infectious_t0["infectious_countryB_vac0_virus1_t0"]
        + infectious_t0["infectious_countryC_vac0_virus1_t0"]
        + infectious_t0["infectious_countryD_vac0_virus1_t0"]
    )
    recovered = (
        recovered_t0["recovered_countryA_vac0_virus1_t0"]
        + recovered_t0["recovered_countryB_vac0_virus1_t0"]
        + recovered_t0["recovered_countryC_vac0_virus1_t0"]
        + recovered_t0["recovered_countryD_vac0_virus1_t0"]
    )

    sums = {"susceptible": susceptible, "infectious": infected, "recovered": recovered}
    dicts = [start_population, infectious_t0, recovered_t0]
    category_names = ["Belgium", "France", "Germany", "Uk"]
    states = ["susceptible", "infectious", "recovered"]
    results = {}
    for i in range(len(states)):
        d = dicts[i]
        ph_str = ""
        if states[i] != "susceptible":
            ph_str = "_virus1"
        results[states[i].capitalize()] = np.round(
            np.array(
                [
                    d[f"{states[i]}_countryA_vac0{ph_str}_t0"],
                    d[f"{states[i]}_countryB_vac0{ph_str}_t0"],
                    d[f"{states[i]}_countryC_vac0{ph_str}_t0"],
                    d[f"{states[i]}_countryD_vac0{ph_str}_t0"],
                ]
            )
            / sums[states[i]],
            2,
        )

    labels = list(results.keys())
    data = np.array(list(results.values()))
    data_cum = data.cumsum(axis=1)
    category_colors = plt.get_cmap("RdYlGn")(np.linspace(0.15, 0.85, data.shape[1]))
    ax.invert_yaxis()
    ax.xaxis.set_visible(False)
    ax.set_xlim(0, np.sum(data, axis=1).max())
    ax.set_title("Relative initial populations")

    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        widths = data[:, i]
        starts = data_cum[:, i] - widths
        rects = ax.barh(
            labels, widths, left=starts, height=0.5, label=colname, color=color
        )

        r, g, b, _ = color
        text_color = "black" if r * g * b < 0.5 else "darkgrey"
        if colname != "Belgium":
            ax.bar_label(
                rects,
                label_type="center",
                fmt="%.2f%%",
                color=text_color,
                fontsize="small",
                padding=0,
            )
        # ax.legend(ncol=2, fontsize="small")
        ax.legend(ncol=len(category_names), bbox_to_anchor=(0, -0.2), loc="lower left")

    ax.text(
        -0.05,
        letter_y,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=letter_size,
    )

    count_plot += 1
    ax = fig.add_subplot(gs[0, 2])
    ax.plot(
        vaccine_available["t"] / 7,
        vaccine_available["vac1"] / scale,
        color=color_vac1,
        label=label_vac1,
    )
    ax.set_ylabel("Doses per day \nin millions")
    ax.fill_between(
        vaccine_available["t"] / 7,
        vaccine_available["vac1"] / scale,
        color=color_vac1,
        step="pre",
        alpha=0.25,
    )
    ax.plot(
        vaccine_available["t"] / 7,
        vaccine_available["vac2"] / scale,
        color=color_vac2,
        label=label_vac2,
    )
    ax.fill_between(
        vaccine_available["t"] / 7,
        vaccine_available["vac2"] / scale,
        color=color_vac2,
        step="pre",
        alpha=0.25,
    )
    ax.set_xlabel("Weeks")
    ax.set_title(title_vac)
    ax.legend()
    ax.text(
        -0.05,
        letter_y,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=letter_size,
    )

    count_plot += 1
    ax = fig.add_subplot(gs[0, 3])
    barWidth = 0.25

    # set heights of bars
    vaccine1 = [
        delta["delta_vac1_virus1"],
        delta["delta_vac1_virus2"],
        omega["omega_vac1_virus1"],
        omega["omega_vac1_virus2"],
    ]
    vaccine2 = [
        delta["delta_vac2_virus1"],
        delta["delta_vac2_virus2"],
        omega["omega_vac2_virus1"],
        omega["omega_vac2_virus2"],
    ]

    # Set position of bar on X axis
    r1 = np.array([0, 0.75, 3, 3.75])
    r2 = np.array([x + barWidth for x in r1])

    # Make the plot
    ax.bar(
        r1, vaccine1, color="#727272", width=barWidth, edgecolor="white", label="mRNA"
    )
    ax.bar(
        r2, vaccine2, color="#cd7058", width=barWidth, edgecolor="white", label="vector"
    )

    midpoints = (r1 + r2) / 2
    # Add xticks on the middle of the group bars
    # ax.set_xlabel('group', fontweight='bold')
    ax.set_xticks(midpoints)
    ax.set_xticklabels(["Alpha", "Delta", "Alpha", "Delta"])
    point1 = (midpoints[1] - midpoints[0]) / 2 + midpoints[0]
    point2 = (midpoints[3] - midpoints[2]) / 2 + midpoints[2]
    ax.text(point1, 1.1, "Infection \nprotection", ha="center", va="center")
    ax.text(point2, 1.1, "Death \nprotection", ha="center", va="center")
    ax.set_ylim([0, 1.2])
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.set_ylabel("Reduction in %")
    ax.set_title("Vaccine infection and \ndeath reduction")
    ax.legend(loc="center")
    ax.text(
        -0.05,
        letter_y,
        chr(count_plot),
        horizontalalignment="center",
        verticalalignment="center",
        transform=ax.transAxes,
        weight="bold",
        size=letter_size,
    )

    for j in range(len(countries)):
        ax = fig.add_subplot(gs[1, j])
        count_plot += 1
        ax.text(
            -0.05,
            letter_y,
            chr(count_plot),
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            weight="bold",
            size=letter_size,
        )
        if j == 0:
            ax.set_ylabel("Degree of NPIs")
        country = areas[j]
        array = np.array([par_R[x] for x in par_R.keys() if country in x])
        spline_R = get_spline(
            array,
            periods=number_xx_R - 1,
            length=length / number_xx_R,
            total_length=length,
            grid_points=total_grid,
            transform=True,
        )

        ax.set_xlabel("Weeks")
        ax.plot(grid_data / 7, 1 - spline_R[0 : (len(grid_data))], color="C2")
        ax.plot(
            grid_sim / 7,
            1 - spline_R[(len(grid_data) + 1) : (total_grid)],
            color="C2",
            linestyle="dotted",
        )
        ax.text(text_x, text_y, text_str)
        ax.text(text_lockdown_x, text_lockdown_y, text_lockdown_str)
        ax.set_ylim(ylim)
        ax.set_title(countries[j].capitalize())
        ax.fill_between(
            grid_data / 7,
            1 - spline_R[0 : (len(grid_data))],
            step="pre",
            alpha=0.4,
            color="C2",
        )
        ax.fill_between(
            grid_sim / 7,
            1 - spline_R[(len(grid_data) + 1) : (total_grid)],
            step="pre",
            alpha=0.25,
            color="C2",
        )
        ax.axvline(
            list(df_inf_true.index)[-1] / 7,
            0,
            1,
            color="firebrick",
            linestyle="dashed",
            label="Last optimized \nspline point",
            linewidth="0.7",
        )
        if j == 0:
            ax.legend(loc="upper right")

    for j in range(len(countries)):
        ax = fig.add_subplot(gs[2, j])
        count_plot += 1
        ax.text(
            -0.05,
            letter_y,
            chr(count_plot),
            horizontalalignment="center",
            verticalalignment="center",
            transform=ax.transAxes,
            weight="bold",
            size=letter_size,
        )
        if j == 0:
            ax.set_ylabel("Active cases \nin millions")

        ax.plot(
            grid_data / 7,
            df_infected.loc[0 : (len(grid_data) - 1), areas[j]] / scale,
            label="Simulated",
            color="C0",
        )
        ax.plot(
            grid_sim / 7,
            df_infected.loc[len(grid_data) : (total_grid), areas[j]] / scale,
            color="C0",
            linestyle="dotted",
        )

        ax.plot(
            np.array(df_inf_true.index) / 7,
            df_inf_true[countries[j]] / scale,
            label="Data",
            color="C1",
        )
        ax.set_title(countries[j].capitalize())
        ax.set_xlabel("Weeks")
        ax.axvline(
            list(df_inf_true.index)[-1] / 7,
            0,
            1,
            color="firebrick",
            linestyle="dashed",
            label="Last optimized \nspline point",
            linewidth="0.7",
        )
        if j == 0:
            ax.legend()

    fig.savefig(
        "/home/manuel/Documents/VaccinationDistribution/paper/images/infected_compare",
        bbox_inches="tight",
    )
    return vaccine_available

def stacked_bar(results, category_names, ax, ylabel=True, legend=True, map_col = 'RdYlGn'):
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
    
    category_colors1 = plt.get_cmap("Set1")(
        np.linspace(0.15, 0.85, 5))[1]
    category_colors2 = plt.get_cmap("Set2")(
        np.linspace(0.15, 0.85, 5))[0]
    category_colors=[category_colors1, category_colors2]
    
    
    ax.invert_yaxis()
    ax.xaxis.set_visible(False)
    ax.set_xlim(0, np.sum(data, axis=1).max())

    for i, (colname, color) in enumerate(zip(category_names, category_colors)):
        widths = data[:, i]
        starts = data_cum[:, i] - widths
        rects = ax.barh(labels, widths, left=starts, height=0.5,
                        label=colname, color=color)

        #r, g, b, _ = color
        text_color = 'black' 
        ax.bar_label(rects, label_type='center', color=text_color)
    if ylabel is False:
        ax.yaxis.set_visible(False)
    if legend is True:
        ax.legend(ncol=len(category_names), bbox_to_anchor=(0.5, 0.48),
                  loc='center', fontsize='small')

def plot_bars_vac(ax, vaccine1, vaccine2, title, barWidth=0.25):
    r1 = np.array([0, 0.75, 3, 3.75])
    r2 = np.array([x + barWidth for x in r1])
    
        # Make the plot
    ax.bar(
            r1, vaccine1, color="#727272", width=barWidth, edgecolor="white", label="Vaccine \none"
        )
    if not(vaccine2 is None):
        ax.bar(
                r2, vaccine2, color="#cd7058", width=barWidth, edgecolor="white", label="VAccine \ntwo"
            )
    
    midpoints = (r1 + r2) / 2
        # Add xticks on the middle of the group bars
        # ax.set_xlabel('group', fontweight='bold')
    ax.set_xticks(midpoints)
    ax.set_xticklabels(["Wild \ntype", "Mutant", "Wild \ntype", "Mutant"])
    point1 = (midpoints[1] - midpoints[0]) / 2 + midpoints[0]
    point2 = (midpoints[3] - midpoints[2]) / 2 + midpoints[2]
    ax.text(point1, 1.2, "Infection \nprotection", ha="center", va="center")
    ax.text(point2, 1.2, "Death \nprotection", ha="center", va="center")
    ax.set_ylim([0, 1.2])
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.set_ylabel(title)
    #ax.set_title("Vaccine infection and \ndeath reduction")
    ax.legend(loc="center")

#---------------------------------------------------------------------------------------------------------------
def plot_horizontal_bars(results, category_names, ax):
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
        category_colors = ["C0", "C1"]
    
        ax.invert_yaxis()
        ax.xaxis.set_visible(False)
        ax.set_xlim(0, np.sum(data, axis=1).max())
    
        for i, (colname, color) in enumerate(zip(category_names, category_colors)):
            widths = data[:, i]
            starts = data_cum[:, i] - widths
            rects = ax.barh(labels, widths, left=starts, height=0.5,
                            label=colname, color=color, alpha=0.9)
    
            #r, g, b, _ = color
            text_color = 'black' # 'darkgrey'
            ax.bar_label(rects, label_type='center', color=text_color)
        ax.legend(ncol=len(category_names), bbox_to_anchor=(0.6, 0),
                  fontsize='small')
    
        return ax


def plot_horizontal_bars_annotated(ax,
                                   dict_use,
                                   scale = 10**5,
                                   color_vline = "black",
                                   linestyle_vline = "dashed", title=""):
    all_results = dict_use["all_strategies"]
    argmin_global = np.argmin(all_results["fval"])
    global_optimum = all_results.iloc[argmin_global]
    
    pareto_results = dict_use["pareto_improvements"]
    argmin_Pareto = np.argmin(pareto_results["fval"])
    pareto_optimum = all_results.iloc[argmin_Pareto]
    
    optimal = np.round([global_optimum["countryA"]/scale, global_optimum["countryB"]/scale],2)
    population = np.round(np.array(dict_use["population_based"])/scale,2)
    pareto = np.round([pareto_optimum["countryA"]/scale, pareto_optimum["countryB"]/scale],2)
    
    category_names = ["Country A", "Country B"]
    results = {
        'Optimal \nStrategy': optimal,
        'Population \nStrategy':  population,
        'Pareto Optimal \nStrategy': pareto,
    }
    
    

    plot_horizontal_bars(results, category_names, ax)
    ax.axvline(x=np.sum(population), color = color_vline, linestyle =linestyle_vline)
    y_ticks = ax.get_yticks()
    #ax.annotate(f"{np.round((np.sum(optimal) / np.sum(population) - 1)*100, 2)}",
    #            xy=(np.sum(optimal), y_ticks[0]), xycoords='data',
    #            xytext=(np.sum(population), y_ticks[0]), #textcoords='offset points',
    #            ha="center", fontsize = 4, va = "center",
    #            arrowprops=dict(arrowstyle="->"),)
    par_opt_list = [optimal, pareto]
    for index in range(len(par_opt_list)):
        if index == 1:
            tick = 2
        else:
            tick = index
        ax.arrow(np.sum(population), y_ticks[tick], np.sum(par_opt_list[index]) - np.sum(population), 0,
                 length_includes_head=True, head_width=0.1, head_length=0.05)     
        ax.text(np.sum(population), y_ticks[tick], f"{np.round((np.sum(par_opt_list[index])/np.sum(population) - 1)*100,2)}%",
                va = "center")    
        ax.set_title(title)
    
def compute_incidences(trajectories,
                       viruses = ["virus1", "virus2"],
                       countries = ["countryA", "countryB"],
                       time = np.linspace(0, 140, 6000),
                       lambda1 = 0.1,
                       habitant_scale = 0.1,):
    
    infected = {}
    for index in range(len(viruses)):
        df_help = pd.DataFrame(np.nan, index=range(trajectories.shape[0]), columns=countries)
        for index_country in range(len(countries)):
            cols = [x for x in trajectories.columns if "infectious" in x and viruses[index] in x and countries[index_country] in x]
            df_help[countries[index_country]] = trajectories[cols].sum(axis=1)
        infected[viruses[index]] = df_help
    
    incidences = {}
    for index in range(len(viruses)):
        df_incidence = pd.DataFrame(np.nan, index=range(trajectories.shape[0]), columns=countries)
        for index_country in range(len(countries)):
            time_course = infected[viruses[index]][countries[index_country]]
            for index_time in range(1, len(time_course)-1):
                delta_t = time[index_time+1] - time[index_time]
                newly_infected = time_course[index_time + 1] - (1 - lambda1*delta_t) * time_course[index_time] 
                df_incidence.loc[index_time, countries[index_country]] = newly_infected
        incidences[viruses[index]] = df_incidence
        
    seven_day_incidences = {}
    for index in range(len(viruses)):
        df_7_day_incidence = pd.DataFrame(np.nan, index=range(trajectories.shape[0]), columns=countries)
        for index_country in range(len(countries)):
            time_course_incidence = incidences[viruses[index]]
            for index_time in range(300, len(time_course)-1):
                df_7_day_incidence.loc[index_time, countries[index_country]] = time_course_incidence.loc[(index_time - 300):index_time, countries[index_country]].sum()
        seven_day_incidences[viruses[index]] = df_7_day_incidence * habitant_scale
    
    return seven_day_incidences

def plot_incidences(incidence, time,
                    countries = ["countryA", "countryB"],
                    ax =None,
                    label_countries = ["Country A", "Country B"],
                    colors = ["C0", "C1"],
                    alpha = 0.3):

    viruses = list(incidence.keys())
    for index_country in range(len(countries)):
    
        sum_infections = 0
        for index in range(len(viruses)):
            sum_infections += incidence[viruses[index]][countries[index_country]] 
        ax.plot(time/7, sum_infections, label = label_countries[index_country],
                color=colors[index_country])
        ax.fill_between(time/7, np.repeat(0, len(sum_infections)), sum_infections,
                                          color=colors[index_country], alpha = alpha)

def plot_incidences_country(ax, trajectories, time,  incidences, viruses = ["virus1", "virus2"],
                            index_country = "countryA",
                            label_type = ["Optimal", "Pareto optimal", "Population\nbased"],
                            colors = ["C0", "C1", "black"], alpha=0.3):
    
    for key in range(len(incidences.keys())):
        incidence = incidences[trajectories[key]]
        sum_infections = 0
        for index in range(len(viruses)):
            sum_infections += incidence[viruses[index]][index_country]
        
        if "pop" in trajectories[key]:
            ax.plot(time/7, sum_infections, label = label_type[key],
                        color=colors[key], linestyle = "dashed")
        else:
            ax.plot(time/7, sum_infections, label = label_type[key],
                        color=colors[key])
            ax.fill_between(time/7, np.repeat(0, len(sum_infections)), sum_infections,
                                                      color=colors[key], alpha = alpha)

def compute_incidences_countries(dicts, time,  trajectories, name = "initalUnequal_vacUnequal_nvacc_60000",
                                 viruses=["virus1", "virus2"],
                                 countries = ["countryA", "countryB"],
                                 lambda1 = 0.1, 
                                 habitant_scale = 0.01):

    incidences = {}
    for index in range(len(trajectories)):
        incidences[trajectories[index]] = compute_incidences(trajectories = dicts[name][trajectories[index]],
                                                             viruses = viruses,
                                                             countries = countries,
                                                             time =time,
                                                             lambda1 =lambda1,
                                                             habitant_scale = habitant_scale,)
    return incidences



def compute_deceased(results, country):
    pareto_deceased = []
    optimal_deceased = []
    pop_deceased = []
    for index in range(len(results)):
        result = results[index]
        pareto_result = result["trajectories_pareto"]
        optimal_result = result["trajectories_best"]
        population_result = result["trajectories_pop"]
        cols = [x for x in pareto_result.columns if "dead" in x and country in x]
        
        deceased_pareto = list(pareto_result[cols].sum(axis=1))
        deceased_optimal = list(optimal_result[cols].sum(axis=1))
        deceased_pop = list(population_result[cols].sum(axis=1))
        
        pareto_deceased.append(deceased_pareto[-1])
        optimal_deceased.append(deceased_optimal[-1])
        pop_deceased.append(deceased_pop[-1])
        
    out = {"pareto" : pareto_deceased,
           "optimal": optimal_deceased,
           "pop" : pop_deceased,}
    return out
def compute_splines_from_results(type_opti = "pareto_improvements",
                                 vac = "vac2",
                                 periods = 10,
                                 length = 14,
                                 total_length = 140,
                                 grid_points = 6000,
                                 add_additional = {"integers" : [9,10],
                                                     "number" : 0.5},n_vaccs=None,results=None, ):
    fractions = []
    for index in range(len(n_vaccs)):
        result = results[index]
        n_vacc = n_vaccs[index]
        df_optimal = result[type_opti]
        if not(add_additional is None):
            for index in add_additional["integers"]:
                for index2 in ["countryA"]:
                    for index3 in ["vac1", "vac2"]:
                        name = f"yy_{index2}_{index3}_{index}"
                        df_optimal[name] = add_additional["number"]
        
        cols = [x for x in df_optimal.columns if vac in x]
        argmin = np.argmin(df_optimal["fval"])
        
        yy_points = pd.Series(list(df_optimal.iloc[argmin][cols]))
        time = np.linspace(0, total_length, grid_points) / 7
        spline = get_spline(yy_points,
                            periods=periods,
                            length=length,
                            total_length=total_length,
                            grid_points=grid_points,
                        )
        area = trapz(spline*n_vacc/2, dx=(time[1] - time[0]))
        total_area = n_vacc*10
        fraction_country_A = area/total_area
        fractions.append(fraction_country_A)
        
    return fractions