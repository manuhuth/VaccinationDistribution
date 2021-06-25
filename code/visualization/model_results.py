import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update(
    {"text.usetex": True, "font.family": "sans-serif", "font.sans-serif": ["Helvetica"]}
)

from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D


def plot_states(
    results,
    model,
    xlabel="$t$ (d)",
    ylabel="$x_i(t)$",
    title="State trajectories",
    state_ids=None,
    time_name="time",
    legend_location="upper left",
    legend_next_to_plot=False,
    ncol=None,
    colors=None,
    custom_label=None,
    ylim=None,
):
    """Plot trajectories of species. Only a part of species can be addressed
    by using the state_ids parameter.

    Parameters
    ----------
    results : numpy.ReturnDataView
        Created by :func:`model_run`.

    model : libsbml.Model
        SBML model to retrieve the names of the species.

    xlabel : str
        Label for the x-axis.

    ylabel : str
        Label for the y-axis.

    title : str
        Label of title.

    state_ids : list of strings
        List of substrings for which the model should be printed. Every
        trajectory of a species is printed if one of the substrings occur in
        their name.

    Returns
    -------
    fig, ax : figure.Figure, AxesSubplot

    """
    fig, ax = plt.subplots()
    df_trajectories = results["states"]
    time = results["observables"][time_name]

    if state_ids is None:
        states = df_trajectories.columns
    else:
        states = state_ids

    index_cols = 0
    for index in states:

        if custom_label is not None:
            label = custom_label[index_cols]
        else:
            label = index

        if colors is not None:
            ax.plot(time, df_trajectories[index], label=label, color=colors[index_cols])
        else:
            ax.plot(time, df_trajectories[index], label=label)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if ylim is not None:
            ax.set_ylim(ylim)
        if ncol is not None:
            ax.legend(loc=legend_location, ncol=ncol)
        else:
            ax.legend(loc=legend_location)
        ax.set_title(title)
        index_cols += 1

    if legend_next_to_plot is True:
        if ncol is not None:
            ax.legend(bbox_to_anchor=(1, 0.8, 0.3, 0.2), loc=legend_location, ncol=ncol)
        else:
            ax.legend(bbox_to_anchor=(1, 0.8, 0.3, 0.2), loc=legend_location)
    ax.spines["top"].set_alpha(0)
    ax.spines["bottom"].set_alpha(0.3)
    ax.spines["right"].set_alpha(0)
    ax.spines["left"].set_alpha(0.3)
    ax.grid(alpha=0.3)
    ax.xaxis.set_ticks_position("none")
    ax.yaxis.set_ticks_position("none")
    return fig, ax


def get_state_trajectory_data_frame(results, model):
    """Create data frame of state trajectories of species.

    Parameters
    ----------
    results : numpy.ReturnDataView
        Created by :func:`model_run`.

    model : libsbml.Model
        SBML model to retrieve the names of the species.

    Returns
    -------
    df_trajectories : pandas.DataFrame
        Data frame of state trajectories of species

    """

    states = model.getStateIds()
    trajectories = results["x"]
    df_trajectories = pd.DataFrame(trajectories, columns=(states))
    return df_trajectories


def plot_observables(
    results,
    model,
    xlabel="Days",
    ylabel="Observables",
    title="Observables trajectories",
    observable_ids=None,
    set_off_scientific_notation=False,
    decimal_floats=4,
    time_name="time",
    legend_next_to_plot=False,
    legend_location="upper left",
    colors=None,
    custom_label=None,
):
    """Plot trajectories of species. Only a part of species can be addressed
    by using the state_ids parameter.

    Parameters
    ----------
    results : numpy.ReturnDataView
        Created by :func:`model_run`.

    model : libsbml.Model
        SBML model to retrieve the names of the species.

    xlabel : str
        Label for the x-axis.

    ylabel : str
        Label for the y-axis.

    title : str
        Label of title.

    state_ids : list of strings
        List of substrings for which the model should be printed. Every
        trajectory of a species is printed if one of the substrings occur in
        their name.

    set_off_scientific_notation : {True, False}
        If True, the numbers of the y-axis are not displayed in scientific
        format in any case.

    decimal_floats : integer
        Number of decimal floats dispalyed at the y-axis if.
        set_off_scientific_notation` is True.

    Returns
    -------
    fig, ax : figure.Figure, AxesSubplot

    """
    fig, ax = plt.subplots()

    df_trajectories = results["observables"]
    df_trajectories.columns = df_trajectories.columns.str.lstrip("observable_")

    if observable_ids is None:
        observables = df_trajectories.columns
    else:
        observables = observable_ids

    # redefinition prevents graphs from applying scientifix notation
    # if unnecessary.
    df_observables_to_be_plotted = df_trajectories[observables]
    index_cols = 0
    for index in observables:
        if custom_label is not None:
            label = custom_label[index_cols]
        else:
            label = index

        if colors is not None:
            ax.plot(
                df_trajectories[time_name],
                df_observables_to_be_plotted[index],
                label=custom_label[index_cols],
                color=colors[index_cols],
            )
        else:
            ax.plot(
                df_trajectories[time_name],
                df_observables_to_be_plotted[index],
                label=label,
            )
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend()
        ax.set_title(title)
        index_cols += 1

    if legend_next_to_plot is True:
        ax.legend(bbox_to_anchor=(1, 0.8, 0.3, 0.2), loc=legend_location)

    if set_off_scientific_notation is True:
        ax.get_yaxis().get_major_formatter().set_useOffset(False)
        ax.yaxis.set_major_formatter(FormatStrFormatter(f"%.{decimal_floats}f"))

    ax.spines["top"].set_alpha(0)
    ax.spines["bottom"].set_alpha(0.3)
    ax.spines["right"].set_alpha(0)
    ax.spines["left"].set_alpha(0.3)
    ax.grid(alpha=0.3)
    ax.xaxis.set_ticks_position("none")
    ax.yaxis.set_ticks_position("none")
    return fig, ax


def get_observable_trajectory_data_frame(results, model):
    """Create data frame of state trajectories of observables.

    Parameters
    ----------
    results : numpy.ReturnDataView
        Created by :func:`model_run`.

    model : libsbml.Model
        SBML model to retrieve the names of the species.

    Returns
    -------
    df_trajectories : pandas.DataFrame
        Data frame of state trajectories of observables

    """

    observables = model.getObservableNames()
    trajectories = results["y"]
    df_trajectories = pd.DataFrame(trajectories, columns=(observables))
    return df_trajectories


def get_observables_by_name(observables, substrings, include_all=True, key="name"):
    observables_names = []
    for index_obs in observables:
        # value = observables[index_dicts][key]
        if include_all is True:
            if all(x in index_obs for x in substrings):
                observables_names.append(index_obs)
        else:
            if any(x in index_obs for x in substrings):
                observables_names.append(index_obs)
    return observables_names


def get_substates(model, substrings, include_all=True):
    """Get states that contain certain substrings.

    Parameters
    ----------
    model : libsbml.Model
        SBML model to retrieve the names of the species.

    substrings : list of strings
        Substrings that occur in the name of the species that should be
        addressed.

    include_all : {True, False}
        If `True`, only states that include all substrings are returned. If
        `False` all states that include one substring are returned.

    Returns
    -------
    substates : list
        List containing all species' names that have the substrings in their
        names.
    """

    states = model.getStateIds()

    if include_all is False:
        substates = list(
            set(get_states_by_substrings_any(states=states, substrings=substrings))
        )
    else:
        substates = get_states_by_substrings_all(states=states, substrings=substrings)

    return substates


def get_states_by_substrings_all(states, substrings):
    """Get state names that contain all substrings.

    Parameters
    ----------
    states : list of strings
        Names of all species in the model.

    substrings : list of strings
        Substrings that occur in the name of the species that should be
        addressed.

    Returns
    -------
    state : list
        List containing all species' names that have the substrings in their
        names.

    """
    state = []
    for index_states in states:
        if all([x in index_states for x in substrings]):
            state.append(index_states)

    return state


def get_states_by_substrings_any(states, substrings):
    """Get state names that contain any substring.

    Parameters
    ----------
    states : list of strings
        Names of all species in the model.

    substrings : list of strings
        Substrings that occur in the name of the species that should be
        addressed.

    Returns
    -------
    state : list
        List containing all species' names that have the substrings in their
        names.

    """
    state = []
    for index in substrings:
        states_removed = [x for x in states if index in x]
        state += states_removed

    return state


def plot_3D_function(
    function,
    xlabel="$vac_1$",
    ylabel="$vac_2$",
    zlabel="infectious",
    set_off_scientific_notation=False,
    decimal_floats=3,
    start_linspace_x=0,
    start_linspace_y=0,
    end_linspace_x=1,
    end_linspace_y=1,
):
    """Plot the function with respect to two input parameters.

    Parameters
    ----------
    function : function
        Function that is used to map the inputs to the output. Must return a
        dictionary with the desired quantity that has value as key.

    number_parameter : int
        Number of parameter.

    xlabel : str
        Label of x-axis.

    ylabel : str
        Label of y-axis.

    zlabel : str
        Label of z-axis.

    Returns
    -------
    None.

    """

    number_parameter = 2
    x1 = np.linspace(start_linspace_x, end_linspace_x)
    x2 = np.linspace(start_linspace_y, end_linspace_y)
    xvalues = []
    array = []
    for i in range(0, len(x1)):
        for j in range(0, len(x2)):
            array = [x1[i], x2[j]]
            xvalues.append(array)

    yvalues = np.linspace(0, 0, len(xvalues))
    for i, j in zip(xvalues, range(0, len(xvalues))):
        theta_df = pd.DataFrame(
            i, columns=["value"], index=[f"theta_{i}" for i in range(number_parameter)]
        )
        yvalues[j] = function(theta=theta_df)["value"]

    xvalues = np.array(xvalues)
    xvalues1 = xvalues[:, 0]
    xvalues2 = xvalues[:, 1]

    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(xvalues1, xvalues2, yvalues, c=yvalues, cmap="viridis", linewidth=0.05)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    if set_off_scientific_notation is True:
        ax.get_zaxis().get_major_formatter().set_useOffset(False)
        ax.zaxis.set_major_formatter(FormatStrFormatter(f"%.{decimal_floats}f"))
