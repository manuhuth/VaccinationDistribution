import pandas as pd
import matplotlib.pyplot as plt


def plot_states(
    results,
    model,
    xlabel="$t$ (d)",
    ylabel="$x_i(t)$",
    title="State trajectories",
    state_ids=None,
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
    df_trajectories = get_state_trajectory_data_frame(results, model)

    if state_ids is None:
        states = df_trajectories.columns
    else:
        states = state_ids

    for index in states:
        label = index
        ax.plot(results["t"], df_trajectories[index], label=label)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend()
        ax.set_title(title)

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
    xlabel="$t$ (d)",
    ylabel="$x_i(t)$",
    title="Observables trajectories",
    observable_ids=None,
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
    df_trajectories = get_observable_trajectory_data_frame(results, model)

    if observable_ids is None:
        observables = df_trajectories.columns
    else:
        observables = observable_ids

    for index in observables:
        label = index
        ax.plot(results["t"], df_trajectories[index], label=label)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.legend()
        ax.set_title(title)

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


def get_substates(model, substrings, include_all = True):
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
        substates = list(set(get_states_by_substrings_any(states=states, substrings=substrings)))
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
