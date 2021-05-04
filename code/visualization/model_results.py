import pandas as pd
import matplotlib.pyplot as plt


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


def get_substates(model, substrings):
    """Get states that contain certain substrings.

    Parameters
    ----------
    model : libsbml.Model
        SBML model to retrieve the names of the species.

    substrings : list of strings
        Substrings that occur in the name of the species that should be
        addressed.

    Returns
    -------
    substates : list
        List containing all species' names that have the substrings in their
        names.
    """

    states = model.getStateIds()
    substates = get_states_by_substrings(states=states, substrings=substrings)
    return substates


def get_states_by_substrings(states, substrings):
    """Get state names that contain substrings.
    
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
