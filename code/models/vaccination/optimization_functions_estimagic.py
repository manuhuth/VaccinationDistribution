import numpy as np
import pandas as pd

from functions.run_sbml import run_model_once
from functions.run_sbml import set_all_initial_conditions_to_zero
from functions.run_sbml import use_end_values_as_start_values

from visualization.model_results import get_substates


def run_model_stepwise_vaccines_estimagic(
    theta,
    model,
    solver,
    set_start_parameter,
    set_parameter,
    observables_names,
    length_periods,
    periods,
    set_initials_zero,
):
    """
    Run model that allows to stepwise adapt the distribution of vaccines.
    After each period the proportions of vaccinations are set to the respective
    value in theta. The first `periods` values in the column `value` of the
    data frame theta correspond to vaccine one, the rest to vaccine two.

    Parameters
    ----------
    theta : pandas.DataFrame
        Data frame to has a column value. The first `periods` values in the
        column `value` of the data frame theta correspond to vaccine one,
        the rest to vaccine two.

    model : amici.model
        Model to use.

    solver : amici.solver
        Solver to use.

    set_start_parameter : dict
        Dictionary including the names of the start parameters as keys
        and the respective start parameters as values. If not specified
        zero is used.

    set_parameter : dict
        Dictionary including the names of the parameters as keys
        and the respective parameters as values. If a parameter is not specified
        the default from :func:`create_sbml` is used.

    observables_names : list
        List of the observables.

    length_periods : int
        Length of each period over wich the fraction of vaccines is decided on.

    periods : int
        Number of decision periods.

    set_initials_zero : {True, False}
        If true, all parameters are reset after running. Should be set to True.
        Otherwise it can be possible that some results are still in the cache.

    Returns
    -------
    trajectory_dict : dict
        Dictionary of state trajectories.

    """

    # TODO generalize over vaccines and countries-----------------------------
    country_A_vac1 = theta["value"][0:periods]
    country_A_vac2 = theta["value"][(periods) : (2 * periods)]
    # ----------------------------------------------------

    if set_initials_zero is True:
        model = set_all_initial_conditions_to_zero(model)

    timepoints = np.linspace(1, length_periods, length_periods)
    periods_minus_one = periods - 1

    # TODO generalize over vaccines and countries-----------------------------
    proportion_dict_start = {
        "proportion_countryA_vac1": country_A_vac1[0],
        "proportion_countryB_vac1": 1 - country_A_vac1[0],
        "proportion_countryA_vac2": country_A_vac2[0],
        "proportion_countryB_vac2": 1 - country_A_vac2[0],
    }
    # ----------------------------------------------------

    set_parameter_first = {
        **set_start_parameter,
        **set_parameter,
        **proportion_dict_start,
    }

    first_period = run_model_once(
        model=model,
        solver=solver,
        timepoints=timepoints,
        set_parameter=set_parameter_first,
    )

    states = model.getStateIds()
    df_trajectories_states = pd.DataFrame(first_period["x"], columns=(states))
    df_trajectories_observables = pd.DataFrame(
        first_period["y"], columns=(observables_names)
    )

    if periods_minus_one > 0:
        for period_index in range(periods_minus_one):
            start_values_it = use_end_values_as_start_values(
                df_trajectories_states=df_trajectories_states
            )
            # TODO generalize over vaccines and countries---------------------
            proportion_dict_start_it = {
                "proportion_countryA_vac1": country_A_vac1[period_index + 1],
                "proportion_countryB_vac1": 1 - country_A_vac1[period_index + 1],
                "proportion_countryA_vac2": country_A_vac2[period_index + 1],
                "proportion_countryB_vac2": 1 - country_A_vac2[period_index + 1],
            }
            # ----------------------------------------------------
            set_parameter_it = {
                **start_values_it,
                **set_parameter,
                **proportion_dict_start_it,
            }
            model_result = run_model_once(
                model=model,
                solver=solver,
                timepoints=timepoints,
                set_parameter=set_parameter_it,
            )

            states_iteration = model_result["x"]
            observables_iteration = model_result["y"]

            df_trajectories_states_iteration = pd.DataFrame(
                states_iteration, columns=(states)
            )

            df_trajectories_observables_iteration = pd.DataFrame(
                observables_iteration, columns=(observables_names)
            )

            df_trajectories_states = pd.concat(
                [df_trajectories_states, df_trajectories_states_iteration],
                ignore_index=True,
            )

            df_trajectories_observables = pd.concat(
                [df_trajectories_observables, df_trajectories_observables_iteration],
                ignore_index=True,
            )

    trajectory_dict = {
        "states": df_trajectories_states,
        "observables": df_trajectories_observables,
    }

    return trajectory_dict


def get_sum_of_states(model, trajectory_dict, state_type=["dead"], final_amount=True):
    """Get the sum of output states.

    Parameters
    ----------
    model : amici.model
        Model to use.

    trajectory_dict : dict
        Dictionary of state trajectories.

    state_type : list
        List including the substrings of the states
        that should be used to sum over.

    final_amount : {True, False}
        If True only the final amount is considered. If False the sum of all
        periods is used.

    Returns
    -------
    sum_states : float
        Sum of the states as specified by the function.

    """

    df_trajectories_states = trajectory_dict["states"]
    substates = get_substates(model=model, substrings=state_type, include_all=True)
    df_substates = df_trajectories_states[substates]

    if final_amount is True:
        df_sum = df_substates.iloc[[-1]]
    else:
        df_sum = df_substates

    sum_states = df_sum.values.sum()

    return sum_states


def run_model_stepwise_vaccines_sum_estimagic(
    theta,
    model,
    solver,
    set_start_parameter,
    set_parameter,
    observables_names,
    length_periods,
    periods,
    set_initials_zero,
    state_type=["infectious"],
    final_amount=False,
):
    """
    Run model that allows to stepwise adapt the distribution of vaccines and
    return sum of specified states.
    After each period the proportions of vaccinations are set to the respective
    value in theta. The first `periods` values in the column `value` of the
    data frame theta correspond to vaccine one, the rest to vaccine two.

    Parameters
    ----------
    theta : pandas.DataFrame
        Data frame to has a column value. The first `periods` values in the
        column `value` of the data frame theta correspond to vaccine one,
        the rest to vaccine two.

    model : amici.model
        Model to use.

    solver : amici.solver
        Solver to use.

    set_start_parameter : dict
        Dictionary including the names of the start parameters as keys
        and the respective start parameters as values. If not specified
        zero is used.

    set_parameter : dict
        Dictionary including the names of the parameters as keys
        and the respective parameters as values. If a parameter is not specified
        the default from :func:`create_sbml` is used.

    observables_names : list
        List of the observables.

    length_periods : int
        Length of each period over wich the fraction of vaccines is decided on.

    periods : int
        Number of decision periods.

    set_initials_zero : {True, False}
        If true, all parameters are reset after running. Should be set to True.
        Otherwise it can be possible that some results are still in the cache.

    state_type : list
        List including the substrings of the states
        that should be used to sum over.

    final_amount : {True, False}
        If True only the final amount is considered. If False the sum of all
        periods is used.

    estimagic_dict : {True, False}
        If True the output is in the form of an estimagic dictionary. If False
        the output is a simple float.

    Returns
    -------
    out : dict
        Dictionary of sum of states with value as key such that function can
        be used by estimagic.

    """
    trajectory_dict = run_model_stepwise_vaccines_estimagic(
        theta,
        model,
        solver,
        set_start_parameter,
        set_parameter,
        observables_names,
        length_periods,
        periods,
        set_initials_zero,
    )

    sum_states = get_sum_of_states(
        model=model,
        trajectory_dict=trajectory_dict,
        state_type=state_type,
        final_amount=final_amount,
    )

    out = {"value": sum_states}

    return out
