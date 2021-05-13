import numpy as np
import pandas as pd

from functions.run_sbml import run_model_once
from functions.run_sbml import set_all_initial_conditions_to_zero
from functions.run_sbml import use_end_values_as_start_values

from visualization.model_results import get_substates


def run_model_stepwise_vaccines(
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

    # TODO generalize-------------------------------------
    country_A_vac1 = theta["value"][0:periods]
    country_A_vac2 = theta["value"][(periods) : (2 * periods)]
    # ----------------------------------------------------

    if set_initials_zero is True:
        model = set_all_initial_conditions_to_zero(model)

    timepoints = np.linspace(1, length_periods, length_periods)
    periods_minus_one = periods - 1

    # TODO generalize-------------------------------------
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
            # TODO generalize-------------------------------------
            proportion_dict_start_it = {
                "proportion_countryA_vac1": country_A_vac1[period_index + 1],
                "proportion_countryB_vac1": 1 - country_A_vac1[period_index + 1],
                "proportion_countryA_vac2": country_A_vac2[period_index + 1],
                "proportion_countryB_vac2": 1 - country_A_vac2[period_index + 1],
            }
            # ----------------------------------------------------
            set_parameter_it = {**start_values_it, **set_parameter, **proportion_dict_start_it}
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

    df_trajectories_states = trajectory_dict["states"]
    substates = get_substates(model=model, substrings=state_type, include_all=True)
    df_substates = df_trajectories_states[substates]

    if final_amount is True:
        df_sum = df_substates.iloc[[-1]]
    else:
        df_sum = df_substates

    sum_states = df_sum.values.sum()

    return sum_states


def run_model_stepwise_vaccines_sum(
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
    trajectory_dict = run_model_stepwise_vaccines(
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

    out = {
        "value": sum_states,
    }
    return out
