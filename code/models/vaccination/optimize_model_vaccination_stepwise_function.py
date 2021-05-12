import numpy as np

from models.vaccination.create_model_vaccination import model_vaccination_create_sbml

from visualization.model_results import plot_states
from visualization.model_results import plot_observables
from visualization.model_results import get_observables_by_name

from functions.run_sbml import get_model_and_solver_from_sbml
from functions.run_sbml import run_model
from functions.run_sbml import create_observables_vaccination_rates

from functions.vaccine_proportions import (
    create_rules_vaccination_proportion_relative_population,
)
from functions.vaccine_proportions import (
    create_rules_vaccination_proportion_relative_infected_population,
)


# --------------------------Create Model--------------------------------------
path_sbml = "stored_models/vaccination/step_function_optim"
vaccination_states = ["vac0", "vac1", "vac2"]
non_vaccination_state = "vac0"
virus_states = ["virW", "virM"]
areas = ["countryA", "countryB"]
species_comp = ["susceptible", "infectious", "recovered", "dead"]


model_vaccination_create_sbml(
    path=path_sbml,
    areas=areas,
    distances=np.array([[0, 10], [10, 0]]),
    t0_susceptible=0,
    t0_infectious=1000,
)

observables_nu = create_observables_vaccination_rates(
    vaccination_states_removed=["vac1", "vac2"], areas=areas, name_parameter="nu"
)
observables_proportion = create_observables_vaccination_rates(
    vaccination_states_removed=["vac1", "vac2"],
    areas=areas,
    name_parameter="proportion",
)
observables = {**observables_nu, **observables_proportion}

model_and_solver = get_model_and_solver_from_sbml(
    path_sbml=path_sbml,
    model_name="vaccination",
    model_directory="stored_models/vaccination/vaccination_step_function_optim",
    observables=observables,
)

# -----------------------Run Model---------------------------------------------

import pandas as pd
from scipy.optimize import minimize
from scipy.optimize import Bounds
from functions.run_sbml import run_model_once
from functions.run_sbml import set_all_initial_conditions_to_zero
from functions.run_sbml import use_end_values_as_start_values

from visualization.model_results import get_substates
from functools import partial

model = model_and_solver["model"]
solver = model_and_solver["solver"]
set_start_parameter = {
    "susceptible_countryA_vac0_t0": 60000,
    "susceptible_countryB_vac0_t0": 40000,
    "infectious_countryA_vac0_virW_t0": 1000,
    "infectious_countryA_vac0_virM_t0": 1000,
    "infectious_countryB_vac0_virW_t0": 1000,
    "infectious_countryB_vac0_virM_t0": 1000,
}
set_parameter = {
    "beta": 2,
    "lambda1": 0.01,
    "p": 0.3,
    "number_vac1": 10,
    "number_vac2": 20,
}
observables_names = observables.keys()

length_periods = 7
periods = 2
observables_names = observables.keys()
set_initials_zero=True

theta = np.concatenate((np.repeat(0.4, periods), np.repeat(0.7, periods)), axis=None) 

def run_model_stepwise_vaccines(theta, model, solver, set_start_parameter, set_parameter,
                                observables_names, length_periods, periods, set_initials_zero):
    
    #TODO generalize-------------------------------------
    country_A_vac1 = theta[0:periods]
    country_A_vac2 = theta[(periods):(2*periods)]
    #----------------------------------------------------
    
    if set_initials_zero is True:
        model = set_all_initial_conditions_to_zero(model)
    
    timepoints = np.linspace(1, length_periods, length_periods)
    periods_minus_one = periods - 1
    
    #TODO generalize-------------------------------------
    proportion_dict_start = { 'proportion_countryA_vac1' : country_A_vac1[0],
                        'proportion_countryB_vac1' : 1- country_A_vac1[0],
                        'proportion_countryA_vac2' : country_A_vac2[0],
                        'proportion_countryB_vac2' : 1- country_A_vac2[0]}
    #----------------------------------------------------
    
    set_parameter_first = {**set_start_parameter, **set_parameter, **proportion_dict_start}
    
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
            #TODO generalize-------------------------------------
            proportion_dict_start = { 'proportion_countryA_vac1' : country_A_vac1[period_index + 1],
                        'proportion_countryB_vac1' : 1- country_A_vac1[period_index + 1],
                        'proportion_countryA_vac2' : country_A_vac2[period_index + 1],
                        'proportion_countryB_vac2' : 1- country_A_vac2[period_index + 1]}
            #----------------------------------------------------
            set_parameter_it = {**start_values_it, **set_parameter}
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

def get_sum_of_states(model, trajectory_dict, state_type=['dead'], final_amount=True):
    
    df_trajectories_states = trajectory_dict['states']
    substates = get_substates(
    model=model, substrings=state_type, include_all=True
    )
    df_substates = df_trajectories_states[substates]
    
    if final_amount is True:
        df_sum = df_substates.iloc[[-1]]
    else:
        df_sum = df_substates
    
    sum_states = df_sum.values.sum()
    
    return sum_states

def run_model_stepwise_vaccines_sum(theta, model, solver, set_start_parameter, set_parameter,
                      observables_names, length_periods, periods, set_initials_zero,
                      state_type=['infectious'], final_amount=False):
    trajectory_dict = run_model_stepwise_vaccines(theta, model, solver, set_start_parameter, set_parameter,
                                observables_names, length_periods, periods, set_initials_zero)
    
    sum_states = get_sum_of_states(model=model, trajectory_dict=trajectory_dict,
                                   state_type=state_type, final_amount=final_amount)
    return sum_states

run_model_stepwise_vaccines_sum_partial = partial(run_model_stepwise_vaccines_sum, model=model, solver=solver,
                                            set_start_parameter=set_start_parameter, set_parameter=set_parameter,
                                            observables_names=observables.keys(), length_periods=length_periods,
                                            periods=periods, set_initials_zero=set_initials_zero)
run_model_stepwise_vaccines_sum_partial(theta=theta)

start = np.repeat(0.2, 2*periods)
bounds = Bounds(np.repeat(0, 2*periods), np.repeat(1, 2*periods))
fit = minimize(fun=run_model_stepwise_vaccines_sum_partial, x0=start,
                        method='L-BFGS-B', bounds=bounds)
fit['jac']