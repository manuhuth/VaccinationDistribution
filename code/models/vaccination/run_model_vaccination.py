import numpy as np

from models.vaccination.create_model_vaccination import model_vaccination_create_sbml

from visualization.model_results import plot_states
from visualization.model_results import plot_observables
from visualization.model_results import get_substates
from visualization.model_results import get_observables_by_name

from functions.run_sbml import get_model_and_solver_from_sbml
from functions.run_sbml import run_model
from functions.run_sbml import create_observables_vaccination_rates

# from functions.vaccine_proportions import (
#    create_rules_vaccination_proportion_relative_population,
# )
# from functions.vaccine_proportions import (
#    create_rules_vaccination_proportion_relative_infected_population,
# )

from functions.vaccine_proportions import create_rules_vaccination_proportion_piecewise
from functions.vaccine_proportions import create_parameters_piecewise


# --------------------------Create Model--------------------------------------
path_sbml = "stored_models/vaccination/vaccination"
vaccination_states = ["vac0", "vac1", "vac2"]
non_vaccination_state = "vac0"
virus_states = ["virW", "virM"]
areas = ["countryA", "countryB"]

species_comp = ["susceptible", "infectious", "recovered", "dead"]

length_decision_period = 3
number_decision_periods = 4
first_vaccination_period = 0


# ----------------------------Create Proportion rules---------------------------
rules_proportions = create_rules_vaccination_proportion_piecewise(
    decision_period_length=length_decision_period,
    number_decision_periods=number_decision_periods,
    first_vaccination_period=first_vaccination_period,
    vaccination_states=vaccination_states,
    non_vaccination_state=non_vaccination_state,
    areas=areas,
)

additional_parameters = create_parameters_piecewise(
    decision_period_length=length_decision_period,
    number_decision_periods=number_decision_periods,
    first_vaccination_period=first_vaccination_period,
    vaccination_states=vaccination_states,
    non_vaccination_state=non_vaccination_state,
    areas=areas,
)


model_vaccination_create_sbml(
    path=path_sbml,
    areas=areas,
    parameter_rules=rules_proportions,
    distances=np.array([[0, 10], [10, 0]]),
    t0_susceptible=0,
    t0_infectious=1000,
    additional_parameters=additional_parameters,
)

observables_nu = create_observables_vaccination_rates(
    vaccination_states_removed=["vac1", "vac2"], areas=areas, name_parameter="nu"
)
observables_proportion = create_observables_vaccination_rates(
    vaccination_states_removed=["vac1", "vac2"],
    areas=areas,
    name_parameter="proportion",
)
observables_time = {"observable_time": {"name": "t", "formula": "t"}}

observables = {**observables_nu, **observables_proportion, **observables_time}

model_and_solver = get_model_and_solver_from_sbml(
    path_sbml=path_sbml,
    model_name="vaccination",
    model_directory="stored_models/vaccination/vaccination_dir",
    observables=observables,
)

# -----------------------Run Model---------------------------------------------
model = model_and_solver["model"]
solver = model_and_solver["solver"]
set_start_parameter = {
    "susceptible_countryA_vac0_t0": 80000,
    "susceptible_countryB_vac0_t0": 80000,
    "infectious_countryA_vac0_virW_t0": 1000,
    "infectious_countryA_vac0_virM_t0": 1000,
    "infectious_countryB_vac0_virW_t0": 1000,
    "infectious_countryB_vac0_virM_t0": 1000,
    "infectious_countryA_vac1_virW_t0": 100,
    "infectious_countryA_vac1_virM_t0": 100,
    "infectious_countryB_vac1_virW_t0": 100,
    "infectious_countryB_vac1_virM_t0": 100,
    "infectious_countryA_vac2_virW_t0": 100,
    "infectious_countryA_vac2_virM_t0": 100,
    "infectious_countryB_vac2_virW_t0": 100,
    "infectious_countryB_vac2_virM_t0": 100,
}
set_fixed_parameter = {
    "beta": 2,
    "lambda1": 0.01,
    "p": 0.3,
    "number_vac1": 200,
    "number_vac2": 200,
    "omega_vac1_virW": 0.5,
    "delta_vac1_virW": 0.6,
    "omega_vac2_virW": 0.5,
    "delta_vac2_virW": 0.6,
    "omega_vac1_virM": 0.5,
    "delta_vac1_virM": 0.6,
    "omega_vac2_virM": 0.5,
    "delta_vac2_virM": 0.6,
    "eta_virW": 1,
    "eta_virM": 1.3,
}
observables_names = list(observables.keys())

set_proportions = {
    "proportion_par_countryA_vac1_0": 0.2,
    "proportion_par_countryA_vac1_3": 0.1,
    "proportion_par_countryA_vac1_6": 0.2,
    "proportion_par_countryA_vac1_9": 0.1,
    "proportion_par_countryA_vac2_0": 0.2,
    "proportion_par_countryA_vac2_3": 0.1,
    "proportion_par_countryA_vac2_6": 0.2,
    "proportion_par_countryA_vac2_9": 0.1,
}

set_parameter = {**set_fixed_parameter, **set_proportions}

model_results = run_model(
    model=model,
    solver=solver,
    periods=1,
    length_periods=12,
    set_start_parameter=set_start_parameter,
    set_parameter=set_parameter,
    observables_names=observables.keys(),
)

#-------------------------optimize--------------------------------------------
from models.vaccination.optimization_functions_estimagic import get_sum_of_states
import pandas as pd
import numpy as np
from estimagic import minimize

def to_optimize(theta):
    model = model_and_solver["model"]
    solver = model_and_solver["solver"]
    set_start_parameter = {
        "susceptible_countryA_vac0_t0": 80000,
        "susceptible_countryB_vac0_t0": 80000,
        "infectious_countryA_vac0_virW_t0": 1000,
        "infectious_countryA_vac0_virM_t0": 1000,
        "infectious_countryB_vac0_virW_t0": 1000,
        "infectious_countryB_vac0_virM_t0": 1000,
        "infectious_countryA_vac1_virW_t0": 100,
        "infectious_countryA_vac1_virM_t0": 100,
        "infectious_countryB_vac1_virW_t0": 100,
        "infectious_countryB_vac1_virM_t0": 100,
        "infectious_countryA_vac2_virW_t0": 100,
        "infectious_countryA_vac2_virM_t0": 100,
        "infectious_countryB_vac2_virW_t0": 100,
        "infectious_countryB_vac2_virM_t0": 100,
    }
    set_fixed_parameter = {
        "beta": 3,
        "lambda1": 0.01,
        "p": 0.3,
        "number_vac1": 800,
        "number_vac2": 800,
        "omega_vac1_virW": 0.5,
        "delta_vac1_virW": 0.6,
        "omega_vac2_virW": 0.5,
        "delta_vac2_virW": 0.6,
        "omega_vac1_virM": 0.5,
        "delta_vac1_virM": 0.6,
        "omega_vac2_virM": 0.5,
        "delta_vac2_virM": 0.6,
        "eta_virW": 1,
        "eta_virM": 1.3,
    }

    
    para = theta["value"]
    set_proportions = {
        "proportion_par_countryA_vac1_0": para[0],
        "proportion_par_countryA_vac1_3": para[1],
        "proportion_par_countryA_vac1_6": para[2],
        "proportion_par_countryA_vac1_9": para[3],
        "proportion_par_countryA_vac2_0": para[4],
        "proportion_par_countryA_vac2_3": para[5],
        "proportion_par_countryA_vac2_6": para[6],
        "proportion_par_countryA_vac2_9": para[7],
    }
    
    set_parameter = {**set_fixed_parameter, **set_proportions}
    
    model_results = run_model(
        model=model,
        solver=solver,
        periods=1,
        length_periods=12,
        set_start_parameter=set_start_parameter,
        set_parameter=set_parameter,
        observables_names=observables.keys(),
    )
    
    trajectory_observables = model_results["observables"]
    trajectory_states = model_results["states"]
    trajectory_dict = {
        "states": trajectory_states,
        "observables": trajectory_observables,
    }
    
    out = {"value" : get_sum_of_states(model, trajectory_dict, state_type=["dead"], final_amount=True)}
    
    return out 

n_multi = 20
periods = 4
cols = (
    [f"start_{i}" for i in range(2*periods)]
    + [f"theta_{i}" for i in range(2*periods)]
    + ["value"]
)
df_multistart = pd.DataFrame(data=None, index=range(n_multi), columns=cols)
for index_multi in range(n_multi):
    print(index_multi)
    start = np.random.uniform(0, 1, 2 * periods)
    start_params = pd.DataFrame(
        data=start,
        columns=["value"],
        index=[f"theta_{i}" for i in range(2 * periods)],
    )
    start_params["lower_bound"] = np.repeat(0, 2 * periods)
    start_params["upper_bound"] = np.repeat(1, 2 * periods)

    res = minimize(
        criterion=to_optimize,
        params=start_params,
        algorithm="nag_pybobyqa",
    )

    df_multistart.iloc[index_multi] = np.append(
        start, np.append(res["solution_x"], res["solution_criterion"])
    )

df_multistart[df_multistart.value == df_multistart.value.min()]

# -----------------------Plot Model--------------------------------------------
trajectory_states = model_results["states"]
trajectory_observables = model_results["observables"]

# visualization
observables_names_nu = get_observables_by_name(
    observables, substrings=["nu"], include_all=True
)
observables_names_proportion = get_observables_by_name(
    observables, substrings=["proportion"], include_all=True
)
# TODO check if values make sense
plot_observables(
    results=model_results, model=model, observable_ids=observables_names_nu
)
plot_observables(
    results=model_results,
    model=model,
    observable_ids=observables_names_proportion,
    set_off_scientific_notation=True,
    decimal_floats=4,
)
# fig, ax = plot_states(results=model_result, model=model)

# dir(model)
# par_values = model.getParameters()
# par_names = model.getParameterNames()
# model.getFixedParameterNames()

substates = get_substates(
    model=model, substrings=["vac0", "countryA"], include_all=True
)
fig, ax = plot_states(results=model_results, model=model, state_ids=substates)
