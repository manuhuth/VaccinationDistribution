import numpy as np
import pandas as pd

from functools import partial
from estimagic import minimize

from models.vaccination.create_model_vaccination import model_vaccination_create_sbml

from functions.run_sbml import get_model_and_solver_from_sbml
from functions.run_sbml import create_observables_vaccination_rates
from models.vaccination.optimization_functions_estimagic import (
    run_model_stepwise_vaccines_sum,
)
from visualization.model_results import plot_3D_function


# --------------------------Create Model--------------------------------------
path_sbml = "stored_models/vaccination/step_function_optim"
vaccination_states = ["vac0", "vac1"]
non_vaccination_state = "vac0"
virus_states = ["virW", "virM"]
areas = ["countryA", "countryB"]
species_comp = ["susceptible", "infectious", "recovered", "dead"]
vaccination_states_removed = [
    x for x in vaccination_states if x != non_vaccination_state
]

model_vaccination_create_sbml(
    path=path_sbml,
    areas=areas,
    distances=np.array([[0, 10], [10, 0]]),
)

observables_nu = create_observables_vaccination_rates(
    vaccination_states_removed=vaccination_states_removed,
    areas=areas,
    name_parameter="nu",
)
observables_proportion = create_observables_vaccination_rates(
    vaccination_states_removed=vaccination_states_removed,
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
set_parameter = {
    "beta": 2,
    "lambda1": 0.01,
    "p": 0.3,
    "number_vac1": 100,
    "number_vac2": 100,
    "omega_vac1_virW": 0.9,
    "delta_vac1_virW": 0.9,
    "omega_vac2_virW": 0.6,
    "delta_vac2_virW": 0.6,
    "omega_vac1_virM": 0.6,
    "delta_vac1_virM": 0.6,
    "omega_vac2_virM": 0.9,
    "delta_vac2_virM": 0.9,
    "eta_virW": 1,
    "eta_virM": 1.3,
}
observables_names = observables.keys()

length_periods = 2
periods = 6
observables_names = observables.keys()
set_initials_zero = True

run_model_stepwise_vaccines_sum_partial = partial(
    run_model_stepwise_vaccines_sum,
    model=model,
    solver=solver,
    set_start_parameter=set_start_parameter,
    set_parameter=set_parameter,
    observables_names=observables.keys(),
    length_periods=length_periods,
    periods=periods,
    set_initials_zero=set_initials_zero,
    state_type=["infectious"],
    final_amount=False,
)

# plot only works for 2 parameter example (1 period)
# plot_3D_function(
#    function=run_model_stepwise_vaccines_sum_partial, number_parameter=2 * periods
# )

# ------------Optimize Estimagic-----------------------------------------------------

start = np.repeat(0.3, 2 * periods)
start_params = pd.DataFrame(
    data=start,
    columns=["value"],
    index=[f"theta_{i}" for i in range(2 * periods)],
)
start_params["lower_bound"] = np.repeat(0, 2 * periods)
start_params["upper_bound"] = np.repeat(1, 2 * periods)

res = minimize(
    criterion=run_model_stepwise_vaccines_sum_partial,
    params=start_params,
    algorithm="nag_pybobyqa",
)
