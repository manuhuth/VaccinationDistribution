import numpy as np
import pandas as pd
import pypesto
import pypesto.visualize as visualize
import scipy as sp

from pypesto import optimize

from functools import partial
from estimagic import minimize

from models.vaccination.create_model_vaccination import model_vaccination_create_sbml

from functions.run_sbml import get_model_and_solver_from_sbml
from functions.run_sbml import create_observables_vaccination_rates
from models.vaccination.optimization_functions_estimagic import (
    run_model_stepwise_vaccines_sum_estimagic,
)
from models.vaccination.optimization_functions_pyPesto import (
    run_model_stepwise_vaccines_sum_pyPesto,
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
    distances=np.array([[0, 3], [3, 0]]),
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
    "infectious_countryA_vac0_virW_t0": 10,
    "infectious_countryA_vac0_virM_t0": 10,
    "infectious_countryB_vac0_virW_t0": 10,
    "infectious_countryB_vac0_virM_t0": 10,
    "infectious_countryA_vac1_virW_t0": 1,
    "infectious_countryA_vac1_virM_t0": 1,
    "infectious_countryB_vac1_virW_t0": 1,
    "infectious_countryB_vac1_virM_t0": 1,
    "infectious_countryA_vac2_virW_t0": 1,
    "infectious_countryA_vac2_virM_t0": 1,
    "infectious_countryB_vac2_virW_t0": 1,
    "infectious_countryB_vac2_virM_t0": 1,
}
set_parameter = {
    "beta": 3, #infection rate
    "lambda1": 0.1, #rate of transition from infectious to recovered/dead
    "p": 0.3, #probability of dying p*lambda is fraction that dies in vac0 
    "number_vac1": 100, #inflow of vaccine 1
    "number_vac2": 100, #inflow of vaccine 2
    "omega_vac1_virW": 0.5, #if 1 vaccine1 gives full death protection against wild type if infected
    "delta_vac1_virW": 0.6, #if 1 vaccine1 gives full infection protection against wild type
    "omega_vac2_virW": 0.5, #if 1 vaccine2 gives full death protection against wild type if infected
    "delta_vac2_virW": 0.6, #if 1 vaccine2 gives full infection protection against wild type
    "omega_vac1_virM": 0.5, #if 1 vaccine1 gives full death protection against mutant type if infected
    "delta_vac1_virM": 0.6, #if 1 vaccine1 gives full infection protection against mutant type
    "omega_vac2_virM": 0.5, #if 1 vaccine1 gives full death protection against mutant type if infected
    "delta_vac2_virM": 0.6, #if 1 vaccine2 gives full infection protection against mutant type
    "eta_virW": 1, #factor with which beta is scaled for the wild type
    "eta_virM": 1.3, #factor with which beta is scaled for the mutant type
}
observables_names = observables.keys()

length_periods = 12
periods = 1
observables_names = observables.keys()
set_initials_zero = True


# ------------Optimize Estimagic-----------------------------------------------
quantity = "dead"
run_model_stepwise_vaccines_sum_partial_estimagic = partial(
    run_model_stepwise_vaccines_sum_estimagic,
    model=model,
    solver=solver,
    set_start_parameter=set_start_parameter,
    set_parameter=set_parameter,
    observables_names=observables.keys(),
    length_periods=length_periods,
    periods=periods,
    set_initials_zero=set_initials_zero,
    state_type=[quantity],
    final_amount=True,
)

# plot only works for 2 parameter example (1 period)
plot_3D_function(
    function=run_model_stepwise_vaccines_sum_partial_estimagic,
    number_parameter=2 * periods, set_off_scientific_notation = True,
    decimal_floats = 3, zlabel=quantity,
)


n_multi = 20
cols = (
    [f"start_{i}" for i in range(2 * periods)]
    + [f"theta_{i}" for i in range(2 * periods)]
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
        criterion=run_model_stepwise_vaccines_sum_partial_estimagic,
        params=start_params,
        algorithm="nag_pybobyqa",
    )

    df_multistart.iloc[index_multi] = np.append(
        start, np.append(res["solution_x"], res["solution_criterion"])
    )

df_multistart[df_multistart.value == df_multistart.value.min()]
# -----------------------Plot Model--------------------------------------------
from visualization.model_results import plot_states
from visualization.model_results import plot_observables
from visualization.model_results import get_substates
from visualization.model_results import get_observables_by_name
from functions.run_sbml import run_model

av1 = 0.5
av2 = 0.7
para ={**set_parameter,  'proportion_countryA_vac1' : av1,
 'proportion_countryB_vac1' : 1- av1 ,
 'proportion_countryA_vac2' : av2,
 'proportion_countryB_vac2' : 1-av2}

model_results = run_model(
    model=model,
    solver=solver,
    periods=1,
    length_periods=12,
    set_start_parameter=set_start_parameter,
    set_parameter=para,
    observables_names=observables.keys(),
)
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















# ------------Optimize pyPesto-----------------------------------------------
###pyPesto has problem with partial functions
def run_model_stepwise_vaccines_sum_partial_pesto(theta):
    out = run_model_stepwise_vaccines_sum_pyPesto(
        theta,
        model=model,
        solver=solver,
        set_start_parameter=set_start_parameter,
        set_parameter=set_parameter,
        observables_names=observables.keys(),
        length_periods=length_periods,
        periods=periods,
        set_initials_zero=set_initials_zero,
        state_type=["dead"],
        final_amount=True,
    )
    return out


lower = np.repeat(0, 2 * periods)
upper = np.repeat(1, 2 * periods)

pp_objective = pypesto.Objective(fun=run_model_stepwise_vaccines_sum_partial_pesto)
pp_problem = pypesto.Problem(objective=pp_objective, lb=lower, ub=upper)

optimizer_cmaes = pypesto.optimize.CmaesOptimizer()

optimizer_cobyla = optimize.ScipyOptimizer(method="COBYLA")

n_starts = 20

# save optimizer trace
history_options = pypesto.HistoryOptions(trace_record=True)

# Run optimizaitons for different optimzers
result_cmaes = optimize.minimize(
    problem=pp_problem,
    optimizer=optimizer_cmaes,
    n_starts=n_starts,
    history_options=history_options,
)

result_cobyla = optimize.minimize(
    problem=pp_problem,
    optimizer=optimizer_cobyla,
    n_starts=n_starts,
    history_options=history_options,
)
visualize.waterfall(result_cmaes, size=(15, 6))
