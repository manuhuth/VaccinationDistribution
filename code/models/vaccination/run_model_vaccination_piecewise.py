import numpy as np
import pandas as pd
from functools import partial

from functions.run_sbml import run_model
from models.vaccination.parameter import general_set_up
from models.vaccination.parameter import fixed_parameter
from models.vaccination.parameter import parameters_vaccine_two
from models.vaccination.parameter import start_parameter_two
from models.vaccination.create_model_vaccination import model_vaccination_create_sbml
from models.vaccination.create_model_vaccination import b_distance_function
from models.vaccination.optimization_functions_estimagic import get_sum_of_states

from visualization.model_results import plot_states
from visualization.model_results import plot_observables
from visualization.model_results import get_substates
from visualization.model_results import get_observables_by_name

from functions.vaccine_proportions import create_rules_vaccination_proportion_piecewise
from functions.vaccine_proportions import create_parameters_piecewise
from functions.run_sbml import create_observables_vaccination_rates
from functions.run_sbml import get_model_and_solver_from_sbml

from estimagic import minimize

first_time = True
# --------------------------Create Model--------------------------------------
model_name = "vaccination_piecewise"
path_sbml = "stored_models/vaccination_piecewise/" + model_name
vaccination_states = ["vac0", "vac1", "vac2"]
non_vaccination_state = general_set_up()["non_vaccination_state"]
virus_states = general_set_up()["virus_states"]
areas = general_set_up()["areas"]
distances = general_set_up()["distances"]
species_comp = general_set_up()["species_compartments"]

length_decision_period = 14
number_decision_periods = 11
first_vaccination_period = 0

vaccination_states_removed = [
    x for x in vaccination_states if x != non_vaccination_state
]
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
    distances=distances,
    additional_parameters=additional_parameters,
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
observables_time = {"observable_time": {"name": "t", "formula": "t"}}

observables = {**observables_nu, **observables_proportion, **observables_time}

model_directory = "stored_models/" + model_name + "/vaccination_dir"

model_and_solver = get_model_and_solver_from_sbml(
    path_sbml=path_sbml,
    model_name=model_name,
    model_directory=model_directory,
    observables=observables,
    only_import=not first_time,
)

created_model = {
    "model": model_and_solver["model"],
    "solver": model_and_solver["solver"],
    "observables": observables,
}
# -----------------------Run Model---------------------------------------------
length = length_decision_period * number_decision_periods
model = created_model["model"]
solver = created_model["solver"]
observables = created_model["observables"]

distance_parameters = {}
distances_df = pd.DataFrame(distances, columns=areas, index=areas)
b_distance_matrix = b_distance_function(distances_df)
for index_areas_row in areas:
    for index_areas_col in areas:
        name = f"distance_{index_areas_row}_{index_areas_col}"
        distance_parameters[name] = b_distance_matrix.loc[
            index_areas_row, index_areas_col
        ]

set_start_parameter = start_parameter_two()
set_fixed_parameter = {
    **fixed_parameter(),
    **parameters_vaccine_two(),
    **distance_parameters,
}
# -------------------------optimize--------------------------------------------


xx = np.linspace(0, 140, 11)


def criterion(
    theta, model, solver, set_fixed_parameter, xx, vaccination_states_removed
):

    para = theta["value"]

    set_proportions = {}
    number_index = 0
    for index_vac in vaccination_states_removed:
        for index_xx in xx:
            name = f"proportion_par_countryA_{index_vac}_{index_xx}"
            set_proportions[name] = para[number_index]
            number_index += 1

    set_parameter = {**set_fixed_parameter, **set_proportions}

    model_results = run_model(
        model=model,
        solver=solver,
        periods=1,
        length_periods=20,
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

    out = {
        "value": get_sum_of_states(
            model, trajectory_dict, state_type=["dead"], final_amount=True
        )
    }

    return out


to_optimize = partial(
    criterion,
    model=model,
    solver=solver,
    set_fixed_parameter=set_fixed_parameter,
    xx=xx,
    vaccination_states_removed=vaccination_states_removed,
)

n_multi = 20
number_vaccines = len(vaccination_states) - 1
periods = number_decision_periods
cols = (
    [f"start_{i}" for i in range(number_vaccines * periods)]
    + [f"theta_{i}" for i in range(number_vaccines * periods)]
    + ["value"]
)
np.random.seed(12345)
df_multistart = pd.DataFrame(data=None, index=range(n_multi), columns=cols)
for index_multi in range(n_multi):
    print(index_multi)
    start = np.random.uniform(0, 1, number_vaccines * periods)
    start_params = pd.DataFrame(
        data=start,
        columns=["value"],
        index=[f"theta_{i}" for i in range(number_vaccines * periods)],
    )
    start_params["lower_bound"] = np.repeat(0, number_vaccines * periods)
    start_params["upper_bound"] = np.repeat(1, number_vaccines * periods)

    res = minimize(
        criterion=to_optimize,
        params=start_params,
        algorithm="nag_pybobyqa",
    )

    df_multistart.iloc[index_multi] = np.append(
        start, np.append(res["solution_x"], res["solution_criterion"])
    )

optimal_column = df_multistart[df_multistart.value == df_multistart.value.min()]

# ------------------------Run Optimum-----------------------------------------
yy_names_str = []
for index_vac in vaccination_states_removed:
    for index_xx in xx:
        yy_names_str.append(f"proportion_par_countryA_{index_vac}_{index_xx}")

theta_names = [col for col in optimal_column if col.startswith("theta")]
theta_optimal = list(optimal_column.loc[optimal_column.index[0], theta_names])
yy_optimal = dict(zip(yy_names_str, theta_optimal))

set_parameter = {**set_fixed_parameter, **yy_optimal}

model_results_optimal = run_model(
    model=model,
    solver=solver,
    periods=1,
    length_periods=length,
    set_start_parameter=set_start_parameter,
    set_parameter=set_parameter,
)
