import numpy as np
import pandas as pd

from estimagic import minimize
from functools import partial

from models.vaccination.optimization_functions_estimagic import get_sum_of_states
from models.vaccination.parameter import general_set_up
from models.vaccination.parameter import fixed_parameter
from models.vaccination.parameter import parameters_vaccine_two
from models.vaccination.parameter import start_parameter_two
from models.vaccination.create_model_vaccination import model_vaccination_create_sbml

from visualization.model_results import plot_states
from visualization.model_results import plot_observables
from visualization.model_results import get_substates
from visualization.model_results import get_observables_by_name

from functions.run_sbml import run_model
from functions.run_sbml import create_observables_vaccination_rates
from functions.run_sbml import get_model_and_solver_from_sbml

from functions.vaccine_proportions import create_splines
from functions.vaccine_proportions import create_spline_parameters
from functions.vaccine_proportions import create_spline_transformation_rules
from functions.vaccine_proportions import create_one_minus_rules_splines

first_time = False
# --------------------------Create Model--------------------------------------
model_name = "vaccination"
path_sbml = "stored_models/vaccination/" + model_name
vaccination_states = ["vac0", "vac1", "vac2"]
non_vaccination_state = general_set_up()["non_vaccination_state"]
virus_states = general_set_up()["virus_states"]
areas = general_set_up()["areas"]
distances = general_set_up()["distances"]
species_comp = general_set_up()["species_compartments"]

length_decision_period = 80
number_decision_periods = 1
first_vaccination_period = 0
number_yy = 20

vaccination_states_removed = [
    x for x in vaccination_states if x != non_vaccination_state
]

areas_without_last = [x for x in areas if x != areas[-1]]

spline_parameters = create_spline_parameters(
    areas=areas, vaccination_states_removed=vaccination_states_removed
)
spline_transformation_rule = create_spline_transformation_rules(
    areas=areas, vaccination_states_removed=vaccination_states_removed
)

minus_one_rules = create_one_minus_rules_splines(
    areas=areas, vaccination_states_removed=vaccination_states_removed, leave_out="last"
)

proportion_rules = {**spline_transformation_rule, **minus_one_rules}

splines = create_splines(
    areas_without=areas_without_last,
    vaccination_states_removed=vaccination_states_removed,
    stop=length_decision_period,
    length=number_yy,
)

model_vaccination_create_sbml(
    path=path_sbml,
    areas=areas,
    parameter_rules=proportion_rules,
    distances=distances,
    additional_parameters=spline_parameters,
    splines=splines,
)

observables = {"observable_time": {"name": "t", "formula": "t"}}

for index in ["nu", "proportion", "spline"]:
    new = create_observables_vaccination_rates(
        vaccination_states_removed=vaccination_states_removed,
        areas=areas,
        name_parameter=index,
    )
    observables = {**observables, **new}

model_directory = "stored_models/" + model_name + "/vaccination_dir"

model_and_solver = get_model_and_solver_from_sbml(
    path_sbml=path_sbml,
    model_name=model_name,
    model_directory=model_directory,
    only_import=not first_time,
)

created_model = {
    "model": model_and_solver["model"],
    "solver": model_and_solver["solver"],
    "observables": observables,
}


# -----------------------Run Model---------------------------------------------
length = length_decision_period * number_decision_periods + first_vaccination_period


model = created_model["model"]
# dir(model)
# par_values = model.getParameters()
# par_names = model.getParameterNames()
# model.getFixedParameterNames()
solver = created_model["solver"]
observables = created_model["observables"]

set_start_parameter = start_parameter_two()
set_fixed_parameter = {**fixed_parameter(), **parameters_vaccine_two()}

# -------------------------optimize--------------------------------------------
def criterion(theta, model, solver, set_fixed_parameter, set_start_parameter, yy_names):

    para = theta["value"]
    para_np = np.array(para)
    transformed_para = np.log(para_np / (1 - para_np))

    set_proportions = dict(zip(yy_names, transformed_para))
    set_fixed_parameter = {**fixed_parameter(), **parameters_vaccine_two()}
    set_parameter = {**set_fixed_parameter, **set_proportions}

    model_results = run_model(
        model=model,
        solver=solver,
        periods=1,
        length_periods=length,
        set_start_parameter=set_start_parameter,
        set_parameter=set_parameter,
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


yy_names_str = []
for index_areas in areas_without_last:
    for index_vaccines in vaccination_states_removed:
        for index in range(number_yy):
            yy_names_str.append(f"yy_{index_areas}_{index_vaccines}_{index}")


to_optimize = partial(
    criterion,
    model=model,
    solver=solver,
    set_fixed_parameter=set_fixed_parameter,
    set_start_parameter=start_parameter_two(),
    yy_names=yy_names_str,
)

n_multi = 2
number_vaccines = len(vaccination_states) - 1

number_parameters = len(yy_names_str)
lower_bound = 0.000001
upper_bound = 0.999999


cols = (
    ["start_" + yy_parameter for yy_parameter in yy_names_str]
    + ["estimated_" + yy_parameter for yy_parameter in yy_names_str]
    + ["value"]
)
np.random.seed(12345)
df_multistart = pd.DataFrame(data=None, index=range(n_multi), columns=cols)
for index_multi in range(n_multi):
    print(index_multi)
    start = np.random.uniform(lower_bound, upper_bound, number_parameters)
    start_params = pd.DataFrame(
        data=start,
        columns=["value"],
        index=yy_names_str,
    )
    start_params["lower_bound"] = np.repeat(lower_bound, number_parameters)
    start_params["upper_bound"] = np.repeat(upper_bound, number_parameters)

    res = minimize(
        criterion=to_optimize,
        params=start_params,
        algorithm="nag_pybobyqa",
    )
    print(res["success"])

    df_multistart.iloc[index_multi] = np.append(
        start, np.append(res["solution_x"], res["solution_criterion"])
    )

optimal_column = df_multistart[df_multistart.value == df_multistart.value.min()]

# ------------------------Run Optimum-----------------------------------------
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

# -----------------------Plot Optimum--------------------------------------------
trajectory_states_optimal = model_results_optimal["states"]
trajectory_observables_optimal = model_results_optimal["observables"]

observables_names_nu_optimal = get_observables_by_name(
    observables, substrings=["nu"], include_all=True
)

observables_names_proportion_optimal = get_observables_by_name(
    observables, substrings=["proportion", "countryA"], include_all=True
)

plot_observables(
    results=model_results_optimal,
    model=model,
    observable_ids=observables_names_nu_optimal,
    time_name="t",
)

plot_observables(
    results=model_results_optimal,
    model=model,
    observable_ids=["pline_countryA_vac1", "pline_countryA_vac2"],
    set_off_scientific_notation=True,
    decimal_floats=4,
    time_name="t",
)

plot_observables(
    results=model_results_optimal,
    model=model,
    observable_ids=observables_names_proportion_optimal,
    set_off_scientific_notation=True,
    decimal_floats=4,
    time_name="t",
)

substates_optimal = get_substates(
    model=model, substrings=["vac0", "countryA"], include_all=True
)

fig, ax = plot_states(
    results=model_results_optimal,
    model=model,
    state_ids=substates_optimal,
    time_name="t",
)
