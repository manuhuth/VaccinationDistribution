import numpy as np
import pandas as pd
import pickle

import pypesto
import pypesto.optimize as optimize
import pypesto.visualize as visualize

from models.vaccination.optimization_functions_estimagic import get_sum_of_states
from models.vaccination.parameter import general_set_up
from models.vaccination.parameter import fixed_parameter
from models.vaccination.parameter import parameters_vaccine_two
from models.vaccination.parameter import start_parameter_two
from models.vaccination.create_model_vaccination import model_vaccination_create_sbml
from models.vaccination.create_model_vaccination import b_distance_function

from functions.run_sbml import run_model
from functions.run_sbml import create_observables_vaccination_rates
from functions.run_sbml import get_model_and_solver_from_sbml

from functions.vaccine_proportions import create_splines
from functions.vaccine_proportions import create_spline_parameters
from functions.vaccine_proportions import create_spline_transformation_rules
from functions.vaccine_proportions import create_one_minus_rules_splines
from functions.vaccine_proportions import create_parameters_piecewise_vaccine_supply
from functions.vaccine_proportions import create_vaccine_supply_rules


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

length_decision_period = 140
number_decision_periods = 1
first_vaccination_period = 0
number_yy = 11
n_intervals = 2000

xx_list = ["xx" + str(i) for i in range(number_yy)]
parameter_vacc_supply = create_parameters_piecewise_vaccine_supply(
    xx=xx_list,
    vaccination_states=vaccination_states,
    non_vaccination_state=non_vaccination_state,
)

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

vaccine_supply_rules = create_vaccine_supply_rules(
    vaccination_states, non_vaccination_state, xx_list
)

model_vaccination_create_sbml(
    path=path_sbml,
    areas=areas,
    parameter_rules={**proportion_rules, **vaccine_supply_rules},
    distances=distances,
    additional_parameters={**spline_parameters, **parameter_vacc_supply},
    splines=splines,
)

observables_vaccinated_vacc = {}
for index_area in areas:
    for index_vacc in vaccination_states_removed:
        name = "quantity_" + index_area + "_" + index_vacc
        proportion = f"proportion_{index_area}_{index_vacc}"
        number_vacc = f"number_{index_vacc}"
        formula = f"{proportion} * {number_vacc}"
        observables_vaccinated_vacc[name] = {"name": name, "formula": formula}

observables = {
    "observable_time": {"name": "t", "formula": "t"},
    **observables_vaccinated_vacc,
}

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
vaccine_supply_parameter_strings = [
    x for x in model.getParameterNames() if "vaccine_supply" in x
]
vaccine_supply_parameter_values = np.repeat(
    float(40000), len(vaccine_supply_parameter_strings)
)
set_vaccine_supply_parameter = dict(
    zip(vaccine_supply_parameter_strings, vaccine_supply_parameter_values)
)
set_fixed_parameter = {
    **fixed_parameter(),
    **parameters_vaccine_two(),
    **distance_parameters,
    **set_vaccine_supply_parameter,
}
# ------------------------Run Current case-----------------------------------------
yy_names_str = []
for index_areas in areas_without_last:
    for index_vaccines in vaccination_states_removed:
        for index in range(number_yy):
            yy_names_str.append(f"yy_{index_areas}_{index_vaccines}_{index}")

yy_current_splines = np.repeat(0.0, len(yy_names_str))
yy_current = dict(zip(yy_names_str, yy_current_splines))
set_parameter_current = {**set_fixed_parameter, **yy_current}

model_results_current = run_model(
    model=model,
    solver=solver,
    periods=1,
    length_periods=length,
    set_start_parameter=set_start_parameter,
    set_parameter=set_parameter_current,
    number_intervals=n_intervals,
)

df_current = model_results_current
states_current = df_current["states"]
trajectory_current = {
    "states": states_current,
}

df_current = model_results_current
states_current = df_current["states"]
trajectory_current = {
    "states": states_current,
}
var_interest = "dead"
current_deceased_A = get_sum_of_states(
    model, trajectory_current, state_type=[var_interest, "countryA"], final_amount=True
)
current_deceased_B = get_sum_of_states(
    model, trajectory_current, state_type=[var_interest, "countryB"], final_amount=True
)
# -------------------------find pareto optimal solution------------------------
def criterion(
    theta,
    model,
    solver,
    set_fixed_parameter,
    set_start_parameter,
    yy_names,
    current_deceased_A,
    current_deceased_B,
):

    para = theta
    para_np = np.array(para)
    transformed_para = np.log(para_np / (1 - para_np))

    set_proportions = dict(zip(yy_names, transformed_para))
    set_fixed_parameter = {**set_fixed_parameter, **parameters_vaccine_two()}
    set_parameter = {**set_fixed_parameter, **set_proportions}

    model_results = run_model(
        model=model,
        solver=solver,
        periods=1,
        length_periods=length,
        set_start_parameter=set_start_parameter,
        set_parameter=set_parameter,
        number_intervals=n_intervals,
    )

    trajectory_observables = model_results["observables"]
    trajectory_states = model_results["states"]
    trajectory_dict = {
        "states": trajectory_states,
        "observables": trajectory_observables,
    }

    deceased_A = get_sum_of_states(
        model, trajectory_dict, state_type=["dead", "countryA"], final_amount=True
    )

    deceased_B = get_sum_of_states(
        model, trajectory_dict, state_type=["dead", "countryB"], final_amount=True
    )

    return deceased_A, deceased_B


large_value = 10 ** 15


def func_to_optimize_max(theta):
    deceased_A, deceased_B = criterion(
        theta=theta,
        model=model,
        solver=solver,
        set_fixed_parameter=set_fixed_parameter,
        set_start_parameter=start_parameter_two(),
        yy_names=yy_names_str,
        current_deceased_A=current_deceased_A,
        current_deceased_B=current_deceased_B,
    )

    percentage_A = deceased_A / current_deceased_A
    percentage_B = deceased_B / current_deceased_B

    if (percentage_A > 1) or (percentage_B > 1):
        value = large_value
    else:
        value = deceased_A + deceased_B

    # value = np.max([percentage_A, percentage_B])

    return value


np.random.seed(12345)
lower_bound = 0.00000001
upper_bound = 0.99999999
n_multi = 100
start_matrix = np.array(np.repeat(np.nan, len(areas) * number_yy))
success_values = []
k = 0
n_rows = 0
while n_rows < n_multi and k < 10000:
    candidate = np.random.uniform(lower_bound, upper_bound, len(areas) * number_yy)
    cand_value = func_to_optimize_max(candidate)

    if cand_value < large_value:
        start_matrix = np.vstack((start_matrix, candidate))
        success_values.append(cand_value)
        print("success")
    print(k)
    k += 1
    n_rows = start_matrix.shape[0] - 1

start_matrix = np.delete(start_matrix, (0), axis=0)

objective = pypesto.Objective(fun=func_to_optimize_max)
number_vaccines = len(vaccination_states) - 1
number_parameters = len(yy_names_str)
lb = np.repeat(lower_bound, number_parameters)
ub = np.repeat(upper_bound, number_parameters)

problem = pypesto.Problem(objective=objective, lb=lb, ub=ub, x_guesses=start_matrix)

# optimizer = optimize.FidesOptimizer()
optimizer = optimize.ScipyOptimizer("L-BFGS-B")
history_options = pypesto.HistoryOptions(trace_record=True)

results_max = optimize.minimize(
    problem=problem,
    optimizer=optimizer,
    n_starts=n_multi,
    history_options=history_options,
)

visualize.waterfall(results_max)

df_results_max = results_max.optimize_result.as_dataframe()
theta_optimal_max = df_results_max.iloc[0]["x"]

# ------------------------Run Optimum-----------------------------------------
transformed_theta_optimal = np.log(theta_optimal_max / (1 - theta_optimal_max))
yy_optimal = dict(zip(yy_names_str, transformed_theta_optimal))

set_parameter = {**set_fixed_parameter, **yy_optimal}

model_results_optimal = run_model(
    model=model,
    solver=solver,
    periods=1,
    length_periods=length,
    set_start_parameter=set_start_parameter,
    set_parameter=set_parameter,
)
trajectory_states_optimal = model_results_optimal["states"]
trajectory_optimal = {
    "states": trajectory_states_optimal,
}

# -----------run seperated-----------------------------------------------------
yy_seperated_splines_wild = np.repeat(20.0, len(yy_names_str) / 2)
yy_seperated_splines_mutant = np.repeat(-20.0, len(yy_names_str) / 2)
yy_seperated_splines = np.concatenate(
    (yy_seperated_splines_wild, yy_seperated_splines_mutant)
)
yy_current = dict(zip(yy_names_str, yy_seperated_splines))
set_parameter_seperated = {**set_fixed_parameter, **yy_current}
model_results_seperated = run_model(
    model=model,
    solver=solver,
    periods=1,
    length_periods=length,
    set_start_parameter=set_start_parameter,
    set_parameter=set_parameter_seperated,
)
df_seperated = model_results_seperated
states_seperated = df_seperated["states"]
trajectory_seperated = {
    "states": states_seperated,
}


# ------------Compare specifications-------------------------------------------
seperated_deceased_A = get_sum_of_states(
    model,
    trajectory_seperated,
    state_type=[var_interest, "countryA"],
    final_amount=True,
)
seperated_deceased_B = get_sum_of_states(
    model,
    trajectory_seperated,
    state_type=[var_interest, "countryB"],
    final_amount=True,
)
seperated_deceased = seperated_deceased_A + seperated_deceased_B


current_deceased_A = get_sum_of_states(
    model, trajectory_current, state_type=[var_interest, "countryA"], final_amount=True
)
current_deceased_B = get_sum_of_states(
    model, trajectory_current, state_type=[var_interest, "countryB"], final_amount=True
)
current_deceased = current_deceased_A + current_deceased_B

dead_A_current = get_sum_of_states(
    model, trajectory_current, state_type=[var_interest, "countryA"], final_amount=True
)
dead_B_current = get_sum_of_states(
    model, trajectory_current, state_type=[var_interest, "countryB"], final_amount=True
)

dead_A_optimal = get_sum_of_states(
    model, trajectory_optimal, state_type=[var_interest, "countryA"], final_amount=True
)
dead_B_optimal = get_sum_of_states(
    model, trajectory_optimal, state_type=[var_interest, "countryB"], final_amount=True
)

optimal_deceased = get_sum_of_states(
    model, trajectory_optimal, state_type=[var_interest], final_amount=True
)

print_compare_A = f"A: Current strategy: {dead_A_current}; Optimal strategy: {dead_A_optimal}; Seperated strategy: {seperated_deceased_A}"
print_compare_B = f"B: Current strategy: {dead_B_current}; Optimal strategy: {dead_B_optimal}; Seperated strategy: {seperated_deceased_B}"
print_compare_total = f"Total: Current strategy: {current_deceased}; Optimal strategy: {optimal_deceased}; Seperated strategy: {seperated_deceased}"
print(print_compare_total)
print(print_compare_A)
print(print_compare_B)

vac1_A_optimal = get_sum_of_states(
    model, trajectory_optimal, state_type=["vac1", "countryA"], final_amount=True
)
vac2_A_optimal = get_sum_of_states(
    model, trajectory_optimal, state_type=["vac2", "countryA"], final_amount=True
)
vac1_B_optimal = get_sum_of_states(
    model, trajectory_optimal, state_type=["vac1", "countryB"], final_amount=True
)
vac2_B_optimal = get_sum_of_states(
    model, trajectory_optimal, state_type=["vac2", "countryB"], final_amount=True
)


dict_out = {
    "df_optimal_results": df_results_max,
    "model_optimal": model_results_optimal,
    "model_seperated": model_results_seperated,
    "model_current": model_results_current,
    "compare_total": print_compare_total,
    "compare_A": print_compare_A,
    "compare_B": print_compare_B,
    "model_directory": model_directory,
    "model_name": model_name,
    "path_sbml": path_sbml,
}


with open(
    "/home/manuel/Documents/VaccinationDistribution/code/objects/output_splines.pkl",
    "wb",
) as output:
    out = dict_out
    pickle.dump(out, output, pickle.HIGHEST_PROTOCOL)
