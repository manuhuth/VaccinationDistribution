import numpy as np
import pandas as pd
import pickle

from functions.run_sbml import run_model
from models.vaccination.parameter import general_set_up
from models.vaccination.parameter import fixed_parameter
from models.vaccination.parameter import parameters_vaccine_two
from models.vaccination.parameter import start_parameter_two
from models.vaccination.create_model_vaccination import model_vaccination_create_sbml
from models.vaccination.create_model_vaccination import b_distance_function
from models.vaccination.optimization_functions_estimagic import get_sum_of_states

from functions.vaccine_proportions import create_rules_vaccination_proportion_piecewise
from functions.vaccine_proportions import create_parameters_piecewise
from functions.run_sbml import create_observables_vaccination_rates
from functions.run_sbml import get_model_and_solver_from_sbml
from functions.vaccine_proportions import create_vaccine_supply_rules
from functions.create_vaccine_doses_inflow import create_inflow_from_data
from functions.vaccine_proportions import create_parameters_piecewise_vaccine_supply
import pypesto
import pypesto.optimize as optimize
import pypesto.visualize as visualize


first_time = False
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
number_decision_periods = 10
first_vaccination_period = 0
n_intervals = 6000

xx_list = ["xx" + str(i) for i in (range(number_decision_periods))]

xx_params = {}
period = 0
for index in xx_list:
    xx_params[index] = {"value" : period, "constant" : True }
    period += 14
    
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

vaccine_supply_rules = create_vaccine_supply_rules(
    vaccination_states, non_vaccination_state, xx_list
)


parameter_vacc_supply = create_parameters_piecewise_vaccine_supply(
    xx=xx_list,
    vaccination_states=vaccination_states,
    non_vaccination_state=non_vaccination_state,
)

additional_parameters_pw = create_parameters_piecewise(
    decision_period_length=length_decision_period,
    number_decision_periods=number_decision_periods,
    first_vaccination_period=first_vaccination_period,
    vaccination_states=vaccination_states,
    non_vaccination_state=non_vaccination_state,
    areas=areas,
)



additional_parameters = {**additional_parameters_pw, **xx_params, **parameter_vacc_supply}

model_vaccination_create_sbml(
    path=path_sbml,
    areas=areas,
    parameter_rules={**rules_proportions, **vaccine_supply_rules},
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

observables_vaccinated_vacc = {}
for index_area in areas:
    for index_vacc in vaccination_states_removed:
        name = "observable_quantity_" + index_area + "_" + index_vacc
        proportion = f"proportion_{index_area}_{index_vacc}"
        number_vacc = f"number_{index_vacc}"
        formula = f"{proportion} * {number_vacc}"
        observables_vaccinated_vacc[name] = {"name": name, "formula": formula}


observables = {**observables_nu, **observables_proportion, **observables_vaccinated_vacc }

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

vaccine_supply_parameter_values = create_inflow_from_data(number_decision_periods=11)
vaccine_supply_parameter_strings = [
    x for x in model.getParameterNames() if "vaccine_supply" in x
]
set_vaccine_supply_parameter = dict(
    zip(vaccine_supply_parameter_strings, vaccine_supply_parameter_values)
)

set_start_parameter = start_parameter_two()
set_fixed_parameter = {
    **fixed_parameter(),
    **parameters_vaccine_two(),
    **distance_parameters,
    **set_vaccine_supply_parameter,
}
#------------------------------run current------------------------------------
areas_without_last = [x for x in areas if x != areas[-1]]
yy_names_str = []
for index_vaccines in vaccination_states_removed:
    for index in range(number_decision_periods):
         yy_names_str.append(f"proportion_par_countryA_{index_vaccines}_{index*14}")

yy_current_par = np.repeat(0.5, len(yy_names_str))
yy_current = dict(zip(yy_names_str, yy_current_par))
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
var_interest = "dead"
current_deceased_A = get_sum_of_states(
    model, trajectory_current, state_type=[var_interest, "countryA"], final_amount=True
)
current_deceased_B = get_sum_of_states(
    model, trajectory_current, state_type=[var_interest, "countryB"], final_amount=True
)


# -------------------------optimize--------------------------------------------
xx = np.linspace(0, length_decision_period, number_decision_periods)
large_value = 10 ** 15

def criterion(
    theta, model, solver, set_fixed_parameter, set_start_parameter, yy_names, current_deceased_A,
    current_deceased_B,
):

    
    set_proportions = dict(zip(yy_names, theta))

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
        #value = deceased_A + deceased_B
    else:
        value = deceased_A + deceased_B

    # value = np.max([percentage_A, percentage_B])

    return float(value)

def func_to_optimize_unrestricted(theta):
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

    value = deceased_A + deceased_B

    return float(value)

np.random.seed(12345)
lower_bound = 0.00000001
upper_bound = 0.99999999
n_multi = 50
start_matrix = np.array(np.repeat(np.nan, len(areas) * number_decision_periods))
success_values = []
k = 0
n_rows = 0
while n_rows < n_multi and k < 10000:
    candidate = np.random.uniform(lower_bound, upper_bound, len(areas) * number_decision_periods)
    cand_value = func_to_optimize_max(candidate)

    if cand_value < large_value:
        start_matrix = np.vstack((start_matrix, candidate))
        success_values.append(cand_value)
        n_rows = start_matrix.shape[0] - 1
        print("success")
    print(k)
    k += 1
    

start_matrix = np.delete(start_matrix, (0), axis=0)
from pypesto import Objective, FD
objective = pypesto.Objective(fun=func_to_optimize_max)

def obj_func(theta):
    
    grad = pypesto.FD(objective, method="central").get_grad(x=np.repeat(0.4, len(yy_names_str)))
    hess = pypesto.FD(objective, method="forward").get_hess(x=np.repeat(0.4, len(yy_names_str)))
    
    return [func_to_optimize_max(theta=theta), grad, hess]
    #return np.array([func_to_optimize_max(theta=theta)])


objective2 = pypesto.Objective(fun=func_to_optimize_max, grad=True, hess=True)
objective3 = pypesto.Objective(fun=obj_func, grad=True, hess=True)

number_vaccines = len(vaccination_states) - 1

lb = np.repeat(lower_bound, 2*number_decision_periods)
ub = np.repeat(upper_bound, 2*number_decision_periods)


history_options = pypesto.HistoryOptions(trace_record=True)


#problem = pypesto.Problem(objective=objective3, lb=lb, ub=ub, x_guesses=start_matrix)
#optimizer = optimize.FidesOptimizer(options={"fatol" : 1})
problem = pypesto.Problem(objective=objective, lb=lb, ub=ub, x_guesses=start_matrix)
optimizer = optimize.ScipyOptimizer("L-BFGS-B")

results_max = optimize.minimize(
    problem=problem,
    optimizer=optimizer,
    n_starts=n_multi,
    history_options=history_options,
)

visualize.waterfall(results_max)

df_results_max = results_max.optimize_result.as_dataframe()
theta_optimal_max = df_results_max.iloc[0]["x"]

np.random.seed(123456)
objective_unrestricted = pypesto.Objective(fun=func_to_optimize_unrestricted)
problem_unrestricted = pypesto.Problem(objective=objective_unrestricted, lb=lb, ub=ub)

results_unrestricted = optimize.minimize(
    problem=problem_unrestricted,
    optimizer=optimizer,
    n_starts=n_multi,
    history_options=history_options,
    #engine=engine,
)
df_results_unrestricted = results_unrestricted.optimize_result.as_dataframe()
theta_optimal_unrestricted = df_results_unrestricted.iloc[0]["x"]

# ------------------------Run Optimum-----------------------------------------
yy_optimal = dict(zip(yy_names_str, theta_optimal_max))

set_parameter = {**set_fixed_parameter, **yy_optimal}

model_results_optimal = run_model(
    model=model,
    solver=solver,
    periods=1,
    length_periods=length,
    set_start_parameter=set_start_parameter,
    set_parameter=set_parameter,
    number_intervals=n_intervals,
)

trajectory_states_optimal = model_results_optimal["states"]
trajectory_optimal = {
    "states": trajectory_states_optimal,
}

# -----------run seperated-----------------------------------------------------
yy_unrestricted = dict(zip(yy_names_str, theta_optimal_unrestricted))
set_parameter_seperated = {**set_fixed_parameter, **yy_unrestricted}
model_results_unrestricted = run_model(
    model=model,
    solver=solver,
    periods=1,
    length_periods=length,
    set_parameter=set_parameter_seperated,
    set_start_parameter=set_start_parameter,
    number_intervals=n_intervals,
)

df_unrestricted= model_results_unrestricted
states_unrestricted = df_unrestricted["states"]
trajectory_unrestricted= {
    "states": states_unrestricted,
}


# ------------Compare specifications-------------------------------------------
seperated_deceased_A = get_sum_of_states(
    model,
    trajectory_unrestricted,
    state_type=[var_interest, "countryA"],
    final_amount=True,
)
seperated_deceased_B = get_sum_of_states(
    model,
    trajectory_unrestricted,
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
    "df_unrestricted_results" : df_results_unrestricted,
    "model_optimal": model_results_optimal,
    "model_seperated": model_results_unrestricted,
    "model_current": model_results_current,
    "compare_total": print_compare_total,
    "compare_A": print_compare_A,
    "compare_B": print_compare_B,
    "model_directory": model_directory,
    "model_name": model_name,
    "path_sbml": path_sbml,
    "time_points" : model.getTimepoints(),
    "type" : "piecewise",
}

model_out = {"state_names" : model.getStateNames(),
             "initial_states" : model.getInitialStates(),
             "par_magnitude" : model.getParameters(),
             "par_names" : model.getParameterNames()
             }


with open(
    "/home/manuel/Documents/VaccinationDistribution/code/objects/output_piecewise.pkl",
    "wb",
) as output:
    out = dict_out
    pickle.dump(out, output, pickle.HIGHEST_PROTOCOL)
    
with open(
    "/home/manuel/Documents/VaccinationDistribution/code/objects/output_model_props_piecewise.pkl",
    "wb",
) as output:
    out = model_out
    pickle.dump(out, output, pickle.HIGHEST_PROTOCOL)