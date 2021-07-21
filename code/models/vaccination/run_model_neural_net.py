import numpy as np
import pandas as pd
import pickle

from functions.run_sbml import run_model
from models.vaccination.parameter import general_set_up
from models.vaccination.parameter import fixed_parameter
from models.vaccination.parameter import parameters_vaccine_two
from models.vaccination.parameter import start_parameter_two
from models.vaccination.create_model_vaccination import create_species_model
from models.vaccination.create_model_vaccination import create_parameters_model
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


first_time = True
# --------------------------Create Model--------------------------------------
model_name = "vaccination_neural_network"
path_sbml = "stored_models/vaccination_neural_network/" + model_name
vaccination_states = ["vac0", "vac1", "vac2"]
non_vaccination_state = general_set_up()["non_vaccination_state"]
virus_states = general_set_up()["virus_states"]
areas = general_set_up()["areas"]
distances = general_set_up()["distances"]
species_comp = general_set_up()["species_compartments"]

n_intervals = 6000


species = create_species_model(
    vaccination_states,
    non_vaccination_state,
    virus_states,
    areas,
    species_comp,
)

parameters = {**parameters_vaccine_two(), **start_parameter_two(), **fixed_parameter()}


number_vacc_supply_decision_periods = 1
xx_list = ["xx" + str(i) for i in (range(number_vacc_supply_decision_periods - 1))]

xx_params = {}
period = 0
for index in xx_list:
    xx_params[index] = {"value": period, "constant": True}
    period += 14

vaccination_states_removed = [
    x for x in vaccination_states if x != non_vaccination_state
]
vaccine_supply_rules = create_vaccine_supply_rules(
    vaccination_states, non_vaccination_state, xx_list
)

parameter_vacc_supply = create_parameters_piecewise_vaccine_supply(
    xx=xx_list,
    vaccination_states=vaccination_states,
    non_vaccination_state=non_vaccination_state,
)


vaccination_parameters = ["proportion_countryA_vac1", "proportion_countryA_vac2"]


# inputs = list({**species, **parameters}.keys())
inputs = choose_paras = [
    "infectious_countryA_vac0_virW",
    "infectious_countryA_vac0_virM",
    "infectious_countryA_vac0_virW",
    "infectious_countryA_vac0_virM",
    "susceptible_countryA_vac0",
    "susceptible_countryB_vac0",
    "beta",
    "eta_virM",
    "omega_vac1_virW",
    "omega_vac2_virW",
    "delta_vac1_virW",
    "delta_vac2_virW",
    "omega_vac1_virM",
    "omega_vac2_virM",
    "delta_vac1_virM",
    "delta_vac2_virM",
    "prob_deceasing",
]

activation_function = "tanh"
depth = 5
width = np.array(np.repeat(6, depth))


def get_parameter_and_rules_NN(
    vaccination_parameters, inputs, activation_function, depth, width
):

    width = np.append(width, len(vaccination_parameters))
    par_layer = {}
    for index_depth in range(depth + 1):
        for index_width in range(width[index_depth]):
            name = f"L{index_depth}{index_width}"
            par_layer[name] = {"value": 0, "constant": False}

    par_inner_prod = {}
    for index_depth in range(depth):
        for index_width in range(width[index_depth]):
            for index_width_next in range(width[index_depth + 1]):
                current_layer = f"L{index_depth}{index_width}"
                next_layer = f"L{index_depth+1}{index_width_next}"
                name = f"w_{current_layer}_{next_layer}"
                par_inner_prod[name] = {"value": 0, "constant": True}

    rules = {}
    for index_depth in range(1, depth + 1):
        for index_width in range(width[index_depth]):
            formula = f"{activation_function}("
            for index_width_prev in range(width[index_depth - 1]):
                current_layer = f"L{index_depth}{index_width}"
                prev_layer = f"L{index_depth-1}{index_width_prev}"
                w = f"w_{prev_layer}_{current_layer}"
                formula += f"+ {prev_layer} * {w} "
            formula += ")"
            rule_name = f"rule_{current_layer}"

            rules[rule_name] = {"parameter_id": current_layer, "formula": formula}

    par_w_first = {}
    for index_width in range(width[0]):
        for index_inputs in inputs:
            layer = f"L0{index_width}"
            name = f"w_{index_inputs}_{layer}"
            par_w_first[name] = {"value": 0, "constant": True}

    rules_first = {}
    for index_width in range(width[0]):
        formula = f"{activation_function}("
        for index_inputs in inputs:
            layer = f"L0{index_width}"
            w = f"w_{index_inputs}_{layer}"
            formula += f"+ {index_inputs}*{w}"
        formula += ")"
        rule_name = f"rule_{index_width}"
        rules_first[rule_name] = {"parameter_id": layer, "formula": formula}

    rules_vaccines = {}

    for index_width in range(width[-1]):
        fraction = vaccination_parameters[index_width]
        formula = f"L{depth}{index_width}"
        name = f"rules_{fraction}"
        rules_vaccines[name] = {"parameter_id": fraction, "formula": formula}

    rules_vaccines["rules_proportion_countryB_vac1"] = {
        "parameter_id": "proportion_countryB_vac1",
        "formula": "1 -proportion_countryA_vac1",
    }
    rules_vaccines["rules_proportion_countryB_vac2"] = {
        "parameter_id": "proportion_countryB_vac2",
        "formula": "1 -proportion_countryA_vac2",
    }

    all_parameter = {**par_layer, **par_inner_prod, **par_w_first}
    all_rules = {**rules, **rules_first, **rules_vaccines}

    return {"parameters": all_parameter, "rules": all_rules}


pars_rules = get_parameter_and_rules_NN(
    vaccination_parameters, inputs, activation_function, depth, width
)


additional_parameters = {
    **xx_params,
    **parameter_vacc_supply,
    **pars_rules["parameters"],
}
additional_rules = {**pars_rules["rules"], **vaccine_supply_rules}


model_vaccination_create_sbml(
    path=path_sbml,
    areas=areas,
    parameter_rules=additional_rules,
    distances=distances,
    additional_parameters=additional_parameters,
)


observables = {"observable_time": {"name": "t", "formula": "t"}}

for index in ["nu", "proportion"]:
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
    observables=observables,
    only_import=not first_time,
)


created_model = {
    "model": model_and_solver["model"],
    "solver": model_and_solver["solver"],
    "observables": observables,
}


# -----------------------Run Model---------------------------------------------
length = 140


model = created_model["model"]
# dir(model)
# par_values = model.getParameters()
# par_names = model.getParameterNames()
# model.getFixedParameterNames()
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

parameter_optimization = [x for x in model.getParameterNames() if "w" in x]


dict_try = dict(
    zip(parameter_optimization, np.random.uniform(0, 1, len(parameter_optimization)))
)

set_parameter_current = {**set_fixed_parameter, **dict_try}

model_results_current = run_model(
    model=model,
    solver=solver,
    periods=1,
    length_periods=length,
    set_start_parameter=set_start_parameter,
    set_parameter=set_parameter_current,
)

current_deceased_A = 1205808.9050508724
current_deceased_B = 1369784.521203611

# -----------------optimize--------------------------------------------------------
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

    set_proportions = dict(zip(yy_names, theta))
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
        yy_names=parameter_optimization,
        current_deceased_A=current_deceased_A,
        current_deceased_B=current_deceased_B,
    )

    percentage_A = deceased_A / current_deceased_A
    percentage_B = deceased_B / current_deceased_B

    if (percentage_A > 1) or (percentage_B > 1):
        value = large_value
        value = deceased_A + deceased_B
    else:
        value = deceased_A + deceased_B

    # value = np.max([percentage_A, percentage_B])

    return float(value)


np.random.seed(12345)
lower_bound = 0
upper_bound = 1
n_multi = 2
start_matrix = np.array(np.repeat(np.nan, len(parameter_optimization)))
success_values = []
k = 0
n_rows = 0
while n_rows < n_multi and k < 10000:
    candidate = np.random.uniform(lower_bound, upper_bound, len(parameter_optimization))
    cand_value = func_to_optimize_max(candidate)
    print(cand_value)

    if cand_value < large_value:
        start_matrix = np.vstack((start_matrix, candidate))
        success_values.append(cand_value)
        n_rows = start_matrix.shape[0] - 1
        print("success")
    print(k)
    k += 1

test2 = dict(zip(model.getParameterNames(), model.getParameters()))

start_matrix = np.delete(start_matrix, (0), axis=0)

objective = pypesto.Objective(fun=func_to_optimize_max)


def obj_func(theta):

    grad = pypesto.FD(objective, method="central").get_grad(
        x=np.repeat(0.4, len(len(parameter_optimization)))
    )
    hess = pypesto.FD(objective, method="forward").get_hess(
        x=np.repeat(0.4, len(len(parameter_optimization)))
    )

    return [func_to_optimize_max(theta=theta), grad, hess]
    # return np.array([func_to_optimize_max(theta=theta)])


objective2 = pypesto.Objective(fun=func_to_optimize_max, grad=True, hess=True)
objective3 = pypesto.Objective(fun=obj_func, grad=True, hess=True)

number_vaccines = len(vaccination_states) - 1
number_parameters = len(len(parameter_optimization))
lb = np.repeat(lower_bound, number_parameters)
ub = np.repeat(upper_bound, number_parameters)


history_options = pypesto.HistoryOptions(trace_record=True)
engine = pypesto.engine.MultiProcessEngine()
engine.n_procs = 8

# problem = pypesto.Problem(objective=objective3, lb=lb, ub=ub, x_guesses=start_matrix)
# optimizer = optimize.FidesOptimizer(options={"fatol" : 1})
problem = pypesto.Problem(
    objective=objective, lb=lb, ub=ub
)  # , x_guesses=start_matrix)
optimizer = optimize.ScipyOptimizer("L-BFGS-B")

results_max = optimize.minimize(
    problem=problem,
    optimizer=optimizer,
    n_starts=n_multi,
    history_options=history_options,
    # engine=engine,
)

visualize.waterfall(results_max)

df_results_max = results_max.optimize_result.as_dataframe()
theta_optimal_max = df_results_max.iloc[0]["x"]
