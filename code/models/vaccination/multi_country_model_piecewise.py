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


first_time = True
# --------------------------Create Model--------------------------------------
model_name = "vaccination_multi_piecewise"
path_sbml = "stored_models/vaccination_multi_piecewise/" + model_name
vaccination_states = ["vac0", "vac1", "vac2", "vac3"]
non_vaccination_state = general_set_up()["non_vaccination_state"]
virus_states = ["virW", "virM", "virMM"]
areas = ["countryA", "countryB", "countryC"]

distances = pd.DataFrame(
    np.matrix([[0, 4999, 3999], [4999, 0, 2999], [3999, 2999, 0]]),
    index=areas,
    columns=areas,
)
species_comp = general_set_up()["species_compartments"]

length_decision_period = 14
number_decision_periods = 10
first_vaccination_period = 0
n_intervals = 6000

xx_list = ["xx" + str(i) for i in (range(number_decision_periods))]

xx_params = {}
period = 0
for index in xx_list:
    xx_params[index] = {"value": period, "constant": True}
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


additional_parameters = {
    **additional_parameters_pw,
    **xx_params,
    **parameter_vacc_supply,
}

model_vaccination_create_sbml(
    path=path_sbml,
    areas=areas,
    vaccination_states=vaccination_states,
    non_vaccination_state=non_vaccination_state,
    virus_states=virus_states,
    omega_matrix=np.array([[0.5, 0.6, 0.4], [0.6, 0.2, 0.4], [0.6, 0.2, 0.4]]),
    delta_matrix=np.array([[0.5, 0.6, 0.4], [0.6, 0.2, 0.4], [0.6, 0.2, 0.4]]),
    eta_vector=np.array([[1, 1.2, 1.3]]),
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


observables = {
    **observables_nu,
    **observables_proportion,
    **observables_vaccinated_vacc,
}

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
