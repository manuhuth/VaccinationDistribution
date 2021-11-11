import numpy as np
import amici
import pandas as pd
import pypesto

from models.vaccination.create_model_vaccination import create_model_splines
from functions.tools import get_name_parameter_to_optimize
from functions.tools import get_sum_of_states
from functions.tools import run_optimization

create_model = True
model_name = "vaccination_multi"
path_sbml = f"stored_models/{model_name}/" + model_name
model_directory = "stored_models/" + model_name + "/vaccination_dir"

max_T = 140

created_model = create_model_splines(
    create_model=create_model,
    number_areas=2,
    number_vaccines=2,
    vaccinated_compartments=["susceptible", "infectious"],
    number_viruses=2,
    length_decision_period=max_T,
    number_yy=int(max_T / 14 + 1),
    model_name=model_name,
    path_sbml=path_sbml,
    model_directory=model_directory,
)

model = created_model["model"]
solver = created_model["solver"]
solver.setSensitivityOrder(amici.SensitivityOrder_first)
solver.setSensitivityMethod(amici.SensitivityMethod_forward)


# set Parameters
timepoints = np.linspace(0, max_T, 6000)
model.setTimepoints(timepoints)


vaccine_supply_parameter_strings = [
    x for x in model.getParameterNames() if "vaccine_supply" in x
]
vaccine_supply_parameter_values = np.repeat(
    float(100), 2 * (created_model["information"]["number_yy"] - 1)
)

vaccine_supply_parameter = dict(
    zip(vaccine_supply_parameter_strings, vaccine_supply_parameter_values)
)

general_parameters = {
    "lambda1": 0.1,
    "prob_deceasing": 0.025,
    "gamma": 0.5,
    "beta": 0.3,
}

omega = {
    "omega_vac1_virus1": 0.99,
    "omega_vac1_virus2": 0.99,
    "omega_vac2_virus1": 0.99,
    "omega_vac2_virus2": 0.99,
}

delta = {
    "delta_vac1_virus1": 0.95,
    "delta_vac1_virus2": 0.9,
    "delta_vac2_virus1": 0.9,
    "delta_vac2_virus2": 0.95,
}

eta = {"eta_virus1": 1.2}

susceptible_t0 = {
    "susceptible_countryA_vac0_t0": 10 ** 6,
    "susceptible_countryB_vac0_t0": 10 ** 6,
}

infectious_t0 = {
    "infectious_countryA_vac0_virus1_t0": 10,
    "infectious_countryA_vac0_virus2_t0": 10,
    "infectious_countryB_vac0_virus1_t0": 10,
    "infectious_countryB_vac0_virus2_t0": 10,
}

distances = {
    "distance_countryA_countryB": 4999,
    "distance_countryB_countryA": 4999,
}

parameters = {
    **general_parameters,
    **omega,
    **delta,
    **eta,
    **susceptible_t0,
    **infectious_t0,
    **distances,
    **vaccine_supply_parameter,
}

for keys in parameters.keys():
    model.setFixedParameterByName(keys, parameters[keys])


# run current case
rdata = amici.runAmiciSimulation(model, solver)
states = pd.DataFrame(rdata["x"], columns=model.getStateIds())
observables = pd.DataFrame(rdata["y"], columns=model.getObservableIds())
countryA = get_sum_of_states(
    model, states, state_type=["countryA", "dead"], final_amount=True
)
countryB = get_sum_of_states(
    model, states, state_type=["countryB", "dead"], final_amount=True
)

infectious = get_sum_of_states(
    model, states, state_type=["infectious"], final_amount=True
)

# get parameter to optimize over
parameter_to_optimize = get_name_parameter_to_optimize(
    observables, created_model, timepoints
)

# last_parameter = int(parameter_to_optimize[-1].rsplit("_", 1)[-1])
# par_optimization = [
#     x
#     for x in model.getParameterNames()
#     if "vaccine_supply" in x and float(x.rsplit("_", 1)[-1][-1]) > last_parameter
#    ]
# for i in par_optimization:
#    model.setParameterByName(i, 0)


# optimize
results_pypesto = run_optimization(
    model,
    solver,
    parameter_to_optimize,
    optimizer=pypesto.optimize.ScipyOptimizer("L-BFGS-B"),
    n_starts=50,
    lower_bound=10 ** (-15),
    trace_record=True,
    pareto_values=None,
)
results_pypesto_df = results_pypesto.optimize_result.as_dataframe()

results_pyboyqa = run_optimization(
    model,
    solver,
    parameter_to_optimize,
    optimizer="pybobyqa",
    n_starts=20,
    lower_bound=10 ** (-15),
    trace_record=True,
    pareto_values=None,
)
