import numpy as np

from models.vaccination.create_model_vaccination import model_vaccination_create_sbml
from visualization.model_results import plot_states
from visualization.model_results import plot_observables
from visualization.model_results import get_substates
from functions.run_sbml import get_model_and_solver_from_sbml
from functions.run_sbml import model_run
from functions.run_sbml import create_observables_vaccination_rates


path_sbml = "stored_models/vaccination/vaccination"

#TODO: write function for proportions
areas = ["countryA", "countryB"]
model_vaccination_create_sbml(
    path=path_sbml,
    areas=areas,
    distances=np.array([[0, 10], [10, 0]]),
    t0_susceptible=0,
    t0_infectious=1000,
)

observables = create_observables_vaccination_rates(vaccination_states_removed=['vac1', 'vac2'], areas=areas)

model_and_solver = get_model_and_solver_from_sbml(
    path_sbml=path_sbml,
    model_name="vaccination",
    model_directory="stored_models/vaccination/vaccination_dir",
    observables=observables,
)

model = model_and_solver["model"]
solver = model_and_solver["solver"]
timepoints = np.linspace(0, 20, 30)

model_result = model_run(
    model=model,
    solver=solver,
    timepoints=timepoints,
    set_parameter={
        "beta": 4,
        "susceptible_countryA_vac0_t0": 400000,
        "susceptible_countryB_vac0_t0": 40000,
        "lambda1": 0.01,
        "p": 0.3,
        "number_vac1": 100,
        "number_vac2": 200,
        'proportion_countryA_vac1' : 0.7,
        'proportion_countryA_vac2' : 0.5,
        'proportion_countryB_vac1' : 0.3,
        'proportion_countryB_vac2' : 0.5,
    },
)

#TODO rewrite plot function
plot_observables(results=model_result, model = model)
fig, ax = plot_states(results=model_result, model=model)

# dir(model)
par_values = model.getParameters()
par_names = model.getParameterNames()
model.getFixedParameterNames()

substates = get_substates(model=model, substrings=['dead', 'countryA'], include_all=True)
fig, ax = plot_states(results=model_result, model=model, state_ids=substates)


