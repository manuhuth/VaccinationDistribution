import numpy as np
from models.vaccination.create_model_vaccination import model_vaccination_create_sbml
from visualization.model_results import plot_states
from visualization.model_results import get_substates
from functions.run_sbml import get_model_and_solver_from_sbml
from functions.run_sbml import model_run

path_sbml = "stored_models/vaccination/vaccination"

model_vaccination_create_sbml(
    path=path_sbml,
    areas=["countryA"],
    distances=np.array([[0]]),
    t0_susceptible=30000,
    t0_infectious=1000,
)
model_and_solver = get_model_and_solver_from_sbml(
    path_sbml=path_sbml,
    model_name="vaccination",
    model_directory="stored_models/vaccination/vaccination_dir",
)

model = model_and_solver["model"]
solver = model_and_solver["solver"]
timepoints = np.linspace(0, 20, 20)

model_result = model_run(
    model=model,
    solver=solver,
    timepoints=timepoints,
    set_parameter={"beta": 4, "susceptible_countryA_vac0_t0": 40000, 'lambda1' : 0.4, 'p' : 0.8 },
)


fig, ax = plot_states(results=model_result, model=model)

# dir(model)
# model.getParameters()
# model.getParameterNames()

substates = get_substates(model=model, substrings=["dead"])
fig, ax = plot_states(results=model_result, model=model, state_ids=substates)

# TODO Find out how to code time dependent rules