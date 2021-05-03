import numpy as np
from models.one_country.create_model_one_country import model_one_country_create_sbml
from visualization.model_results import plot_states
from visualization.model_results import get_substates
from functions.run_sbml import get_model_and_solver_from_sbml
from functions.run_sbml import model_run

path_sbml = "stored_models/one_country/one_country"

model_one_country_create_sbml(
    path=path_sbml, check_error=False, t0_susceptible=30000, t0_infectious=1000
)
model_and_solver = get_model_and_solver_from_sbml(
    path_sbml=path_sbml,
    model_name="one_country",
    model_directory="stored_models/one_country/one_country_dir",
)

model = model_and_solver["model"]
solver = model_and_solver["solver"]
timepoints = np.linspace(0, 100, 100)

model_result = model_run(
    model=model,
    solver=solver,
    timepoints=timepoints,
    set_parameter={"beta": 4, "susceptible_countryA_vac0_t0": 40000},
)


fig, ax = plot_states(results=model_result, model=model)

#
substates = get_substates(model=model, substrings=["infectious"])
fig, ax = plot_states(results=model_result, model=model, state_ids=substates)
