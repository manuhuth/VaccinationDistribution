import numpy as np
from amici.plotting import plotStateTrajectories
from models.one_country.create_model_one_country import model_one_country_create_sbml
from functions.run_sbml import get_model_and_solver_from_sbml
from functions.run_sbml import model_run

path_sbml = "stored_models/one_country/one_country"

model_one_country_create_sbml(path=path_sbml, check_error=False)
model_and_solver = get_model_and_solver_from_sbml(
    path_sbml=path_sbml,
    model_name="one_country",
    model_directory="stored_models/one_country/one_country_dir",
)

model = model_and_solver["model"]
solver = model_and_solver["solver"]
timepoints = np.linspace(0, 100, 100)

model_result = model_run(
    model=model, solver=solver, timepoints=timepoints, set_parameter=None
)

plotStateTrajectories(model_result, model=model)
