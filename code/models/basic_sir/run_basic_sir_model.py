import numpy as np
from amici.plotting import plotStateTrajectories
from models.create_basic_sir_model import model_basic_sir_create_sbml
from functions.run_sbml import get_model_and_solver_from_sbml
from functions.run_sbml import model_run

path_sbml = 'stored_models/basic_sir/basic_sir'
S_0 = 1000
I_0 = 0.000001

model_basic_sir_create_sbml(path=path_sbml, S_0=S_0, I_0=I_0, R_0=0)
model_and_solver = get_model_and_solver_from_sbml(path_sbml=path_sbml, model_name='basic_sir',
                                                  model_directory='stored_models/basic_sir/basic_sir_dir')

model = model_and_solver['model']
solver = model_and_solver['solver']
timepoints = np.linspace(0,100,100)

model_result = model_run(model=model, solver=solver, timepoints=timepoints, set_parameter=None)

plotStateTrajectories(model_result, model = model)