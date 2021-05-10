import numpy as np

from models.vaccination.create_model_vaccination import model_vaccination_create_sbml

from visualization.model_results import plot_states
from visualization.model_results import plot_observables
from visualization.model_results import get_substates
from visualization.model_results import get_observables_by_name

from functions.run_sbml import get_model_and_solver_from_sbml
from functions.run_sbml import run_model
from functions.run_sbml import create_observables_vaccination_rates

from functions.vaccine_proportions import (
    create_rules_vaccination_proportion_relative_population,
)
from functions.vaccine_proportions import (
    create_rules_vaccination_proportion_relative_infected_population,
)


# --------------------------Create Model--------------------------------------
path_sbml = "stored_models/vaccination/vaccination"
vaccination_states = ["vac0", "vac1", "vac2"]
non_vaccination_state = "vac0"
virus_states = ["virW", "virM"]
areas = ["countryA", "countryB"]
species_comp = ["susceptible", "infectious", "recovered", "dead"]

# TODO: write more function for proportions
rules_proportions = create_rules_vaccination_proportion_relative_infected_population(
    species_comp, vaccination_states, non_vaccination_state, virus_states, areas
)

model_vaccination_create_sbml(
    path=path_sbml,
    areas=areas,
    parameter_rules=rules_proportions,
    distances=np.array([[0, 10], [10, 0]]),
    t0_susceptible=0,
    t0_infectious=1000,
)

observables_nu = create_observables_vaccination_rates(
    vaccination_states_removed=["vac1", "vac2"], areas=areas, name_parameter="nu"
)
observables_proportion = create_observables_vaccination_rates(
    vaccination_states_removed=["vac1", "vac2"],
    areas=areas,
    name_parameter="proportion",
)
observables = {**observables_nu, **observables_proportion}

model_and_solver = get_model_and_solver_from_sbml(
    path_sbml=path_sbml,
    model_name="vaccination",
    model_directory="stored_models/vaccination/vaccination_dir",
    observables=observables,
)

# -----------------------Run Model---------------------------------------------
model = model_and_solver["model"]
solver = model_and_solver["solver"]
set_start_parameter = {
    "susceptible_countryA_vac0_t0": 60000,
    "susceptible_countryB_vac0_t0": 40000,
    "infectious_countryA_vac0_virW_t0": 1000,
    "infectious_countryA_vac0_virM_t0": 1000,
    "infectious_countryB_vac0_virW_t0": 1000,
    "infectious_countryB_vac0_virM_t0": 1000,
}
set_parameter = {
    "beta": 2,
    "lambda1": 0.01,
    "p": 0.3,
    "number_vac1": 10,
    "number_vac2": 20,
}
observables_names = observables.keys()

model_results = run_model(
    model=model,
    solver=solver,
    periods=1,
    length_periods=12,
    set_start_parameter=set_start_parameter,
    set_parameter=set_parameter,
    observables_names=observables.keys(),
)

trajectory_states = model_results["states"]
trajectory_observables = model_results["observables"]

# visualization
observables_names_nu = get_observables_by_name(
    observables, substrings=["nu"], include_all=True
)
observables_names_proportion = get_observables_by_name(
    observables, substrings=["proportion"], include_all=True
)
# TODO check if values make sense
plot_observables(
    results=model_results, model=model, observable_ids=observables_names_nu
)
plot_observables(
    results=model_results,
    model=model,
    observable_ids=observables_names_proportion,
    set_off_scientific_notation=True,
    decimal_floats=4,
)
# fig, ax = plot_states(results=model_result, model=model)

# dir(model)
# par_values = model.getParameters()
# par_names = model.getParameterNames()
# model.getFixedParameterNames()

substates = get_substates(
    model=model, substrings=["infectious", "vac0"], include_all=True
)
fig, ax = plot_states(results=model_results, model=model, state_ids=substates)
