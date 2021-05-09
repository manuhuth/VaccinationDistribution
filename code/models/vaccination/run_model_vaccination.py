import numpy as np

from models.vaccination.create_model_vaccination import model_vaccination_create_sbml
from visualization.model_results import plot_states
from visualization.model_results import plot_observables
from visualization.model_results import get_substates
from visualization.model_results import get_observables_by_name
from functions.run_sbml import get_model_and_solver_from_sbml
from functions.run_sbml import model_run
from functions.run_sbml import create_observables_vaccination_rates
from functions.run_sbml import create_rules_vaccination_proportion_relative_population



path_sbml = "stored_models/vaccination/vaccination"
vaccination_states = ["vac0", "vac1", "vac2"]
non_vaccination_state = "vac0"
virus_states = ["virW", "virM"]
areas = ["countryA", "countryB"]
species_comp = ["susceptible", "infectious", "recovered", "dead"]

# create rules for vaccination proportions
# TODO: write more function for proportions
rules_proportions = create_rules_vaccination_proportion_relative_population(
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
        #        'proportion_countryA_vac1' : 0.7,
        #        'proportion_countryA_vac2' : 0.5,
        #        'proportion_countryB_vac1' : 0.3,
        #        'proportion_countryB_vac2' : 0.5,
    },
)


#visualization
observables_names_nu = get_observables_by_name(
    observables, substrings=["nu"], include_all=True
)
observables_names_proportion = get_observables_by_name(
    observables, substrings=["proportion"], include_all=True
)
plot_observables(results=model_result, model=model, observable_ids=observables_names_nu)
plot_observables(
    results=model_result, model=model, observable_ids=observables_names_proportion
)
fig, ax = plot_states(results=model_result, model=model)

# dir(model)
par_values = model.getParameters()
par_names = model.getParameterNames()
model.getFixedParameterNames()

substates = get_substates(
    model=model, substrings=["dead", "countryA"], include_all=True
)
fig, ax = plot_states(results=model_result, model=model, state_ids=substates)
