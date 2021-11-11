import numpy as np
import pickle

from models.vaccination.create_model_vaccination import create_model_splines
from functions.tools import run_intervention

create_model = False
model_name = "vaccination_multi_test"
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


# type (used for name while saving pickle object)
appears_week = 4
wor = 8

country = "countryB"  # countryB
reduction_R = 0.4  # 0.4


specification = f"R_value_{wor}_{country}_{reduction_R}"
path = (
    f"/home/manuel/Documents/VaccinationDistribution/code/objects/{specification}.pkl"
)


# par_to_optimize = [x for x in model.getParameterNames() if "vac1" in x]
par_to_optimize = list(model.getParameterNames())
vaccine_supply_parameter_values = np.concatenate(
    (
        np.repeat(  # change here
            float(60000), (created_model["information"]["number_yy"] - 1)
        ),
        np.repeat(float(60000), (created_model["information"]["number_yy"] - 1)),
    )
)

# par_to_optimize = [x for x in par_to_optimize if float(x.rsplit('_', 1)[-1]) <= 8]

# areas
areas = created_model["information"]["areas"]

# set Parameters
timepoints = np.linspace(0, max_T, 6000)
model.setTimepoints(timepoints)


vaccine_supply_parameter_strings = [
    x for x in model.getFixedParameterNames() if "vaccine_supply" in x
]

vaccine_supply_parameter = dict(
    zip(vaccine_supply_parameter_strings, vaccine_supply_parameter_values)
)

general_parameters = {
    "lambda1": 0.1,
    "prob_deceasing": 0.04,
    "gamma": 0.5,
    "beta": 0.3,
}

R0s = {"R0_modifier_countryA": 1, "R0_modifier_countryB": 1}
R0s[f"R0_modifier_{country}"] = reduction_R

eta = {"eta_virus2": 1.3}

susceptible_t0 = {
    "susceptible_countryA_vac0_t0": 10 ** 7,
    "susceptible_countryB_vac0_t0": 10 ** 7,
}

delta = {
    "delta_vac1_virus1": 0.95,
    "delta_vac1_virus2": 0.6,
    "delta_vac2_virus1": 0.6,
    "delta_vac2_virus2": 0.95,
}

omega = {
    "omega_vac1_virus1": 0.99,
    "omega_vac1_virus2": 0.99,
    "omega_vac2_virus1": 0.99,
    "omega_vac2_virus2": 0.99,
}


infectious_t0 = {
    "infectious_countryA_vac0_virus1_t0": 10,
    "infectious_countryB_vac0_virus1_t0": 10,
}

mutant_appears = {
    "virus2_countryB_appears_time": appears_week * 7,
    "virus2_countryB_appears_quantity": 50,
}


distances = {
    "distance_countryA_countryB": 1 / 5000,
    "distance_countryB_countryA": 1 / 5000,
    "distance_countryA_countryA": 1,
    "distance_countryB_countryB": 1,
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
    **R0s,
    **mutant_appears,
}


# compute pareto optimal strategy for intervention case

# run_model_intervention:
dict_out = run_intervention(
    model,
    solver,
    par_to_optimize,
    areas,
    parameters,
    week_of_reaction=wor,
    T=140,
    n_grid=6000,
    number_yy=11,
    lb=10 ** -15,
    starts_bobyqa=30,
    pop_size=40,
    number_generations=400,
    p_seed=12345,
    pareto_seed=1234,
    crossover_eta=15,
    crossover_prob=0.9,
    mutation_eta=20,
)

with open(
    path,
    "wb",
) as output:
    out = dict_out
    pickle.dump(out, output, pickle.HIGHEST_PROTOCOL)
