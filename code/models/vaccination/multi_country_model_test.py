import numpy as np
import pickle

from models.vaccination.create_model_vaccination import create_model_splines
from functions.tools import get_model_output

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
    number_xx_R=3,
)

model = created_model["model"]
solver = created_model["solver"]

vac = "Equal"
initial_states = "Equal"  # make equal


# type (used for name while saving pickle object)
specification = f"inital{initial_states}_vac{vac}"
path = (
    f"/home/manuel/Documents/VaccinationDistribution/code/objects/{specification}.pkl"
)


if vac == "Equal":
    par_to_optimize = [x for x in model.getParameterNames() if "vac1" in x]
    vaccine_supply_parameter_values = np.concatenate(
        (
            np.repeat(  # change here
                float(120000), (created_model["information"]["number_yy"] - 1)
            ),
            np.repeat(float(0), (created_model["information"]["number_yy"] - 1)),
        )
    )

    delta = {
        "delta_vac1_virus1": 0.95,
        "delta_vac1_virus2": 0.6,
        "delta_vac2_virus1": 0,
        "delta_vac2_virus2": 0,
    }

    omega = {
        "omega_vac1_virus1": 0.9,
        "omega_vac1_virus2": 0.9,
        "omega_vac2_virus1": 0,
        "omega_vac2_virus2": 0,
    }

elif vac == "Unequal":
    par_to_optimize = model.getParameterNames()
    vaccine_supply_parameter_values = np.repeat(
        float(60000), 2 * (created_model["information"]["number_yy"] - 1)
    )
    delta = {
        "delta_vac1_virus1": 0.95,
        "delta_vac1_virus2": 0.6,
        "delta_vac2_virus1": 0.6,
        "delta_vac2_virus2": 0.95,
    }

    omega = {
        "omega_vac1_virus1": 0.9,
        "omega_vac1_virus2": 0.9,
        "omega_vac2_virus1": 0.9,
        "omega_vac2_virus2": 0.9,
    }

if initial_states == "Equal":
    infectious_t0 = {
        "infectious_countryA_vac0_virus1_t0": 5,
        "infectious_countryA_vac0_virus2_t0": 5,
        "infectious_countryB_vac0_virus1_t0": 5,
        "infectious_countryB_vac0_virus2_t0": 5,
    }
elif initial_states == "Unequal":
    infectious_t0 = {
        "infectious_countryA_vac0_virus1_t0": 10,
        "infectious_countryA_vac0_virus2_t0": 0,
        "infectious_countryB_vac0_virus1_t0": 0,
        "infectious_countryB_vac0_virus2_t0": 10,
    }


par_to_optimize = [x for x in par_to_optimize if float(x.rsplit("_", 1)[-1]) <= 8]

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
    "prob_deceasing": 0.02,
    "gamma": 0.5,
    "beta": 0.24,
}

R0s = {
    "yR0_countryA_0": 1,
    "yR0_countryA_1": 1,
    "yR0_countryA_2": 1,
    "yR0_countryB_0": 1,
    "yR0_countryB_1": 1,
    "yR0_countryB_2": 1,
}

eta = {"eta_virus2": 2}

susceptible_t0 = {
    "susceptible_countryA_vac0_t0": 10 ** 7,
    "susceptible_countryB_vac0_t0": 10 ** 7,
}


distances = {
    "distance_countryA_countryB": 1 / 1000,
    "distance_countryB_countryA": 1 / 1000,
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
}

solver.setAbsoluteTolerance(1e-01)
solver.setRelativeTolerance(1e-01)
dict_out = get_model_output(
    model,
    solver,
    parameters,
    areas,
    par_to_optimize,
    n_starts_pb=50,
    n_starts_pareto=500,
    number_generations_pareto=50,
)

with open(
    path,
    "wb",
) as output:
    out = dict_out
    pickle.dump(out, output, pickle.HIGHEST_PROTOCOL)
