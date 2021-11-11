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
)

model = created_model["model"]
solver = created_model["solver"]

vac = "Equal"  # 1
initial_states = "Equal"  # 2


# type (used for name while saving pickle object)
specification = f"inital{initial_states}_vac{vac}"
path = f"/home/manuel/Documents/VaccinationDistribution/code/objects/{specification}_distance.pkl"


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
        "omega_vac1_virus1": 0.99,
        "omega_vac1_virus2": 0.99,
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
        "omega_vac1_virus1": 0.99,
        "omega_vac1_virus2": 0.99,
        "omega_vac2_virus1": 0.99,
        "omega_vac2_virus2": 0.99,
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
    "prob_deceasing": 0.04,
    "gamma": 0.5,
    "beta": 0.3,
}

R0s = {"R0_modifier_countryA": 1, "R0_modifier_countryB": 1}

eta = {"eta_virus2": 1.3}

susceptible_t0 = {
    "susceptible_countryA_vac0_t0": 10 ** 7,
    "susceptible_countryB_vac0_t0": 10 ** 7,
}


import pandas as pd

distances = (
    list(np.linspace(10 ** -15, 4 * 10 ** -4, 10))
    + list(np.linspace(5 * 10 ** -4, 10 ** -3, 4))
    + list(np.linspace(2 * 10 ** -3, 1, 10))
)
cols = areas + ["fval"]
save_pop = pd.DataFrame(
    np.nan, columns=cols + ["distance"], index=range(len(distances))
)
save_pareto = pd.DataFrame(
    np.nan, columns=cols + ["distance"], index=range(len(distances))
)
save_optimal = pd.DataFrame(
    np.nan, columns=cols + ["distance"], index=range(len(distances))
)


def run_multi_distances(dist):

    distances_dict = {
        "distance_countryA_countryB": dist,
        "distance_countryB_countryA": dist,
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
        **distances_dict,
        **vaccine_supply_parameter,
        **R0s,
    }

    dict_out = get_model_output(
        model,
        solver,
        parameters,
        areas,
        par_to_optimize,
        n_starts_pb=1,
        n_starts_pareto=100,
        number_generations_pareto=100,
    )

    pop_based = pd.Series(
        dict_out["population_based"] + [np.sum(dict_out["population_based"])],
        index=cols,
    )

    pareto_dfs = dict_out["pareto_improvements"].reset_index()
    pareto_optimal = pareto_dfs.iloc[np.argmin(pareto_dfs["fval"])][cols]

    optimal_df = dict_out["all_strategies"].reset_index()
    unrestricted_optimal = optimal_df.iloc[np.argmin(optimal_df["fval"])][cols]

    return {
        "pop": pop_based,
        "pareto": pareto_optimal,
        "optimal": unrestricted_optimal,
        "distance": dist,
    }


import multiprocessing as mp
import tqdm

with mp.Pool() as p:
    results = list(
        tqdm.tqdm(
            p.imap_unordered(run_multi_distances, distances), total=len(distances)
        )
    )


for i in range(len(results)):
    di = results[i]
    dist = di["distance"]
    pop_based = di["pop"]
    pareto_optimal = di["pareto"]
    unrestricted_optimal = di["optimal"]
    save_pop.iloc[i] = pop_based.append(pd.Series(dist, index=["distance"]))
    save_pareto.iloc[i] = pareto_optimal.append(pd.Series(dist, index=["distance"]))
    save_optimal.iloc[i] = unrestricted_optimal.append(
        pd.Series(dist, index=["distance"])
    )


dict_end = {"pop_based": save_pop, "pareto": save_pareto, "optimal": save_optimal}


with open(
    path,
    "wb",
) as output:
    out = dict_end
    pickle.dump(out, output, pickle.HIGHEST_PROTOCOL)
