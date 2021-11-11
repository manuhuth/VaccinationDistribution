import numpy as np
import pickle

from models.vaccination.create_model_vaccination import create_model_splines
from functions.tools import run_intervention
from functions.tools import simulate_intervention, get_sum_of_states

create_model = False
model_name = "vaccination_multi_test"
path_sbml = f"stored_models/{model_name}/" + model_name
model_directory = "stored_models/" + model_name + "/vaccination_dir"

max_T = 200

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


# type (used for name while saving pickle object)
appears_week = 4

for wor in [12]:

    # wor = 12
    specification = f"reaction_{wor}"
    path = f"/home/manuel/Documents/VaccinationDistribution/code/objects/{specification}.pkl"

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
    timepoints = np.linspace(0, max_T, 1000)
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

    eta = {"eta_virus2": 2.0}

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
        "omega_vac1_virus1": 0.9,
        "omega_vac1_virus2": 0.9,
        "omega_vac2_virus1": 0.9,
        "omega_vac2_virus2": 0.9,
    }

    infectious_t0 = {
        "infectious_countryA_vac0_virus1_t0": 10,
        "infectious_countryB_vac0_virus1_t0": 10,
    }

    mutant_appears = {
        "virus2_countryB_appears_time": appears_week * 7,
        "virus2_countryB_appears_quantity": 20,
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
        **mutant_appears,
    }

    # compute pareto optimal strategy for intervention case
    solver.setAbsoluteTolerance(1e-01)
    solver.setRelativeTolerance(1e-01)
    # run_model_intervention:
    dict_out = run_intervention(
        model,
        solver,
        par_to_optimize,
        areas,
        parameters,
        week_of_reaction=wor,
        T=200,
        n_grid=1000,
        number_yy=11,
        lb=10 ** -8,
        starts_bobyqa=50,
        pop_size=50,
        number_generations=50,
        p_seed=123,
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

# ---------------------------------
o = simulate_intervention(
    model, solver, parameters, 200, 1000, ["countryA", "countryB"], par_optim_dict=None
)
o["trajectories"]["t"] = np.linspace(0, 200, 1000)
import matplotlib.pyplot as plt

cols = [x for x in o["trajectories"].columns if "infectious" in x and "virus2" in x]
y = o["trajectories"][cols].sum(axis=1).reset_index(drop=True)
plt.plot(o["trajectories"]["t"] / 7, y)
plt.scatter(np.linspace(0, 200, 11) / 7, y.iloc[np.linspace(0, 200, 11) / 200 * 999])
