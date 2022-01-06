import numpy as np
import pickle
import multiprocessing as mp
import tqdm

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

vac = "Unequal" #only 
initial_states = "Unequal"

dict_R_cases = {}
case = 0
for index1 in range(3, 11):
    for index2 in range(3,11):
        dict_R_cases[f"case{case}"] = np.array([index1, index2])/10
        case += 1

def run_n_vacc(case_R):
    
    
    
    initial = 0
    n_vacc = 60000
    # type (used for name while saving pickle object)
    model = created_model["model"]
    solver = created_model["solver"]
    
    vac = "Unequal" #only 
    initial_states = "Unequal"  
    specification = f"inital{initial_states}_vac{vac}"
    path = (
        f"/home/manuel/Documents/VaccinationDistribution/code/objects/degree_R_{case_R}.pkl"
    )
    
    
    if vac == "Equal":
        par_to_optimize = [x for x in model.getParameterNames() if "vac1" in x]
        vaccine_supply_parameter_values = np.concatenate(
            (
                np.repeat(  # change here
                    float(n_vacc), (created_model["information"]["number_yy"] - 1)
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
            float(n_vacc), 2 * (created_model["information"]["number_yy"] - 1)
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
            "infectious_countryA_vac0_virus1_t0": 10 - initial,
            "infectious_countryA_vac0_virus2_t0": 0 + initial,
            "infectious_countryB_vac0_virus1_t0": 0 + initial,
            "infectious_countryB_vac0_virus2_t0": 10 - initial,
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
    
    case = dict_R_cases[f"case{case_R}"]
    A = case[0]
    B = case[1]
    
    R0s = {
        "yR0_countryA_0": A,
        "yR0_countryA_1": A,
        "yR0_countryA_2": A,
        "yR0_countryB_0": B,
        "yR0_countryB_1": B,
        "yR0_countryB_2": B,
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
        n_starts_pb=5,
        n_starts_pareto=500,
        number_generations_pareto=50,
    )
    

    with open(
        path,
        "wb",
    ) as output:
        out = dict_out
        pickle.dump(out, output, pickle.HIGHEST_PROTOCOL)
    
    return dict_out

with mp.Pool() as p:
    results = list(
        tqdm.tqdm(
            p.imap_unordered(
                run_n_vacc, range(len(dict_R_cases.keys()))
            ),
            total=len(dict_R_cases.keys()),
        )
    )

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

    with open(
        "/home/manuel/Documents/VaccinationDistribution/code/objects/R_values_parallel.pkl",
        "wb",
    ) as output:
        out = results
        pickle.dump(out, output, pickle.HIGHEST_PROTOCOL)