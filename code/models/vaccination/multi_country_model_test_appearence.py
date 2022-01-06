import numpy as np
import pickle
import multiprocessing as mp
import tqdm

from models.vaccination.create_model_vaccination import create_model_splines
from functions.tools import run_intervention
from functions.tools import simulate_intervention, get_sum_of_states

create_model = False
model_name = "vaccination_multi_test"
path_sbml = f"stored_models/{model_name}/" + model_name
model_directory = "stored_models/" + model_name + "/vaccination_dir"

max_T = 294

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
appears_week = 12

def run_week_appearence(appears_week):
    for wor in [0]:
    
        #wor = 0
        specification = f"appears_week_{appears_week}"
        path = f"/home/manuel/Documents/VaccinationDistribution/code/objects/{specification}.pkl"
    
        # par_to_optimize = [x for x in model.getParameterNames() if "vac1" in x]
        par_to_optimize = list(model.getParameterNames())
        
        vac_times_par = [f"xx{i}" for i in range(15)]
        vac_times = np.linspace(0, max_T, 15)
        vacs = dict(zip(vac_times_par, vac_times))
        
        vac_two_weeks = 60000
        vaccine_supply_parameter_values = np.concatenate(
            (
                np.repeat(  # change here
                    float(vac_two_weeks), (created_model["information"]["number_yy"] - 1)
                ),
                np.repeat(float(vac_two_weeks), (created_model["information"]["number_yy"] - 1)),
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
        
        spline_R_val = {"xx_R0_0":0,
                        "xx_R0_1":1,
                        "xx_R0_2":max_T}
        
    
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
            **vacs,
            **spline_R_val,
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
            number_yy=15,
            lb=10 ** -8,
            starts_bobyqa=5,
            pop_size=500,
            number_generations=50,
            p_seed=1235,
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
        
    return dict_out

#----------------------------------------------------------------------------
with mp.Pool() as p:
    results = list(
        tqdm.tqdm(
            p.imap_unordered(
                run_week_appearence, [0.05, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5]
            ),
            total=13,
        )
    )


# ---------------------------------
o = simulate_intervention(
    model, solver, parameters, max_T, 1001, ["countryA", "countryB"], par_optim_dict=None
)
o["trajectories"]["t"] = np.linspace(0, max_T, 1000)
import matplotlib.pyplot as plt

cols = [x for x in o["trajectories"].columns if "infectious" in x]# and "virus2" in x]
y = o["trajectories"][cols].sum(axis=1).reset_index(drop=True)
plt.plot(o["trajectories"]["t"] / 7, y)
plt.scatter(np.linspace(0, max_T, 11) / 7, y.iloc[np.linspace(0, max_T, 11) / max_T * 999])
