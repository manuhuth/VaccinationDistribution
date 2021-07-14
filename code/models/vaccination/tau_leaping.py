import pickle
import pandas as pd
import numpy as np
import scipy

from models.vaccination.create_model_vaccination import create_reactions_model
from models.vaccination.parameter import general_set_up
from models.vaccination.create_model_vaccination import b_distance_function
from models.vaccination.create_model_vaccination import create_rules_vaccination_rate
from models.vaccination.create_model_vaccination import get_all_individuals_to_be_vaccinated

with open(
    "/home/manuel/Documents/VaccinationDistribution/code/objects/output_model_props.pkl",
    "rb",
) as input:
    model_out = pickle.load(input)
    
with open(
    "/home/manuel/Documents/VaccinationDistribution/code/objects/output_splines.pkl",
    "rb",
) as input:
    dict_out = pickle.load(input)
df_optimal_results = dict_out["df_optimal_results"]
model_optimal = dict_out["model_optimal"]
model_seperated = dict_out["model_seperated"]
model_current = dict_out["model_current"]
time_points = dict_out["time_points"]

vaccination_states = ["vac0", "vac1", "vac2"]
non_vaccination_state = general_set_up()["non_vaccination_state"]
virus_states = general_set_up()["virus_states"]
areas = general_set_up()["areas"]
distances = general_set_up()["distances"]
species_comp = general_set_up()["species_compartments"]

vaccination_states_removed = [
    x for x in vaccination_states if x != non_vaccination_state
]

species = model_out["state_names"]
vaccinated_compartments = ["susceptible", "infectious", "recovered"]

par_dict = dict(zip(model_out["par_names"], model_out["par_magnitude"]))

distance_parameters = {}
distances_df = pd.DataFrame(distances, columns=areas, index=areas)
b_distance_matrix = b_distance_function(distances_df)
D = pd.DataFrame(np.nan, index=areas, columns=areas)

for index_areas_row in areas:
    for index_areas_col in areas:
            name = f"distance_{index_areas_row}_{index_areas_col}"
            D.loc[index_areas_row, index_areas_col] = name
            distance_parameters[name] = {
                "value": b_distance_matrix.loc[index_areas_row, index_areas_col]
            }

reactions = create_reactions_model(
        areas=areas,
        vaccination_states=vaccination_states,
        vaccination_states_removed=vaccination_states_removed,
        virus_states=virus_states,
        species=dict(zip(species, np.repeat(0,len(species)))),
        distances=D,
        vaccinated_compartments=vaccinated_compartments,
        non_vaccination_state=non_vaccination_state,
    )

#supply_1 = dict(filter(lambda item: "vaccine_supply_par_vac1" in item[0], par_dict.items()))
#supply_2 = dict(filter(lambda item: "vaccine_supply_par_vac2" in item[0], par_dict.items()))
#vaccine_supply_df = pd.DataFrame({"supply_vac1" : supply_1.values(), "supply_vac2" : supply_2.values(), "countryA_vac1" : np.repeat(0.5, 10),
#                                  "countryA_vac2" : np.repeat(0.5, 10), "countryB_vac1" : np.repeat(0.5, 10),
#                                  "countryB_vac2" : np.repeat(0.5, 10),
#                                  "time" : np.linspace(0, (len(supply_1.values())-1)*14, len(supply_1.values()))})

vaccine_supply_df = model_optimal["observables"][["observable_quantity_countryA_vac1", "observable_quantity_countryA_vac2", "observable_quantity_countryB_vac1", "observable_quantity_countryB_vac2"]]
vaccine_supply_df["time"] = time_points

def round_to_base(x, base=14):
    return base * round(x/base)

def get_nu_parameter(species, vaccination_states_removed, areas,
                     vaccinated_compartments, non_vaccination_state, virus_states, name_df = "Y"):
    nu_parameter = {}
    for index_vaccination in vaccination_states_removed:
          for index_areas in areas:
              id_nu = f"nu_{index_areas}_{index_vaccination}"
              to_be_vaccinated = get_all_individuals_to_be_vaccinated(
                  vaccinated_compartments=vaccinated_compartments,
                  non_vaccination_state=non_vaccination_state,
                  virus_states=virus_states,
                  areas=[index_areas],
              )
              
              for index_species in species:
                  new = f"{name_df}[str('{index_species}')][period]"
                  to_be_vaccinated = to_be_vaccinated.replace(index_species, new)
                  
              formula = (
                  f"vaccine_supply_df.loc[vaccine_supply_df[str('time')] == tau, str('observable_quantity_{index_areas}_{index_vaccination}')] / ({to_be_vaccinated})"
              )
              nu_parameter[id_nu] = formula
    return nu_parameter

nu_parameter = get_nu_parameter(species, vaccination_states_removed, areas,
                     vaccinated_compartments, non_vaccination_state, virus_states, name_df = "Y")

def create_formulas_eval(dictionary, species, par_dict, name_dict = "par_dict", name_df = "Y"):

    for index_species in species:
        for index_reactions in dictionary.keys():
            formula = dictionary[index_reactions]["formula"]
            if index_species in formula:
                new = f"{name_df}[str('{index_species}')][period]"
                dictionary[index_reactions]["formula"] = formula.replace(index_species, new)
    
    for index_parameters in par_dict.keys():
        for index_reactions in dictionary.keys():
            formula = dictionary[index_reactions]["formula"]
            if index_parameters in formula:
                new = f"{name_dict}[str('{index_parameters}')]"
                dictionary[index_reactions]["formula"] = formula.replace(index_parameters, new)
    
    
    formulas_dict = {}
    for index_reactions in dictionary.keys():
            formula = dictionary[index_reactions]["formula"]
            formulas_dict[index_reactions] = formula
    return formulas_dict

formulas = create_formulas_eval(reactions, species, par_dict, name_dict = "par_dict", name_df = "Y")



#-----------------tau-leaping-------------------------------------------------

#create stoichiometric matrix


def create_stoichimetric_matrix(reactions, species, init=0):
    S = pd.DataFrame(init, columns=reactions.keys(), index = species)
    
    for index in reactions.keys():
        products = list(reactions[index]["products"].keys())
        reactants = list(reactions[index]["reactants"].keys())
        for index_prod in products:
            product_stoichiometric_coef = reactions[index]["products"][index_prod]
            S[index][index_prod] += product_stoichiometric_coef
        
        for  index_reactant in reactants:
            reactant_stoichiometric_coef = reactions[index]["reactants"][index_reactant]
            S[index][index_reactant] += -reactant_stoichiometric_coef
    
    return S


S = create_stoichimetric_matrix(reactions, species)
null_space = scipy.linalg.null_space(S)
tau = time_points[1] - time_points[0] 



def run_tau_leaping(seed, S, model_out, nu_parameter, time_points, formulas):
    np.random.seed(seed)
    Y = pd.DataFrame([model_out["initial_states"]], columns = list(species))
    for period in range(len(time_points)-1):
        nu_countryA_vac1 = float(eval(nu_parameter["nu_countryA_vac1"]))
        nu_countryA_vac2 = float(eval(nu_parameter["nu_countryA_vac2"]))
        nu_countryB_vac1 = float(eval(nu_parameter["nu_countryB_vac1"]))
        nu_countryB_vac2 = float(eval(nu_parameter["nu_countryB_vac2"]))
        
        propensities = np.repeat(np.nan, len(formulas.keys()))
        num = 0
        for index in formulas.keys():
            propensities[num] = eval(formulas[index]) * tau
            num += 1
            
        numb_reactions = np.random.poisson(propensities)
        
        Y_new = Y.iloc[period] + S @ numb_reactions
        Y.loc[period + 1] = Y_new
        
    return Y

def run_tau_leaping_parallel(seed):
    return run_tau_leaping(seed,S, model_out, nu_parameter, time_points, formulas)

#np.random.seed(123)
#results2 = {}
#for run in range(2):
#    results2[f"run{run}"] = run_tau_leaping(S, model_out, nu_parameter, time_points, formulas)


import multiprocessing as mp
processes = mp.cpu_count()  # Specify number of processes here 
p = mp.Pool(processes)
results = p.map(run_tau_leaping_parallel, range(100))
p.close()


from visualization.model_results import plot_states
from visualization.model_results import get_substates
substates_optimal_A = get_substates(
    model=model, substrings=["vac0", "infectious", "countryA"], include_all=True
)
ylim_infectious = [0, 2.5 * 10 ** 7]
colors = ["C2", "C3"]
fig, ax = plot_states(
    results=results[2],
    state_ids=substates_optimal_A,
    time = time_points,
    time_name="t",
    title="Unvaccinated infectious individuals (A) (Pareto optimal)",
    colors=["C2", "C3"],
    xlabel="Days",
    ylabel="Number of individuals",
    custom_label=["Virus wild type", "Virus mutant"],
    ylim=ylim_infectious,
    legend_next_to_plot=False,
)










