import numpy as np
import pandas as pd
from functions.create_sbml import write_entities_to_dict
from functions.create_sbml import write_dict_to_sbml_file


def model_one_country_create_sbml(path, S_0, I_0, R_0=0):
    compartments = { 'countryA' : {"units" : 'dimensionless', 'constant' : True, 'volume' : 1}
    }
    
    
    species = { 'susceptible' : {'compartment' : 'A', 'initial_amount' : S_0},
                'infectious'  : {'compartment' : 'A', 'initial_amount' : I_0},
                'recovered'   : {'compartment' : 'A', 'initial_amount' : R_0}
               }
    
    parameters = { 'beta' : {'value' : 1 },
                   'gamma' : {'value' : 0.13 }
                 }
    
    infection_reaction =  {'reactants' : {'susceptible' : 1, 'infectious' : 1},
                           'products'  : {'infectious' : 2},
                           'modifiers' : {'recovered'},
                           'formula' : "beta * susceptible * infectious / (susceptible + infectious + recovered)"
                           } 
    
    recover_reaction = {'reactants' : {'infectious' : 1},
                           'products'  : {'recovered' : 1},
                           'formula' : "gamma * infectious"
                           } 
    
    reactions = { 'infection' : infection_reaction,
                   'recover'  : recover_reaction} 
    
    model_input_dictionary = write_entities_to_dict(compartments=compartments, species=species, parameters=parameters,
                           reactions=reactions)
    
    write_dict_to_sbml_file(dictionary=model_input_dictionary, path=path)
    


#needed to create species
vaccination_states = ['vac0', 'vac1', 'vac2']
virus_states = ['virW', 'virM']
areas = ['countryA']
species_comp = ['susceptible', 'infectious', 'recovered', 'dead']

#needed to define parameetrs
omega_matrix = np.array([[0.9, 0.6], [0.7, 0.9]])
delta_matrix = np.array([[0.9, 0.6], [0.6, 0.9]])
eta_vector = np.array([[1, 1.3]]) 

vaccination_states_removed = [x for x in vaccination_states if x != 'vac0']
omega_df = pd.DataFrame(omega_matrix, columns=vaccination_states_removed, index=virus_states)
delta_df = pd.DataFrame(delta_matrix, columns=vaccination_states_removed, index=virus_states)
eta_df = pd.DataFrame(eta_vector, columns=virus_states)    

#create species
species = {}
for index_compartments in species_comp:
    for index_areas in areas:
        for index_vaccination in vaccination_states:
            for index_virus in virus_states:
                if index_compartments == 'susceptible':
                    species_combined = f'{index_compartments}_{index_areas}_{index_vaccination}'
                    amount_t0 = 100000
                else:
                    species_combined = f'{index_compartments}_{index_areas}_{index_vaccination}_{index_virus}'
                    
                    if index_compartments in ['recovered', 'dead']:
                        amount_t0 = 0
                    elif index_compartments == 'infectious':
                        amount_t0 = 10**(-13)
                if index_vaccination != 'vac0':
                    amount_t0 = 0
                species[species_combined] = {'compartment' : index_areas, 'initial_amount' : amount_t0}         
                    
#create parameters       
omega_delta_dict = {}
eta_dict = {}
for index_virus in virus_states:
    for index_vaccination in vaccination_states_removed:

        key_omega = f'omega_{index_vaccination}_{index_virus}'
        key_delta = f'delta_{index_vaccination}_{index_virus}'
        omega_delta_dict[key_omega] = {'value' : omega_df[index_vaccination][index_virus]}
        omega_delta_dict[key_delta] = {'value' : delta_df[index_vaccination][index_virus]}
    
    key_eta = f'eta_{index_virus}'
    eta_dict[key_eta] = {'value' : eta_df[index_virus][0]}



parameters_fixed = { 'lambda' :   {'value' : 0.01 },
               'p_lambda' : {'value' : 0.01 },
               'gamma' : {'value' : 0.5}, 
               'beta' : {'value' : 2},
               'nu_vac1' : {'value' : 0.01},
               'nu_vac2' : {'value' : 0.01},
             }

parameters = {**parameters_fixed, **omega_delta_dict, **eta_dict}

#create reactions

#Dead/Recovery
dead_recover_reactions = {}
for index_areas in areas:
    for index_vaccination in vaccination_states:
        for index_virus in virus_states:
            numb_infected = f'infectious_{index_areas}_{index_vaccination}_{index_virus}'
            numb_recovered = f'recovered_{index_areas}_{index_vaccination}_{index_virus}'
            numb_dead = f'dead_{index_areas}_{index_vaccination}_{index_virus}'
            key_recover = f'{numb_infected}_to_{numb_recovered}'
            key_death = f'{numb_infected}_to_{numb_dead}'
            
            if vaccination_states == 'vac0':
                omega_term = '0'
            else:
                omega_term = f'omega_{index_vaccination}_{index_virus}'
                
            dead_recover_reactions[key_recover] = {'reactants' : {f'{numb_infected}' : 1},
                                      'products'  : {f'{numb_recovered}' : 1},
                                      'formula' : f"(1 - (1 - {omega_term}) * p) * lambda * {numb_infected}"
                                      } 
            dead_recover_reactions[key_death]= {'reactants' : {f'{numb_infected}' : 1},
                                    'products'  : {f'{numb_dead}' : 1},
                                    'formula' : f" (1 - {omega_term}) * p * lambda * {numb_infected}"
                                    }
    
  

#infection
total_population_alive = [ x for x in species.keys() if ("dead" not in x) ]
numb_total_population_alive = '+'.join(total_population_alive)
for index_areas_infected in areas:
    for index_vaccination_infected in vaccination_states:
        for index_virus_infected in virus_states:
            for index_areas_susceptible in areas:
                for index_vaccination_susceptible in vaccination_states:
                    
                    numb_susceptible = f'susceptible__{index_areas}_{index_vaccination}'
                    numb_infectious = f'infectious_{index_areas}_{index_vaccination}_{index_virus}'
                    susceptible_area_population_alive = [ x for x in species.keys() if ((f"{index_areas_susceptible}" in x) and ("dead" not in x)) ]
                    numb_susceptible_area_population_alive = '+'.join(susceptible_area_population_alive)
                    infected_area_population_alive = [ x for x in species.keys() if ((f"{index_areas_infected}" in x) and ("dead" not in x)) ]
                    numb_infected_area_population_alive = '+'.join(infected_area_population_alive)
                    
                    eta_term = f'eta_{index_virus_infected}'
                    
                    if index_vaccination_susceptible != 'vac0':
                        delta_term = f'delta_{index_vaccination_susceptible}_{index_virus_infected}'
                    else:
                        delta_term = '0'
                        
                    if index_vaccination_infected != 'vac0':
                        gamma_term = 'gamma'
                    else:
                        gamma_term = '0'

#TODO: give key proper name, define mll (define b matrix first and with that the M matrix) and define modifiers (list of all used) 
                    infection_reaction[key] =  {'reactants' : {f'{numb_susceptible}' : 1, '{numb_infectious}' : 1},
                                                'products'  : {'{numb_infectious}' : 2},
                                                'modifiers' : {'recovered'},
                                                'formula' : "beta * {eta_Term} * (1-{delta_term}) * (1-{gamma_term})*{numb_susceptible}*{numb_infectious}*mll*/{numb_susceptible_area_population_alive}"
                                               } 
        
                    
                    
                    
                    
                    
                    
                    
                    
                    