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
    



vaccination_states = ['vac0', 'vac1', 'vac2']
virus_states = ['virW', 'virM']
areas = ['countryA']
species_comp = ['sus', 'inf', 'rec', 'dead']
species = {}

for index_compartments in species_comp:
    for index_areas in areas:
        for index_vaccination in vaccination_states:
            for index_virus in virus_states:
                if index_compartments == 'sus':
                    species_combined = f'{index_compartments}_{index_areas}_{index_vaccination}'
                else:
                    species_combined = f'{index_compartments}_{index_areas}_{index_vaccination}_{index_virus}'
                    
                    if index_compartments in ['rec', 'dead']:
                        amount_t0 = 0
                    elif index_compartments == 'inf':
                        amount_t0 = 10^-13
                    else: #susceptible case
                        if vaccination_states == 'vac0':
                            amount_t0 = 100000
                species[species_combined] = {'compartment' : index_areas, 'initial_amount' : amount_t0}         
                    
                    
#next: parameters and reactions                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    