from functions.create_sbml import write_entities_to_dict
from functions.create_sbml import write_dict_to_sbml_file


def model_basic_sir_create_sbml(path, S_0, I_0, R_0=0):
    compartments = { 'A' : {"units" : 'dimensionless', 'constant' : True, 'volume' : 1}
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

    
    

    