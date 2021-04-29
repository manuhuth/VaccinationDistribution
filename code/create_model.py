from sbml_definitions import create_sbml_file
from sbml_definitions import write_entities_to_dict
from sbml_definitions import write_dict_to_sbml_file


def model_basic_sir_create(path, S_0, I_0, R_0=0):
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
    
    model_input_dictionary = write_entities_to_dict(compartments, species, parameters,
                           reactions)
    
    write_dict_to_sbml_file(dictionary=model_input_dictionary, path=path)


def model_basic_sir_run():
    
    
    
    sbml_importer = amici.SbmlImporter('test.xml')

    model_name = 'test'
    model_dir = 'model_dir'
    sbml_importer.sbml2amici(model_name, model_dir)
    
    # load the model module
    model_module = amici.import_model_module(model_name, model_dir)
    # instantiate model
    model = model_module.getModel()
    # instantiate solver
    solver = model.getSolver()
    
    
    model.setParameterByName('beta' , 0.4)
    model.setParameterByName('gamma', 0.13)
    
    
    # set timepoints
    model.setTimepoints(np.linspace(0, 100, 100))
    rdata = amici.runAmiciSimulation(model, solver)
    
    
    
    
    
    

    