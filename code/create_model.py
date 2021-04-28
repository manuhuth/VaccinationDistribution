from libsbml import SBMLDocument
from libsbml import writeSBMLToFile
from sbml_definitions import create_species
from sbml_definitions import create_parameter
from sbml_definitions import create_reaction
from sbml_definitions import create_compartment


compartments = { 'A' : {"units" : 'dimensionless', 'constant' : True, 'volume' : 1}, 
                 'B' : {"units" : 'dimensionless', 'constant' : True, 'volume' : 1} 
}

species = { 'susceptible' : {'compartment' : 'A'},
            'infectious'  : {'compartment' : 'A'},
            'recovered'   : {'compartment' : 'A'}
           }

parameters = { 'beta' : {'value' : 1 },
               'gamma' : {'value' : 0.13 }
             }

infection_reaction =  {'reactants' : {'susceptible' : 1, 'infectious' : 1},
                       'products'  : {'infectious' : 2},
                       'formula' : "beta * susceptible * infectious / (susceptible + infectious + recovered)"
                       } 

recover_reaction = {'reactants' : {'infectious' : 1},
                       'products'  : {'recovered' : 1},
                       'formula' : "gamma * infectious"
                       } 

reactions = { 'infection' : infection_reaction,
               'recover'  : recover_reaction} 

def create_sbml_file(compartments, species, parameters, reactions):
    

    ## Create new SBML
    document = SBMLDocument()
    model = document.createModel()
    
    
    compartment_identifier = compartments.keys()
    for keys in compartment_identifier:
        create_compartment(model = model, compartment_id = keys, **compartments[keys])
    
    
    species_identifier = species.keys()
    for keys in species_identifier:
        create_species(model=model, species_id=keys, **species[keys])
        
    parameters_identifier = parameters.keys()
    for keys in parameters_identifier:
        create_parameter(model=model, parameter_id=keys, **parameters[keys])
        
    reactions_identifier = reactions.keys()
    for keys in reactions_identifier:
        create_reaction(model=model, reaction_id=keys, **reactions[keys])
    
    return document
    
    
        
def check_model_for_errors(document):
    for i in range(document.getNumErrors()):
        print(document.getError(i).getMessage())
    

doc = create_sbml_file(compartments=compartments, species=species, parameters=parameters, reactions=reactions)    
writeSBMLToFile(doc, 'test.xml')
    
    
    
    
    
    

    