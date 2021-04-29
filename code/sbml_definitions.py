import libsbml
from libsbml import parseL3Formula
from libsbml import SBMLDocument
from libsbml import writeSBMLToFile


def create_sbml_file(compartments, species, parameters, reactions, assignments = None):
    

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

    if assignments is not None:
        assignments_identifier = assignments.keys()
        for keys in assignments_identifier:
            create_initial_assignment(model=model, assignment_id=keys, **assignments[keys])
        
    reactions_identifier = reactions.keys()
    for keys in reactions_identifier:
        create_reaction(model=model, reaction_id=keys, **reactions[keys])
    
    return document

        
def check_model_for_errors(document):
    for i in range(document.getNumErrors()):
        print(document.getError(i).getMessage())



def create_species(model: libsbml.Model, species_id, compartment,
                   constant = False, initial_amount = 0.0,
                   substance_units = 'Mole',
                   boundary_condition = False,
                   has_only_substance_units = False):
    """
    Creates new species for libsbml.Model
    :param model: libsbml.Model for which species will be created
    :param s_id: species id for the new species
    :param compartment: compartment to which new species will be assigned
    :param constant: True if the new species should be constant
    :param initial_amount: start value for species
    :param substance_units: substance unit for species
    :param boundary_condition: True if there is a boundary condition
    :param has_only_substance_units: True if species has only substance units
    """
    species = model.createSpecies()
    species.setId(species_id)
    species.setName(species_id)
    species.setCompartment(compartment)
    species.setConstant(constant)
    species.setInitialAmount(initial_amount)
    species.setSubstanceUnits(substance_units)
    species.setBoundaryCondition(boundary_condition)
    species.setHasOnlySubstanceUnits(has_only_substance_units)

    return species


def create_parameter(model: libsbml.Model, parameter_id, value, 
                     constant=True, units="dimensionless"):
    """
    Creates new parameter for libsbml.Model
    :param model: libsbml.Model to which new parameter will be added
    :param p_id: Id for the new parameter
    :param constant: True if the new parameter should be constant
    :param value: Initial value of the parameter
    :param units: Units for the parameter
    """
    k = model.createParameter()
    k.setId(parameter_id)
    k.setName(parameter_id)
    k.setConstant(constant)
    k.setValue(value)
    k.setUnits(units)

    return k



def create_reaction(model: libsbml.Model, reaction_id: str,
                    reactants,
                    products,
                    formula: str, modifiers = None,
                    reversible: bool = False, fast: bool = False):
    """
    Creates new reaction for libsbml.Model
    :param model: libsbml.Model to which new reaction will be added
    :param m_id: Id for the new reaction
    :param reactants: List or reactants of the reaction
    :param products: List of products of the reaction
    :param formula: Formula for the rate of the reaction
    :param modifiers: List of modifiers of the reaction. Used but no reactant
                        or product
    :param reversible: True if the reaction should be reversible
    :param fast: True if the reaction should be fast
    """
    if modifiers is None:
        modifiers = []
    reaction = model.createReaction()
    reaction.setId(reaction_id)
    reaction.setReversible(reversible)
    reaction.setFast(fast)


    for keys in reactants.keys():
        species_ref = reaction.createReactant()
        species_ref.setSpecies(keys)
        species_ref.setConstant(True)  # TODO ?
        species_ref.setStoichiometry(reactants[keys])

    for keys in products.keys():
        species_ref = reaction.createProduct()
        species_ref.setSpecies(keys)
        species_ref.setConstant(True)  # TODO ?
        species_ref.setStoichiometry(products[keys])

    for name in modifiers:
        species_ref = reaction.createModifier()
        species_ref.setSpecies(name)
        # species_ref.setConstant(True)  # TODO ?

    math_ast = parseL3Formula(formula)
    kinetic_law = reaction.createKineticLaw()
    kinetic_law.setMath(math_ast)

    return reaction

def create_compartment(model: libsbml.Model, compartment_id, units="dimensionless",
                       constant = False, volume = 1):

    compartment = model.createCompartment()
    compartment.setName(compartment_id)
    compartment.setId(compartment_id)
    compartment.setUnits(units)
    compartment.setConstant(constant)
    compartment.setVolume(volume)
    
    return

def create_initial_assignment(model: libsbml.Model,
                             species_id,
                             formula,
                             assignment_id,
                             assignment_name = None):
    """Create SBML InitialAssignment
    Arguments:
        sbml_model: Model to add output to
        species_id: Target of assignment
        formula: Formula string for model output
        rule_id: SBML id for created rule
        rule_name: SBML name for created rule
    Returns:
        The created ``InitialAssignment``
    """

    if assignment_name is None:
        assignment_name = assignment_id

    assignment = model.createInitialAssignment()
    assignment.setId(assignment_id)
    assignment.setName(assignment_name)
    assignment.setSymbol(species_id)
    math_ast = libsbml.parseL3Formula(formula)
    assignment.setMath(math_ast)

    return assignment

def write_entities_to_dict(compartments, species, parameters,
                           reactions, assignments=None):
    if assignments is None:
        output = {'compartments' : compartments,
                  'species' : species,
                  'parameters' : parameters,
                  'reactions' : reactions}
    else:
        output = {'compartments' : compartments,
                  'species' : species,
                  'parameters' : parameters,
                  'reactions' : reactions,
                  'assignments' : assignments}
        
    return output


def write_dict_to_sbml_file(dictionary, path):
    
    doc = create_sbml_file(**dictionary)    
    
    output_path = path + '.xml'
    writeSBMLToFile(doc, output_path)