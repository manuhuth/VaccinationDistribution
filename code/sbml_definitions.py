import libsbml
from typing import List, Tuple
from libsbml import parseL3Formula


def create_species(model: libsbml.Model, species_id, compartment,
                   constant = False, initial_amount = 0.0,
                   substance_units = 'nanoMole',
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
                    formula: str, modifiers: List[str] = None,
                    reversible: bool = False, fast: bool = False) -> None:
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