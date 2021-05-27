import libsbml
from libsbml import parseL3Formula
from libsbml import SBMLDocument
from libsbml import writeSBMLToFile

from amici.splines import CubicHermiteSpline



def create_sbml_file(
    compartments,
    species,
    parameters,
    reactions,
    assignments=None,
    parameter_rules=None,
    rate_rules=None,
    splines=None,
):
    """Create SBML object for given inputs. Function does not save the file.
    Saving is done by :func:`write_dict_to_sbml_file`.

    Parameters
    ----------
    compartment : dict
        Dictionary containing the information for the compartments. Keys are
        names of compartments. Values can be dictionaries as well specifying
        all other arguments passed for the compartment.

    species : dict
        Dictionary containing the information for the species. Keys are
        names of species. Values can be dictionaries as well specifying
        all other arguments passed for the species.

    parameters : dict
        Dictionary containing the information for the parameters. Keys are
        names of parameters. Values can be dictionaries as well specifying
        all other arguments passed for the parameters.

    reactions : dict
        Dictionary containing the information for the reactions. Keys are
        names ofreactions. Values can be dictionaries as well specifying
        all other arguments passed for the reactions.

    assignments : dict
        Dictionary containing the information for the assignments. Keys are
        names of assignments. Values can be dictionaries as well specifying
        all other arguments passed for the assignments.

    parameter_rules : dict
        Dictionary containing the information for the parameter rules. Keys are
        names of rules. Values can be dictionaries as well specifying
        all other arguments passed for the rules.

    Returns
    -------
    output : dict
        Dictionary containing all information for compartments, species,
        parameters, reactions, assignments and parameter rules that can be
        passed to :func:`write_dict_to_sbml_file`.
    """

    document = SBMLDocument()
    model = document.createModel()

    compartment_identifier = compartments.keys()
    for keys in compartment_identifier:
        create_compartment(model=model, compartment_id=keys, **compartments[keys])

    species_identifier = species.keys()
    for keys in species_identifier:
        create_species(model=model, species_id=keys, **species[keys])

    parameters_identifier = parameters.keys()
    for keys in parameters_identifier:
        create_parameter(model=model, parameter_id=keys, **parameters[keys])

    if assignments is not None:
        assignments_identifier = assignments.keys()
        for keys in assignments_identifier:
            create_initial_assignment(
                model=model, assignment_id=keys, **assignments[keys]
            )

    reactions_identifier = reactions.keys()
    for keys in reactions_identifier:
        create_reaction(model=model, reaction_id=keys, **reactions[keys])

    if parameter_rules is not None:
        parameter_rule_identifier = parameter_rules.keys()
        for keys in parameter_rule_identifier:
            create_parameter_rule(model=model, rule_id=keys, **parameter_rules[keys])

   
    if splines is not None:
        spline_identifier = splines.keys()
        for keys in spline_identifier:
            create_cubic_hermite_spline(model=model, parameter_name=keys, **splines[keys])

    if rate_rules is not None:
        rate_rule_identifier = rate_rules.keys()
        for keys in rate_rule_identifier:
            create_parameter_rate_rule(model=model, rule_id=keys, **rate_rules[keys])
        
    return document

def create_cubic_hermite_spline(model, parameter_name, time_symbol, xx, yy, xx_names,
                                yy_names):
    spline = CubicHermiteSpline(
            sbmlId = parameter_name,
            x = time_symbol,
            xx = xx, #make to xx_names if absolute value problem is solved
            yy = yy_names,
        )
    
    spline.addToSbmlModel(model=model, auto_add = True,
        x_nominal = xx,
        y_nominal = yy,
        #y_units: Optional[str] = None,
        #x_units: Optional[str] = None)
        )
    
    return spline

def write_dict_to_sbml_file(dictionary, path, error_check=False):
    """
    Creates SBML document from input dictionary and saves it on the path.

    Parameters
    ----------
    dictionary: dict
        Dictionary containing all information for compartments, species,
        parameters, reactions, assignments and parameter rules.

    path : str
        Path at which SBML model should be stored.

    error_check : {True, False}
        If True the function :func:`check_model_for_errors` is run within the
        function execution.

    """

    doc = create_sbml_file(**dictionary)
    if error_check is True:
        check_model_for_errors(doc)

    output_path = path + ".xml"
    writeSBMLToFile(doc, output_path)


def write_entities_to_dict(
    compartments,
    species,
    parameters,
    reactions,
    assignments=None,
    parameter_rules=None,
    rate_rules=None,
    splines=None,
):
    """
    Transforms inputs to dictionaries with the entities as keys and the
    objects as containing all the information as values.

    Parameters
    ----------
    compartment : dict
        Dictionary containing the information for the compartments. Keys are
        names of compartments. Values can be dictionaries as well specifying
        all other arguments passed for the compartment.

    species : dict
        Dictionary containing the information for the species. Keys are
        names of species. Values can be dictionaries as well specifying
        all other arguments passed for the species.

    parameters : dict
        Dictionary containing the information for the parameters. Keys are
        names of parameters. Values can be dictionaries as well specifying
        all other arguments passed for the parameters.

    reactions : dict
        Dictionary containing the information for the reactions. Keys are
        names ofreactions. Values can be dictionaries as well specifying
        all other arguments passed for the reactions.

    assignments : dict
        Dictionary containing the information for the assignments. Keys are
        names of assignments. Values can be dictionaries as well specifying
        all other arguments passed for the assignments.

    parameter_rules : dict
        Dictionary containing the information for the parameter rules. Keys are
        names of rules. Values can be dictionaries as well specifying
        all other arguments passed for the rules.

    Returns
    -------
    output : dict
        Dictionary containing all information for compartments, species,
        parameters, reactions, assignments and parameter rules that can be
        passed to :func:`write_dict_to_sbml_file`.
    """

    output = {
        "compartments": compartments,
        "species": species,
        "parameters": parameters,
        "reactions": reactions,
    }

    if assignments is not None:
        output["assignments"] = assignments

    if parameter_rules is not None:
        output["parameter_rules"] = parameter_rules

    if rate_rules is not None:
        output["rate_rules"] = rate_rules
    
    if splines is not None:
        output["splines"] = splines

    return output


def check_model_for_errors(document):
    """Checks SBML model for errors.

    Parameters
    ----------
    document : sbml.Document
       Document for which potential errors should be checked.

    Returns
    -------
    None.

    """

    for i in range(document.getNumErrors()):
        print(document.getError(i).getMessage())


def create_species(
    model,
    species_id,
    compartment,
    constant=False,
    initial_amount=0.0,
    substance_units="dimensionless",
    boundary_condition=False,
    has_only_substance_units=False,
):
    """
    Creates new species for libsbml.Model. It is not necessary to assign the
    function to a value. The species is automatically added to the input model
    outside of the function.

    Parameters
    ----------
    model : libsbml.Model
        Model for which species will be created.

    species_id : str
        Species id for the new species.

    compartment : str
        Compartment to which new species will be assigned.

    constant: {True, False}
        True if the new parameter can only be constant. And False if not.
        False does not mean that the parameter has to change but rather that
        it is allowed to.

    initial_amount : float
        Initial value of the species at time :math:`t=0`.

    substance_units : str
        Substance unit for the species.

    boundary_condition : {True, False}
        True if there is a boundary condition.

    has_only_substance_units : {True, False}
        If True species is only allwoed to have substance units.

    Returns
    -------
    species : libsbml.Species
        Species defined for the model.

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


def create_parameter(model, parameter_id, value, constant=True, units="dimensionless"):
    """
    Creates new parameter for libsbml.Model. It is not necessary to assign the
    function to a value. The parameter is automatically added to the input model
    outside of the function.

    Parameters
    ----------
    model : libsbml.Model
        Model for which species will be created.

    parameter_id : str
        Id for the new parameter.

    constant : {True, False}
        True if the new parameter can only be constant. And False if not.
        False does not mean that the parameter has to change but rather that
        it is allowed to.

    value : float
        Value of the parameter.

    units : str
        Units of the parameter.

    Returns
    -------
    parameter : libsbml.Parameter
        Parameter defined for the model.
    """

    parameter = model.createParameter()
    parameter.setId(parameter_id)
    parameter.setName(parameter_id)
    parameter.setConstant(constant)
    parameter.setValue(value)
    parameter.setUnits(units)

    return parameter


def create_reaction(
    model,
    reaction_id: str,
    reactants,
    products,
    formula: str,
    modifiers=None,
    reversible: bool = False,
    fast: bool = False,
):

    """
    Creates new reaction for libsbml.Model. It is not necessary to assign the
    function to a value. The reaction is automatically added to the input model
    outside of the function.

    Parameters
    ----------
    model : libsbml.Model
        Model for which species will be created.

    reaction_id : str
        Reaction id for the new species.

    reactants : dict
        Dictionary of reactants with name of reactants as keys and stoichometric
        coefficients as values.

    products : dict
        Dictionary of products with name of products as keys and stoichometric
        coefficients as values.

    formula: str
        Formula for the rate of the reaction.

    modifiers : list of strings
        List of modifiers of the reaction. Modifiers are not explicitly
        reactants or products but appear in the formula.

    param_reversible : {True, False}
        True if the reaction is reversible.

    fast: {True, False}
        True if the reaction is fast.

    Returns
    -------
    reaction: libsbml.Reaction
        Reaction defined for the model.

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
        species_ref.setConstant(True)
        species_ref.setStoichiometry(reactants[keys])

    for keys in products.keys():
        species_ref = reaction.createProduct()
        species_ref.setSpecies(keys)
        species_ref.setConstant(True)
        species_ref.setStoichiometry(products[keys])

    for name in modifiers:
        species_ref = reaction.createModifier()
        species_ref.setSpecies(name)

    math_ast = parseL3Formula(formula)
    kinetic_law = reaction.createKineticLaw()
    kinetic_law.setMath(math_ast)

    return reaction


def create_compartment(
    model,
    compartment_id,
    units="dimensionless",
    constant=False,
    volume=1,
):
    """
    Creates new comaprtment for libsbml.Model. It is not necessary to assign the
    function to a value. The compartment is automatically added to the input model
    outside of the function. Note that a compartment in the language of sbml
    is an area in our model.

    Parameters
    ----------
    model : libsbml.Model
        Model for which species will be created.

    compartment_id : str
        Compartment id for the new species.

    units : str
        Units of the entities inside the compartment.

    constant : {True, False}
        True if the new compartment can only be constant. And False if not.
        False does not mean that the compartment has to change but rather that
        it is allowed to.


    volume : {1,2,3}
        Integer describing the type of compartment volume.

    Returns
    -------
    reaction: libsbml.Compartment
        Compartment defined for the model.
    """

    compartment = model.createCompartment()
    compartment.setName(compartment_id)
    compartment.setId(compartment_id)
    compartment.setUnits(units)
    compartment.setConstant(constant)
    compartment.setVolume(volume)

    return


def create_initial_assignment(
    model, species_id, formula, assignment_id, assignment_name=None
):
    """
    Creates new initial assignment for libsbml.Model. It is not necessary
    to assign the function to a value. The uinitial assignment is
    automatically added to the input model outside of the function.
    Initial assignments overwrite initial amounts of species but can be
    changed after the model has been created by
    :func:`get_model_and_solver_from_sbml` and should thus be prefered
    over initial amounts.

    Parameters
    ----------
    model : libsbml.Model
        Model for which species will be created.

    species_id : str
        Species for which the initial assignment is created.

    formula : str
        Formula that describes the type of initial assignment.

    assignment_id: str
        Assignment id for created initial asignment.

    assignment_name: str
        SBML name for created initial asignment.

    Returns
    -------
    assignment: libsbml.InitialAssignment
        Initial assignment defined for the species.
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


def create_parameter_rule(model, parameter_id, formula, rule_id, rule_name=None):
    """
    Creates new rule for libsbml.Model. It is not necessary
    to assign the function to a value. The rule is automatically
    added to the input model outside of the function.

    Parameters
    ----------
    model : libsbml.Model
        Model for which species will be created.

    parameter_id : str
        Parameter for which the rule is created.

    formula : str
        Formula that describes the type of rule.

    rule_id: str
        Rule id for created rule.

    rule_name: str
        SBML name for created rule.

    Returns
    -------
    parameter_rule : libsbml.AssignmentRule
        Rule defined for the parameter.
    """

    if rule_name is None:
        rule_name = rule_id
    parameter_rule = model.createAssignmentRule()
    parameter_rule.setId(rule_id)
    parameter_rule.setName(rule_name)
    parameter_rule.setVariable(parameter_id)
    math_ast = libsbml.parseL3Formula(formula)
    parameter_rule.setMath(math_ast)

    return parameter_rule


def create_parameter_rate_rule(model, parameter_id, formula, rule_id, rule_name=None):
    """
    Creates new rule for libsbml.Model. It is not necessary
    to assign the function to a value. The rule is automatically
    added to the input model outside of the function.

    Parameters
    ----------
    model : libsbml.Model
        Model for which species will be created.

    parameter_id : str
        Parameter for which the rule is created.

    formula : str
        Formula that describes the type of rule.

    rule_id: str
        Rule id for created rule.

    rule_name: str
        SBML name for created rule.

    Returns
    -------
    parameter_rule : libsbml.AssignmentRule
        Rule defined for the parameter.
    """

    if rule_name is None:
        rule_name = rule_id
    parameter_rule = model.createRateRule()
    parameter_rule.setId(rule_id)
    parameter_rule.setName(rule_name)
    parameter_rule.setVariable(parameter_id)
    math_ast = libsbml.parseL3Formula(formula)
    parameter_rule.setMath(math_ast)

    return parameter_rule
