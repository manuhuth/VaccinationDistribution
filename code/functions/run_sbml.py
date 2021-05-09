import amici

from models.vaccination.create_model_vaccination import get_all_species_alive
from models.vaccination.create_model_vaccination import get_all_species_alive_area
from models.vaccination.create_model_vaccination import create_species_model


def create_rules_vaccination_proportion_relative_population(
    species_comp, vaccination_states, non_vaccination_state, virus_states, areas
):
    """Create rules for the vaccination proportions as relative country sizes.

    Parameters
    ----------
    species_comp : list of strings
        List containing the names of the compartments.

    vaccination_states : list of strings
        List containing the names of the vaccinations containing a state for
        non-vaccinated individuals.

    non_vaccination_state : str
        Name of state indicates non-vaccinated individuals.

    virus_states : list of strings
        List containing the names of the virus types.

    areas : list of strings
        List containing the names of the areas.

    Returns
    -------
    rules_proportions_vaccination : dict
        Dictionary including all rules for the proportion parameters.
    """
    vaccination_states_removed = vaccination_states_removed = [
        x for x in vaccination_states if x != non_vaccination_state
    ]
    species = create_species_model(
        species_comp=species_comp,
        areas=areas,
        vaccination_states=vaccination_states,
        non_vaccination_state=non_vaccination_state,
        virus_states=virus_states,
    )

    population_alive_all = get_all_species_alive(species_list=species)

    rules_proportions_vaccination = {}
    for index_vaccination in vaccination_states_removed:
        for index_areas in areas:
            population_alive_area = get_all_species_alive_area(
                species_list=species, area=index_areas
            )
            id_rule = f"rule_proportion_{index_areas}_{index_vaccination}"
            id_proportion = f"proportion_{index_areas}_{index_vaccination}"

            formula = f"({population_alive_area}) / ({population_alive_all})"
            rules_proportions_vaccination[id_rule] = {
                "parameter_id": id_proportion,
                "formula": formula,
            }

    return rules_proportions_vaccination


def create_observables_vaccination_rates(
    vaccination_states_removed, areas, name_parameter="nu"
):
    """Create dictionary with observables that is used by
    :func:`get_model_and_solver_from_sbml`.

    Parameters
    ----------
    vaccination_states : list of strings
        List containing the names of the vaccinations containing a state for
        non-vaccinated individuals.

    areas : list of strings
        List containing the names of the areas.

    name_parameter : str
        String that is used as name for the vaccination parameters.

    Returns
    -------
    observables : dict
        dictionary with observables that is used by
        :func:`get_model_and_solver_from_sbml`.
    """
    observables = {}
    for index_vaccinations in vaccination_states_removed:
        for index_areas in areas:
            observable_id = (
                f"observable_{name_parameter}_{index_areas}_{index_vaccinations}"
            )
            formula = f"{name_parameter}_{index_areas}_{index_vaccinations}"

            observables[observable_id] = {"name": formula, "formula": formula}

    return observables


def get_model_and_solver_from_sbml(
    path_sbml, model_name, model_directory, observables=None
):
    """Get model and solver objects from a SBML file.

     Parameters
     ----------
     path_sbml : str
         Path to sbml file that contains the model.

     model_name : str
         Name of the model.

     model_directory : str
         Path and name of where model directory should be stored.

     Returns
     -------
    model_solver : dict
        Dictionary containing model and solver as keys and as values
        the respective objects from the sbml file.

    """

    filename = path_sbml + ".xml"
    sbml_importer = amici.SbmlImporter(filename)
    sbml_importer.sbml2amici(model_name, model_directory, observables=observables)

    model_module = amici.import_model_module(model_name, model_directory)
    model = model_module.getModel()
    solver = model.getSolver()

    model_solver = {"model": model, "solver": solver}

    return model_solver


def model_run(model, solver, timepoints, set_parameter=None):
    """Run model and return outputs for specification.

    Parameters
    ----------
    model : amici.Model
        Model to run.

    solver : amici.Solver
        Solver for amici model.

    timepoints : array
        Defining the grid for which model should be computed.

    set_parameter : dict
        Allows to change parameter values by passing a dictionary. Keys must
        be the names of the parameters and the values the magnitudes of the
        parameters.

    Returns
    -------
    rdata : numpy.ReturnDataView
        Output data.

    """

    if set_parameter is not None:
        for keys in set_parameter.keys():
            model.setParameterByName(keys, set_parameter[keys])

    model.setTimepoints(timepoints)
    rdata = amici.runAmiciSimulation(model, solver)

    return rdata
