from models.vaccination.create_model_vaccination import get_all_species_alive
from models.vaccination.create_model_vaccination import get_all_species_alive_area
from models.vaccination.create_model_vaccination import create_species_model
from visualization.model_results import get_substates


def create_rules_vaccination_proportion_relative_infected_population(
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
    vaccination_states_removed = [
        x for x in vaccination_states if x != non_vaccination_state
    ]
    species = create_species_model(
        species_comp=species_comp,
        areas=areas,
        vaccination_states=vaccination_states,
        non_vaccination_state=non_vaccination_state,
        virus_states=virus_states,
    )

    numb_total_infectious = get_all_infectious(species)

    rules_proportions_vaccination = {}
    for index_vaccination in vaccination_states_removed:
        for index_areas in areas:
            population_infectious_area = get_infectious_area(
                species=species, area=index_areas
            )
            id_rule = f"rule_proportion_{index_areas}_{index_vaccination}"
            id_proportion = f"proportion_{index_areas}_{index_vaccination}"

            formula = f"({population_infectious_area}) / ({numb_total_infectious})"
            rules_proportions_vaccination[id_rule] = {
                "parameter_id": id_proportion,
                "formula": formula,
            }

    return rules_proportions_vaccination


def get_all_infectious(species):
    """Get sum of all names of states that are infectiious.

    Parameters
    ----------
    species : dict of species
        List of all species.

    Returns
    -------
    numb_total_infectious : list of strings
        List of all infectious states.
    """

    species_all = species.keys()
    infectious_all = [x for x in species_all if "infectious" in x]
    numb_total_infectious = "+".join(infectious_all)
    return numb_total_infectious


def get_infectious_area(species, area):
    """Get sum of all names of states that are infectiious from a speciffic
    area.

    Parameters
    ----------
    species : dict of species
        List of all species.

    Returns
    -------
    numb_infectious : list of strings
        List of all infectious states from country `area`.
    """

    species_all = species.keys()
    infectious = [x for x in species_all if ("infectious" in x) and (area in x)]
    numb_infectious = "+".join(infectious)
    return numb_infectious


def create_rules_vaccination_proportion_relative_population(
    species_comp, vaccination_states, non_vaccination_state, virus_states, areas
):
    """Create rules for the vaccination proportions as relative population
    sizes.

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
