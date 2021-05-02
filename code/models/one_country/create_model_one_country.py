import numpy as np
import pandas as pd
from functions.create_sbml import write_entities_to_dict
from functions.create_sbml import write_dict_to_sbml_file


def model_one_country_create_sbml(
    path,
    vaccination_states=["vac0", "vac1", "vac2"],
    non_vaccination_state = 'vac0',
    virus_states=["virW", "virM"],
    areas=["countryA"],
    species_comp=["susceptible", "infectious", "recovered", "dead"],
    vaccinated_compartments = ["susceptible", "recovered"],
    omega_matrix=np.array([[0.9, 0.6], [0.7, 0.9]]),
    delta_matrix=np.array([[0.9, 0.6], [0.6, 0.9]]),
    eta_vector=np.array([[1, 1.3]]),
    t0_susceptible = 10000,
    t0_infectious = 1,
    distances =  np.array([[0]]),
    check_error=False,
):

    vaccination_states_removed = [x for x in vaccination_states if x != non_vaccination_state]
    omega_df = pd.DataFrame(
        omega_matrix, columns=vaccination_states_removed, index=virus_states
    )
    delta_df = pd.DataFrame(
        delta_matrix, columns=vaccination_states_removed, index=virus_states
    )
    eta_df = pd.DataFrame(eta_vector, columns=virus_states)

    # create compartments
    compartments = create_comaprtments_model(areas)

    # create species
    species = create_species_model(
        species_comp=species_comp,
        areas=areas,
        vaccination_states=vaccination_states, non_vaccination_state=non_vaccination_state,
        virus_states=virus_states, t0_susceptible=t0_susceptible, t0_infectious=t0_infectious,
    )

    # create parameters
    parameters = create_parameters_model(
        virus_states=virus_states,
        vaccination_states_removed=vaccination_states_removed,
        omega_df=omega_df,
        delta_df=delta_df,
        eta_df=eta_df,
    )

    # create reactions
    reactions = create_reactions_model(
        areas=areas,
        vaccination_states=vaccination_states,
        vaccination_states_removed=vaccination_states_removed,
        virus_states=virus_states,
        species=species,
        distances=distances,
        vaccinated_compartments=vaccinated_compartments, non_vaccination_state=non_vaccination_state,
    )

    model_input_dictionary = write_entities_to_dict(
        compartments=compartments,
        species=species,
        parameters=parameters,
        reactions=reactions,
    )

    write_dict_to_sbml_file(dictionary=model_input_dictionary, path=path, error_check=check_error)


def b_distance_function(x):
    return 1 / (1 + x)


def get_all_species_alive(species_list):
    total_population_alive = [x for x in species_list.keys() if ("dead" not in x)]
    numb_total_population_alive = "+".join(total_population_alive)
    return numb_total_population_alive


def get_all_species_alive_area(species_list, area):
    total_population_alive = [
        x for x in species_list.keys() if ("dead" not in x) and (area in x)
    ]
    numb_total_population_alive = "+".join(total_population_alive)
    return numb_total_population_alive


def get_M_matrix(areas, species, distance_matrix):
    distances_df = pd.DataFrame(distance_matrix, columns=areas, index=areas)
    b_distance_matrix = b_distance_function(distances_df)
    numb_total_population_alive = get_all_species_alive(species_list=species)
    M = pd.DataFrame(np.nan, index=areas, columns=areas)
    for index_areas_row in areas:
        for index_areas_col in areas:
            if index_areas_row == index_areas_col:
                M.loc[index_areas_row, index_areas_col] = ""
            else:
                M.loc[
                    index_areas_row, index_areas_col
                ] = f"({get_all_species_alive_area(species_list=species, area=index_areas_row )})/({numb_total_population_alive })*{b_distance_matrix[index_areas_col][index_areas_row]}"
    if len(areas) == 1:
        add_term = "1"
    else:
        add_term = "1-"
    diagonal = add_term + M.sum(0)
    np.fill_diagonal(M.values, diagonal)

    return M


def create_comaprtments_model(areas):
    compartments = {}
    for index_area in areas:
        compartments[index_area] = {
            "units": "dimensionless",
            "constant": True,
            "volume": 1,
        }
    return compartments


def create_species_model(species_comp, areas, vaccination_states, virus_states, t0_susceptible,
                         t0_infectious, non_vaccination_state):
    species = {}
    for index_compartments in species_comp:
        for index_areas in areas:
            for index_vaccination in vaccination_states:
                for index_virus in virus_states:
                    if index_compartments == "susceptible":
                        species_combined = (
                            f"{index_compartments}_{index_areas}_{index_vaccination}"
                        )
                        amount_t0 = t0_susceptible
                    else:
                        species_combined = f"{index_compartments}_{index_areas}_{index_vaccination}_{index_virus}"

                        if index_compartments in ["recovered", "dead"]:
                            amount_t0 = 0
                        elif index_compartments == "infectious":
                            amount_t0 = t0_infectious
                    if index_vaccination != non_vaccination_state:
                        amount_t0 = 0
                    species[species_combined] = {
                        "compartment": index_areas,
                        "initial_amount": amount_t0,
                    }
    return species


def create_parameters_model(
    virus_states, vaccination_states_removed, omega_df, delta_df, eta_df
):
    omega_delta_dict = {}
    eta_dict = {}
    for index_virus in virus_states:
        for index_vaccination in vaccination_states_removed:

            key_omega = f"omega_{index_vaccination}_{index_virus}"
            key_delta = f"delta_{index_vaccination}_{index_virus}"
            omega_delta_dict[key_omega] = {
                "value": omega_df[index_vaccination][index_virus]
            }
            omega_delta_dict[key_delta] = {
                "value": delta_df[index_vaccination][index_virus]
            }

        key_eta = f"eta_{index_virus}"
        eta_dict[key_eta] = {"value": eta_df[index_virus][0]}

    parameters_fixed = {
        "lambda1": {"value": 0.01},
        "p": {"value": 0.01},
        "gamma": {"value": 0.5},
        "beta": {"value": 2},
        "nu_vac1": {"value": 0.01},
        "nu_vac2": {"value": 0.01},
    }

    parameters = {**parameters_fixed, **omega_delta_dict, **eta_dict}
    return parameters


def create_dead_recover_reactions_model(areas, vaccination_states, virus_states):
    dead_recover_reactions = {}
    for index_areas in areas:
        for index_vaccination in vaccination_states:
            for index_virus in virus_states:
                numb_infected = (
                    f"infectious_{index_areas}_{index_vaccination}_{index_virus}"
                )
                numb_recovered = (
                    f"recovered_{index_areas}_{index_vaccination}_{index_virus}"
                )
                numb_dead = f"dead_{index_areas}_{index_vaccination}_{index_virus}"
                key_recover = f"{numb_infected}_to_{numb_recovered}"
                key_death = f"{numb_infected}_to_{numb_dead}"

                if index_vaccination == "vac0":
                    omega_term = "0"
                else:
                    omega_term = f"omega_{index_vaccination}_{index_virus}"

                dead_recover_reactions[key_recover] = {
                    "reactants": {f"{numb_infected}": 1},
                    "products": {f"{numb_recovered}": 1},
                    "formula": f" (1 - (1 - {omega_term}) * p) * lambda1 + {numb_infected}",
                }
                dead_recover_reactions[key_death] = {
                    "reactants": {f"{numb_infected}": 1},
                    "products": {f"{numb_dead}": 1},
                    "formula": f"(1 - {omega_term}) * p * lambda1 * {numb_infected}",
                }
    return dead_recover_reactions


def create_infection_reactions_model(
    species, areas, vaccination_states, virus_states, distance_M
):
    infection_reactions = {}
    for index_areas_infected in areas:
        for index_vaccination_infected in vaccination_states:
            for index_virus_infected in virus_states:
                for index_areas_susceptible in areas:
                    for index_vaccination_susceptible in vaccination_states:

                        numb_susceptible = f"susceptible_{index_areas_susceptible}_{index_vaccination_susceptible}"
                        numb_infectious = f"infectious_{index_areas_infected}_{index_vaccination_infected}_{index_virus_infected}"
                        numb_susceptible_area_population_alive = (
                            get_all_species_alive_area(
                                species_list=species, area=index_areas_susceptible
                            )
                        )
                        eta_term = f"eta_{index_virus_infected}"

                        if index_vaccination_susceptible != "vac0":
                            delta_term = f"delta_{index_vaccination_susceptible}_{index_virus_infected}"
                        else:
                            delta_term = "0"

                        if index_vaccination_infected != "vac0":
                            gamma_term = "gamma"
                        else:
                            gamma_term = "0"
                        modifier = list(species.keys()) #",".join(species.keys())
                        keys = f"{numb_susceptible}_to_{numb_infectious}"
                        infection_reactions[keys] = {
                            "reactants": {
                                f"{numb_susceptible}": 1,
                                f"{numb_infectious}": 1,
                            },
                            "products": {f"{numb_infectious}": 2},
                            "modifiers": modifier,
                            "formula": f"beta * {eta_term} * (1-{delta_term}) * (1-{gamma_term})*{numb_susceptible}*{numb_infectious}*({distance_M.loc[index_areas_susceptible, index_areas_infected]})/({numb_susceptible_area_population_alive})",
                        }
    return infection_reactions


def create_vaccination_reactions_model(
    species_comp, areas, virus_states, vaccination_states, vaccination_states_removed, non_vaccination_state
):
    vaccination_reactions = {}
    for index_compartments in species_comp:
        for index_areas in areas:
            for index_virus in virus_states:
                if index_compartments == 'susceptible':
                    numb_to_be_vaccinated = (
                        f"{index_compartments}_{index_areas}_{non_vaccination_state}"
                    )
                else:
                    numb_to_be_vaccinated = (
                        f"{index_compartments}_{index_areas}_{non_vaccination_state}_{index_virus}"
                    )                    
                for index_vaccination in vaccination_states_removed:
                    if index_compartments == 'susceptible':
                        numb_vaccinated = f"{index_compartments}_{index_areas}_{index_vaccination}"
                    else:
                        numb_vaccinated = f"{index_compartments}_{index_areas}_{index_vaccination}_{index_virus}"
                    vaccination_reactions[
                        f"{numb_to_be_vaccinated}_to_{numb_vaccinated}"
                    ] = {
                        "reactants": {f"{numb_to_be_vaccinated}": 1},
                        "products": {f"{numb_vaccinated}": 1},
                        "formula": f" nu_{index_vaccination}* {numb_to_be_vaccinated}",
                    }
    return vaccination_reactions


def create_reactions_model(
    areas,
    vaccination_states,
    vaccination_states_removed,
    virus_states,
    species,
    distances,
    vaccinated_compartments,
    non_vaccination_state,
):
    # Dead/Recovery
    dead_recover_reactions = create_dead_recover_reactions_model(
        areas=areas, vaccination_states=vaccination_states, virus_states=virus_states
    )

    # infection
    M = get_M_matrix(areas=areas, species=species, distance_matrix=distances)
    infection_reactions = create_infection_reactions_model(
        species=species,
        areas=areas,
        vaccination_states=vaccination_states,
        virus_states=virus_states,
        distance_M=M,
    )

    # vaccination reactions
    vaccination_reactions = create_vaccination_reactions_model(
        species_comp=vaccinated_compartments,
        vaccination_states=vaccination_states,
        areas=areas,
        virus_states=virus_states, vaccination_states_removed=vaccination_states_removed, non_vaccination_state=non_vaccination_state,
    )

    # concatenate reactions
    reactions = {
        **dead_recover_reactions,
        **infection_reactions,
        **vaccination_reactions,
    }
    return reactions

