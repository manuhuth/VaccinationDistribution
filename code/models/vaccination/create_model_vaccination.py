import numpy as np
import pandas as pd

from functions.create_sbml import write_entities_to_dict
from functions.create_sbml import write_dict_to_sbml_file
from functions.run_sbml import create_observables_vaccination_rates
from functions.run_sbml import get_model_and_solver_from_sbml
from functions.vaccine_proportions import create_splines_R
from functions.vaccine_proportions import create_splines
from functions.vaccine_proportions import create_spline_parameters
from functions.vaccine_proportions import create_spline_transformation_rules
from functions.vaccine_proportions import create_one_minus_rules_splines
from functions.vaccine_proportions import create_parameters_piecewise_vaccine_supply
from functions.vaccine_proportions import create_vaccine_supply_rules
from functions.tools import get_all_species_alive
from functions.tools import get_all_species_alive_area
from functions.tools import get_all_individuals_to_be_vaccinated
from functions.tools import get_all_individuals_to_be_vaccinated_by_area
from functions.tools import get_states_by_substrings_any
from functions.tools import create_species_model


def create_model_splines(
    create_model,
    number_areas,
    number_vaccines,
    vaccinated_compartments,
    number_viruses,
    length_decision_period,
    number_yy,
    number_xx_R,
    model_name,
    path_sbml,
    model_directory,
    # R0_periods,
):
    """Create raw AMICI model with splines as functional form of the vaccination
    pathways. All paremeters set to zero. Parameters can subsequently be
    modified.

    Parameters
    ----------
    create_model : bool
        If True, the whole model including c++ code is build into the
        specified directoys. If False, an already created model is loaded
        from the specified directorys.

    number_areas : int
        Number of distinct areas.

    number_vaccines : int
        Number of distinct vaccines.

    number_viruses : int
        Number of distinct viruses.

    length_decision_period : int
        Length of the ODE simulation interval.

    number_yy : int
        Number of spline parameters. number_yy - 1 yields the number of
        polynomial intervals.

    model_name : str
        Name under which the model is stored

    path_sbml : str
        Path where sbml file should is stored.

    model_directory : str
        Path where the c++ code is stored.
    """

    species_comp = ["susceptible", "infectious", "recovered", "dead"]
    virus_states = ["virus" + str(i) for i in range(1, number_viruses + 1)]
    vaccination_states = ["vac" + str(i) for i in range(number_vaccines + 1)]
    non_vaccination_state = "vac0"

    areas = ["country" + chr(i) for i in range(65, 65 + number_areas)]

    distances = pd.DataFrame(
        np.zeros([number_areas, number_areas]), index=areas, columns=areas
    )

    xx_list = ["xx" + str(i) for i in (range(number_yy - 1))]
    parameter_vacc_supply = create_parameters_piecewise_vaccine_supply(
        xx=xx_list,
        vaccination_states=vaccination_states,
        non_vaccination_state=non_vaccination_state,
    )

    vaccination_states_removed = [
        x for x in vaccination_states if x != non_vaccination_state
    ]

    areas_without_last = [x for x in areas if x != areas[-1]]

    spline_parameters = create_spline_parameters(
        areas=areas, vaccination_states_removed=vaccination_states_removed
    )
    spline_transformation_rule = create_spline_transformation_rules(
        areas=areas, vaccination_states_removed=vaccination_states_removed
    )

    minus_one_rules = create_one_minus_rules_splines(
        areas=areas,
        vaccination_states_removed=vaccination_states_removed,
        leave_out="last",
    )

    proportion_rules = {**spline_transformation_rule, **minus_one_rules}

    splines = create_splines(
        areas_without=areas_without_last,
        vaccination_states_removed=vaccination_states_removed,
        stop=length_decision_period,
        length=number_yy,
    )

    splines_R_values = create_splines_R(
        areas=areas,
        stop=length_decision_period,
        length=number_xx_R,
    )

    vaccine_supply_rules = create_vaccine_supply_rules(
        vaccination_states, non_vaccination_state, xx_list
    )

    # R0_rules = create_R0_rules(areas, R0_periods)

    rules = {**proportion_rules, **vaccine_supply_rules}  # , **R0_rules}

    omega = np.zeros([number_viruses, number_vaccines], dtype="float")
    delta = omega
    eta = np.repeat(1.0, number_viruses)

    parameters = model_vaccination_create_sbml(
        path=path_sbml,
        xx=xx_list,
        species_comp=species_comp,
        vaccinated_compartments=vaccinated_compartments,
        areas=areas,
        vaccination_states=vaccination_states,
        non_vaccination_state=non_vaccination_state,
        virus_states=virus_states,
        omega_matrix=omega,
        delta_matrix=delta,
        eta_vector=eta,
        parameter_rules=rules,
        distances=distances,
        additional_parameters={**spline_parameters, **parameter_vacc_supply},
        splines={**splines, **splines_R_values},
        # R0_periods=R0_periods,
    )

    observables_vaccinated_vacc = {}
    for index_area in areas:
        for index_vacc in vaccination_states_removed:
            name = "observable_quantity_" + index_area + "_" + index_vacc
            proportion = f"proportion_{index_area}_{index_vacc}"
            number_vacc = f"number_{index_vacc}"
            formula = f"{proportion} * {number_vacc}"
            observables_vaccinated_vacc[name] = {"name": name, "formula": formula}

    obs = {}
    for index in ["nu", "proportion"]:
        new = create_observables_vaccination_rates(
            vaccination_states_removed=vaccination_states_removed,
            areas=areas,
            name_parameter=index,
        )

        obs = {**obs, **new}

    species = create_species_model(
        vaccination_states,
        non_vaccination_state,
        virus_states,
        areas,
        species_comp,
    )

    proportion_obs = {}
    proportion_obs_sus_inf = {}
    for index_areas in areas:
        name = f"proportion_to_be_vaccinated_{index_areas}"
        name_sus_inf = f"proportion_to_be_vaccinated_sus_inf_{index_areas}"

        to_be_vaccinated = get_all_individuals_to_be_vaccinated_by_area(
            vaccinated_compartments, non_vaccination_state, virus_states, index_areas
        )
        to_be_vaccinated_sus_inf = get_all_individuals_to_be_vaccinated_by_area(
            ["susceptible", "infectious"],
            non_vaccination_state,
            virus_states,
            index_areas,
        )
        whole_population = get_all_species_alive_area(species, index_areas)

        dictionary_save = {
            "name": name,
            "formula": f"({to_be_vaccinated}) / ({whole_population})",
        }
        dictionary_save_sus_inf = {
            "name": name_sus_inf,
            "formula": f"({to_be_vaccinated_sus_inf}) / ({whole_population})",
        }
        save_name = f"observable_to_be_vaccinated_{index_areas}"
        save_name_sus_inf = f"observable_to_be_vaccinated_sus_inf_{index_areas}"
        proportion_obs[save_name] = dictionary_save
        proportion_obs_sus_inf[save_name_sus_inf] = dictionary_save_sus_inf

    dead_states = get_states_by_substrings_any(species, ["dead"])
    objective = "+".join(dead_states)
    objective_dict = {
        "observable_objective": {"name": "objective", "formula": objective}
    }

    observables = {
        # **obs,
        # **observables_vaccinated_vacc,
        # **proportion_obs,
        # **proportion_obs_sus_inf,
        **objective_dict,
    }

    yR_splines = []
    for i in areas:
        yR_spline_country = ["yR0_" + i + f"_{j}" for j in (range(number_xx_R))]
        yR_splines = yR_splines + yR_spline_country

    fixed_parameters = (
        [x for x in parameters if not ("yy_" in x)]
        + ["xx" + str(i) for i in (range(number_yy))]
        + yR_splines
        + ["xx_R0_" + str(i) for i in (range(number_xx_R))]
    )

    model_and_solver = get_model_and_solver_from_sbml(
        path_sbml=path_sbml,
        model_name=model_name,
        model_directory=model_directory,
        only_import=not create_model,
        fixed_parameters=fixed_parameters,
        observables=observables,
    )

    model = model_and_solver["model"]
    solver = model_and_solver["solver"]

    model_information = {
        "vaccine_states": vaccination_states,
        "areas": areas,
        "number_yy": number_yy,
        "T": length_decision_period,
    }
    dict_out = {"model": model, "solver": solver, "information": model_information}

    return dict_out


def model_vaccination_create_sbml(
    path,
    xx,
    vaccination_states=["vac0", "vac1", "vac2"],
    non_vaccination_state="vac0",
    virus_states=["virW", "virM"],
    areas=["countryA", "CountryB"],
    species_comp=["susceptible", "infectious", "recovered", "dead"],
    vaccinated_compartments=["susceptible", "infectious", "recovered"],
    omega_matrix=np.array([[0.5, 0.6], [0.5, 0.6]]),
    delta_matrix=np.array([[0.5, 0.6], [0.5, 0.6]]),
    eta_vector=np.array([[1, 1.3]]),
    single_parameter={
        "lambda": 0.01,
        "prob_deceasing": 0.1,
        "gamma": 0.5,
        "beta": 2,
    },
    parameters_constant=False,
    parameter_rules=None,
    additional_parameters=None,
    splines=None,
    t0_susceptible=0.0,
    t0_infectious=0.0,
    distances=np.array([[0, 3], [3, 0]]),
    check_error=False,
    # R0_periods=10,
):
    """Create model specified by the input parameters. Parameters can
    subsequently be modified.

    Parameters
    ----------
    path : str
        Path where sbml file should be stored.

    vaccination_states : list of strings
        List containing the names of the vaccinations containing a state for
        non-vaccinated individuals.

    non_vaccination_state : str
        Name of state indicates non-vaccinated individuals.

    virus_states : list of strings
        List containing the names of the virus types.

    areas : list of strings
        List containing the names of the areas.

    species_comp : list of strings
        List containing the names of the compartments. Should not be changed.

    vaccinated_compartments : list of strings
        List of compartments from which individuals are vaccinated.

    omega_matrix : numpy.array
        Array indicating the omega coefficients regarding the combinations
        of vaccines and virus types. Each cell indicates that an
        infected individual has a lower probability of :math:`1-\omega` to die.
        Rows indicate the vaccinations and columns the virus types. Number of
        rows must equal number of vaccination states (without non-vaccination
        state) and the number of columns must equal the number of virus types.

    delta_matrix : numpy.array
        Array indicating the delta coefficients regarding the combinations
        of vaccines and virus types. Each cell indicates that a
        vaccinated individual has a lower probability of :math:`1-\delta` to
        become infected while meeting an infectious person.
        Rows indicate the vaccinations and columns the virus types. Number of
        rows must equal number of vaccination states (without non-vaccination
        state) and the number of columns must equal the number of virus types.

    eta_vector : np.array
        Array indicating the degree of infectiouness of the virus types relative
        to the wild type. Number of elements must be the same as the number
        of virus types. The first entry must be the degree of the wild type and
        should be set to one.

    single_parameter : dict
        Dictionary containing the values for :math:`\lambda, p, \gamma, \beta`.
        the keys are the respective names and the values are the values that
        are desired.

    parameters_constant : {True, False}
        If True, all parameters are constant. If False, they are allowed to be
        non-constant. Rules can eb set by the `parameter_rules' argument.

    parameter_rules : dict, None
        Dictionary that defines parameter rules. Keys must be valid inputs for
        :func:`functions.create_sbml.create_parameter_rule`.

    additional_parameters : dict, None
        Allows to specify new additional parameters in a dict. The keys must be
        the names and the values the respective magnitudes. Can be used to
        define parameters that are used for `parameter_rules`.

    t0_susceptible : float
        Number of initial assignment for susceptible states.

    t0_infectious : float
        Number of initial assignment for infectious states.

    distances : numpy.array
        array indicating the distances between countries. Does not have to be
        symmetric and can therefore be used to model seasonal behavior or the
        like.

    check_error : {True, False}
        If True, the SBML model is checked for errors.

    """

    vaccination_states_removed = [
        x for x in vaccination_states if x != non_vaccination_state
    ]
    omega_df = pd.DataFrame(
        omega_matrix, columns=vaccination_states_removed, index=virus_states
    )
    delta_df = pd.DataFrame(
        delta_matrix, columns=vaccination_states_removed, index=virus_states
    )
    eta_df = pd.DataFrame(eta_vector, virus_states).transpose()

    # create compartments
    compartments = create_comaprtments_model(areas)

    # create species
    species = create_species_model(
        species_comp=species_comp,
        areas=areas,
        vaccination_states=vaccination_states,
        non_vaccination_state=non_vaccination_state,
        virus_states=virus_states,
    )

    distance_parameters = {}
    distances_df = pd.DataFrame(distances, columns=areas, index=areas)
    b_distance_matrix = b_distance_function(distances_df)
    D = pd.DataFrame(np.nan, index=areas, columns=areas)

    for index_areas_row in areas:
        for index_areas_col in areas:
            name = f"distance_{index_areas_row}_{index_areas_col}"
            D.loc[index_areas_row, index_areas_col] = name
            distance_parameters[name] = {
                "value": b_distance_matrix.loc[index_areas_row, index_areas_col]
            }

    # create parameters
    parameters = create_parameters_model(
        species=species,
        areas=areas,
        virus_states=virus_states,
        vaccination_states_removed=vaccination_states_removed,
        non_vaccination_state=non_vaccination_state,
        omega_df=omega_df,
        delta_df=delta_df,
        eta_df=eta_df,
        distance_parameters=distance_parameters,
        single_parameter=single_parameter,
        parameters_constant=parameters_constant,
        additional_parameters=additional_parameters,
        t0_susceptible=t0_susceptible,
        t0_infectious=t0_infectious,
        # R0_periods=R0_periods,
    )

    # create reactions
    reactions = create_reactions_model(
        areas=areas,
        vaccination_states=vaccination_states,
        vaccination_states_removed=vaccination_states_removed,
        virus_states=virus_states,
        species=species,
        distances=D,
        vaccinated_compartments=vaccinated_compartments,
        non_vaccination_state=non_vaccination_state,
    )

    assignments = create_initial_assignments_model(species=species)

    nu_rules = create_rules_vaccination_rate(
        vaccinated_compartments=vaccinated_compartments,
        vaccination_states_removed=vaccination_states_removed,
        non_vaccination_state=non_vaccination_state,
        virus_states=virus_states,
        areas=areas,
        xx=xx,
    )

    # events = create_events_model(areas, virus_states) currently not supported

    # rules must be specified in the right order. Since nu needs parameter_rules
    # parameter rules must appear in the dictionary before nu rules!
    if parameter_rules is None:
        parameter_rules_all = nu_rules
    else:
        parameter_rules_all = {**parameter_rules, **nu_rules}

    exclude = []
    for index in parameter_rules_all.keys():
        exclude.append(parameter_rules_all[index]["parameter_id"])

    spline_par = []
    for index1 in areas[0 : (len(areas) - 1)]:
        for index2 in vaccination_states_removed:
            spline_par.append(f"spline_{index1}_{index2}")

    spline_R_par = []
    for index in areas:
        par = f"R0_modifier_{index}"
        spline_R_par.append(par)

    exclude_all = exclude + spline_par + spline_R_par
    parameter_out = [x for x in parameters.keys() if not (x in exclude_all)]

    model_input_dictionary = write_entities_to_dict(
        compartments=compartments,
        species=species,
        parameters=parameters,
        reactions=reactions,
        assignments=assignments,
        parameter_rules=parameter_rules_all,
        rate_rules=None,
        splines=splines,
        # events=events, currently not supported
    )

    write_dict_to_sbml_file(
        dictionary=model_input_dictionary, path=path, error_check=check_error
    )

    return parameter_out


def b_distance_function(x):
    """Returns the function :math:`f(x)= \frac{1}{1+x}`. The function is used
    to transform the distance between two countries to a measure between zero
    and one that is used to regulate the probability of meeting an individual
    from a certain country.

    Parameters
    ----------
    x : float
        Distance between two countries.

    Returns
    -------
    float
        Function value evaluated at :math:`x`

    """
    return 1 / (1 + x)


def get_M_matrix(areas, species, distance_matrix):
    """Computes strings of formulas for the matrix that defines the probability
    of meeting another individual.

    Parameters
    ----------
    areas : list of strings
        List containing the names of the areas.

    species : list of strings
        List containing the names of all species.

    distance_matrix : np.array
        Matrix indicating the distances between countries. Does not have to be
        symmetric and can therefore be used to model seasonal behavior or the
        like.

    Returns
    -------
    M : pd.DataFrame
        Matrix containing the formulas for the M matrix as strings.

    """
    # distances_df = pd.DataFrame(distance_matrix, columns=areas, index=areas)
    # b_distance_matrix = b_distance_function(distances_df)
    numb_total_population_alive = get_all_species_alive(species_list=species)
    M = pd.DataFrame(np.nan, index=areas, columns=areas)
    for index_areas_row in areas:
        for index_areas_col in areas:
            if index_areas_row == index_areas_col:
                M.loc[index_areas_row, index_areas_col] = ""
            else:
                M.loc[
                    index_areas_row, index_areas_col
                ] = f"({get_all_species_alive_area(species_list=species, area=index_areas_row )})/({numb_total_population_alive })*{distance_matrix.loc[index_areas_col, index_areas_row]}"

    if len(areas) == 1:
        diagonal = "1"
    else:
        np.fill_diagonal(M.values, np.nan)
        diagonal = "1-" + M.apply(lambda x: "-".join(x.dropna()), axis=1)

    np.fill_diagonal(M.values, diagonal)

    return M


def create_comaprtments_model(areas):
    """Creates compartments of the model

    Parameters
    ----------
    areas : list of strings
        List containing the names of the areas.

    Returns
    -------
    compartments : dict
        Dictionary that contains the compartment names as keys and dicitonaries,
        that specify all parameters for the compartments, as values.

    """

    compartments = {}
    for index_area in areas:
        compartments[index_area] = {
            "units": "dimensionless",
            "constant": True,
            "volume": 1,
        }
    return compartments


def create_parameters_model(
    species,
    vaccination_states_removed,
    non_vaccination_state,
    virus_states,
    areas,
    omega_df,
    delta_df,
    eta_df,
    single_parameter,
    distance_parameters,
    parameters_constant,
    additional_parameters,
    t0_susceptible,
    t0_infectious,
    # R0_periods,
):
    """Create parameters for the model.

    Parameters
    ----------
    species : list of strings
        List containing the names of all species.

    vaccination_states_removed : list of strings
        List containing the names of the vaccinations not containing a state for
        non-vaccinated individuals.

    non_vaccination_state : str
        Name of state indicates non-vaccinated individuals.

    virus_states : list of strings
        List containing the names of the virus types.

    areas : list of strings
        List containing the names of the areas.

    omega_df : pandas.DataFrame
        Data frame indicating the omega coefficients regarding the combinations
        of vaccines and virus types. Each cell indicates that an
        infected individual has a lower probability of :math:`1-\omega` to die.
        Rows indicate the vaccinations and columns the virus types. Number of
        rows must equal number of vaccination states (without non-vaccination
        state) and the number of columns must equal the number of virus types.

    delta_df : pandas.DataFrame
        data frame indicating the delta coefficients regarding the combinations
        of vaccines and virus types. Each cell indicates that a
        vaccinated individual has a lower probability of :math:`1-\delta` to
        become infected while meeting an infectious person.
        Rows indicate the vaccinations and columns the virus types. Number of
        rows must equal number of vaccination states (without non-vaccination
        state) and the number of columns must equal the number of virus types.

    eta_df : pandas.DataFrame
        Data frame indicating the degree of infectiouness of the virus types
        relative to the wild type. Number of elements must be the same as the
        number of virus types. The first entry must be the degree of the wild
        type and should be set to one.

    single_parameter : dict
        Dictionary containing the values for :math:`\lambda, p, \gamma, \beta`.
        the keys are the respective names and the values are the values that
        are desired.

    parameters_constant : {True, False}
        If True, all parameters are constant. If False, they are allowed to be
        non-constant. Rules can eb set by the `parameter_rules' argument.

    additional_parameters : dict, None
        Allows to specify new additional parameters in a dict. The keys must be
        the names and the values the respective magnitudes. Can be used to
        define parameters that are used for `parameter_rules`.

    t0_susceptible : float
        Number of initial assignment for susceptible states.

    t0_infectious : float
        Number of initial assignment for infectious states.

    Returns
    -------
    parameters : dict
        Dictionary that contains the parameter names as keys and dicitonaries,
        that specify all parameters for the parameters, as values.
    """
    omega_delta_dict = {}
    eta_dict = {}
    for index_virus in virus_states:
        for index_vaccination in vaccination_states_removed:

            key_omega = f"omega_{index_vaccination}_{index_virus}"
            key_delta = f"delta_{index_vaccination}_{index_virus}"
            omega_delta_dict[key_omega] = {
                "value": omega_df[index_vaccination][index_virus],
                "constant": parameters_constant,
            }
            omega_delta_dict[key_delta] = {
                "value": delta_df[index_vaccination][index_virus],
                "constant": parameters_constant,
            }

        key_eta = f"eta_{index_virus}"
        eta_dict[key_eta] = {"value": eta_df[index_virus][0]}

    parameters_fixed = {
        "lambda1": {
            "value": single_parameter["lambda"],
            "constant": parameters_constant,
        },
        "prob_deceasing": {
            "value": single_parameter["prob_deceasing"],
            "constant": parameters_constant,
        },
        "gamma": {"value": single_parameter["gamma"], "constant": parameters_constant},
        "beta": {"value": single_parameter["beta"], "constant": parameters_constant},
        # "t": {"value": 0, "constant": False},
    }

    country_modifier = {}
    for index in areas:
        name = f"R0_modifier_{index}"
        country_modifier[name] = {"value": 1, "constant": False}

    # R0_time = {}
    # for i in range(R0_periods):
    #    name = f"R0_modifier_t{i}"
    #    R0_time[name] = {"value":i*14, "constant":False}
    #    for keys in country_modifier.keys():
    #        name_p = keys + f"_t{i}"
    #        R0_time[name_p] =  {"value": 1, "constant": False}

    virus_events = {}
    for index1 in areas:
        for index2 in virus_states:
            name_trigger_time = f"{index2}_{index1}_appears_time"
            name_trigger_quantity = f"{index2}_{index1}_appears_quantity"
            virus_events[name_trigger_time] = {"value": 0, "constant": False}
            virus_events[name_trigger_quantity] = {"value": 0, "constant": False}

    total_numb_vacc = {}

    for index_vaccination in vaccination_states_removed:
        keys = f"number_{index_vaccination}"
        total_numb_vacc[keys] = {"value": 0, "constant": False}

    vaccination_parameters = {}
    vaccination_proportions = {}
    for index_vaccination in vaccination_states_removed:
        for index_areas in areas:
            keys_nu = f"nu_{index_areas}_{index_vaccination}"
            keys_proportion = f"proportion_{index_areas}_{index_vaccination}"
            vaccination_parameters[keys_nu] = {"value": 0, "constant": False}
            vaccination_proportions[keys_proportion] = {"value": 0, "constant": False}

    all_species = species.keys()
    initial_amount = {}

    for index_species in all_species:
        keys = index_species + "_t0"

        if "susceptible" in index_species:
            amount_t0 = t0_susceptible
        else:
            if "infectious" in index_species:
                amount_t0 = t0_infectious
            else:
                amount_t0 = 0
        if non_vaccination_state not in index_species:
            amount_t0 = 0
        initial_amount[keys] = {"value": amount_t0}

    if additional_parameters is not None:
        parameters = {
            **parameters_fixed,
            **omega_delta_dict,
            **eta_dict,
            **initial_amount,
            **vaccination_parameters,
            **total_numb_vacc,
            **vaccination_proportions,
            **additional_parameters,
            **distance_parameters,
            **country_modifier,
            **virus_events,
            # **R0_time,
        }

    else:
        parameters = {
            **parameters_fixed,
            **omega_delta_dict,
            **eta_dict,
            **initial_amount,
            **vaccination_parameters,
            **total_numb_vacc,
            **vaccination_proportions,
            **distance_parameters,
            **country_modifier,
            **virus_events,
            # **R0_time,
        }

    return parameters


def create_piecewise_R0(area, R0_periods):

    formula = "piecewise("
    for index in range(R0_periods - 1):
        domain = (
            f"((time >= R0_modifier_t{index}) && (time < R0_modifier_t{index + 1}) ),"
        )
        parameter = f"R0_modifier_{area}_t{index},"
        formula += parameter + domain
    formula += f"R0_modifier_{area}_t{index+1})"

    return formula


def create_R0_rules(areas, R0_periods):
    rules = {}
    for area in areas:
        id_rule = f"rule_R0_{area}"
        formula = create_piecewise_R0(area, R0_periods)
        parameter_id = f"R0_modifier_{area}"
        rules[id_rule] = {"parameter_id": parameter_id, "formula": formula}
    return rules


def create_rules_vaccination_rate(
    vaccinated_compartments,
    vaccination_states_removed,
    non_vaccination_state,
    virus_states,
    areas,
    xx,
):
    """Create rules for the vaccination rates.

    Parameters
    ----------
    vaccinated_compartments : list of strings
        List of compartments from which individuals are vaccinated.

    non_vaccination_state : str
        Name of state indicates non-vaccinated individuals.

    virus_states : list of strings
        List containing the names of the virus types.

    areas : list of strings
        List containing the names of the areas.

    Returns
    -------
    rules_nu : dict
        Dictionary including all rules for the vaccination parameters.
    """
    rules_nu = {}
    for index_vaccination in vaccination_states_removed:
        id_number_vaccinations = f"number_{index_vaccination}"
        for index_areas in areas:
            id_rule = f"rule_nu_{index_areas}_{index_vaccination}"
            id_nu = f"nu_{index_areas}_{index_vaccination}"
            id_proportion = f"proportion_{index_areas}_{index_vaccination}"
            to_be_vaccinated = get_all_individuals_to_be_vaccinated(
                vaccinated_compartments=vaccinated_compartments,
                non_vaccination_state=non_vaccination_state,
                virus_states=virus_states,
                areas=[index_areas],
            )

            str_if = "( "
            for index_comp in vaccinated_compartments:
                if index_comp == "susceptible":
                    add = f"(susceptible_{index_areas}_vac0 < 0)"
                    str_if += add
                else:
                    for index_virus in virus_states:
                        add = f" || ({index_comp}_{index_areas}_vac0_{index_virus} < 0)"
                        str_if += add
            str_if += " )"

            formula = "piecewise("
            for index in range(len(xx) - 1):
                domain = f"((time >= {xx[index]}) && (time < {xx[index + 1]}) && ({id_proportion}*vaccine_supply_par_{index_vaccination}_{xx[index]} + 500 <= ({to_be_vaccinated}) )  ),"
                parameter = f"{id_proportion}*vaccine_supply_par_{index_vaccination}_{xx[index]} / ({to_be_vaccinated}),"
                formula += parameter + domain
            formula += f"{id_proportion}*vaccine_supply_par_{index_vaccination}_{xx[len(xx)-1]}/({to_be_vaccinated}), ( (time >= {xx[len(xx)-1]}) && ({id_proportion}*vaccine_supply_par_{index_vaccination}_{xx[len(xx)-1]} + 500 <= ({to_be_vaccinated}))), 0  )"

            # formula = f"piecewise(0, {str_if}, {id_proportion} * {id_number_vaccinations} / ({to_be_vaccinated}) )"

            # formula = f"piecewise({id_proportion} * {id_number_vaccinations} / ({to_be_vaccinated}), {id_proportion}  < ({to_be_vaccinated}) , ({to_be_vaccinated}) )"

            rules_nu[id_rule] = {"parameter_id": id_nu, "formula": formula}

    return rules_nu


def create_dead_recover_reactions_model(
    vaccination_states,
    non_vaccination_state,
    virus_states,
    areas,
):
    """Create death and recover reactions for the model.

    Parameters
    ----------
    species : list of strings
        List containing the names of all species.

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
    dead_recover_reactions: dict
        Dictionary that contains the dead/recover reaction names as keys
        and dicitonaries that contain the reactants, products and formulas
        as values.
    """
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

                if index_vaccination == non_vaccination_state:
                    omega_term = "0"
                else:
                    omega_term = f"omega_{index_vaccination}_{index_virus}"

                dead_recover_reactions[key_recover] = {
                    "reactants": {f"{numb_infected}": 1},
                    "products": {f"{numb_recovered}": 1},
                    "formula": f" (1 - (1 - {omega_term}) * prob_deceasing) * lambda1 * {numb_infected}",
                }
                dead_recover_reactions[key_death] = {
                    "reactants": {f"{numb_infected}": 1},
                    "products": {f"{numb_dead}": 1},
                    "formula": f"(1 - {omega_term}) * prob_deceasing * lambda1 * {numb_infected}",
                }
    return dead_recover_reactions


def create_infection_reactions_model(
    species, vaccination_states, non_vaccination_state, virus_states, areas, distance_M
):
    """Create infection reactions for the model.

    Parameters
    ----------
    species : list of strings
        List containing the names of all species.

    vaccination_states: list of strings
        List containing the names of the vaccinations containing a state for
        non-vaccinated individuals.

    non_vaccination_state : str
        Name of state indicates non-vaccinated individuals.

    virus_states : list of strings
        List containing the names of the virus types.

    areas : list of strings
        List containing the names of the areas.

    distance_M : pd.DataFrame
        Data frame containing the formulas for the M matrix as strings.

    Returns
    -------
    infection_reactions : dict
        Dictionary that contains the infection reaction names as keys
        and dicitonaries that contain the reactants, products and formulas
        as values.
    """

    infection_reactions = {}
    for index_areas_infected in areas:
        for index_vaccination_infected in vaccination_states:
            for index_virus_infected in virus_states:
                for index_areas_susceptible in areas:
                    for index_vaccination_susceptible in vaccination_states:

                        numb_susceptible = f"susceptible_{index_areas_susceptible}_{index_vaccination_susceptible}"
                        numb_susceptible_as_infectious = f"infectious_{index_areas_susceptible}_{index_vaccination_susceptible}_{index_virus_infected}"
                        numb_infectious = f"infectious_{index_areas_infected}_{index_vaccination_infected}_{index_virus_infected}"
                        numb_susceptible_area_population_alive = (
                            get_all_species_alive_area(
                                species_list=species, area=index_areas_susceptible
                            )
                        )
                        eta_term = f"eta_{index_virus_infected}"

                        if index_vaccination_susceptible != non_vaccination_state:
                            delta_term = f"delta_{index_vaccination_susceptible}_{index_virus_infected}"
                        else:
                            delta_term = "0"

                        if index_vaccination_infected != non_vaccination_state:
                            gamma_term = "gamma"
                        else:
                            gamma_term = "0"
                        modifier = list(species.keys())
                        keys = f"{numb_susceptible}_to_{numb_infectious}"

                        if f"{numb_susceptible_as_infectious}" == f"{numb_infectious}":
                            prod = {
                                f"{numb_susceptible_as_infectious}": 2,
                            }
                        else:
                            prod = {
                                f"{numb_susceptible_as_infectious}": 1,
                                f"{numb_infectious}": 1,
                            }

                        infection_reactions[keys] = {
                            "reactants": {
                                f"{numb_susceptible}": 1,
                                f"{numb_infectious}": 1,
                            },
                            "products": prod,
                            "modifiers": modifier,
                            "formula": f"beta * (R0_modifier_{index_areas_infected} + R0_modifier_{index_areas_susceptible}) / 2 * {eta_term} * (1-{delta_term}) * (1-{gamma_term})*{numb_susceptible}*{numb_infectious}*({distance_M.loc[index_areas_susceptible, index_areas_infected]}) / ({numb_susceptible_area_population_alive})",
                        }

    return infection_reactions


def create_vaccination_reactions_model(
    species_comp,
    vaccination_states,
    vaccination_states_removed,
    non_vaccination_state,
    virus_states,
    areas,
):
    """Create vaccination reactions for the model.

    Parameters
    ----------
    species_comp : list of strings
        List containing the names of the compartments.

    vaccination_states : list of strings
        List containing the names of the vaccinations containing a state for
        non-vaccinated individuals.

    vaccination_states_removed : list of strings
        List containing the names of the vaccinations not containing a state for
        non-vaccinated individuals.

    non_vaccination_state : str
        Name of state indicates non-vaccinated individuals.

    virus_states : list of strings
        List containing the names of the virus types.

    areas : list of strings
        List containing the names of the areas.

    Returns
    -------
    vaccination_reactions : dict
        Dictionary that contains the vaccination reaction names as keys
        and dicitonaries that contain the reactants, products and formulas
        as values.
    """
    vaccination_reactions = {}
    for index_compartments in species_comp:
        for index_areas in areas:
            for index_virus in virus_states:
                if index_compartments == "susceptible":
                    numb_to_be_vaccinated = (
                        f"{index_compartments}_{index_areas}_{non_vaccination_state}"
                    )
                else:
                    numb_to_be_vaccinated = f"{index_compartments}_{index_areas}_{non_vaccination_state}_{index_virus}"
                for index_vaccination in vaccination_states_removed:
                    if index_compartments == "susceptible":
                        numb_vaccinated = (
                            f"{index_compartments}_{index_areas}_{index_vaccination}"
                        )
                    else:
                        numb_vaccinated = f"{index_compartments}_{index_areas}_{index_vaccination}_{index_virus}"
                    vaccination_reactions[
                        f"{numb_to_be_vaccinated}_to_{numb_vaccinated}"
                    ] = {
                        "reactants": {f"{numb_to_be_vaccinated}": 1},
                        "products": {f"{numb_vaccinated}": 1},
                        "formula": f" nu_{index_areas}_{index_vaccination} * {numb_to_be_vaccinated}",
                    }

    return vaccination_reactions


def create_reactions_model(
    vaccination_states,
    vaccination_states_removed,
    non_vaccination_state,
    virus_states,
    areas,
    species,
    distances,
    vaccinated_compartments,
):
    """Create reactions for the model.

    Parameters
    ----------
    vaccination_states : list of strings
        List containing the names of the vaccinations containing a state for
        non-vaccinated individuals.

    vaccination_states_removed : list of strings
        List containing the names of the vaccinations not containing a state for
        non-vaccinated individuals.

    non_vaccination_state : str
        Name of state indicates non-vaccinated individuals.

    virus_states : list of strings
        List containing the names of the virus types.

    areas : list of strings
        List containing the names of the areas.

    species : list of strings
        List containing the names of all species.

    distances : np.array
        Matrix indicating the distances between countries. Does not have to be
        symmetric and can therefore be used to model seasonal behavior or the
        like.

    vaccinated_compartments : list of strings
        List of compartments from which individuals are vaccinated.

    Returns
    -------
    reactions : dict
        Dictionary that contains the reaction names as keys
        and dicitonaries that contain the reactants, products and formulas
        as values.
    """
    # Dead/Recovery
    dead_recover_reactions = create_dead_recover_reactions_model(
        areas=areas,
        vaccination_states=vaccination_states,
        virus_states=virus_states,
        non_vaccination_state=non_vaccination_state,
    )

    # infection
    M = get_M_matrix(areas=areas, species=species, distance_matrix=distances)

    infection_reactions = create_infection_reactions_model(
        species=species,
        areas=areas,
        vaccination_states=vaccination_states,
        non_vaccination_state=non_vaccination_state,
        virus_states=virus_states,
        distance_M=M,
    )

    # vaccination reactions
    vaccination_reactions = create_vaccination_reactions_model(
        species_comp=vaccinated_compartments,
        vaccination_states=vaccination_states,
        areas=areas,
        virus_states=virus_states,
        vaccination_states_removed=vaccination_states_removed,
        non_vaccination_state=non_vaccination_state,
    )

    # concatenate reactions
    reactions = {
        **dead_recover_reactions,
        **infection_reactions,
        **vaccination_reactions,
    }
    return reactions


def create_initial_assignments_model(species):
    """Create initial assignments for the model. Here only the formal type
    is defined. Assigning values to the parameters is done by
    :func:`create_parameters_model`.

    Parameters
    ----------
    species : list of strings
        List containing the names of all species.

    Returns
    -------
    assignments: dict
        Dictionary that contains the assignment names as keys
        and dicitonaries that contain the species' ids, and formulas
        as values.
    """

    all_species = species.keys()
    assignments = {}
    for index_species in all_species:
        keys = f"assignment_{index_species}"
        assignments[keys] = {
            "species_id": index_species,
            "formula": f"{index_species}_t0",
        }
    return assignments


def create_events_model(areas, virus_states):
    """Create events for the model.

    Parameters
    ----------
    virus_states : list of strings
        List containing the names of all virus variants.

    Returns
    -------
    events: dict
        Dictionary that contains the event names as keys
        and dicitonaries that contain the event ids, and formulas
        as values.
    """

    events = {}
    for index1 in areas:
        for index2 in virus_states:
            keys = f"event_{index1}_{index2}"
            assignee = f"infectious_{index1}_vac0_{index2}"
            trigger_time_par = f"{index2}_{index1}_appears_time"
            trigger_quantity = f"{index2}_{index1}_appears_quantity"

            events[keys] = {
                "trigger_formula": f"geq(time, {trigger_time_par})",
                "assignee_id": assignee,
                "assign_formula": trigger_quantity,
            }
    return events
