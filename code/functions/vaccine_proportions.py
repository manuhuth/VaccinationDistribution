import numpy as np

from sympy import symbols
from amici.splines import UniformGrid

from models.vaccination.create_model_vaccination import get_all_species_alive
from models.vaccination.create_model_vaccination import get_all_species_alive_area
from models.vaccination.create_model_vaccination import create_species_model


def create_spline(area, vaccine, stop, length):
    """Generate dictionary with information for the spline.
    Spline is associated to a vaccine and an area.

    Parameters
    ----------
    area : str
        Name of area for which spline is created.

    vaccine : str
        Name of vaccine for which spline is created.

    stop : int
        Last period for the uniform grid on whoch the spline is computed.

    length : int
        Number of intervalss on which the spline is computed.

    Returns
    -------
    splines : dict
        Dicitionary containing the information of the spline.
    """
    xx = UniformGrid(start=0, stop=stop, length=length)
    yy_names = []
    xx_names = []

    for index in range(length):
        yy_names.append(symbols(f"yy_{area}_{vaccine}_{index}"))
        xx_names.append(symbols(f"xx{index}"))

    yy = np.repeat(0, length)
    splines = {
        f"spline_{area}_{vaccine}": {
            "time_symbol": "time",
            "xx": xx,
            "yy": yy,
            "xx_names": xx_names,
            "yy_names": yy_names,
        }
    }

    return splines


def create_splines(areas_without, vaccination_states_removed, stop, length):
    """Generate dictionary with information for spline that is passed to
    :func: `model_vaccination_create_sbml`. Spline is associated to a vaccine
    and an area.

    Parameters
    ----------
    areas_without : list
        Names of areas for which splines are created.

    vaccination_states_removed : list
        Names of vaccines for which splines are created.

    stop : int
        Last period for the uniform grid on whoch the spline is computed.

    length : int
        Number of intervalss on which the spline is computed.

    Returns
    -------
    splines : dict
        Dicitionary containing the information of that is passed to
        :func: `model_vaccination_create_sbml`.
    """

    splines = {}
    for index_areas in areas_without:
        for index_vaccines in vaccination_states_removed:
            new_spline = create_spline(
                area=index_areas, vaccine=index_vaccines, stop=stop, length=length
            )
            splines = {**splines, **new_spline}

    return splines


def create_spline_parameters(areas, vaccination_states_removed, leave_out="last"):
    """Generate dictionary with information for the actual spline parameter
    that is passed to :func: `model_vaccination_create_sbml`.

    Parameters
    ----------
    areas : list
        Names of areas for which splines are created.

    vaccination_states_removed : list
        Names of vaccines for which splines are created.

    leave_out : {'last', str, None}
        Name of area for which no spline should be computed. If 'last' the last
        area from the list areas is left out. We need to leave on area out
        such that the last area becomes one minus the sum of all other areas.

    Returns
    -------
    spline_parameters : dict
        Dicitionary containing the information of that is passed to
        :func: `model_vaccination_create_sbml`.
    """

    if leave_out == "last":
        areas_to_use = areas[:-1]
    elif leave_out is None:
        areas_to_use = areas
    else:
        areas_to_use = leave_out

    spline_parameters = {}
    for index_areas in areas_to_use:
        for index_vaccine in vaccination_states_removed:
            parameter_name = f"spline_{index_areas}_{index_vaccine}"
            spline_parameters[parameter_name] = {"value": 0, "constant": False}

    return spline_parameters


# create transformation
def create_spline_transformation_rules(
    areas, vaccination_states_removed, leave_out="last"
):
    """Generate dictionary with information for the actual spline transformation
    rules that is passed to :func: `model_vaccination_create_sbml`.

    Parameters
    ----------
    areas : list
        Names of areas for which splines are created.

    vaccination_states_removed : list
        Names of vaccines for which splines are created.

    leave_out : {'last', str, None}
        Name of area for which no spline should be computed. If 'last' the last
        area from the list areas is left out. We need to leave on area out
        such that the last area becomes one minus the sum of all other areas.

    Returns
    -------
    rules_transformation_splines: dict
        Dicitionary containing the information of that is passed to
        :func: `model_vaccination_create_sbml`.
    """
    if leave_out == "last":
        areas_to_use = areas[:-1]
    elif leave_out is None:
        areas_to_use = areas
    else:
        areas_to_use = leave_out
    rules_transformation_splines = {}
    for index_areas in areas_to_use:
        for index_vaccine in vaccination_states_removed:
            parameter_id = f"proportion_{index_areas}_{index_vaccine}"
            spline_id = f"spline_{index_areas}_{index_vaccine}"
            rule_id = "rule_transformation" + spline_id
            formula = f"1 / (1 + 2.718281828^(-({spline_id})))"
            rules_transformation_splines[rule_id] = {
                "parameter_id": parameter_id,
                "formula": formula,
            }
    return rules_transformation_splines


def create_one_minus_rules_splines(areas, vaccination_states_removed, leave_out="last"):
    """Generate dictionary with information for the rules for the splie
    parameter that is passed to :func: `model_vaccination_create_sbml`.

    Parameters
    ----------
    areas : list
        Names of areas for which splines are created.

    vaccination_states_removed : list
        Names of vaccines for which splines are created.

    leave_out : {'last', str, None}
        Name of area for which no spline should be computed. If 'last' the last
        area from the list areas is left out. We need to leave on area out
        such that the last area becomes one minus the sum of all other areas.

    Returns
    -------
    rules_minus_one : dict
        Dicitionary containing the information of that is passed to
        :func: `model_vaccination_create_sbml`.
    """
    if leave_out == "last":
        area = areas[-1]
    elif leave_out is None:
        area = leave_out

    areas_to_use = [x for x in areas if x != area]

    rules_minus_one = {}
    for index_vaccines in vaccination_states_removed:
        formula = get_string_formula_one_minus_splines(
            other_areas=areas_to_use,
            vaccine=index_vaccines,
        )
        parameter_id = f"proportion_{area}_{index_vaccines}"
        rule_id = f"rule_proportion_one_minus_{index_vaccines}"
        rules_minus_one[rule_id] = {"parameter_id": parameter_id, "formula": formula}

    return rules_minus_one


def create_parameters_piecewise_vaccine_supply(
    xx,
    vaccination_states,
    non_vaccination_state,
):
    """Generate dictionary with information for the piecewise parameter
    that is passed to :func: `model_vaccination_create_sbml`.

    Parameters
    ----------
    decision_period_length : int
        Last time point on which the piecewise function is computed.

    number_decision_periods : int
        Number of intervals for piecewise function.

    first_vaccination_period : int
        Start period of avccination.

    vaccination_states : list
        Names of vaccination states including a state for non-vaccinated.

    non_vaccination_state : str
        State that is not vaccinated.

    areas : list
        Names of areas for which parameters are created.

    Returns
    -------
    parameter_proportions_vaccination: dict
        Dicitionary containing the information of that is passed to
        :func: `model_vaccination_create_sbml`.
    """

    vaccination_states_removed = [
        x for x in vaccination_states if x != non_vaccination_state
    ]

    parameter_vaccine_supply = {}
    for index_vaccination in vaccination_states_removed:
        for index_time in xx:
            vaccine_supply_par_id = (
                f"vaccine_supply_par_{index_vaccination}_{index_time}"
            )
            parameter_vaccine_supply[vaccine_supply_par_id] = {
                "value": 0,
                "constant": True,
            }

    return parameter_vaccine_supply


def get_piecewise_formula_xx_period(vaccine, xx):

    formula = "piecewise("
    for index in range(len(xx) - 1):
        domain = f"((time >= {xx[index]}) && (time < {xx[index + 1]})),"
        parameter = f"vaccine_supply_par_{vaccine}_{xx[index]},"
        formula += parameter + domain
    formula += "0)"

    return formula


def create_vaccine_supply_rules(vaccination_states, non_vaccination_state, xx):
    vaccination_states_removed = [
        x for x in vaccination_states if x != non_vaccination_state
    ]

    rules_vaccine_supply = {}
    for index_vaccine in vaccination_states_removed:
        formula = get_piecewise_formula_xx_period(vaccine=index_vaccine, xx=xx)
        id_rule = f"supply_{index_vaccine}_rule"
        id_parameter = f"number_{index_vaccine}"

        rules_vaccine_supply[id_rule] = {
            "parameter_id": id_parameter,
            "formula": formula,
        }

    return rules_vaccine_supply


def create_parameters_piecewise(
    decision_period_length,
    number_decision_periods,
    first_vaccination_period,
    vaccination_states,
    non_vaccination_state,
    areas,
):
    """Generate dictionary with information for the piecewise parameter
    that is passed to :func: `model_vaccination_create_sbml`.

    Parameters
    ----------
    decision_period_length : int
        Last time point on which the piecewise function is computed.

    number_decision_periods : int
        Number of intervals for piecewise function.

    first_vaccination_period : int
        Start period of avccination.

    vaccination_states : list
        Names of vaccination states including a state for non-vaccinated.

    non_vaccination_state : str
        State that is not vaccinated.

    areas : list
        Names of areas for which parameters are created.

    Returns
    -------
    parameter_proportions_vaccination: dict
        Dicitionary containing the information of that is passed to
        :func: `model_vaccination_create_sbml`.
    """

    vaccination_states_removed = [
        x for x in vaccination_states if x != non_vaccination_state
    ]

    number_periods = decision_period_length * number_decision_periods
    decision_periods = (
        np.array(range(0, number_periods, decision_period_length))
        + first_vaccination_period
    )

    areas_no_last = areas[:-1]
    parameter_proportions_vaccination = {}
    for index_vaccination in vaccination_states_removed:
        for index_areas in areas_no_last:
            for index_time in decision_periods:
                proportion_par_id = (
                    f"proportion_par_{index_areas}_{index_vaccination}_{index_time}"
                )
                parameter_proportions_vaccination[proportion_par_id] = {
                    "value": 0,
                    "constant": True,
                }

    return parameter_proportions_vaccination


def get_piecewise_string_formula(
    areas,
    area_index,
    vac_index,
    first_vaccination_period,
    decision_periods,
    decision_period_length,
    minus_others=False,
):
    """
    Generate piecewise formulas for proportions.

    Parameters
    ----------
    area_index : str
        Name of area for which proportion is computed.

    vac_index : str
        Name of vaccine for which proportion is computed.

    first_vaccination_period : int
        First period in which vaccines are used.

    decision_periods : int
        Number of decision periods.

    decision_period_length : int
        Length of decision periods.

    Returns
    -------
    string_formula : str
    Formula that determines the piecewise proportions.
    """
    no_last_area = areas[:-1]
    string_formula = f"piecewise(0, time < {first_vaccination_period},"
    for index_decision_periods in decision_periods:
        if minus_others is False:
            proportion_id = (
                f" proportion_par_{area_index}_{vac_index}_{index_decision_periods},"
            )

        else:
            proportion_id = get_string_formula_one_minus(
                other_areas=no_last_area,
                vaccine=vac_index,
                decision_period=index_decision_periods,
            )

        period_id = f" ((time >= {index_decision_periods}) && (time < {index_decision_periods + decision_period_length})),"

        string_formula = string_formula + proportion_id + period_id

    string_formula = string_formula + " 0)"

    return string_formula


def get_string_formula_one_minus(other_areas, vaccine, decision_period):

    string_formula = "1 "
    for index_areas in other_areas:
        string_formula = (
            string_formula
            + f"- proportion_par_{index_areas}_{vaccine}_{decision_period},"
        )

    return string_formula


def get_string_formula_one_minus_splines(
    other_areas,
    vaccine,
):
    string_formula = "1 "
    for index_areas in other_areas:
        string_formula = string_formula + f"- proportion_{index_areas}_{vaccine}"

    return string_formula


def create_rules_vaccination_proportion_piecewise(
    decision_period_length,
    number_decision_periods,
    first_vaccination_period,
    vaccination_states,
    non_vaccination_state,
    areas,
):
    """Create rules for the vaccination proportions as relative country sizes.

    Parameters
    ----------
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

    number_periods = decision_period_length * number_decision_periods
    decision_periods = (
        np.array(range(0, number_periods, decision_period_length))
        + first_vaccination_period
    )

    last_area = areas[-1]

    rules_proportions_vaccination = {}
    for index_areas in areas:
        for index_vaccination in vaccination_states_removed:
            id_rule = f"rule_proportion_{index_areas}_{index_vaccination}"
            id_proportion = f"proportion_{index_areas}_{index_vaccination}"

            if index_areas != last_area:
                formula_string = get_piecewise_string_formula(
                    areas=areas,
                    area_index=index_areas,
                    vac_index=index_vaccination,
                    first_vaccination_period=first_vaccination_period,
                    decision_periods=decision_periods,
                    decision_period_length=decision_period_length,
                    minus_others=False,
                )
            else:
                formula_string = get_piecewise_string_formula(
                    areas=areas,
                    area_index=index_areas,
                    vac_index=index_vaccination,
                    first_vaccination_period=first_vaccination_period,
                    decision_periods=decision_periods,
                    decision_period_length=decision_period_length,
                    minus_others=True,
                )

            formula = formula_string
            rules_proportions_vaccination[id_rule] = {
                "parameter_id": id_proportion,
                "formula": formula,
            }

    return rules_proportions_vaccination


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
