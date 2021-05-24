import amici
import numpy as np
import pandas as pd


def run_model(
    model,
    solver,
    periods,
    length_periods,
    set_start_parameter,
    set_parameter,
    observables_names,
    set_initials_zero=True,
    number_intervals=1000
):
    """
    Parameters
    ----------
    model : amici.Model
        Model to run.

    solver : amici.Solver
        Solver for amici model.

    periods : int
        Number of start and stop periods to be simulated.

    length_periods : int
        Length of each start and stop period.

    set_start_parameter : dict
        Dictionary specifying the start parameter.

    set_parameter : dict
        Dictionary containing other parameter specifications.

    observables_names : list of strings
        Names of the observables.

    Returns
    -------
    trajectory_dict : dict
        Dictionary with keys `states` and `observables` and the respective
        data frames of trajectories as values.
    """
    # set all first values to zero, to prevent state
    if set_initials_zero is True:
        model = set_all_initial_conditions_to_zero(model)

    timepoints = np.linspace(0, length_periods-1, number_intervals)
    periods_minus_one = periods - 1

    set_parameter_first = {**set_start_parameter, **set_parameter}

    first_period = run_model_once(
        model=model,
        solver=solver,
        timepoints=timepoints,
        set_parameter=set_parameter_first,
    )

    states = model.getStateIds()
    df_trajectories_states = pd.DataFrame(first_period["x"], columns=(states))
    df_trajectories_observables = pd.DataFrame(
        first_period["y"], columns=(observables_names)
    )

    if periods_minus_one > 0:
        for _ in range(periods_minus_one):
            start_values_it = use_end_values_as_start_values(
                df_trajectories_states=df_trajectories_states
            )
            set_parameter_it = {**start_values_it, **set_parameter}
            model_result = run_model_once(
                model=model,
                solver=solver,
                timepoints=timepoints,
                set_parameter=set_parameter_it,
            )

            states_iteration = model_result["x"]
            observables_iteration = model_result["y"]

            df_trajectories_states_iteration = pd.DataFrame(
                states_iteration, columns=(states)
            )

            df_trajectories_observables_iteration = pd.DataFrame(
                observables_iteration, columns=(observables_names)
            )

            df_trajectories_states = pd.concat(
                [df_trajectories_states, df_trajectories_states_iteration],
                ignore_index=True,
            )

            df_trajectories_observables = pd.concat(
                [df_trajectories_observables, df_trajectories_observables_iteration],
                ignore_index=True,
            )

    trajectory_dict = {
        "states": df_trajectories_states,
        "observables": df_trajectories_observables,
    }
    return trajectory_dict


def set_all_initial_conditions_to_zero(model):
    """Set all initial conditions to zero.

    Parameters
    ----------
    model : amici.model
        Amici model.

    Returns
    -------
    model : amici.model
        Amici model with initial parameters parameters zero.
    """
    parameter_names_all = model.getParameterNames()
    parameter_names_start = [x for x in parameter_names_all if "t0" in x]
    for keys in parameter_names_start:
        model.setParameterByName(keys, 0)

    return model


def use_end_values_as_start_values(df_trajectories_states, identifier_start="t0"):
    """Returns the end values of a model run as dictionary that can be passed
    to set new starting values for a new model run.

    Parameters
    ----------
    df_trajectories : pandas.DataFrame
        Data frame containing the trajectories of all states.

    identifier_start : str
        identifier tha is used to identify starting values.

    Returns
    -------
    new_start_values_dict : dict
    Dictionary that contains the names of the starting parameters as keys and
    the last values of the previous trajectory as values.
    """
    df_end_values = df_trajectories_states.iloc[[-1]]

    new_start_values_dict = {}
    for index_states in df_trajectories_states.columns:
        string_with_identifier = f"{index_states}_{identifier_start}"
        new_start_values_dict[string_with_identifier] = float(
            df_end_values[index_states]
        )

    return new_start_values_dict


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


def run_model_once(model, solver, timepoints, set_parameter=None):
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
