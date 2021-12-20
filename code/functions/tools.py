import numpy as np
import pandas as pd
import amici
import pypesto
import pypesto.optimize as optimize
import pymoo.optimize
import requests
import multiprocessing as mp
import tqdm

from pymoo.core.problem import ElementwiseProblem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.factory import get_sampling, get_crossover, get_mutation
from pymoo.factory import get_termination
from pymoo.factory import get_reference_directions
from pymoo.algorithms.moo.ctaea import CTAEA

from estimagic import minimize


def get_R_spline_values(
    areas,
    spline_xx,
    url="/home/manuel/Documents/VaccinationDistribution/mobility.ods",
    sheet_name="Sheet1",
    countries=["Belgium", "France", "Germany", "UK"],
    scale=1,
    end=0.95,
):

    df_raw = pd.read_excel(url, sheet_name=sheet_name, skipfooter=0)

    times = df_raw["Time"]
    times_xx = np.array(list(spline_xx.values()))

    time_mapping = {}
    for j in times_xx:
        key = f"{j}"
        if j <= np.max(times):
            value = times[np.argmin(np.abs(times - j))]
        else:
            value = np.max(times) + 100
        time_mapping[key] = value

    par_R = {}
    for j in range(len(spline_xx.keys())):
        for k in range(len(areas)):
            key = f"yR0_{areas[k]}_{j}"
            if times_xx[j] <= np.max(times):
                time_odf = time_mapping[f"{times_xx[j]}"]
                value = (
                    float(1 + df_raw[countries[k]]
                          [(df_raw["Time"] == time_odf)])
                    * scale
                )
            else:
                value = end

            par_R[key] = value

    return par_R


def get_all_species_alive(species_list):
    """Get all names of species that are alive.

    Parameters
    ----------
    species_list : list of strings
        List of all species.

    Returns
    -------
    numb_total_population_alive : list of strings
        List of all species excluding dead species.
    """

    total_population_alive = [
        x for x in species_list.keys() if ("dead" not in x)]
    numb_total_population_alive = "+".join(total_population_alive)
    return numb_total_population_alive


def get_all_species_alive_area(species_list, area):
    """Get all names of species that are alive from a certain area.

    Parameters
    ----------
    species_list : list of strings
        List of all species.

    area : str
        name of area for which the species should be returned.

    Returns
    -------
    numb_total_population_alive : list of strings
        List of all species excluding dead species from the specified area.
    """

    total_population_alive = [
        x for x in species_list.keys() if ("dead" not in x) and (area in x)
    ]
    numb_total_population_alive = "+".join(total_population_alive)
    return numb_total_population_alive


def get_all_individuals_to_be_vaccinated(
    vaccinated_compartments, non_vaccination_state, virus_states, areas
):
    """Get sum of all names of species that have to be vaccinated.

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
    vaccinated_individuals : str
        Sum of all names of species that have to be vaccinated.
    """

    vaccinated_individuals = ""
    for index_compartments in vaccinated_compartments:
        for index_virus in virus_states:
            for index_areas in areas:
                if index_compartments != "susceptible":
                    state = f"{index_compartments}_{index_areas}_{non_vaccination_state}_{index_virus}"
                    vaccinated_individuals = vaccinated_individuals + "+" + state

    for index_areas in areas:
        state = f"susceptible_{index_areas}_{non_vaccination_state}"
        vaccinated_individuals = vaccinated_individuals + "+" + state

    return vaccinated_individuals


def get_all_individuals_to_be_vaccinated_by_area(
    vaccinated_compartments, non_vaccination_state, virus_states, area
):
    """Get sum of all names of species that have to be vaccinated.

    Parameters
    ----------
    vaccinated_compartments : list of strings
        List of compartments from which individuals are vaccinated.

    non_vaccination_state : str
        Name of state indicates non-vaccinated individuals.

    virus_states : list of strings
        List containing the names of the virus types.

    area : str
        Name of the area.

    Returns
    -------
    vaccinated_individuals : str
        Sum of all names of species that have to be vaccinated.
    """

    vaccinated_individuals = ""
    for index_compartments in vaccinated_compartments:
        for index_virus in virus_states:
            for index_areas in [area]:
                if index_compartments != "susceptible":
                    state = f"{index_compartments}_{index_areas}_{non_vaccination_state}_{index_virus}"
                    vaccinated_individuals = vaccinated_individuals + "+" + state

    for index_areas in [area]:
        state = f"susceptible_{index_areas}_{non_vaccination_state}"
        vaccinated_individuals = vaccinated_individuals + "+" + state

    return vaccinated_individuals


def create_species_model(
    vaccination_states,
    non_vaccination_state,
    virus_states,
    areas,
    species_comp,
):
    """Create species of the model.

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

    species_comp : list of strings
        List containing the names of the compartments. Should not be changed.


    Returns
    -------
    species : dict
        Dictionary that contains the species names as keys and dicitonaries,
        that specify all parameters for the species, as values.
    """

    species = {}
    for index_compartments in species_comp:
        for index_areas in areas:
            for index_vaccination in vaccination_states:
                for index_virus in virus_states:
                    if index_compartments == "susceptible":
                        species_combined = (
                            f"{index_compartments}_{index_areas}_{index_vaccination}"
                        )
                    else:
                        species_combined = f"{index_compartments}_{index_areas}_{index_vaccination}_{index_virus}"

                    amount_t0 = 0
                    species[species_combined] = {
                        "compartment": index_areas,
                        "initial_amount": amount_t0,
                    }
    return species


def get_name_parameter_to_optimize(observables, created_model, timepoints):
    """Get the parameter to optimize..

    Parameters
    ----------
    observables : pd.DataFrame
        Data Frame of the observables.

    created_model : dict
        Dictionary of containing the information of the creatde model ith the key
        'informaton'.

    timepoints : np.array
        Timepoints used for the simulation.

    Returns
    -------
    par_optimization : list
        List of all vaccine parameters that have an impact on the results.

    """
    model = created_model["model"]
    string = "observable_nu"
    nu_columns = [col for col in observables.columns if string in col]
    observables_nu = observables[nu_columns]

    last_vac_sim = (observables_nu == 0).all(axis=1).idxmax() - 1
    last_vac_period = timepoints[last_vac_sim]
    period_length = created_model["information"]["T"] / (
        created_model["information"]["number_yy"] - 1
    )
    number_periods = int(np.ceil(last_vac_period / period_length))

    par_optimization = [
        x
        for x in model.getParameterNames()
        if "yy_" in x and float(x.rsplit("_", 1)[-1]) <= number_periods
    ]

    return par_optimization


def get_sum_of_states(model, states, state_type=["dead"], final_amount=True):
    """Get the sum of output states.

    Parameters
    ----------
    model : amici.model
        Model to use.

    trajectory_dict : dict
        Dictionary of state trajectories.

    state_type : list
        List including the substrings of the states
        that should be used to sum over.

    final_amount : {True, False}
        If True only the final amount is considered. If False the sum of all
        periods is used.

    Returns
    -------
    sum_states : float
        Sum of the states as specified by the function.

    """

    df_trajectories_states = states
    substates = get_substates(
        model=model, substrings=state_type, include_all=True)
    df_substates = df_trajectories_states[substates]

    if final_amount is True:
        df_sum = df_substates.iloc[[-1]]
    else:
        df_sum = df_substates

    sum_states = df_sum.values.sum()

    return sum_states


def get_substates(model, substrings, include_all=True):
    """Get states that contain certain substrings.

    Parameters
    ----------
    model : libsbml.Model
        SBML model to retrieve the names of the species.

    substrings : list of strings
        Substrings that occur in the name of the species that should be
        addressed.

    include_all : {True, False}
        If `True`, only states that include all substrings are returned. If
        `False` all states that include one substring are returned.

    Returns
    -------
    substates : list
        List containing all species' names that have the substrings in their
        names.
    """

    states = model.getStateIds()

    if include_all is False:
        substates = list(
            set(get_states_by_substrings_any(
                states=states, substrings=substrings))
        )
    else:
        substates = get_states_by_substrings_all(
            states=states, substrings=substrings)

    return substates


def get_states_by_substrings_all(states, substrings):
    """Get state names that contain all substrings.

    Parameters
    ----------
    states : list of strings
        Names of all species in the model.

    substrings : list of strings
        Substrings that occur in the name of the species that should be
        addressed.

    Returns
    -------
    state : list
        List containing all species' names that have the substrings in their
        names.

    """
    state = []
    for index_states in states:
        if all([x in index_states for x in substrings]):
            state.append(index_states)

    return state


def get_states_by_substrings_any(states, substrings):
    """Get state names that contain any substring.

    Parameters
    ----------
    states : list of strings
        Names of all species in the model.

    substrings : list of strings
        Substrings that occur in the name of the species that should be
        addressed.

    Returns
    -------
    state : list
        List containing all species' names that have the substrings in their
        names.

    """
    state = []
    for index in substrings:
        states_removed = [x for x in states if index in x]
        state += states_removed

    return state


def optimize_pyPesto(
    model,
    solver,
    parameter_to_optimize,
    optimizer=pypesto.optimize.ScipyOptimizer("L-BFGS-B"),
    n_starts=50,
    lower_bound=10 ** (-5),
    trace_record=True,
    pareto_values=None,
):
    """Run optimization over specified parameters with pyPesto.

    Parameters
    ----------
    model : amici.model
        Amici model which should be evaluated.

    solver : amici.solver
        Solver object of the corresponding amici model.

    parameter_to_optimize : list
        List of parameters you want to optimize over.

    optimizer : pypesto.optimize object.
        pyPesto routine used for optimization.

    n_starts : int
        Number of strats for the multi-start optimization.

    lower_bound : float
        Lower bound for the parameter estimation.

    trace_record : bool
        If True, trace record is set to True for pyPesto optimizers.

    pareto_values : None or array-like
        If an array is passed, its length must equal the number of areas
        within the model. The array is sued to enforce maximum values each areas
        can have in the optimal solution.


    Returns
    -------
    results : pd.DataFrame
        List of results.
    """

    def function_to_optimize_pyPesto(theta):
        theta_trans = np.log(theta / (1 - theta))
        inputs = dict(zip(parameter_to_optimize, theta_trans))
        for i in inputs:
            model.setParameterByName(i, inputs[i])

        rdata = amici.runAmiciSimulation(model, solver)
        states = pd.DataFrame(rdata["x"], columns=model.getStateIds())

        dead = get_sum_of_states(model, states, state_type=[
                                 "dead"], final_amount=True)

        if pareto_values is not None:
            observed_pareto = []
            for index in range(65, 65 + len(pareto_values)):
                country = f"country{chr(index)}"
                observed_pareto.append(
                    get_sum_of_states(
                        model, states, state_type=[country, "dead"], final_amount=True
                    )
                )
            if any(pareto_values > observed_pareto):
                dead = 10 ** 15

        return dead

    objective = pypesto.Objective(fun=function_to_optimize_pyPesto)
    lb = np.repeat(lower_bound, len(parameter_to_optimize))
    ub = 1 - lb
    problem_unrestricted = pypesto.Problem(objective=objective, lb=lb, ub=ub)

    results = pypesto.optimize.minimize(
        problem=problem_unrestricted,
        optimizer=optimizer,
        n_starts=n_starts,
        history_options=pypesto.HistoryOptions(trace_record=trace_record),
    )

    return results


def optimize_estimagic(
    model,
    solver,
    parameter_to_optimize,
    n_starts=50,
    lower_bound=10 ** (-6),
    pareto_values=None,
):
    """Run optimization over specified parameters with estimagic.

    Parameters
    ----------
    model : amici.model
        Amici model which should be evaluated.

    solver : amici.solver
        Solver object of the corresponding amici model.

    parameter_to_optimize : list
        List of parameters you want to optimize over.

    n_starts : int
        Number of strats for the multi-start optimization.

    lower_bound : float
        Lower bound for the parameter estimation.

    pareto_values : None or array-like
        If an array is passed, its length must equal the number of areas
        within the model. The array is sued to enforce maximum values each areas
        can have in the optimal solution.


    Returns
    -------
    results : pd.DataFrame
        List of results.
    """

    def function_to_optimize_estimagic(theta):

        theta_trans = np.log(theta["value"] / (1 - theta["value"]))
        inputs = dict(zip(parameter_to_optimize, theta_trans))
        for i in inputs:
            model.setParameterByName(i, inputs[i])

        rdata = amici.runAmiciSimulation(model, solver)
        states = pd.DataFrame(rdata["x"], columns=model.getStateIds())

        dead = get_sum_of_states(model, states, state_type=[
                                 "dead"], final_amount=True)

        return {"value": dead}

    start_par = np.random.uniform(
        lower_bound, 1 - lower_bound, len(parameter_to_optimize) * n_starts
    ).reshape(n_starts, len(parameter_to_optimize))

    save_df = pd.DataFrame(
        np.NaN,
        columns=list(parameter_to_optimize) + ["fval", "feval"],
        index=range(n_starts),
    )

    for i in range(n_starts):
        print(f"Bobyqa: {i}/{n_starts} start.")

        start_theta = start_par[i]
        parameter_np = np.array(
            [
                start_theta,
                np.repeat(lower_bound, len(parameter_to_optimize)),
                np.repeat(1 - lower_bound, len(parameter_to_optimize)),
            ]
        )
        parameter_df = pd.DataFrame(
            parameter_np.transpose(),
            columns=["value", "lower_bound", "upper_bound"],
            index=parameter_to_optimize,
            dtype="float",
        )
        r = minimize(
            criterion=function_to_optimize_estimagic,
            params=parameter_df,
            algorithm="nag_pybobyqa",
        )

        save_df.iloc[i] = list(r["solution_x"]) + [
            r["solution_criterion"],
            r["n_criterion_evaluations"],
        ]

    return save_df


def optimize_fides(
    model,
    solver,
    parameter_to_optimize,
    n_starts,
    lower_bound,
):
    """Run optimization over specified parameters with pyPesto's fides
    optimizer.

    Parameters
    ----------
    model : amici.model
        Amici model which should be evaluated.

    solver : amici.solver
        Solver object of the corresponding amici model.

    parameter_to_optimize : list
        List of parameters you want to optimize over.

    n_starts : int
        Number of strats for the multi-start optimization.

    lower_bound : float
        Lower bound for the parameter estimation.


    Returns
    -------
    result : amici.returnData
        Data frame of results.
    """

    solver.setSensitivityMethod(amici.SensitivityMethod_forward)
    solver.setSensitivityOrder(amici.SensitivityOrder_second)
    edata = amici.ExpData(
        1, 1, 1, [model.getTimepoints()[-1]], [0], [1], [0], [1])

    def ob(theta_unif):

        theta = np.log(theta_unif / (1 - theta_unif))
        theta_dict = dict(zip(parameter_to_optimize, theta))
        for keys in theta_dict:
            model.setParameterByName(keys, theta_dict[keys])
        rdata = amici.runAmiciSimulation(model, solver, edata)
        func = rdata["y"][-1][-1]
        gradient = rdata["sllh"]
        hessian = np.linalg.inv(rdata["FIM"])

        return func, gradient, hessian

    objective = pypesto.Objective(fun=ob, grad=True, hess=True)
    # objective = pypesto.AmiciObjective(model, solver, [edata], max_sensi_order=2)
    upper_bound = 1 - lower_bound
    problem = pypesto.Problem(
        objective=objective,
        lb=np.repeat(lower_bound, len(parameter_to_optimize)),
        ub=np.repeat(upper_bound, len(parameter_to_optimize)),
    )
    optimizer = pypesto.optimize.FidesOptimizer()
    result = optimize.minimize(
        problem=problem, optimizer=optimizer, n_starts=50)
    return result


def run_optimization(
    model,
    solver,
    parameter_to_optimize,
    optimizer=pypesto.optimize.ScipyOptimizer("L-BFGS-B"),
    n_starts=50,
    lower_bound=10 ** (-5),
    trace_record=True,
    pareto_values=None,
):
    """Run optimization over specified parameters.

    Parameters
    ----------
    model : amici.model
        Amici model which should be evaluated.

    solver : amici.solver
        Solver object of the corresponding amici model.

    parameter_to_optimize : list
        List of parameters you want to optimize over.

    optimizer : string or pypesto.optimize object.
        If "pybobyqa" minimization is conducted using estimagic. Otherwise,
        the passed pypesto routine is used.

    n_starts : int
        Number of strats for the multi-start optimization.

    lower_bound : float
        Lower bound for the parameter estimation.

    trace_record : bool
        If True, trace record is set to True for pyPesto optimizers.

    pareto_values : None or array-like
        If an array is passed, its length must equal the number of areas
        within the model. The array is used to enforce maximum values each areas
        can have in the optimal solution.


    Returns
    -------
    results : pd.DataFrame
        List of results.
    """

    if optimizer == "pybobyqa":
        results = optimize_estimagic(
            model,
            solver,
            parameter_to_optimize,
            n_starts=n_starts,
            lower_bound=lower_bound,
            pareto_values=pareto_values,
        )
    elif optimizer == "fides":
        results = optimize_fides(
            model,
            solver,
            parameter_to_optimize,
            n_starts=n_starts,
            lower_bound=lower_bound,
        )
    else:
        results = optimize_pyPesto(
            model,
            solver,
            parameter_to_optimize,
            optimizer=optimizer,
            n_starts=n_starts,
            lower_bound=lower_bound,
            trace_record=trace_record,
            pareto_values=pareto_values,
        )

    return results


def get_pareto_front(
    model,
    solver,
    lb,
    number_areas,
    parameter_to_optimize,
    pop_size=100,
    number_generations=100,
    seed=1234,
    crossover_eta=15,
    crossover_prob=0.9,
    mutation_eta=20,
):
    """Compute the pareto front with respect to the areas using the NSGA-II
    algorithm.

    Parameters
    ----------
    model : amici.model
        Amici model which should be evaluated.

    solver : amici.solver
        Solver object of the corresponding amici model.

    lb : float
        Lower bound for the parameter estimation.

    parameter_to_optimize : list
        List of parameters you want to optimize over.

    pop_size : int
        The population sized used by the algorithm.

    number_generations : int
        Number of points for which Pareto front is computed.

    seed : int
        Seed for optimization.

    crossover_eta : bool
        See pymoo documentation.

    crossover_prob : float
        See pymoo documentation.

    mutation_eta : int
        See pymoo documentation.


    Returns
    -------
    res : pymoo.
        Results of the pareto optimization
    """

    areas = ["country" + chr(i) for i in range(65, 65 + number_areas)]
    number_par = len(parameter_to_optimize)
    ub = 1 - lb

    def run_model(theta):
        theta_trans = np.log(theta / (1 - theta))

        inputs = dict(zip(parameter_to_optimize, theta_trans))
        for i in inputs:
            model.setParameterByName(i, inputs[i])

        rdata = amici.runAmiciSimulation(model, solver)
        df_rdata = pd.DataFrame(rdata["x"], columns=model.getStateNames())
        death_countries = []
        for index in areas:
            death_countries.append(
                get_sum_of_states(
                    model, df_rdata, state_type=["dead", index], final_amount=True
                )
            )
        return death_countries

    # pareto_values = np.array(run_model(np.repeat(0.5, 22)))
    best_a = np.array(run_model(np.repeat(lb, number_par)))[0]
    best_b = np.array(run_model(np.repeat(ub, number_par)))[1]
    pareto_values = np.array([best_a, best_b])

    class multi_objective_problem(ElementwiseProblem):
        def __init__(self):
            super().__init__(
                n_var=number_par,
                n_obj=number_areas,
                n_constr=number_areas,
                xl=np.repeat(lb, number_par),
                xu=np.repeat(ub, number_par),
            )

        def _evaluate(self, x, out, *args, **kwargs):
            out["F"] = run_model(x)
            out["G"] = out["F"] - np.array(pareto_values)

    problem = multi_objective_problem()

    algorithm = NSGA2(
        pop_size=pop_size,
        sampling=get_sampling("real_random"),
        crossover=get_crossover(
            "real_sbx", prob=crossover_prob, eta=crossover_eta),
        mutation=get_mutation("real_pm", eta=mutation_eta),
        eliminate_duplicates=True,
    )

    termination = get_termination("n_gen", number_generations)

    res = pymoo.optimize.minimize(
        problem, algorithm, termination, seed=seed, save_history=True, verbose=True
    )
    return res


def to_be_vaccinated(model, area):
    """Get all states that can be vaccinated (by area).

    Parameters
    ----------
    model : amici.model
        Amici model which should be evaluated.

    areas : list
        List of area names as strings.

    Returns
    -------
    states : list
        List of states that can be vaccinated.
    """
    if area is None:
        states = [
            x
            for x in model.getStateNames()
            if not ("vac0" in x) and (("susceptible" in x) or ("infectious" in x))
        ]
    else:
        states = [
            x
            for x in model.getStateNames()
            if not ("vac0" in x)
            and (("susceptible" in x) or ("infectious" in x))
            and (area in x)
        ]

    return states


def get_vaccinated_model(model, area=None):
    """Get all states that can be vaccinated or recovered (by area).

    Parameters
    ----------
    model : amici.model
        Amici model which should be evaluated.

    areas : list
        List of area names as strings.

    Returns
    -------
    states : list
        List of states that can be vaccinated.
    """
    if area is None:
        states = [
            x
            for x in model.getStateNames()
            if not ("vac0" in x)
            and (("susceptible" in x) or ("infectious" in x))
            or ("recovered" in x)
        ]
    else:
        states = [
            x
            for x in model.getStateNames()
            if (
                not ("vac0" in x)
                and (("susceptible" in x) or ("infectious" in x))
                or ("recovered" in x)
            )
            and (area in x)
        ]

    return states


def get_alive_model(model, area=None):
    """Get all states that for which individuals are alive (by area).

    Parameters
    ----------
    model : amici.model
        Amici model which should be evaluated.

    areas : list
        List of area names as strings.

    Returns
    -------
    states : list
        List of states for which individuals are alive.
    """
    if area is None:
        states = [
            x
            for x in model.getStateNames()
            if (("susceptible" in x) or ("infectious" in x) or ("recovered" in x))
        ]
    else:
        states = [
            x
            for x in model.getStateNames()
            if (("susceptible" in x) or ("infectious" in x) or ("recovered" in x))
            and (area in x)
        ]

    return states


def get_fraction_vaccinated(model, trajectories, area=None, include_recovered=True):
    """Get fraction of individuals that are vaccinated or immune (by area) by
    state.

    Parameters
    ----------
    model : amici.model
        Amici model which should be evaluated.

    trajectories : pd.DataFrame
        Trajectories of the model simulation.

    areas : list
        List of area names as strings.

    include_recovered : bool
        If True, recovered individuals are counted as well.

    Returns
    -------
    percentage_vaccinated: pd.Series
        Trajectories of the fraction that is vaccinated or immune.
    """
    vaccinated = get_vaccinated_model(model, area=area)
    sus_inf = get_alive_model(model, area=area)

    df_vaccinated = trajectories[vaccinated]
    df_sus_inf = trajectories[sus_inf]

    total_vaccinated = df_vaccinated.sum(axis=1)
    sus_inf = vaccinated = df_sus_inf.sum(axis=1)

    percentage_vaccinated = total_vaccinated / sus_inf

    return percentage_vaccinated


def get_sum_vaccinated(model, trajectories, area=None, include_recovered=True):
    """Get fraction of individuals that are vaccinated or immune (by area).

    Parameters
    ----------
    model : amici.model
        Amici model which should be evaluated.

    trajectories : pd.DataFrame
        Trajectories of the model simulation.

    areas : list
        List of area names as strings.

    include_recovered : bool
        If True, recovered individuals are counted as well.

    Returns
    -------
    percentage_vaccinated: pd.Series
        Trajectory of the fraction that is vaccinated or immune.
    """
    unvaccinated = to_be_vaccinated(model, area)
    df_unvaccinated = trajectories[unvaccinated]

    return df_unvaccinated.sum(axis=1)


def get_start_splines_pop_based(model, areas, number_yy):
    """Modify model such that spline parameters match population size based
    weights.

    Parameters
    ----------
    model : amici.model
        Amici model which should be evaluated.

    areas : list
        List of area names as strings.

    number_yy: int
        Number of parameters per area.

    Returns
    -------
    transformed_fractions : np.array
        Start parameters.
    """
    fixed_parameters = dict(
        zip(model.getFixedParameterNames(), model.getFixedParameters())
    )
    pop_t0 = []
    for area in areas:
        s = 0
        for keys in fixed_parameters.keys():
            if "t0" in keys and area in keys:
                s += fixed_parameters[keys]

        pop_t0.append(s)
    fractions = pop_t0 / np.sum(pop_t0)
    transformed_fractions = np.log(fractions / (1 - fractions))

    dict_start = dict(
        zip(model.getParameterNames(), np.repeat(
            transformed_fractions, 2 * number_yy))
    )
    for keys in dict_start.keys():
        model.setParameterByName(keys, dict_start[keys])

    return dict_start


def transform_pybobyqa_results(model, solver, results_pybobyqa, areas):
    """Modify pybobyqa results such that areawise death numbers are added.

    Parameters
    ----------
    model : amici.model
        Amici model which should be evaluated.

    solver : amici.solver
        Solver object of the corresponding amici model.

    results_pybobyqa : pd.DataFrame
        Results of the oybobyqa run.


    areas : list
        List of area names as strings.

    include_recovered : bool
        If True, recovered individuals are counted as well.

    Returns
    -------
    results : pd.dataFrame
        pybobyqa results with added areawise death numbers.
    """
    df_strategies = results_pybobyqa.drop(columns=["fval", "feval"])
    results = results_pybobyqa.reindex(
        columns=list(results_pybobyqa.columns) + list(areas)
    )
    ncol = df_strategies.shape[0]
    for i in range(ncol):
        pars = df_strategies.iloc[i]
        dict_pars = dict(zip(df_strategies.columns, pars))
        for keys in dict_pars.keys():
            model.setParameterByName(keys, dict_pars[keys])
        rdata = amici.runAmiciSimulation(model, solver)
        states = pd.DataFrame(rdata["x"], columns=model.getStateIds())
        for area in areas:
            dead = get_sum_of_states(
                model, states, state_type=["dead", area], final_amount=True
            )
            results.xs(i)[area] = dead

    return results


def get_fraction_vaccinated_pybobyqa(
    re_pybobyqa, model, solver, par_to_optimize, areas, include_recovered=True
):
    """Get vaccination fractions from optimal bobyqa solution.

    Parameters
    ----------
    re_pybobyqa : pd.DataFrame
        Results of bobyqa simulation.

    model : amici.model
        Amici model which should be evaluated.

    solver : amici.solver
        Solver object of the corresponding amici model.

    par_to_optimize : list
        List of parameters to optimize over.

    areas : list
        List of area names as strings.

    Returns
    -------
    results : pd.dataFrame
        pybobyqa results with added areawise death numbers.
    """

    optimal_strategy_unif = re_pybobyqa[list(par_to_optimize)].iloc[
        np.argmin(re_pybobyqa["fval"])
    ]
    optimal_strategy = np.log(
        optimal_strategy_unif / (1 - optimal_strategy_unif))
    for keys in optimal_strategy.index:
        model.setParameterByName(keys, optimal_strategy[keys])
    rdata = amici.runAmiciSimulation(model, solver)
    states = pd.DataFrame(rdata["x"], columns=model.getStateIds())
    fraction_vaccinated = pd.DataFrame(
        np.nan, columns=areas, index=range(states.shape[0])
    )

    for area in areas:
        fraction_vaccinated[area] = get_fraction_vaccinated(
            model=model,
            trajectories=states,
            area=area,
            include_recovered=include_recovered,
        )
    fraction_vaccinated["total"] = get_fraction_vaccinated(
        model=model, trajectories=states, include_recovered=include_recovered
    )

    return fraction_vaccinated, states


def get_model_output(
    model,
    solver,
    parameters,
    areas,
    par_to_optimize,
    lb=10 ** -15,
    n_starts_pb=50,
    n_starts_pareto=100,
    number_generations_pareto=100,
    seed_pareto=1234,
    crossover_eta=15,
    crossover_prob=0.9,
    mutation_eta=20,
    seed_pb=12345,
):
    """Get model simulations and output.

    Parameters
    ----------
    model : amici.model
        Amici model which should be evaluated.

    solver : amici.solver
        Solver object of the corresponding amici model.

    parameters : dict
        All aprameters as keys and their magnitude as value.

    par_to_optimize : list
        List of parameters to optimize over.

    lb : float
        Lower bound for the parameter estimation.

    n_starts_pb : int
        Number of starts for pybobyqa.

    n_starts_pareto : int
        Number of starts for pareto optimum.

    number_generations_pareto : int
        Number of points for which Pareto front is computed.

    seed_pareto : int
        Seed for Pareto optimization.

    crossover_eta : bool
        See pymoo documentation.

    crossover_prob : float
        See pymoo documentation.

    mutation_eta : int
        See pymoo documentation.

    seed_pb : int
        Seed for pybobyqa optimization.

    Returns
    -------
    results : dict
        Resultstrance
    """

    for keys in parameters.keys():
        model.setFixedParameterByName(keys, parameters[keys])

    number_yy = len([x for x in model.getParameterNames() if "yy_" in x]) / 2
    get_start_splines_pop_based(model, areas=areas, number_yy=number_yy)

    # compute parameters
    rdata = amici.runAmiciSimulation(model, solver)
    trajectories = pd.DataFrame(rdata["x"], columns=model.getStateNames())

    pareto = []
    for index in areas:
        pareto.append(
            get_sum_of_states(
                model, trajectories, state_type=["dead", index], final_amount=True
            )
        )

    results_pareto_front = get_pareto_front(
        model,
        solver,
        parameter_to_optimize=par_to_optimize,
        lb=lb,
        number_areas=len(areas),
        pop_size=n_starts_pareto,
        number_generations=number_generations_pareto,
        seed=seed_pareto,
        crossover_eta=crossover_eta,
        crossover_prob=crossover_prob,
        mutation_eta=mutation_eta,
    )

    np.random.seed(seed_pb)
    results_pybobyqa = run_optimization(
        model,
        solver,
        par_to_optimize,
        optimizer="pybobyqa",
        n_starts=n_starts_pb,
        lower_bound=lb,
        trace_record=True,
        pareto_values=None,
    )
    re_pybobyqa = transform_pybobyqa_results(
        model, solver, results_pybobyqa, areas)

    add_row = dict(
        zip(
            list(par_to_optimize) + areas + ["fval"],
            list(np.repeat(0.5, len(par_to_optimize))) +
            pareto + [np.sum(pareto)],
        )
    )

    if not (results_pareto_front.F is None):
        pareto_strategies = pd.DataFrame(
            results_pareto_front.X, columns=par_to_optimize
        )
        for i in range(len(areas)):
            pareto_strategies[areas[i]] = results_pareto_front.F[:, i]
        pareto_strategies["fval"] = results_pareto_front.F.sum(axis=1)
    else:
        pareto_strategies = pd.DataFrame(
            np.matrix(list(add_row.values())), columns=add_row.keys(), index=[0]
        )

    appended_df = re_pybobyqa.append(pareto_strategies)

    appended_df = appended_df.append(add_row, ignore_index=True)

    pareto_optimal = appended_df[appended_df["countryA"] <= pareto[0]]
    for i in range(1, len(areas)):
        pareto_optimal = pareto_optimal[appended_df[areas[i]] <= pareto[i]]

    vaccinated_best_strategy, states_best = get_fraction_vaccinated_pybobyqa(
        appended_df,
        model,
        solver,
        par_to_optimize,
        areas,
    )

    vaccinated_pareto_strategy, states_pareto = get_fraction_vaccinated_pybobyqa(
        pareto_optimal,
        model,
        solver,
        par_to_optimize,
        areas,
    )

    out = {
        "optimal_strategies": re_pybobyqa,
        "optimal_frac_vaccinated": vaccinated_best_strategy,
        "pareto_frac_vaccinated": vaccinated_pareto_strategy,
        "pareto_frontier": pareto_strategies,
        "population_based": pareto,
        "pareto_improvements": pareto_optimal,
        "all_strategies": appended_df,
        "trajectories_best": states_best,
        "trajectories_pareto": states_pareto,
        "trajectories_pop": trajectories,
    }

    return out


# ------------------------------------------------------------------------------
def set_parameter_by_name(model, parameters, fixed=False):
    if fixed is True:
        for keys in parameters.keys():
            model.setFixedParameterByName(keys, float(parameters[keys]))
    elif fixed is False:
        for keys in parameters.keys():
            model.setParameterByName(keys, float(parameters[keys]))
    else:
        raise ValueError("fixed must be boolean")


def set_start_states_default(model, parameters, T, n_grid):
    states = model.getStateNames()
    dict_states = dict(zip(states, np.repeat(0.0, len(states))))
    for keys in dict_states.keys():
        model.setFixedParameterByName(keys + "_t0", dict_states[keys])

    for keys in parameters.keys():
        model.setFixedParameterByName(keys, parameters[keys])
    model.setT0(0)
    timepoints = np.linspace(0, T, n_grid)
    model.setTimepoints(timepoints)


def set_end_values_to_start(rdata, model, t0, T, n_grid):
    trajectories = pd.DataFrame(rdata["x"], columns=model.getStateNames())
    end_trajectories = trajectories.iloc[-1]

    for i in end_trajectories.index:
        model.setFixedParameterByName(i + "_t0", end_trajectories[i])
    # model.setFixedParameterByName("infectious_countryB_vac0_virus2_t0",
    #                              model.getFixedParameterByName('virus2_countryB_appears_quantity'))
    model.setT0(t0)
    timepoints = np.linspace(t0, T, n_grid)
    model.setTimepoints(timepoints)


def simulate_intervention(
    model, solver, parameters, T, n_grid, areas, par_optim_dict=None
):
    if not (par_optim_dict is None):
        set_parameter_by_name(model, par_optim_dict, fixed=False)

    set_parameter_by_name(model, parameters, fixed=True)
    appears_day = model.getFixedParameterByName("virus2_countryB_appears_time")

    set_start_states_default(
        model, parameters, appears_day, int(appears_day / T * n_grid)
    )
    rdata1 = amici.runAmiciSimulation(model, solver)
    trajectories1 = pd.DataFrame(rdata1["x"], columns=model.getStateNames())

    set_end_values_to_start(
        rdata1, model, appears_day, T, int((T - appears_day) / T * n_grid)
    )

    if len(areas) == 2:
        model.setFixedParameterByName(
            "infectious_countryB_vac0_virus2_t0",
            model.getFixedParameterByName("virus2_countryB_appears_quantity"),
        )

    rdata2 = amici.runAmiciSimulation(model, solver)
    trajectories2 = pd.DataFrame(rdata2["x"], columns=model.getStateNames())

    deaths = [
        get_sum_of_states(model, trajectories2, state_type=[
                          "dead"], final_amount=True)
    ]

    for i in areas:
        deaths.append(
            get_sum_of_states(
                model, trajectories2, state_type=["dead", i], final_amount=True
            )
        )
    trajectory = trajectories1.append(trajectories2)

    return {"deaths": deaths, "trajectories": trajectory}


def optimize_intervention_pybobyqa(
    model,
    solver,
    parameters,
    T,
    n_grid,
    areas,
    week_of_reaction,
    par_to_optimize,
    number_yy,
    lower_bound=10 ** -15,
    n_starts=10,
):

    set_parameter_by_name(model, parameters, fixed=True)

    par_optimize = [
        x for x in par_to_optimize if float(x.rsplit("_", 1)[-1]) > week_of_reaction / 2
    ]
    par_rest = [
        x
        for x in par_to_optimize
        if float(x.rsplit("_", 1)[-1]) <= week_of_reaction / 2
    ]
    transformed_frac_pareto = get_start_splines_pop_based(
        model, areas, number_yy)

    dict_par_rest = {}
    for i in par_rest:
        dict_par_rest[i] = float(transformed_frac_pareto[i])

    def get_optim_function(theta):
        par_optim_dict = dict(
            zip(par_optimize, np.log(theta["value"] / (1 - theta["value"])))
        )
        deaths = simulate_intervention(
            model, solver, parameters, T, n_grid, areas, par_optim_dict=par_optim_dict
        )

        return {"value": deaths["deaths"][0]}

    start_par = np.random.uniform(
        lower_bound, 1 - lower_bound, len(par_optimize) * n_starts
    ).reshape(n_starts, len(par_optimize))

    save_df = pd.DataFrame(
        np.NaN,
        columns=list(dict_par_rest.keys()) +
        list(par_optimize) + ["fval", "feval"],
        index=range(n_starts),
    )
    for i in range(n_starts):

        print(f"{i}/{n_starts} start.")

        start_theta = start_par[i]
        parameter_np = np.array(
            [
                start_theta,
                np.repeat(lower_bound, len(par_optimize)),
                np.repeat(1 - lower_bound, len(par_optimize)),
            ]
        )
        parameter_df = pd.DataFrame(
            parameter_np.transpose(),
            columns=["value", "lower_bound", "upper_bound"],
            index=par_optimize,
            dtype="float",
        )
        r = minimize(
            criterion=get_optim_function,
            params=parameter_df,
            algorithm="nag_pybobyqa",
        )
        save_df.iloc[i] = (
            list(1 / (1 + np.exp(np.array(list(dict_par_rest.values())))))
            + list(r["solution_x"])
            + [
                r["solution_criterion"],
                r["n_criterion_evaluations"],
            ]
        )

    return save_df


def trannsform_pybobyqa_intervention(df, areas, model, solver, parameters, T, n_grid):
    df[areas] = np.repeat(np.nan, len(areas) * df.shape[0]).reshape(
        df.shape[0], len(areas)
    )
    par = [x for x in df.columns if "yy_" in x]
    for i in range(df.shape[0]):
        theta = np.log(df.iloc[i][par] / (1 - df.iloc[i][par]))
        dict_theta = dict(zip(par, theta))
        out = simulate_intervention(
            model, solver, parameters, T, n_grid, areas, par_optim_dict=dict_theta
        )
        deaths = out["deaths"]
        df.loc[i, areas] = deaths[1: len(deaths)]

    return df


def get_pareto_front_intervention(
    model,
    solver,
    lb,
    areas,
    week_of_reaction,
    T,
    n_grid,
    parameters,
    number_yy,
    par_to_optimize,
    pop_size=100,
    number_generations=100,
    seed=1234,
    crossover_eta=15,
    crossover_prob=0.9,
    mutation_eta=20,
):
    """Compute the pareto front with respect to the areas using the NSGA-II
    algorithm.

    Parameters
    ----------
    model : amici.model
        Amici model which should be evaluated.

    solver : amici.solver
        Solver object of the corresponding amici model.

    lb : float
        Lower bound for the parameter estimation.

    parameter_to_optimize : list
        List of parameters you want to optimize over.

    pop_size : int
        The population sized used by the algorithm.

    number_generations : int
        Number of points for which Pareto front is computed.

    seed : int
        Seed for optimization.

    crossover_eta : bool
        See pymoo documentation.

    crossover_prob : float
        See pymoo documentation.

    mutation_eta : int
        See pymoo documentation.


    Returns
    -------
    res : pymoo.
        Results of the pareto optimization
    """

    set_parameter_by_name(model, parameters, fixed=True)

    par_optimize = [
        x for x in par_to_optimize if float(x.rsplit("_", 1)[-1]) > week_of_reaction / 2
    ]
    par_rest = [
        x
        for x in par_to_optimize
        if float(x.rsplit("_", 1)[-1]) <= week_of_reaction / 2
    ]
    transformed_frac_pareto = get_start_splines_pop_based(
        model, areas, number_yy)
    number_par = len(par_optimize)
    number_areas = len(areas)

    dict_par_rest = {}
    for i in par_rest:
        dict_par_rest[i] = float(transformed_frac_pareto[i])

    def get_optim_function_pareto(theta):

        par_optim_dict = dict(zip(par_optimize, np.log(theta / (1 - theta))))
        deaths = simulate_intervention(
            model, solver, parameters, T, n_grid, areas, par_optim_dict=par_optim_dict
        )
        d = deaths["deaths"]
        return d[1: len(d)]

    ub = 1 - lb
    # pareto_values = np.array(run_model(np.repeat(0.5, 22)))
    best_a = np.array(get_optim_function_pareto(np.repeat(lb, number_par)))[0]
    best_b = np.array(get_optim_function_pareto(np.repeat(ub, number_par)))[1]
    pareto_values = np.array([best_a, best_b])

    class multi_objective_problem(ElementwiseProblem):
        def __init__(self):
            super().__init__(
                n_var=number_par,
                n_obj=number_areas,
                n_constr=number_areas,
                xl=np.repeat(lb, number_par),
                xu=np.repeat(ub, number_par),
            )

        def _evaluate(self, x, out, *args, **kwargs):
            out["F"] = get_optim_function_pareto(x)
            out["G"] = out["F"] - np.array(pareto_values)

    problem = multi_objective_problem()

    algorithm = NSGA2(
        pop_size=pop_size,
        sampling=get_sampling("real_random"),
        crossover=get_crossover(
            "real_sbx", prob=crossover_prob, eta=crossover_eta),
        mutation=get_mutation("real_pm", eta=mutation_eta),
        eliminate_duplicates=True,
    )

    termination = get_termination("n_gen", number_generations)

    res = pymoo.optimize.minimize(
        problem, algorithm, termination, seed=seed, save_history=True, verbose=True
    )
    return res


def get_frac_vaccinated_intervention(model, states, areas, include_recovered=True):
    fraction_vaccinated = pd.DataFrame(
        np.nan, columns=areas, index=range(states.shape[0])
    )

    for area in areas:
        fraction_vaccinated[area] = get_fraction_vaccinated(
            model=model,
            trajectories=states,
            area=area,
            include_recovered=include_recovered,
        )
    fraction_vaccinated["total"] = get_fraction_vaccinated(
        model=model, trajectories=states, include_recovered=include_recovered
    )

    return fraction_vaccinated


def run_intervention(
    model,
    solver,
    par_to_optimize,
    areas,
    parameters,
    week_of_reaction=12,
    T=140,
    n_grid=6000,
    number_yy=11,
    lb=10 ** -15,
    starts_bobyqa=10,
    pop_size=2,
    number_generations=2,
    p_seed=12345,
    pareto_seed=1234,
    crossover_eta=15,
    crossover_prob=0.9,
    mutation_eta=20,
):

    np.random.seed(p_seed)
    pybobyqa_results = optimize_intervention_pybobyqa(
        model,
        solver,
        parameters,
        T,
        n_grid,
        areas,
        week_of_reaction,
        par_to_optimize,
        number_yy,
        lower_bound=lb,
        n_starts=starts_bobyqa,
    )

    optimal_results = trannsform_pybobyqa_intervention(
        pybobyqa_results, areas, model, solver, parameters, T, n_grid
    )

    results_pareto_front = get_pareto_front_intervention(
        model,
        solver,
        lb,
        areas,
        week_of_reaction,
        T,
        n_grid,
        parameters,
        number_yy,
        par_to_optimize,
        pop_size=pop_size,
        number_generations=number_generations,
        seed=pareto_seed,
        crossover_eta=crossover_eta,
        crossover_prob=crossover_prob,
        mutation_eta=mutation_eta,
    )

    par_dict = dict(zip(par_to_optimize, np.repeat(0, len(par_to_optimize))))
    dict_pop = simulate_intervention(
        model, solver, parameters, T, n_grid, areas, par_optim_dict=par_dict
    )

    pareto = dict_pop["deaths"][1: len(dict_pop["deaths"])]
    add_row = dict(
        zip(
            list(par_to_optimize) + areas + ["fval"],
            list(np.repeat(0.5, len(par_to_optimize))) +
            pareto + [np.sum(pareto)],
        )
    )

    par_optimize = [
        x for x in par_to_optimize if float(x.rsplit("_", 1)[-1]) > week_of_reaction / 2
    ]
    par_rest = [
        x
        for x in par_to_optimize
        if float(x.rsplit("_", 1)[-1]) <= week_of_reaction / 2
    ]

    if not (results_pareto_front.F is None):
        pareto_strategies = pd.DataFrame(
            results_pareto_front.X, columns=par_optimize)
        pareto_strategies[par_rest] = np.repeat(
            0.5, len(par_rest) * pareto_strategies.shape[0]
        ).reshape(pareto_strategies.shape[0], len(par_rest))
        for i in range(len(areas)):
            pareto_strategies[areas[i]] = results_pareto_front.F[:, i]
            pareto_strategies["fval"] = results_pareto_front.F.sum(axis=1)
        pareto_strategies = pareto_strategies.append(
            add_row, ignore_index=True)
    else:
        pareto_strategies = pd.DataFrame(
            np.matrix(list(add_row.values())), columns=add_row.keys(), index=[0]
        )

    complete_df = optimal_results.append(pareto_strategies).reset_index()

    pareto_optimal = complete_df[complete_df["countryA"] <= pareto[0]]
    for i in range(1, len(areas)):
        pareto_optimal = pareto_optimal[pareto_optimal[areas[i]] <= pareto[i]]

    pareto_optimal = pareto_optimal.reset_index()

    optimal_vac = complete_df.iloc[np.argmin(
        complete_df["fval"])][par_to_optimize]
    pareto_vac = pareto_optimal.iloc[np.argmin(
        pareto_optimal["fval"])][par_to_optimize]

    optimal_out = simulate_intervention(
        model,
        solver,
        parameters,
        T,
        n_grid,
        areas,
        par_optim_dict=dict(zip(par_to_optimize, optimal_vac)),
    )
    pareto_out = simulate_intervention(
        model,
        solver,
        parameters,
        T,
        n_grid,
        areas,
        par_optim_dict=dict(zip(par_to_optimize, pareto_vac)),
    )

    optimal_frac_vac = get_frac_vaccinated_intervention(
        model,
        states=optimal_out["trajectories"].reset_index(),
        areas=areas,
        include_recovered=True,
    )
    pareto_frac_vac = get_frac_vaccinated_intervention(
        model,
        states=pareto_out["trajectories"].reset_index(),
        areas=areas,
        include_recovered=True,
    )

    out = {
        "optimal_strategies": pybobyqa_results,
        "optimal_frac_vaccinated": optimal_frac_vac,
        "pareto_frac_vaccinated": pareto_frac_vac,
        "pareto_frontier": pareto_strategies,
        "population_based": pareto,
        "pareto_improvements": pareto_optimal,
        "all_strategies": complete_df,
    }

    return out


# ------------------------------COVID19 pandemic--------------------------------
def run_interventions_model(
    model,
    solver,
    interventions,
    parameters,
    total_grid,
    length,
    areas,
    number_yy,
    theta=None,
):

    df_trajectories = pd.DataFrame(columns=(model.getStateNames()))
    set_start_states_default(
        model,
        parameters,
        interventions[0]["t"],
        int(total_grid * interventions[0]["t"] / length),
    )
    if theta is None:
        get_start_splines_pop_based(
            model=model, areas=areas, number_yy=number_yy)
    for i in range(len(interventions)):

        rdata = amici.runAmiciSimulation(model, solver)
        trajectories = pd.DataFrame(rdata["x"], columns=model.getStateNames())
        df_trajectories = df_trajectories.append(trajectories)
        t0 = interventions[i]["t"]

        if i != (len(interventions) - 1):
            T = interventions[i + 1]["t"]
        elif i == (len(interventions) - 1):
            T = length
        n_grid = int((T - t0) * total_grid / length)
        set_end_values_to_start(rdata, model, t0, T, n_grid)
        set_parameter_by_name(model, interventions[i]["parameter"], fixed=True)

    rdata = amici.runAmiciSimulation(model, solver)
    trajectories = pd.DataFrame(rdata["x"], columns=model.getStateNames())
    df_trajectories = df_trajectories.append(trajectories)

    deaths = [
        get_sum_of_states(
            model, df_trajectories, state_type=["dead"], final_amount=True
        )
    ]

    for i in areas:
        deaths.append(
            get_sum_of_states(
                model, df_trajectories, state_type=["dead", i], final_amount=True
            )
        )
    set_start_states_default(
        model,
        parameters,
        interventions[0]["t"],
        int(total_grid * interventions[0]["t"] / length),
    )

    return {"deaths": deaths, "trajectories": df_trajectories}


def run_interventions_pybobyqa(
    model,
    solver,
    interventions,
    parameters,
    total_grid,
    length,
    areas,
    number_yy,
    par_optimize,
    lower_bound=10 ** -15,
    n_starts=10,
):

    set_start_states_default(
        model,
        parameters,
        interventions[0]["t"],
        int(total_grid * interventions[0]["t"] / length),
    )
    get_start_splines_pop_based(model=model, areas=areas, number_yy=number_yy)

    def get_optim_function(theta):
        theta_trans = np.log(theta["value"] / (1 - theta["value"]))
        dict_theta = dict(zip(par_optimize, theta_trans))

        set_parameter_by_name(model, dict_theta, fixed=False)

        deaths = run_interventions_model(
            model,
            solver,
            interventions,
            parameters,
            total_grid,
            length,
            areas,
            number_yy,
            theta_trans,
        )

        return {"value": deaths["deaths"][0]}

    start_par = np.random.uniform(
        lower_bound, 1 - lower_bound, number_yy * n_starts
    ).reshape(number_yy, n_starts)

    save_df = pd.DataFrame(
        np.NaN,
        columns=list(par_optimize) + ["fval", "feval"],
        index=range(n_starts),
    )
    for i in range(n_starts):

        print(f"{i}/{n_starts} start.")

        start_theta = start_par[i]
        parameter_np = np.array(
            [
                start_theta,
                np.repeat(lower_bound, len(par_optimize)),
                np.repeat(1 - lower_bound, len(par_optimize)),
            ]
        )
        parameter_df = pd.DataFrame(
            parameter_np.transpose(),
            columns=["value", "lower_bound", "upper_bound"],
            index=par_optimize,
            dtype="float",
        )
        r = minimize(
            criterion=get_optim_function,
            params=parameter_df,
            algorithm="nag_pybobyqa",
        )
        save_df.iloc[i] = list(r["solution_x"]) + [
            r["solution_criterion"],
            r["n_criterion_evaluations"],
        ]

    return save_df


def get_infected_true():
    belgium_infected = [
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        1,
        7,
        14,
        22,
        49,
        108,
        168,
        199,
        238,
        266,
        312,
        394,
        553,
        680,
        872,
        1025,
        1182,
        1281,
        1529,
        1899,
        2421,
        2829,
        3064,
        3426,
        3890,
        4898,
        5640,
        7066,
        8459,
        9204,
        9474,
        10092,
        10867,
        11695,
        12616,
        13371,
        14077,
        14676,
        15125,
        15980,
        16975,
        17547,
        18647,
        19163,
        19177,
        20915,
        21472,
        22168,
        22599,
        23566,
        24756,
        25092,
        25455,
        25860,
        26903,
        27410,
        27983,
        28317,
        28489,
        28611,
        28848,
        28955,
        29248,
        29493,
        29697,
        29554,
        29495,
        29226,
        28608,
        27613,
        26414,
        25431,
        24132,
        23392,
        23169,
        21071,
        20180,
        19142,
        18376,
        17295,
        16000,
        15279,
        14622,
        14013,
        12799,
        12017,
        11321,
        10905,
        10515,
        10202,
        9667,
        9349,
        9000,
        8709,
        8418,
        8258,
        8126,
        7652,
        7215,
        6752,
        6356,
        6120,
        5932,
        5838,
        5630,
        5385,
        5111,
        4875,
        4685,
        4557,
        4493,
        4315,
        4039,
        3740,
        3718,
        3556,
        3552,
        3514,
        3360,
        3234,
        3175,
        3046,
        2992,
        2983,
        3042,
        3071,
        3002,
        2944,
        2832,
        2775,
        2773,
        2788,
        2758,
        2787,
        2789,
        2752,
        2772,
        2884,
        2994,
        3151,
        3230,
        3343,
        3544,
        3708,
        3817,
        3949,
        4192,
        4621,
        4817,
        5133,
        5301,
        5908,
        6497,
        7153,
        7675,
        8011,
        8405,
        8632,
        9100,
        9893,
        10574,
        11044,
        11683,
        12014,
        12301,
        12866,
        13319,
        14074,
        14631,
        14824,
        14828,
        15004,
        15385,
        15920,
        16267,
        16621,
        16737,
        16365,
        16421,
        16602,
        16838,
        16617,
        16593,
        16291,
        15834,
        15638,
        15597,
        15896,
        16016,
        15809,
        15583,
        15368,
        14989,
        15071,
        15560,
        15890,
        16287,
        16342,
        16437,
        16472,
        17414,
        19079,
        20170,
        20570,
        21401,
        21924,
        23290,
        24795,
        26321,
        27946,
        29303,
        30227,
        30754,
        32073,
        33216,
        35572,
        38323,
        41079,
        43061,
        44378,
        45924,
        49099,
        54455,
        59290,
        66363,
        70721,
        73402,
        79785,
        87205,
        97164,
        106203,
        115139,
        122604,
        129732,
        137864,
        149994,
        164906,
        180813,
        196641,
        210055,
        220915,
        233110,
        252984,
        275143,
        293862,
        308170,
        316784,
        319732,
        323075,
        336010,
        345050,
        350176,
        350572,
        351808,
        346251,
        345217,
        349511,
        347153,
        343501,
        339321,
        333788,
        324756,
        318352,
        315307,
        310105,
        300294,
        286772,
        271079,
        254493,
        240772,
        231249,
        220243,
        202492,
        182190,
        164285,
        148116,
        138194,
        134897,
        131547,
        119226,
        110601,
        103712,
        98447,
        93783,
        94461,
        92922,
        87768,
        85573,
        83125,
        77931,
        75483,
        77187,
        77599,
        75252,
        73496,
        72250,
        70007,
        70330,
        72419,
        73101,
        72279,
        70498,
        68048,
        65286,
        64812,
        66897,
        67284,
        66037,
        64312,
        62572,
        60945,
        60979,
        63117,
        64083,
        63060,
        62297,
        61104,
        59257,
        59185,
        61038,
        61464,
        60031,
        59005,
        57800,
        56011,
        55320,
        57345,
        57718,
        56950,
        56950,
        56420,
        56828,
        57765,
        59928,
        61026,
        60878,
        61341,
        61864,
        62248,
        63260,
        65508,
        66112,
        65700,
        65215,
        64907,
        63789,
        64071,
        65728,
        65778,
        64971,
        64229,
        63542,
        62279,
        62366,
        64111,
        65052,
        65483,
        65552,
        65202,
        63740,
        64030,
        65951,
        67243,
        66992,
        66977,
        66391,
        64457,
        64191,
        66191,
        66980,
        66615,
        66962,
        66459,
        65138,
        65561,
        68122,
        69774,
        70772,
        72221,
        73063,
        72746,
        74395,
        78654,
        81596,
        84059,
        86449,
        87267,
        86859,
        88302,
        93496,
        96912,
        99922,
        102178,
        104103,
        103522,
        105331,
        110159,
        113576,
        114977,
        117223,
        117986,
        117155,
        116782,
        119204,
        121886,
        123312,
        124379,
        123753,
        121908,
        122320,
        125418,
        126250,
        125026,
        124162,
        122026,
        123650,
        117998,
        120898,
        125246,
        125287,
        122915,
        113653,
        109984,
        108238,
        109987,
        113605,
        107731,
        110908,
        102944,
        99391,
        98586,
        102666,
        102553,
        106237,
        106111,
        103512,
        104747,
        97819,
        101049,
        95855,
        97467,
        94733,
        96562,
        85895,
        87688,
        86869,
        90035,
        92926,
        86945,
        84336,
        77770,
        78605,
        76548,
        79426,
        81778,
        76793,
        74559,
        71757,
        67924,
        69769,
        70048,
        64641,
        66104,
        59801,
        60262,
        57814,
        56739,
        56321,
        57098,
        49286,
        48227,
        48463,
        48965,
        45824,
        45349,
        43084,
        40723,
        37705,
        37925,
        32904,
        33925,
        32083,
        31791,
        29705,
        27041,
        24486,
        23348,
        23127,
        23874,
        21145,
        20680,
        20199,
        19952,
        20690,
        20291,
        21206,
        20340,
        21685,
        22768,
        21257,
        22558,
        24037,
        25489,
        26832,
        27967,
        29295,
        30623,
        29260,
        30993,
        32122,
        32401,
        33922,
        35288,
        34397,
        35728,
        37786,
        39644,
        39578,
        40975,
        42517,
        43875,
        43060,
        45229,
        45796,
        47611,
        48315,
        47531,
        47673,
        49416,
        52022,
        50823,
        52835,
        52595,
        52252,
        51454,
        52072,
        53482,
        52876,
        55036,
        55001,
        54901,
        55116,
        56726,
        55405,
        57866,
        58873,
        58163,
        57750,
        56684,
        57004,
        58486,
        60680,
        61245,
        61489,
        61355,
        60370,
        60277,
        60948,
        61575,
        62104,
        61595,
        60222,
        60468,
        60413,
        61646,
        62586,
        60445,
        60345,
        60110,
        59540,
        61330,
        60625,
        62902,
        63499,
        63620,
        61840,
        60435,
        60910,
        61855,
        62137,
        59681,
        59920,
        59787,
        59829,
        58091,
        60686,
        62016,
        60036,
        60147,
        59688,
        59320,
        61742,
        64090,
    ]
    france_infected = [
        7,
        7,
        6,
        4,
        4,
        4,
        4,
        1,
        1,
        0,
        2,
        5,
        24,
        42,
        80,
        107,
        162,
        181,
        248,
        371,
        578,
        840,
        1073,
        1247,
        1582,
        2019,
        2546,
        3242,
        3991,
        4795,
        5874,
        6254,
        7442,
        8332,
        9431,
        10997,
        11689,
        14991,
        15894,
        17706,
        19857,
        22268,
        26140,
        26708,
        29542,
        34414,
        36833,
        35909,
        37965,
        39360,
        39801,
        41457,
        41386,
        42455,
        43059,
        44293,
        45029,
        45139,
        46461,
        49609,
        48415,
        48206,
        46206,
        46337,
        46061,
        46548,
        46669,
        46311,
        45867,
        45684,
        45611,
        45479,
        47853,
        47494,
        47401,
        46925,
        46534,
        46758,
        46688,
        46332,
        45666,
        48050,
        47403,
        46943,
        43060,
        38913,
        36224,
        35270,
        33070,
        28153,
        26044,
        23827,
        23716,
        21683,
        21051,
        19389,
        17018,
        15557,
        14137,
        12471,
        11473,
        7832,
        7424,
        9924,
        10077,
        11585,
        11399,
        11055,
        10865,
        10615,
        10811,
        11037,
        11468,
        11331,
        11040,
        11029,
        11244,
        11277,
        11725,
        12112,
        12160,
        12046,
        12191,
        9523,
        9405,
        8463,
        8826,
        8778,
        8845,
        9000,
        8356,
        7783,
        8730,
        8898,
        9190,
        9075,
        9072,
        9534,
        9471,
        9523,
        9667,
        10124,
        9966,
        9982,
        10166,
        9988,
        10003,
        10364,
        10641,
        10426,
        10830,
        11220,
        10231,
        10526,
        10848,
        11397,
        11217,
        10904,
        11223,
        11674,
        12208,
        12564,
        13354,
        13390,
        13448,
        14171,
        14846,
        15481,
        14854,
        17230,
        17272,
        17857,
        18947,
        19667,
        21002,
        22100,
        23681,
        23870,
        24244,
        25616,
        27061,
        28776,
        30927,
        33274,
        33057,
        33851,
        36103,
        39317,
        43621,
        44355,
        48429,
        49289,
        50799,
        54389,
        57977,
        62992,
        66200,
        70543,
        72124,
        74430,
        78510,
        82555,
        87872,
        93067,
        99240,
        101084,
        103681,
        107253,
        112186,
        117633,
        122948,
        127854,
        130532,
        132806,
        136253,
        139269,
        146554,
        154141,
        161167,
        161464,
        164271,
        169822,
        176504,
        183305,
        190194,
        196688,
        194366,
        193872,
        196689,
        200972,
        202461,
        211647,
        217660,
        215062,
        215723,
        223374,
        227985,
        234405,
        249727,
        259865,
        258454,
        258380,
        264475,
        278386,
        288403,
        308396,
        332577,
        337450,
        344604,
        356528,
        384188,
        407706,
        438540,
        482577,
        497857,
        511625,
        528806,
        554424,
        575369,
        593706,
        630983,
        669894,
        684606,
        695750,
        728500,
        756650,
        811974,
        837605,
        839129,
        836727,
        833155,
        826658,
        808179,
        791308,
        793555,
        772840,
        752276,
        734207,
        707870,
        695896,
        667418,
        628661,
        596927,
        565638,
        524626,
        478773,
        407141,
        380888,
        369338,
        350464,
        322561,
        302811,
        290578,
        269171,
        253928,
        253627,
        243198,
        229431,
        223280,
        214711,
        210512,
        211253,
        217898,
        212181,
        207723,
        211526,
        216964,
        219867,
        227169,
        235422,
        233277,
        231148,
        233244,
        243016,
        249904,
        242463,
        247541,
        237450,
        234447,
        246371,
        252490,
        257559,
        249989,
        258835,
        251787,
        254484,
        261171,
        266829,
        268939,
        275863,
        285386,
        277679,
        282205,
        284287,
        285193,
        302252,
        314061,
        326901,
        319713,
        317040,
        323475,
        326760,
        345365,
        356096,
        369623,
        354374,
        351283,
        356175,
        359888,
        362404,
        370331,
        385020,
        370563,
        370079,
        374897,
        376940,
        377628,
        381329,
        396324,
        378220,
        370793,
        373175,
        371083,
        368059,
        370682,
        382230,
        365611,
        358736,
        359907,
        359572,
        359313,
        362256,
        378866,
        361326,
        355416,
        362990,
        366053,
        370389,
        374407,
        389080,
        375785,
        373411,
        378784,
        383080,
        385216,
        391560,
        407935,
        394550,
        392940,
        400262,
        403125,
        405807,
        413045,
        433406,
        420650,
        419201,
        431493,
        440682,
        451074,
        465519,
        489805,
        503614,
        501825,
        509436,
        530207,
        547628,
        567141,
        596878,
        601853,
        599995,
        613828,
        637693,
        653570,
        647420,
        698519,
        681965,
        700075,
        696121,
        678102,
        683637,
        695558,
        693046,
        677707,
        683078,
        681074,
        677485,
        671689,
        670606,
        671301,
        651568,
        652123,
        637397,
        625798,
        637576,
        611028,
        622411,
        573756,
        573317,
        588023,
        574223,
        556406,
        547749,
        549017,
        515841,
        497885,
        486582,
        472758,
        457052,
        448982,
        450701,
        413908,
        399666,
        387636,
        375550,
        351519,
        343279,
        350782,
        325475,
        312028,
        304996,
        297470,
        286210,
        289560,
        295748,
        275108,
        253740,
        246061,
        242040,
        233893,
        236040,
        241514,
        224068,
        213743,
        204254,
        205874,
        198070,
        191662,
        193619,
        178594,
        166726,
        156875,
        148545,
        139805,
        134139,
        134786,
        132241,
        122813,
        111968,
        102726,
        95244,
        89361,
        89973,
        80586,
        74036,
        68213,
        63256,
        58599,
        55674,
        56077,
        50539,
        47319,
        45309,
        44098,
        42814,
        42972,
        44840,
        42374,
        42903,
        44946,
        46210,
        48173,
        51053,
        54837,
        53869,
        58473,
        65389,
        67010,
        75784,
        85161,
        97195,
        99017,
        114728,
        133607,
        152851,
        169391,
        192467,
        206931,
        208609,
        231413,
        254889,
        275507,
        295083,
        314311,
        332674,
        330870,
        348812,
        373977,
        389525,
        403639,
        416887,
        433185,
        420734,
        427764,
        436783,
        445745,
        446569,
        455772,
        471639,
        450534,
        450688,
        453890,
        453526,
        452409,
        455439,
        467587,
        445849,
        441857,
        439154,
        433743,
        426451,
        423601,
        431469,
        406629,
        395096,
        384213,
        373701,
        362750,
        354924,
        359539,
        334414,
        320493,
        309400,
        298047,
        285400,
        277735,
        280265,
        257417,
        243977,
        233506,
        223399,
        213181,
        206989,
        209022,
        191072,
        181217,
        172195,
        164954,
        157510,
        153132,
        154811,
        141508,
        135477,
        130361,
        125585,
        120929,
        118214,
        119917,
        110687,
        107076,
        103935,
        101177,
        98229,
        97172,
        99655,
        92902,
        91950,
        91344,
        90637,
    ]
    germany_infected = [
        13,
        13,
        9,
        7,
        7,
        3,
        2,
        2,
        2,
        2,
        3,
        11,
        32,
        58,
        63,
        114,
        149,
        187,
        246,
        528,
        652,
        782,
        1022,
        1204,
        1545,
        1938,
        2714,
        3621,
        4544,
        5754,
        7188,
        9274,
        12194,
        15161,
        19600,
        22071,
        24513,
        28480,
        29542,
        33570,
        37998,
        43862,
        48781,
        52683,
        52740,
        54933,
        58349,
        61245,
        65306,
        68244,
        69834,
        72859,
        69559,
        64639,
        63212,
        65513,
        65171,
        64521,
        62566,
        60502,
        58335,
        56631,
        53915,
        53769,
        53083,
        50685,
        48148,
        45913,
        44233,
        42417,
        40814,
        39772,
        38109,
        36174,
        34647,
        32860,
        30415,
        29129,
        28172,
        26433,
        24888,
        23164,
        22111,
        21351,
        20448,
        19883,
        19271,
        18206,
        17510,
        16720,
        15971,
        15590,
        15175,
        14539,
        13907,
        13334,
        12685,
        12334,
        11693,
        11630,
        11134,
        10763,
        10535,
        10655,
        10298,
        9767,
        9662,
        9220,
        8990,
        8399,
        8360,
        8124,
        8000,
        7966,
        7795,
        7458,
        6939,
        6717,
        6761,
        6629,
        6574,
        6532,
        6345,
        6950,
        7053,
        7273,
        7528,
        7686,
        7823,
        8065,
        7924,
        7946,
        8246,
        8136,
        8108,
        8224,
        7653,
        7436,
        7326,
        6900,
        7010,
        6745,
        6738,
        6525,
        6923,
        6446,
        6431,
        6151,
        6189,
        6170,
        6095,
        6077,
        6252,
        5958,
        5883,
        5855,
        6487,
        6583,
        6661,
        6928,
        6332,
        6503,
        6911,
        6747,
        6717,
        7572,
        8405,
        9114,
        8224,
        8609,
        8361,
        9121,
        10132,
        9731,
        10834,
        10209,
        10594,
        11308,
        11335,
        11647,
        12161,
        11908,
        12261,
        12780,
        14463,
        15873,
        16459,
        17133,
        17866,
        15549,
        15530,
        17154,
        18600,
        15388,
        15684,
        15666,
        15951,
        16062,
        15494,
        14793,
        14871,
        14788,
        16253,
        15420,
        16088,
        17193,
        14920,
        15012,
        15361,
        16975,
        15792,
        16143,
        16208,
        16272,
        18289,
        18258,
        18753,
        19315,
        19980,
        19743,
        19758,
        20170,
        22299,
        24649,
        25966,
        25977,
        26646,
        26683,
        27313,
        28017,
        29240,
        29604,
        30042,
        31314,
        31857,
        33734,
        36320,
        38964,
        40235,
        41862,
        44446,
        46812,
        50044,
        54379,
        59329,
        61853,
        65188,
        69005,
        72616,
        79229,
        87703,
        96933,
        103570,
        110372,
        118549,
        126429,
        136335,
        148691,
        161470,
        175480,
        177797,
        186625,
        194720,
        204902,
        216396,
        227367,
        235015,
        241571,
        247783,
        252547,
        259260,
        269527,
        278583,
        283144,
        287816,
        289398,
        290346,
        294502,
        301681,
        308442,
        310890,
        313926,
        313222,
        311297,
        311804,
        313393,
        315004,
        314543,
        316723,
        313849,
        309448,
        308364,
        310424,
        312932,
        315359,
        319224,
        317104,
        316199,
        319382,
        327030,
        334785,
        340849,
        348104,
        349979,
        351448,
        357965,
        365589,
        374401,
        381958,
        390529,
        391327,
        391574,
        398582,
        400144,
        396567,
        389209,
        388016,
        385562,
        381431,
        383871,
        383005,
        377086,
        370460,
        367033,
        358909,
        352567,
        351995,
        355932,
        360314,
        362054,
        362475,
        353332,
        344238,
        339893,
        338158,
        335799,
        333484,
        330776,
        319338,
        306017,
        298020,
        295362,
        292856,
        289566,
        287307,
        277260,
        265484,
        258123,
        254473,
        251337,
        247716,
        246002,
        238111,
        226369,
        218919,
        213390,
        206288,
        201878,
        199959,
        192204,
        181252,
        173987,
        168869,
        164151,
        159254,
        156803,
        151883,
        144590,
        140468,
        138518,
        136595,
        134727,
        135294,
        132480,
        128099,
        128476,
        130064,
        130031,
        130133,
        131020,
        128889,
        125878,
        126526,
        128639,
        129229,
        129990,
        131561,
        129907,
        127764,
        130614,
        135684,
        139639,
        143076,
        146609,
        145902,
        145775,
        152222,
        169880,
        168141,
        173783,
        179258,
        179171,
        179423,
        188386,
        199538,
        209822,
        218513,
        225000,
        224724,
        224499,
        231217,
        241613,
        246477,
        247837,
        248393,
        243625,
        237392,
        234615,
        241672,
        247536,
        252159,
        259091,
        257767,
        251940,
        265579,
        278484,
        284271,
        291988,
        294946,
        296179,
        295829,
        302894,
        310764,
        314530,
        321394,
        321765,
        315026,
        312380,
        314155,
        320547,
        313187,
        316267,
        316868,
        310852,
        302091,
        299848,
        294043,
        289847,
        286588,
        282659,
        273963,
        261347,
        251527,
        248755,
        240096,
        233388,
        229493,
        220406,
        208131,
        199654,
        192546,
        183828,
        176939,
        168808,
        162831,
        150387,
        139441,
        134558,
        129252,
        123541,
        119741,
        113747,
        104968,
        98569,
        93160,
        87627,
        82884,
        80603,
        77351,
        72286,
        67758,
        64209,
        61088,
        57914,
        55699,
        52749,
        48340,
        44758,
        41494,
        39246,
        37079,
        35824,
        34067,
        31780,
        29430,
        26809,
        25236,
        24108,
        23434,
        22351,
        21057,
        20239,
        19385,
        18913,
        18332,
        18142,
        17623,
        17100,
        17051,
        16925,
        16970,
        17054,
        17221,
        17068,
        17263,
        17960,
        18779,
        19519,
        20215,
        20698,
        20869,
        21565,
        22130,
        23491,
        24621,
        25389,
        25740,
        25944,
        26379,
        28010,
        28951,
        30519,
        30925,
        31279,
        30970,
        31493,
        33516,
        35661,
        37530,
        38888,
        39830,
        40032,
        41095,
        44353,
        47943,
        51353,
        54032,
        56065,
        57097,
        59760,
        63033,
        68337,
        77261,
        80978,
        84621,
        86657,
        90310,
        98467,
        106623,
        113135,
        117413,
        120543,
        122120,
        124363,
        132351,
        140164,
        145862,
        148720,
        152220,
        151846,
        152775,
        158399,
        164075,
        168087,
        168803,
        170156,
        169996,
        168772,
        170641,
        172359,
        171480,
        169548,
        169000,
        165356,
        161048,
        159962,
        158861,
        157091,
        154699,
        153953,
        150011,
        145836,
        146475,
        147031,
        146316,
        142415,
        141573,
        137915,
        135306,
        134930,
        147247,
        147211,
        146315,
        146649,
        143267,
        140695,
        142822,
        145809,
    ]
    uk_infected = [
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        5,
        5,
        5,
        5,
        8,
        12,
        15,
        27,
        30,
        41,
        74,
        89,
        132,
        172,
        234,
        272,
        324,
        394,
        509,
        697,
        988,
        1200,
        1283,
        1623,
        2202,
        2742,
        3356,
        4210,
        4788,
        5554,
        6701,
        7836,
        9534,
        11859,
        13867,
        15852,
        17845,
        20182,
        23428,
        26623,
        29911,
        32527,
        37218,
        40087,
        42333,
        46261,
        49078,
        55560,
        59450,
        63545,
        66724,
        70427,
        73746,
        76888,
        81013,
        84897,
        89690,
        93359,
        96077,
        99277,
        102700,
        106564,
        110173,
        113803,
        117377,
        120086,
        123002,
        127797,
        132693,
        136449,
        140080,
        143404,
        146710,
        151622,
        156204,
        159816,
        163049,
        166354,
        169665,
        172177,
        174665,
        177403,
        180293,
        182984,
        186078,
        188378,
        190085,
        191144,
        193224,
        192292,
        191414,
        190404,
        186712,
        185233,
        184009,
        181007,
        179188,
        173217,
        170442,
        167216,
        164733,
        161980,
        159746,
        157177,
        153687,
        150246,
        146089,
        143389,
        140713,
        138052,
        135443,
        132056,
        129170,
        126357,
        123599,
        121005,
        118763,
        114781,
        110529,
        107394,
        104388,
        101615,
        98363,
        93985,
        90208,
        87014,
        84420,
        81642,
        78749,
        76505,
        74263,
        71695,
        69193,
        66669,
        63842,
        61945,
        60448,
        59793,
        57977,
        55881,
        53916,
        52292,
        51185,
        49970,
        48851,
        47902,
        46898,
        46199,
        45045,
        44063,
        43159,
        42322,
        41461,
        40765,
        40126,
        39625,
        39103,
        38306,
        38296,
        38041,
        37436,
        36899,
        36461,
        36243,
        36006,
        35952,
        35729,
        35290,
        35181,
        34884,
        35158,
        35336,
        35871,
        36292,
        36145,
        36387,
        36288,
        36629,
        36821,
        37278,
        37731,
        38495,
        39034,
        39368,
        40182,
        40647,
        41535,
        42173,
        42766,
        43657,
        44411,
        45187,
        46297,
        47497,
        48799,
        49920,
        52076,
        54294,
        56139,
        58351,
        60706,
        63469,
        66184,
        68753,
        70618,
        73008,
        76405,
        79025,
        82476,
        86084,
        89453,
        93127,
        97232,
        103411,
        109831,
        116360,
        123027,
        129960,
        133538,
        141749,
        150696,
        160658,
        171204,
        176775,
        183637,
        195149,
        208878,
        221876,
        238508,
        251092,
        265137,
        276654,
        289400,
        305610,
        324004,
        341777,
        355774,
        370516,
        386318,
        403301,
        422949,
        448158,
        467727,
        486311,
        507224,
        525055,
        542841,
        562377,
        584329,
        604474,
        625708,
        643796,
        663419,
        678925,
        695931,
        717594,
        737447,
        757067,
        777367,
        793369,
        810386,
        825861,
        843324,
        869445,
        889138,
        908179,
        925509,
        939046,
        953890,
        964173,
        976687,
        985614,
        993622,
        1005097,
        1012486,
        1010691,
        1014020,
        1016084,
        1014383,
        1016228,
        1013387,
        1012764,
        1011781,
        1010559,
        1005857,
        1003205,
        1003031,
        1004318,
        1002010,
        995050,
        990318,
        984554,
        984999,
        985904,
        981485,
        982070,
        979278,
        981484,
        992191,
        997527,
        1000052,
        1014075,
        1024050,
        1041273,
        1060191,
        1073967,
        1082265,
        1095629,
        1101134,
        1121607,
        1152936,
        1181977,
        1214258,
        1233757,
        1263925,
        1291896,
        1325169,
        1363874,
        1405353,
        1437596,
        1504043,
        1520552,
        1574705,
        1620155,
        1664254,
        1710018,
        1757249,
        1714809,
        1739234,
        1761371,
        1786312,
        1805992,
        1829898,
        1851154,
        1875930,
        1892311,
        1906492,
        1911184,
        1915292,
        1927318,
        1954691,
        1982341,
        2004267,
        2024634,
        1948678,
        1944223,
        1944771,
        1940763,
        1924604,
        1914694,
        1904311,
        1882737,
        1861205,
        1837941,
        1813257,
        1789802,
        1771084,
        1746546,
        1726129,
        1695093,
        1655491,
        1618721,
        1575487,
        1533344,
        1486378,
        1442370,
        1391879,
        1341844,
        1290443,
        1246864,
        1186925,
        1133676,
        1084433,
        1044612,
        1050649,
        966099,
        923778,
        874631,
        879705,
        805268,
        773470,
        746623,
        714876,
        683723,
        652592,
        620634,
        595749,
        578793,
        564738,
        545937,
        522168,
        498796,
        480966,
        456794,
        386118,
        331129,
        282348,
        280122,
        275005,
        268200,
        260085,
        252025,
        244007,
        238021,
        231554,
        224307,
        218101,
        210908,
        203292,
        197508,
        193066,
        190131,
        187233,
        182559,
        179686,
        175560,
        172074,
        168640,
        166146,
        163585,
        159663,
        156625,
        152420,
        148232,
        145425,
        143471,
        140459,
        136911,
        133242,
        129637,
        127008,
        123873,
        120947,
        117538,
        114141,
        110223,
        105983,
        102520,
        100404,
        99010,
        96410,
        94152,
        92451,
        90448,
        89327,
        88615,
        88546,
        87854,
        87439,
        86698,
        86129,
        85724,
        86044,
        87166,
        86345,
        86142,
        86125,
        86006,
        86492,
        87876,
        90216,
        90716,
        91488,
        92529,
        93027,
        94738,
        97994,
        102553,
        106284,
        109105,
        112780,
        116519,
        121816,
        127563,
        134198,
        140489,
        146263,
        152091,
        157347,
        164065,
        173308,
        182258,
        190483,
        197607,
        206220,
        215327,
        229304,
        244031,
        257131,
        273165,
        285321,
        305108,
        322247,
        345053,
        369901,
        394408,
        416485,
        437843,
        461487,
        486151,
        513958,
        542598,
        574571,
        603091,
        631203,
        660797,
        691519,
        726912,
        769045,
        815008,
        863516,
        905109,
        936934,
        975327,
        1010699,
        1042252,
        1070624,
        1094186,
        1115450,
        1131093,
        1143173,
        1160079,
        1182373,
        1202404,
        1217642,
        1230209,
        1235854,
        1240550,
        1254479,
        1266353,
        1283272,
        1289127,
        1296305,
        1295532,
        1291030,
        1293576,
        1301806,
        1310331,
        1312773,
        1311072,
        1307298,
        1301671,
        1300127,
        1304496,
        1310277,
        1308195,
        1304210,
        1294224,
        1276692,
        1260942,
        1244738,
        1234922,
        1227683,
        1214856,
        1197612,
        1190310,
        1189580,
        1195957,
        1208604,
        1221100,
        1234562,
        1248027,
        1253975,
        1262970,
        1274425,
        1287289,
        1294712,
        1302124,
        1303449,
        1299538,
        1298088,
        1296145,
        1301059,
        1305875,
        1311849,
        1317933,
        1316000,
        1317267,
        1323995,
        1333849,
        1336250,
        1341459,
        1345036,
        1342026,
        1340835,
        1344625,
        1347377,
        1345244,
        1344308,
        1342891,
        1337558,
        1338541,
        1345985,
        1353767,
        1361626,
        1363335,
        1367403,
        1367482,
        1368664,
        1376530,
    ]

    df_infected_true = pd.DataFrame(
        np.array([belgium_infected, france_infected,
                 germany_infected, uk_infected]).T,
        columns=["belgium", "france", "germany", "uk"],
    )
    df_inf_true = df_infected_true.loc[
        list(range(321, df_infected_true.shape[0]))
    ].reset_index(drop=True)

    return df_inf_true


# --------------------------------------------------------------------------------------------
def create_R_rule_points(areas, number_xx_R, spline_xx_R):
    measures_sd = pd.DataFrame(np.nan, index=range(604 - 332), columns=areas)
    measures_masks = pd.DataFrame(
        np.nan, index=range(604 - 332), columns=areas)
    ids = {"countryA": 76, "countryB": 80, "countryC": 81, "countryD": 95}

    for i in ids.keys():
        url_sd = f"https://covid19.healthdata.org/api/data/hospitalization?location={ids[i]}&measure=22&scenario%5B%5D=6&scenario%5B%5D=1&scenario%5B%5D=2&scenario%5B%5D=3"
        url_mask = f"https://covid19.healthdata.org/api/data/hospitalization?location={ids[i]}&measure=25&scenario%5B%5D=6&scenario%5B%5D=1&scenario%5B%5D=3"
        data_sd = requests.get(url_sd).json()
        data_mask = requests.get(url_mask).json()
        values_sd = data_sd["values"][332:604]
        values_mask = data_mask["values"][332:604]
        val_sd = [values_sd[x][3] / 100 + 1 for x in range(len(values_sd))]
        val_masks = [values_mask[x][3] / 100 for x in range(len(values_mask))]

        measures_sd[i] = val_sd
        measures_masks[i] = val_masks

    mask_reduction = 0.33

    measures = (1 - mask_reduction * measures_masks) * measures_sd
    measures["Time"] = np.linspace(0, measures.shape[0] - 1, measures.shape[0])

    par_R = {}
    for i in range(len(areas)):
        for j in range(number_xx_R):
            key = f"yR0_{areas[i]}_{j}"
            time_key = f"xx_R0_{j}"
            time_value = spline_xx_R[time_key]
            closest_val = measures.loc[
                np.argmin(np.abs(measures["Time"] - time_value)), areas[i]
            ]
            if i == 1 and time_value / 7 >= 20:
                closest_val = closest_val * 0.93

            if i == 0 and time_value / 7 >= 20:
                closest_val = closest_val * 1.08

            if i == 2 and time_value / 7 >= 20:
                closest_val = closest_val * 0.73

            if i == 3 and time_value / 7 >= 20:
                closest_val = closest_val * 1.1

            if i in [0, 1] and time_value / 7 < 20:
                closest_val = closest_val * 1.16

            if i == 2 and time_value / 7 < 20:
                closest_val = closest_val * 1.05

            if i == 3 and time_value / 7 >= 10 and time_value / 7 < 10:
                closest_val = closest_val * 0.95

            if i == 3 and time_value / 7 < 10:
                closest_val = closest_val * 1
            par_R[key] = closest_val

    return par_R


def get_spline_xx(leave_out, length, end_data, number_yy, key_str):
    spline_xx = {}
    for i in range(number_yy):
        key = f"{key_str}{i}"
        if i <= number_yy - leave_out:
            spline_xx[key] = end_data / ((number_yy - leave_out) - 1) * (i)
        else:
            spline_xx[key] = ((length - end_data) / leave_out) * (
                leave_out + 1 - number_yy + i
            ) + end_data
    return spline_xx


# ------------------------------------------------------------------------------
def run_pareto_front(
    par_optimize,
    areas,
    model,
    number_yy,
    interventions,
    total_grid,
    length,
    parameters,
    solver,
    relative_population,
    rim,
    pop_size=100,
    number_generations=500,
    seed=1234,
    crossover_eta=15,
    crossover_prob=0.9,
    mutation_eta=20,
    lb=10 ** -15,
):

    set_parameter_by_name(model, parameters, fixed=True)

    transformed_frac_pareto = get_start_splines_pop_based(
        model, areas, number_yy)
    number_par = len(par_optimize)
    number_areas = len(areas)

    def get_optim_function_pareto(theta):

        dict_theta = dict(zip(par_optimize, np.log(theta / (1 - theta))))
        set_parameter_by_name(model, dict_theta, fixed=False)
        deaths = run_interventions_model(
            model,
            solver,
            interventions,
            parameters,
            total_grid,
            length,
            areas,
            number_yy,
            theta,
        )
        d = deaths["deaths"]
        return d[1: len(d)]

    ub = 1 - lb
    # pareto_values = np.array(run_model(np.repeat(0.5, 22)))
    pareto_values = []
    help_dict = dict(zip(par_optimize, np.repeat(0, len(par_optimize))))
    for i in range(len(areas)):
        area = areas[i]

        par_area = [x for x in par_optimize if area in x]
        par_rest = [x for x in par_optimize if not (area in x)]

        relative = {}
        for k in areas:
            if area != k:
                relative[k] = relative_population[k]

        for j in par_area:
            help_dict[j] = lb
        for j in par_rest:
            for l in relative.keys():
                if l in j:
                    help_dict[j] = relative[l] / \
                        np.sum(list(relative.values()))

        theta_par = np.array(list(help_dict.values()))
        pareto_values.append(get_optim_function_pareto(theta_par)[i])

    # problem with linear constraints
    class multi_objective_problem(ElementwiseProblem):
        def __init__(self):
            super().__init__(
                n_var=number_par,
                n_obj=number_areas,
                n_constr=number_areas,
                xl=np.repeat(lb, number_par),
                xu=np.repeat(ub, number_par),
            )

        def _evaluate(self, x, out, *args, **kwargs):
            out["F"] = get_optim_function_pareto(x)

            dict_para = dict(zip(par_optimize, x))
            vacs = ["vac1", "vac2"]
            names = []
            for i in vacs:
                for j in range(number_yy):
                    names.append(f"{i}_{j}")

            sums = []
            for i in names:
                s = 0
                for j in dict_para.keys():
                    if (i in j) and (
                        float(i.rsplit("_", 1)[-1]
                              ) == float(j.rsplit("_", 1)[-1])
                    ):
                        s += dict_para[j]
                sums.append(s)
            # out["G"] = np.concatenate([out["F"] - np.array(rim["deaths"][1:len(rim["deaths"])]), np.array(sums) - 1])
            out["G"] = np.array(sums) - 1

    problem = multi_objective_problem()
    ref_dirs = get_reference_directions("energy", 4, 100, seed=1234)

    # , sampling=np.array(di["x"].loc[23:24,:]))
    algorithm = CTAEA(ref_dirs=ref_dirs)
    termination = get_termination("n_gen", number_generations)

    res = pymoo.optimize.minimize(
        problem, algorithm, termination, seed=seed, save_history=True, verbose=True
    )
    f = pd.DataFrame(res.F)
    x = pd.DataFrame(res.X)
    p = np.array(rim["deaths"][1:5] * x.shape[0]).reshape(x.shape[0], 4)

    dict_out = {"f": f, "x": x, "p": p}

    return dict_out


def create_constraints(
    vaccines,
    par_optimize,
    areas,
):
    number_par = len(par_optimize)
    index_one = []
    index_two = []
    count = 0
    for j in vaccines[1: (len(vaccines))]:
        for i in range(
            int((len(par_optimize) / (len(vaccines) - 1) / (len(areas) - 1)))
        ):
            name = f"{j}_par{i}"
            list_name = [name]
            index_one += list_name
            count += 1

    A_constraints = pd.DataFrame(0, index=index_one, columns=par_optimize)

    for i in A_constraints.index:
        for j in A_constraints.columns:
            vac_par = i.split("_")
            last_char = j.rsplit("_", 1)[-1]
            if (vac_par[0] in j) and (vac_par[1] == f"par{last_char}"):
                A_constraints.loc[i, j] = 1

    A_bounds = np.diag(np.repeat(1, number_par))

    A = A_constraints.append(pd.DataFrame(A_bounds, columns=par_optimize)).reset_index(
        drop=True
    )

    return {"bounds": A_bounds, "lin_constraints": A_constraints, "both": A}
