import numpy as np

from models.vaccination.parameter import general_set_up
from models.vaccination.parameter import fixed_parameter
from models.vaccination.parameter import parameters_vaccine_one
from models.vaccination.parameter import start_parameter_one
from models.vaccination.create_model_vaccination import model_vaccination_create_sbml

from visualization.model_results import plot_states
from visualization.model_results import plot_observables
from visualization.model_results import get_substates
from visualization.model_results import get_observables_by_name

from functions.run_sbml import run_model
from functions.run_sbml import create_observables_vaccination_rates
from functions.run_sbml import get_model_and_solver_from_sbml

from functions.vaccine_proportions import create_spline_parameters
from functions.vaccine_proportions import create_spline_transformation_rules
from functions.vaccine_proportions import create_one_minus_rules_splines

from amici.sbml_utils import amici_time_symbol
from amici.splines import UniformGrid

# --------------------------Create Model--------------------------------------
model_name = "vaccination"
path_sbml = "stored_models/vaccination/" + model_name
vaccination_states = ["vac0", "vac1"]
non_vaccination_state = general_set_up()["non_vaccination_state"]
virus_states = general_set_up()["virus_states"]
areas = general_set_up()["areas"]
distances = general_set_up()["distances"]
species_comp = general_set_up()["species_compartments"]

length_decision_period = 40
number_decision_periods = 1
first_vaccination_period = 0

vaccination_states_removed = [
    x for x in vaccination_states if x != non_vaccination_state
]

areas_without_last = [
    x for x in areas if x != areas[-1]
]

spline_parameters = create_spline_parameters(areas=areas,
                                             vaccination_states_removed=vaccination_states_removed)
spline_transformation_rule = create_spline_transformation_rules(areas=areas,
                                             vaccination_states_removed=vaccination_states_removed)

minus_one_rules = create_one_minus_rules_splines(areas=areas, vaccination_states_removed=vaccination_states_removed,
                       leave_out = "last")

proportion_rules = {**spline_transformation_rule, **minus_one_rules}

#TODO Define better grid space xx, fix issue with absolut values
xx = np.array([0, 0.2, 0.6]) #UniformGrid(start=0, stop=40, length=s) 
yy = np.round(np.random.uniform(size=len(xx)))
from sympy import symbols
xx_names = [symbols("xx1"), symbols("xx2"), symbols("xx3"), symbols("xx4"), symbols("xx5")] #x is not a feasible name!
yy_names = [symbols("y1"), symbols("y2"), symbols("y3")]

#TODO extend for multiple vaccines and countries; first see how to optimize over yy
splines = { 'spline_countryA_vac1' : {"time_symbol" : "t",
                                          "xx" : xx,
                                          "yy" : yy,
                                          "xx_names" : xx_names, #not used currently
                                          "yy_names" : yy_names}}

model_vaccination_create_sbml(
    path=path_sbml,
    areas=areas,
    parameter_rules=proportion_rules,
    distances=distances,
    additional_parameters=spline_parameters,
    splines=splines
)

observables_nu = create_observables_vaccination_rates(
    vaccination_states_removed=vaccination_states_removed,
    areas=areas,
    name_parameter="nu",
)

observables_proportion = create_observables_vaccination_rates(
    vaccination_states_removed=vaccination_states_removed,
    areas=areas,
    name_parameter="proportion",
)

observables_splines = create_observables_vaccination_rates(
    vaccination_states_removed=vaccination_states_removed,
    areas=areas_without_last ,
    name_parameter="spline",
)
observables_time = {"observable_time": {"name": "t", "formula": "t"}}

observables = {**observables_nu, **observables_proportion, **observables_splines ,**observables_time}

model_directory = "stored_models/" + model_name + "/vaccination_dir"

model_and_solver = get_model_and_solver_from_sbml(
    path_sbml=path_sbml,
    model_name=model_name,
    model_directory=model_directory,
)

created_model = {
    "model": model_and_solver["model"],
    "solver": model_and_solver["solver"],
    "observables": observables,
}
# -----------------------Run Model---------------------------------------------
length = min(
    length_decision_period * number_decision_periods + first_vaccination_period, 370
)

model = created_model["model"]
# dir(model)
# par_values = model.getParameters()
# par_names = model.getParameterNames()
# model.getFixedParameterNames()
solver = created_model["solver"]
observables = created_model["observables"]

set_start_parameter = start_parameter_one()
set_fixed_parameter = {**fixed_parameter(), **parameters_vaccine_one()}

set_parameter = {**set_fixed_parameter}

model_results = run_model(
    model=model,
    solver=solver,
    periods=1,
    length_periods=length,
    set_start_parameter=set_start_parameter,
    set_parameter=set_parameter,
)

# -----------------------Plot Model--------------------------------------------
trajectory_states = model_results["states"]
trajectory_observables = model_results["observables"]

observables_names_nu = get_observables_by_name(
    observables, substrings=["nu"], include_all=True
)

observables_names_proportion = get_observables_by_name(
    observables, substrings=["proportion"], include_all=True
)

plot_observables(
    results=model_results, model=model, observable_ids=observables_names_nu, time_name="t"
)

plot_observables(
    results=model_results,
    model=model,
    observable_ids=["pline_countryA_vac1"],
    set_off_scientific_notation=True,
    decimal_floats=4,
    time_name="t"
)

plot_observables(
    results=model_results,
    model=model,
    observable_ids= observables_names_proportion,
    set_off_scientific_notation=True,
    decimal_floats=4,
    time_name="t"
)

substates = get_substates(
    model=model, substrings=["vac0", "countryA"], include_all=True
)

fig, ax = plot_states(results=model_results, model=model, state_ids=substates, time_name="t")

# -------------------------optimize--------------------------------------------
from estimagic import minimize


def criterion(theta, model, solver, set_fixed_parameter):

    para = theta["value"]
    set_proportions = {
        "proportion_par_countryA_vac1_0": para[0],
        "proportion_par_countryA_vac1_200": para[1],
    }

    set_parameter = {**set_fixed_parameter, **set_proportions}

    model_results = run_model(
        model=model,
        solver=solver,
        periods=1,
        length_periods=20,
        set_start_parameter=set_start_parameter,
        set_parameter=set_parameter,
        observables_names=observables.keys(),
    )

    trajectory_observables = model_results["observables"]
    trajectory_states = model_results["states"]
    trajectory_dict = {
        "states": trajectory_states,
        "observables": trajectory_observables,
    }

    out = {
        "value": get_sum_of_states(
            model, trajectory_dict, state_type=["dead"], final_amount=True
        )
    }

    return out


to_optimize = partial(
    criterion, model=model, solver=solver, set_fixed_parameter=set_fixed_parameter
)

n_multi = 20
number_vaccines = len(vaccination_states) - 1
periods = number_decision_periods
cols = (
    [f"start_{i}" for i in range(number_vaccines * periods)]
    + [f"theta_{i}" for i in range(number_vaccines * periods)]
    + ["value"]
)
df_multistart = pd.DataFrame(data=None, index=range(n_multi), columns=cols)
for index_multi in range(n_multi):
    print(index_multi)
    start = np.random.uniform(0, 1, number_vaccines * periods)
    start_params = pd.DataFrame(
        data=start,
        columns=["value"],
        index=[f"theta_{i}" for i in range(number_vaccines * periods)],
    )
    start_params["lower_bound"] = np.repeat(0, number_vaccines * periods)
    start_params["upper_bound"] = np.repeat(1, number_vaccines * periods)

    res = minimize(
        criterion=to_optimize,
        params=start_params,
        algorithm="nag_pybobyqa",
    )

    df_multistart.iloc[index_multi] = np.append(
        start, np.append(res["solution_x"], res["solution_criterion"])
    )

df_multistart[df_multistart.value == df_multistart.value.min()]
