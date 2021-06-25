import numpy as np
import pandas as pd
import pickle

from estimagic import minimize
from functools import partial

from models.vaccination.optimization_functions_estimagic import get_sum_of_states
from models.vaccination.parameter import general_set_up
from models.vaccination.parameter import fixed_parameter
from models.vaccination.parameter import parameters_vaccine_two
from models.vaccination.parameter import start_parameter_two
from models.vaccination.create_model_vaccination import model_vaccination_create_sbml
from models.vaccination.create_model_vaccination import b_distance_function

from visualization.model_results import plot_states
from visualization.model_results import plot_observables
from visualization.model_results import get_substates
from visualization.model_results import get_observables_by_name

from functions.run_sbml import run_model
from functions.run_sbml import create_observables_vaccination_rates
from functions.run_sbml import get_model_and_solver_from_sbml

from functions.vaccine_proportions import create_splines
from functions.vaccine_proportions import create_spline_parameters
from functions.vaccine_proportions import create_spline_transformation_rules
from functions.vaccine_proportions import create_one_minus_rules_splines

import matplotlib.pyplot as plt

first_time = False
# --------------------------Create Model--------------------------------------
model_name = "vaccination"
path_sbml = "stored_models/vaccination/" + model_name
vaccination_states = ["vac0", "vac1", "vac2"]
non_vaccination_state = general_set_up()["non_vaccination_state"]
virus_states = general_set_up()["virus_states"]
areas = general_set_up()["areas"]
distances = general_set_up()["distances"]
species_comp = general_set_up()["species_compartments"]

length_decision_period = 140
number_decision_periods = 1
first_vaccination_period = 0
number_yy = 11

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
    areas=areas, vaccination_states_removed=vaccination_states_removed, leave_out="last"
)

proportion_rules = {**spline_transformation_rule, **minus_one_rules}

splines = create_splines(
    areas_without=areas_without_last,
    vaccination_states_removed=vaccination_states_removed,
    stop=length_decision_period,
    length=number_yy,
)

model_vaccination_create_sbml(
    path=path_sbml,
    areas=areas,
    parameter_rules=proportion_rules,
    distances=distances,
    additional_parameters=spline_parameters,
    splines=splines,
)

observables = {"observable_time": {"name": "t", "formula": "t"}}

for index in ["nu", "proportion", "spline"]:
    new = create_observables_vaccination_rates(
        vaccination_states_removed=vaccination_states_removed,
        areas=areas,
        name_parameter=index,
    )
    observables = {**observables, **new}

model_directory = "stored_models/" + model_name + "/vaccination_dir"

model_and_solver = get_model_and_solver_from_sbml(
    path_sbml=path_sbml,
    model_name=model_name,
    model_directory=model_directory,
    only_import=not first_time,
)

created_model = {
    "model": model_and_solver["model"],
    "solver": model_and_solver["solver"],
    "observables": observables,
}


# -----------------------Run Model---------------------------------------------
length = length_decision_period * number_decision_periods + first_vaccination_period


model = created_model["model"]
# dir(model)
# par_values = model.getParameters()
# par_names = model.getParameterNames()
# model.getFixedParameterNames()
solver = created_model["solver"]
observables = created_model["observables"]

distance_parameters = {}
distances_df = pd.DataFrame(distances, columns=areas, index=areas)
b_distance_matrix = b_distance_function(distances_df)
for index_areas_row in areas:
    for index_areas_col in areas:
        name = f"distance_{index_areas_row}_{index_areas_col}"
        distance_parameters[name] = b_distance_matrix.loc[
            index_areas_row, index_areas_col
        ]

set_start_parameter = start_parameter_two()
set_fixed_parameter = {
    **fixed_parameter(),
    **parameters_vaccine_two(),
    **distance_parameters,
}
# ------------------------Run Current case-----------------------------------------
yy_names_str = []
for index_areas in areas_without_last:
    for index_vaccines in vaccination_states_removed:
        for index in range(number_yy):
            yy_names_str.append(f"yy_{index_areas}_{index_vaccines}_{index}")

yy_current_splines = np.repeat(0.0, len(yy_names_str))

yy_current = dict(zip(yy_names_str, yy_current_splines))

set_parameter_current = {**set_fixed_parameter, **yy_current}

model_results_current = run_model(
    model=model,
    solver=solver,
    periods=1,
    length_periods=length,
    set_start_parameter=set_start_parameter,
    set_parameter=set_parameter_current,
)

df_current = model_results_current
states_current = df_current["states"]
trajectory_current = {
    "states": states_current,
}
current_deceased_A = get_sum_of_states(
    model, trajectory_current, state_type=["dead", "countryA"], final_amount=True
)
current_deceased_B = get_sum_of_states(
    model, trajectory_current, state_type=["dead", "countryB"], final_amount=True
)
current_deceased = current_deceased_A + current_deceased_B
# -------------------------optimize--------------------------------------------
def criterion(
    theta,
    model,
    solver,
    set_fixed_parameter,
    set_start_parameter,
    yy_names,
    current_deceased_A,
    current_deceased_B,
):

    para = theta["value"]
    para_np = np.array(para)
    transformed_para = np.log(para_np / (1 - para_np))

    set_proportions = dict(zip(yy_names, transformed_para))
    set_fixed_parameter = {**fixed_parameter(), **parameters_vaccine_two()}
    set_parameter = {**set_fixed_parameter, **set_proportions}

    model_results = run_model(
        model=model,
        solver=solver,
        periods=1,
        length_periods=length,
        set_start_parameter=set_start_parameter,
        set_parameter=set_parameter,
    )

    trajectory_observables = model_results["observables"]
    trajectory_states = model_results["states"]
    trajectory_dict = {
        "states": trajectory_states,
        "observables": trajectory_observables,
    }

    deceased_A = get_sum_of_states(
        model, trajectory_dict, state_type=["dead", "countryA"], final_amount=True
    )

    deceased_B = get_sum_of_states(
        model, trajectory_dict, state_type=["dead", "countryB"], final_amount=True
    )

    if (deceased_A > current_deceased_A) or (deceased_B > current_deceased_B):
        value = 10 ** 15
    else:
        value = deceased_A + deceased_B

    out = {"value": value}

    return out


to_optimize = partial(
    criterion,
    model=model,
    solver=solver,
    set_fixed_parameter=set_fixed_parameter,
    set_start_parameter=start_parameter_two(),
    yy_names=yy_names_str,
    current_deceased_A=current_deceased_A,
    current_deceased_B=current_deceased_B,
)

n_multi = 100
number_vaccines = len(vaccination_states) - 1

number_parameters = len(yy_names_str)
lower_bound = 0.000001
upper_bound = 0.999999


cols = (
    ["start_" + yy_parameter for yy_parameter in yy_names_str]
    + ["theta_" + yy_parameter for yy_parameter in yy_names_str]
    + ["value"]
)

np.random.seed(12345)
res_optimal = {"solution_criterion": np.inf}
df_multistart = pd.DataFrame(data=None, index=range(n_multi), columns=cols)
for index_multi in range(n_multi):
    print(index_multi)
    start = np.random.uniform(lower_bound, upper_bound, number_parameters)
    start_params = pd.DataFrame(
        data=start,
        columns=["value"],
        index=yy_names_str,
    )
    start_params["lower_bound"] = np.repeat(lower_bound, number_parameters)
    start_params["upper_bound"] = np.repeat(upper_bound, number_parameters)

    res = minimize(
        criterion=to_optimize,
        params=start_params,
        algorithm="nag_pybobyqa",
    )
    if res["solution_criterion"] < res_optimal["solution_criterion"]:
        res_optimal = res
        print(f"new minimum at {index_multi} iteration.")

    print(res["success"])

    df_multistart.iloc[index_multi] = np.append(
        start, np.append(res["solution_x"], res["solution_criterion"])
    )

with open(
    "/home/manuel/Documents/VaccinationDistribution/code/objects/optimal_df.pkl", "wb"
) as output:
    optimal = res_optimal
    pickle.dump(optimal, output, pickle.HIGHEST_PROTOCOL)

with open(
    "/home/manuel/Documents/VaccinationDistribution/code/objects/optimal_df.pkl", "rb"
) as input:
    b = pickle.load(input)
# ------------------------Run Optimum-----------------------------------------
theta_optimal = np.array(res_optimal["solution_x"])
transformed_theta_optimal = np.log(theta_optimal / (1 - theta_optimal))
yy_optimal = dict(zip(yy_names_str, transformed_theta_optimal))

set_parameter = {**set_fixed_parameter, **yy_optimal}

model_results_optimal = run_model(
    model=model,
    solver=solver,
    periods=1,
    length_periods=length,
    set_start_parameter=set_start_parameter,
    set_parameter=set_parameter,
)

# -----------------------Plot Optimum--------------------------------------------
plot_path = "/home/manuel/Documents/VaccinationDistribution/paper/images/distinguished_"
trajectory_states_optimal = model_results_optimal["states"]
trajectory_observables_optimal = model_results_optimal["observables"]

A_peak_vacM = trajectory_states_optimal.iloc[
    trajectory_states_optimal["infectious_countryA_vac0_virM"].idxmax()
]["amici_t"]
A_peak_vacW = trajectory_states_optimal.iloc[
    trajectory_states_optimal["infectious_countryA_vac0_virW"].idxmax()
]["amici_t"]
B_peak_vacM = trajectory_states_optimal.iloc[
    trajectory_states_optimal["infectious_countryB_vac0_virM"].idxmax()
]["amici_t"]
B_peak_vacW = trajectory_states_optimal.iloc[
    trajectory_states_optimal["infectious_countryB_vac0_virW"].idxmax()
]["amici_t"]

observables_names_nu_optimal = get_observables_by_name(
    observables, substrings=["nu"], include_all=True
)

observables_names_proportion_optimal = get_observables_by_name(
    observables, substrings=["proportion", "countryA"], include_all=True
)

plot_observables(
    results=model_results_optimal,
    model=model,
    observable_ids=observables_names_nu_optimal,
    time_name="t",
)

plot_observables(
    results=model_results_optimal,
    model=model,
    observable_ids=["pline_countryA_vac1", "pline_countryA_vac2"],
    set_off_scientific_notation=True,
    decimal_floats=4,
    time_name="t",
)

# Proportion_observables
fig_obs, ax_obs = plot_observables(
    results=model_results_optimal,
    model=model,
    observable_ids=observables_names_proportion_optimal,
    set_off_scientific_notation=True,
    decimal_floats=1,
    time_name="t",
    title="Vaccines assigned to country A",
    ylabel="Fraction",
    xlabel="Days",
    colors=["C2", "C3"],
    custom_label=["Vaccine wild type", "Vaccine mutant"],
    legend_next_to_plot=False,
)
ax_obs.spines["right"].set_visible(False)
ax_obs.spines["top"].set_visible(False)
ax_obs.get_legend().remove()
# legend_location="upper left"
# ax_obs.axvline(x=A_peak_vacW, label='max_infectious_countryA_virW', color="green", linestyle='--')
# ax_obs.axvline(x=A_peak_vacM, label='max_infectious_countryA_virM', color="red", linestyle='--')
# ax_obs.axvline(x=B_peak_vacW, label='max_infectious_countryB_virW', color="black", linestyle='--')
# ax_obs.axvline(x=B_peak_vacM, label='max_infectious_countryB_virM', color="purple", linestyle='--')
# ax_obs.legend(bbox_to_anchor=(1, 0.8, 0.3, 0.2), loc=legend_location)
fig_obs.savefig(plot_path + "observables", bbox_inches="tight")

# ------- Infectious A
substates_optimal_A = get_substates(
    model=model, substrings=["vac0", "infectious", "countryA"], include_all=True
)

fig, ax = plot_states(
    results=model_results_optimal,
    model=model,
    state_ids=substates_optimal_A,
    time_name="t",
    title="Unvaccinated infectious individuals (A)",
    colors=["C2", "C3"],
    xlabel="Days",
    ylabel="Number of individuals",
    custom_label=["Virus wild type", "Virus mutant"],
    ylim=[0, 2.3 * 10 ** 6],
    legend_next_to_plot=False,
)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.get_legend().remove()
fig.savefig(plot_path + "infectious_A", bbox_inches="tight")
# -------

# ------- Infectious B
substates_optimal_B = get_substates(
    model=model, substrings=["vac0", "infectious", "countryB"], include_all=True
)

fig, ax = plot_states(
    results=model_results_optimal,
    model=model,
    state_ids=substates_optimal_B,
    time_name="t",
    title="Unvaccinated infectious individuals (B)",
    colors=["C2", "C3"],
    xlabel="Days",
    ylabel="Number of individuals",
    custom_label=["Virus wild type", "Virus mutant"],
    ylim=[0, 2.3 * 10 ** 6],
    legend_next_to_plot=False,
)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.get_legend().remove()
fig.savefig(plot_path + "infectious_B", bbox_inches="tight")
# --------
# ------- Infectious A
substates_optimal_A = get_substates(
    model=model, substrings=["vac0", "countryA"], include_all=True
)

fig, ax = plot_states(
    results=model_results_optimal,
    model=model,
    state_ids=substates_optimal_A,
    time_name="t",
    title="Unvaccinated infectious individuals (A)",
    colors=["C0", "C2", "C3", "C1", "C4", "C5", "C6"],
    xlabel="Days",
    ylabel="Number of individuals",
    custom_label=[
        "Susceptible",
        "Infectious (W)",
        "Infectious (M)",
        "Recovered (W)",
        "Recovered (M)",
        "Deceased (W)",
        "Deceased (M)",
    ],
    # ylim = [0, 2.3*10**6],
    legend_next_to_plot=True,
)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.get_legend().remove()
fig.savefig(plot_path + "unvaccinated_A", bbox_inches="tight")


# ------- Infectious B
substates_optimal_B = get_substates(
    model=model, substrings=["vac0", "countryB"], include_all=True
)

fig, ax = plot_states(
    results=model_results_optimal,
    model=model,
    state_ids=substates_optimal_B,
    time_name="t",
    title="Unvaccinated infectious individuals (B)",
    colors=["C0", "C2", "C3", "C1", "C4", "C5", "C6"],
    xlabel="Days",
    ylabel="Number of individuals",
    custom_label=[
        "Susceptible",
        "Infectious (W)",
        "Infectious (M)",
        "Recovered (W)",
        "Recovered (M)",
        "Deceased (W)",
        "Deceased (M)",
    ],
    # ylim = [0, 2.3*10**6],
    legend_next_to_plot=True,
)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.get_legend().remove()
fig.savefig(plot_path + "unvaccinated_B", bbox_inches="tight")
# --------

# -----------run seperated-----------------------------------------------------
yy_seperated_splines_wild = np.repeat(20.0, len(yy_names_str) / 2)
yy_seperated_splines_mutant = np.repeat(-20.0, len(yy_names_str) / 2)
yy_seperated_splines = np.concatenate(
    (yy_seperated_splines_wild, yy_seperated_splines_mutant)
)
yy_current = dict(zip(yy_names_str, yy_seperated_splines))
set_parameter_seperated = {**set_fixed_parameter, **yy_current}
model_results_seperated = run_model(
    model=model,
    solver=solver,
    periods=1,
    length_periods=length,
    set_start_parameter=set_start_parameter,
    set_parameter=set_parameter_seperated,
)
df_seperated = model_results_seperated
states_seperated = df_seperated["states"]
trajectory_seperated = {
    "states": states_seperated,
}
seperated_deceased_A = get_sum_of_states(
    model, trajectory_seperated, state_type=["dead", "countryA"], final_amount=True
)
seperated_deceased_B = get_sum_of_states(
    model, trajectory_seperated, state_type=["dead", "countryB"], final_amount=True
)
seperated_deceased = seperated_deceased_A + seperated_deceased_B


df_current = model_results_current
states_current = df_current["states"]
trajectory_current = {
    "states": states_current,
}
trajectory_optimal = {
    "states": trajectory_states_optimal,
}

current_deceased_A = get_sum_of_states(
    model, trajectory_current, state_type=["dead", "countryA"], final_amount=True
)
current_deceased_B = get_sum_of_states(
    model, trajectory_current, state_type=["dead", "countryB"], final_amount=True
)
current_deceased = current_deceased_A + current_deceased_B

dead_A_current = get_sum_of_states(
    model, trajectory_current, state_type=["dead", "countryA"], final_amount=True
)
dead_B_current = get_sum_of_states(
    model, trajectory_current, state_type=["dead", "countryB"], final_amount=True
)

dead_A_optimal = get_sum_of_states(
    model, trajectory_optimal, state_type=["dead", "countryA"], final_amount=True
)
dead_B_optimal = get_sum_of_states(
    model, trajectory_optimal, state_type=["dead", "countryB"], final_amount=True
)

optimal_deceased = get_sum_of_states(
    model, trajectory_optimal, state_type=["dead"], final_amount=True
)

print_compare_A = f"A: Current strategy: {dead_A_current}; Optimal strategy: {dead_A_optimal}; Seperated strategy: {seperated_deceased_A}"
print_compare_B = f"B: Current strategy: {dead_B_current}; Optimal strategy: {dead_B_optimal}; Seperated strategy: {seperated_deceased_B}"
print_compare_total = f"Total: Current strategy: {current_deceased}; Optimal strategy: {optimal_deceased}; Seperated strategy: {seperated_deceased}"
print(print_compare_total)
print(print_compare_A)
print(print_compare_B)

vac1_A_optimal = get_sum_of_states(
    model, trajectory_optimal, state_type=["vac1", "countryA"], final_amount=True
)
vac2_A_optimal = get_sum_of_states(
    model, trajectory_optimal, state_type=["vac2", "countryA"], final_amount=True
)
vac1_B_optimal = get_sum_of_states(
    model, trajectory_optimal, state_type=["vac1", "countryB"], final_amount=True
)
vac2_B_optimal = get_sum_of_states(
    model, trajectory_optimal, state_type=["vac2", "countryB"], final_amount=True
)
# data = [[415344, 206562, 208782],
# [406495, 204607, 201888],
# [404909, 218045, 186863]]
# X = np.arange(3)
# fig = plt.figure()
# ax = fig.add_axes([0,0,1,1])
# ax.bar(X + 0.00, data[0], color = 'C0', width = 0.15)
# ax.bar(X + 0.15, data[1], color = 'C1', width = 0.15)
# ax.bar(X + 0.30, data[2], color = 'C2', width = 0.15)

fig_obs, ax_obs = plot_observables(
    results=model_results_current,
    model=model,
    observable_ids=observables_names_proportion_optimal,
    set_off_scientific_notation=True,
    decimal_floats=4,
    time_name="t",
    title="Equiallocated fraction of vaccines assigned to country A",
    ylabel="$f_{j,A}(t)$",
    xlabel="t (in days)",
    colors=["C2", "C3"],
    custom_label=["Vaccine wild type", "Vaccine mutant"],
)
fig.savefig(plot_path + "current_allocation", bbox_inches="tight")

# plot infectious in A and B current vs. optimal
# ------- Infectious A
substates_optimal_A = get_substates(
    model=model, substrings=["vac0", "infectious", "countryA"], include_all=True
)

fig, ax = plot_states(
    results=model_results_current,
    model=model,
    state_ids=substates_optimal_A,
    time_name="t",
    title="Unvaccinated infectious individuals (A)",
    colors=["C2", "C3"],
    xlabel="Days",
    ylabel="Number of individuals",
    custom_label=["Virus wild type", "Virus mutant"],
    ylim=[0, 2.3 * 10 ** 6],
    legend_next_to_plot=False,
)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.get_legend().remove()
fig.savefig(plot_path + "infectious_A_current", bbox_inches="tight")
# -------
# ------- Infectious B
substates_optimal_B = get_substates(
    model=model, substrings=["vac0", "infectious", "countryB"], include_all=True
)

fig, ax = plot_states(
    results=model_results_current,
    model=model,
    state_ids=substates_optimal_B,
    time_name="t",
    title="Unvaccinated infectious individuals (B)",
    colors=["C2", "C3"],
    xlabel="Days",
    ylabel="Number of individuals",
    custom_label=["Virus wild type", "Virus mutant"],
    ylim=[0, 2.3 * 10 ** 6],
    legend_next_to_plot=False,
)
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.get_legend().remove()
fig.savefig(plot_path + "infectious_B_current", bbox_inches="tight")

# -----------------------------------------------------------------------------
# plot the top 4 results
minus_df = -df_multistart
minus_df["value"] = pd.to_numeric(minus_df["value"])
four = -minus_df.nlargest(20, "value")
four_thetas = four.loc[:, four.columns.str.startswith("theta")]

indices = four.index
c = 0
for index in range(4):
    ind = indices[index]
    theta_sim = np.array(four_thetas.iloc[index], dtype="float")
    transformed_theta_sim = np.log(theta_sim / (1 - theta_sim))
    yy_sim = dict(zip(yy_names_str, transformed_theta_sim))

    set_parameter_sim = {**set_fixed_parameter, **yy_sim}

    model_results_sim = run_model(
        model=model,
        solver=solver,
        periods=1,
        length_periods=length,
        set_start_parameter=set_start_parameter,
        set_parameter=set_parameter_sim,
    )

    fig_obs, ax_obs = plot_observables(
        results=model_results_sim,
        model=model,
        observable_ids=observables_names_proportion_optimal,
        set_off_scientific_notation=True,
        decimal_floats=1,
        time_name="t",
        title="Vaccines assigned to country A",
        ylabel="Fraction",
        xlabel="Days",
        colors=["C2", "C3"],
        custom_label=["Vaccine wild type", "Vaccine mutant"],
        legend_next_to_plot=False,
    )
    ax_obs.spines["right"].set_visible(False)
    ax_obs.spines["top"].set_visible(False)
    # ax_obs.get_legend().remove()
    c += 1
    save = f"sim_{c}"
    fig_obs.savefig(plot_path + save, bbox_inches="tight")
