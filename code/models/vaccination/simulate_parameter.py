import numpy as np
import pandas as pd
import pickle

import pypesto
import pypesto.optimize as optimize
import pypesto.visualize as visualize

from models.vaccination.optimization_functions_estimagic import get_sum_of_states
from models.vaccination.parameter import general_set_up
from models.vaccination.parameter import fixed_parameter
from models.vaccination.parameter import parameters_vaccine_two
from models.vaccination.parameter import start_parameter_two
from models.vaccination.create_model_vaccination import model_vaccination_create_sbml
from models.vaccination.create_model_vaccination import b_distance_function

from functions.run_sbml import run_model
from functions.run_sbml import create_observables_vaccination_rates
from functions.run_sbml import get_model_and_solver_from_sbml

from functions.vaccine_proportions import create_splines
from functions.vaccine_proportions import create_spline_parameters
from functions.vaccine_proportions import create_spline_transformation_rules
from functions.vaccine_proportions import create_one_minus_rules_splines
from functions.vaccine_proportions import create_parameters_piecewise_vaccine_supply
from functions.vaccine_proportions import create_vaccine_supply_rules
from functions.create_vaccine_doses_inflow import create_inflow_from_data
from visualization.model_results import get_substates


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
n_intervals = 6000

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
    areas=areas, vaccination_states_removed=vaccination_states_removed, leave_out="last"
)

proportion_rules = {**spline_transformation_rule, **minus_one_rules}

splines = create_splines(
    areas_without=areas_without_last,
    vaccination_states_removed=vaccination_states_removed,
    stop=length_decision_period,
    length=number_yy,
)


vaccine_supply_rules = create_vaccine_supply_rules(
    vaccination_states, non_vaccination_state, xx_list
)

model_vaccination_create_sbml(
    path=path_sbml,
    areas=areas,
    parameter_rules={**proportion_rules, **vaccine_supply_rules},
    distances=distances,
    additional_parameters={**spline_parameters, **parameter_vacc_supply},
    splines=splines,
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
for index in ["nu", "proportion", "spline"]:
    new = create_observables_vaccination_rates(
        vaccination_states_removed=vaccination_states_removed,
        areas=areas,
        name_parameter=index,
    )

    obs = {**obs, **new}

observables = {
    **obs,
    **observables_vaccinated_vacc,
}

model_directory = "stored_models/" + model_name + "/vaccination_dir"

model_and_solver = get_model_and_solver_from_sbml(
    path_sbml=path_sbml,
    model_name=model_name,
    model_directory=model_directory,
    only_import=not first_time,
    observables=observables,
)

created_model = {
    "model": model_and_solver["model"],
    "solver": model_and_solver["solver"],
    "observables": observables,
}


# -----------------------Run Model---------------------------------------------
length = length_decision_period * number_decision_periods + first_vaccination_period


model = created_model["model"]
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
vaccine_supply_parameter_strings = [
    x for x in model.getParameterNames() if "vaccine_supply" in x
]


vaccine_supply_parameter_values = create_inflow_from_data(
    number_decision_periods=number_yy - 1
)

# comment out if true vaccine inflow
# vaccine_supply_parameter_values = np.repeat(
#    float(300000), len(vaccine_supply_parameter_strings) )

set_vaccine_supply_parameter = dict(
    zip(vaccine_supply_parameter_strings, vaccine_supply_parameter_values)
)
set_fixed_parameter = {
    **fixed_parameter(),
    **parameters_vaccine_two(),
    **distance_parameters,
    **set_vaccine_supply_parameter,
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
    number_intervals=n_intervals,
)


df_current = model_results_current
states_current = df_current["states"]
trajectory_current = {
    "states": states_current,
}
var_interest = "dead"
current_deceased_A = get_sum_of_states(
    model, trajectory_current, state_type=[var_interest, "countryA"], final_amount=True
)
current_deceased_B = get_sum_of_states(
    model, trajectory_current, state_type=[var_interest, "countryB"], final_amount=True
)
total_deceased = current_deceased_A + current_deceased_B
# -------------------------find pareto optimal solution------------------------
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

    para = theta
    para_np = np.array(para)
    transformed_para = np.log(para_np / (1 - para_np))

    set_proportions = dict(zip(yy_names, transformed_para))
    set_fixed_parameter = {**set_fixed_parameter, **parameters_vaccine_two()}
    set_parameter = {**set_fixed_parameter, **set_proportions}

    model_results = run_model(
        model=model,
        solver=solver,
        periods=1,
        length_periods=length,
        set_start_parameter=set_start_parameter,
        set_parameter=set_parameter,
        number_intervals=n_intervals,
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

    return deceased_A, deceased_B


def criterion_sim(
    para1,
    para2,
    name1,
    name2,
    model,
    solver,
    set_fixed_parameter,
    set_start_parameter,
    yy_names_str,
):

    par = np.array([para1, para2])
    transformed_para = np.log(par / (1 - par))

    set_proportions = dict(zip(yy_names_str, np.repeat(float(0), len(yy_names_str))))
    set_proportions[name1] = transformed_para[0]
    set_proportions[name2] = transformed_para[1]

    set_fixed_parameter = {**set_fixed_parameter, **parameters_vaccine_two()}
    set_parameter = {**set_fixed_parameter, **set_proportions}

    model_results = run_model(
        model=model,
        solver=solver,
        periods=1,
        length_periods=length,
        set_start_parameter=set_start_parameter,
        set_parameter=set_parameter,
        number_intervals=n_intervals,
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

    return deceased_A, deceased_B


def unrestricted_sim(
    para1,
    para2,
    name1,
    name2,
):

    A, B = criterion_sim(
        para1=para1,
        para2=para2,
        name1=name1,
        name2=name2,
        model=model,
        solver=solver,
        set_fixed_parameter=set_fixed_parameter,
        set_start_parameter=start_parameter_two(),
        yy_names_str=yy_names_str,
    )

    total = A + B

    return total


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def plot_simulation_para(
    para_numb,
    length_interval=15,
    total_deceased=total_deceased,
    three_d=False,
    path="/home/manuel/Documents/VaccinationDistribution/paper/images/simulation_para",
):

    interval = np.linspace(10 ** (-7), 1 - 10 ** (-7), length_interval)
    save_sim = pd.DataFrame(
        columns=["par1", "par2", "value"], index=range(len(interval) ** 2)
    )
    index = 0
    for index1 in range(len(interval)):
        for index2 in range(len(interval)):
            para1 = interval[index1]
            para2 = interval[index2]
            deceased = unrestricted_sim(
                para1=para1,
                para2=para2,
                name1=f"yy_countryA_vac1_{para_numb}",
                name2=f"yy_countryA_vac2_{para_numb}",
            )

            vec_save = np.array([para1, para2, deceased])
            save_sim.iloc[index] = vec_save
            index += 1

    if three_d is False:
        fig, ax = plt.subplots()

        save_pivot = save_sim.pivot(index="par2", columns="par1", values="value")
        image = ax.pcolormesh(
            interval,
            interval,
            np.array(save_pivot, dtype="float") / 10 ** 6,
            cmap="viridis",
        )
        fig.colorbar(image, ax=ax)
        save_path = f"{path}_mesh_{para_numb}"
        ax.scatter([0.5], [0.5], c="firebrick", marker="x", linewidth=2)

        maximal_row = save_sim[save_sim["value"] == np.min(save_sim["value"])]
        max_par1 = float(maximal_row["par1"])
        max_par2 = float(maximal_row["par2"])

        adjust1 = 0.01
        if max_par1 > 0.5:
            adjust_1 = -0.01
        adjust2 = -0.01
        if max_par2 > 0.5:
            adjust2 = 0.01

        ax.arrow(
            0.5,
            0.5,
            max_par1 - (0.5 + adjust1),
            max_par2 - (0.5 + adjust2),
            color="firebrick",
            linestyle="--",
        )

    else:
        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(
            save_sim["par1"],
            save_sim["par2"],
            save_sim["value"],
            c=save_sim["value"],
            cmap="viridis",
            linewidth=0.05,
        )
        ax.scatter([0.5], [0.5], [total_deceased], c="r", marker="x", linewidth=5)
        ax.set_zlabel("Deceased")
        save_path = f"{path}_{para_numb}"

    xlab = f"mRNA week {para_numb}"
    ylab = f"Vector week {para_numb}"
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)

    fig.savefig(save_path, bbox_inches="tight")


for index in range(10):
    plot_simulation_para(index, three_d=False)
    # plot_simulation_para(index, three_d=True)
