from pymoo.visualization.scatter import Scatter
from pymoo.factory import get_reference_directions
from pymoo.factory import get_termination
from pymoo.factory import get_sampling, get_crossover, get_mutation
from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.algorithms.moo.ctaea import CTAEA
from pymoo.core.problem import ElementwiseProblem
import pymoo.optimize
import pandas as pd
import amici
from estimagic import minimize
from functions.tools import set_parameter_by_name
from functions.tools import get_sum_of_states
from functions.tools import get_start_splines_pop_based
from functions.tools import set_end_values_to_start
from functions.tools import set_start_states_default
from functions.tools import run_interventions_model
from functions.tools import run_interventions_pybobyqa
from functions.tools import get_R_spline_values
from functions.tools import get_infected_true
import numpy as np
import pickle

from models.vaccination.create_model_vaccination import create_model_splines
from functions.create_vaccine_doses_inflow import create_inflow_from_data
from functions.tools import get_model_output

create_model = True
model_name = "vaccination_multi_six_recovered_test"
path_sbml = f"stored_models/{model_name}/" + model_name
model_directory = "stored_models/" + model_name + "/vaccination_dir"

max_T = 28  # 112 works
number_yy = int(max_T / 14 + 1)
created_model = create_model_splines(
    create_model=create_model,
    number_areas=2,
    number_vaccines=2,
    vaccinated_compartments=["susceptible", "infectious"],
    number_viruses=1,
    length_decision_period=max_T,
    number_yy=int(max_T / 14 + 1),
    model_name=model_name,
    path_sbml=path_sbml,
    model_directory=model_directory,
    # R0_periods=4,
)

model = created_model["model"]
solver = created_model["solver"]

areas = created_model["information"]["areas"]

# Countries:
# - Belgium, France, Germany, United Kingdom
# https://ec.europa.eu/eurostat/documents/2995521/11081093/3-10072020-AP-EN.pdf/d2f799bf-4412-05cc-a357-7b49b93615f1
infectious_t0 = {  # *1.2 due to undetected cases: https://www.cdc.gov/coronavirus/2019-ncov/cases-updates/burden.html
    # https://covid19.who.int/region/euro/country/be
    "infectious_countryA_vac0_virus1_t0": 67000,
    # https://covid19.who.int/region/euro/country/fr
    "infectious_countryB_vac0_virus1_t0": 0,
    # https://covid19.who.int/region/euro/country/de
    # "infectious_countryC_vac0_virus1_t0": 377000,
    # https://covid19.who.int/region/euro/country/gb
    # "infectious_countryD_vac0_virus1_t0": 1233000,
    # mid april; https://www.nature.com/articles/d41586-021-01696-3 ; https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/
}

recovered_t0 = {
    "recovered_countryA_vac0_virus1_t0": 651225 * 1.2
    - 19834,  # https://covid19.who.int/region/euro/country/be
    # https://covid19.who.int/region/euro/country/fr
    "recovered_countryB_vac0_virus1_t0": 2608686 * 1.2 - 64543,
}
# https://covid19.who.int/region/euro/country/de
# "recovered_countryC_vac0_virus1_t0": 1765666*1.2 - 34272,
# "recovered_countryD_vac0_virus1_t0": 2688021*1.2 - 74570, }  # https://covid19.who.int/region/euro/country/gb


start_population = {
    "susceptible_countryA_vac0_t0": 11455.5 * 1000
    - infectious_t0["infectious_countryA_vac0_virus1_t0"]
    - recovered_t0["recovered_countryA_vac0_virus1_t0"],
    "susceptible_countryB_vac0_t0": 67012.9 * 1000
    - infectious_t0["infectious_countryB_vac0_virus1_t0"]
    - recovered_t0["recovered_countryB_vac0_virus1_t0"],
}
# 'susceptible_countryC_vac0_t0': 83019.2 * 1000 - infectious_t0["infectious_countryC_vac0_virus1_t0"] - recovered_t0["recovered_countryC_vac0_virus1_t0"],
#'susceptible_countryD_vac0_t0': 66647.1*1000 - infectious_t0["infectious_countryD_vac0_virus1_t0"] - recovered_t0["recovered_countryD_vac0_virus1_t0"], }

pop_sizes = np.array(list(infectious_t0.values())) + np.array(
    list(start_population.values())
)
relative_population = dict(zip(areas, pop_sizes / np.sum(pop_sizes)))

length = 420
grid_length = 6000
total_grid = int(length / max_T * grid_length)
time_grid = np.linspace(0, length, int(length / max_T * grid_length))

vaccine_supply_parameter_values = create_inflow_from_data(
    column_interest="FirstDose",
    number_decision_periods=number_yy - 1,
    weeks_decision_period=int(length / 7 / (number_yy - 1)),
    model_population=np.sum(pop_sizes),
)

vaccine_supply_parameter_strings = [
    x for x in model.getFixedParameterNames() if "vaccine_supply" in x
]

vaccine_inflow = dict(
    zip(vaccine_supply_parameter_strings, vaccine_supply_parameter_values * 100)
)

general_parameters = {
    "lambda1": 0.1,
    "prob_deceasing": 0.02,
    "gamma": 0.5,
    "beta": 0.24,
}

# http://www.healthdata.org/covid/covid-19-vaccine-efficacy-summary
delta = {
    "delta_vac1_virus1": 0.87,
    # "delta_vac1_virus2": 0.78,
    "delta_vac2_virus1": 0.62,
    # "delta_vac2_virus2": 0.53,
}

omega = {
    "omega_vac1_virus1": 0.94,
    # "omega_vac1_virus2": 0.89,
    "omega_vac2_virus1": 0.72,
    # "omega_vac2_virus2": 0.50,
}


# https://www.cdc.gov/coronavirus/2019-ncov/variants/delta-variant.html
# eta = {"eta_virus2": 2}


# R-values country

d = 1 / 1000
r = 1 / 1500
distances = {  # should be redefined
    "distance_countryA_countryB": d,
    "distance_countryB_countryA": d,
    # "distance_countryA_countryC": d,
    # "distance_countryC_countryA": d,
    # "distance_countryA_countryD": r,
    # "distance_countryD_countryA": d,
    # "distance_countryB_countryC": d,
    # "distance_countryC_countryB": d,
    # "distance_countryB_countryD": d,
    # "distance_countryD_countryB": d,
    # "distance_countryC_countryD": r,
    # "distance_countryD_countryC": r,
}


spline_xx = {
    "xx0": 0,
    "xx1": length / (number_yy - 1),
    "xx2": length / (number_yy - 1) * 2,
    #'xx3': length / (number_yy - 1) * 3,
    #'xx4': length / (number_yy - 1) * 4,
    #'xx5': length / (number_yy - 1) * 5,
    #'xx6': length / (number_yy - 1) * 6,
    #'xx7': length / (number_yy - 1) * 7,
    #'xx8': length / (number_yy - 1) * 8,
}

par_R = get_R_spline_values(
    areas,
    spline_xx,
    url="/home/manuel/Documents/VaccinationDistribution/mobility.ods",
    sheet_name="Sheet1",
    countries=["Belgium", "France", "Germany", "UK"],
    scale=0.85,
    end=0.45,
)
# par_R["yR0_countryC_6"] = 0.2
# par_R["yR0_countryC_7"] = 0.2
# par_R["yR0_countryC_8"] = 0.2
# par_R["yR0_countryB_6"] = 0.2
# par_R["yR0_countryB_7"] = 0.2
# par_R["yR0_countryB_7"] = 0.2
# par_R["yR0_countryA_6"] = 0.3
# par_R["yR0_countryA_7"] = 0.3
# par_R["yR0_countryA_8"] = 0.1
# par_R["yR0_countryD_7"] = 0.2
# par_R["yR0_countryD_8"] = 0.2


parameters = {
    **spline_xx,
    **infectious_t0,
    **start_population,
    **general_parameters,
    # **eta,
    **distances,
    **vaccine_inflow,
    **par_R,
}


interventionA = {
    "t": 20,  # first detected at day 53
    "parameter": {
        "infectious_countryB_vac0_virus1_t0": 20,
    },
}
interventions = [interventionA]


# results_pybobyqa = run_interventions_pybobyqa(model, solver, interventions, parameters, total_grid, length,
#                                              areas, number_yy, par_optimize, lower_bound = 10**-15, n_starts=1)

rim = run_interventions_model(
    model,
    solver,
    interventions,
    parameters,
    total_grid,
    length,
    areas,
    number_yy,
    theta=None,
)
rim["trajectories"]["t"] = time_grid[0 : (len(time_grid) - 1)]
import matplotlib.pyplot as plt

# beautify plot!
df_infected = pd.DataFrame(time_grid[0 : (len(time_grid) - 1)], columns=["Time"])
for i in areas:
    df_infected[i] = (
        rim["trajectories"][
            [x for x in rim["trajectories"].columns if ("infectious" in x) & (i in x)]
        ]
        .sum(axis=1)
        .reset_index(drop=True)
    )

df_inf_true = get_infected_true()
countries = list(df_inf_true.columns)
fig, axs = plt.subplots(2, 2)
for j in range(len(countries)):
    if j < 2:
        j_m = 0
    else:
        j_m = 1
    j_n = j % 2
    if j_n == 0:
        axs[j_m][j_n].set_ylabel("Active cases")

    axs[j_m][j_n].plot(
        df_infected["Time"] / 7, df_infected[areas[j]], label="Simulated Infectious"
    )
    axs[j_m][j_n].plot(
        np.array(list(df_inf_true.index)) / 7,
        df_inf_true[countries[j]],
        label="True infectious",
    )
    axs[j_m][j_n].set_title(countries[j].capitalize())
    axs[j_m][j_n].set_xlabel("Weeks")
handles, labels = axs[0][0].get_legend_handles_labels()
fig.legend(
    handles=handles,
    labels=labels,
    loc="lower center",
    ncol=4,
    bbox_to_anchor=[0.5, -0.1],
)
plt.tight_layout()
fig.savefig(
    "/home/manuel/Documents/VaccinationDistribution/paper/images/infected_compare",
    bbox_inches="tight",
)


par_optimize = model.getParameterNames()
pop_size = 10
number_generations = 100
seed = 1234
crossover_eta = 15
crossover_prob = 0.9
mutation_eta = 20
lb = 10 ** -15


set_parameter_by_name(model, parameters, fixed=True)


transformed_frac_pareto = get_start_splines_pop_based(model, areas, number_yy)
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
    return d[1 : len(d)]


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
                help_dict[j] = relative[l] / np.sum(list(relative.values()))

    theta_par = np.array(list(help_dict.values()))
    pareto_values.append(get_optim_function_pareto(theta_par)[i])


# estimagic
index_one = []
index_two = []
vaccines = created_model["information"]["vaccine_states"]

for area in range(len(areas) - 1):
    count = 0
    for j in vaccines[1 : (len(vaccines))]:
        for i in range(int((number_par / (len(vaccines) - 1) / (len(areas) - 1)))):
            name = f"{j}_par{i}"
            list_name = [name]
            index_one += list_name
            count += 1

    index_two += [areas[area]] * count

tuples = list(zip(index_one, index_two))
m_index = pd.MultiIndex.from_tuples(tuples, names=["first", "second"])
start = np.repeat(0.1, number_par)
lower = np.repeat(lb, number_par)
upper = np.repeat(1 - lb, number_par)


start_params = pd.DataFrame(
    np.array([start, lower]).T, index=m_index, columns=["value", "lower_bound"]
)


def get_optim_function_estimagic(theta):
    return {"value": np.sum(get_optim_function_pareto(theta["value"]))}


set_index_one = list(set(index_one))
constraints = []
for i in set_index_one:
    dictionary = {"loc": i, "type": "linear", "weights": 1, "upper": 1}
    constraints.append(dictionary)


res = minimize(
    criterion=get_optim_function_estimagic,
    params=start_params,
    algorithm="nelder-mead",
    constraints=constraints,
)


# scipy
def get_optim_function_scipy(theta):
    return np.sum(get_optim_function_pareto(theta))


index_one = []
index_two = []
vaccines = created_model["information"]["vaccine_states"]
count = 0
for j in vaccines[1 : (len(vaccines))]:
    for i in range(int((number_par / (len(vaccines) - 1) / (len(areas) - 1)))):
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
lb = 0
ub = 1

from scipy.optimize import LinearConstraint
from scipy.optimize import minimize

constraints = LinearConstraint(A, lb, ub)

scipy_min = minimize(
    fun=get_optim_function_scipy,
    x0=np.repeat(0.3, number_par),
    method="nelder-mead",
    constraints=constraints,
)


# create linear constraints
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
                    float(i.rsplit("_", 1)[-1]) == float(j.rsplit("_", 1)[-1])
                ):
                    s += dict_para[j]
            sums.append(s)
        # out["G"] = np.concatenate([out["F"] - np.array(pareto_values) - 10**6, np.array(sums) - 1])
        out["G"] = np.array(sums) - 1


problem = multi_objective_problem()

algorithm = NSGA2(
    pop_size=pop_size,
    sampling=get_sampling("real_random"),
    crossover=get_crossover("real_sbx", prob=crossover_prob, eta=crossover_eta),
    mutation=get_mutation("real_pm", eta=mutation_eta),
    eliminate_duplicates=True,
)
ref_dirs = get_reference_directions("das-dennis", 4, n_partitions=8)
algorithm = NSGA3(
    ref_dirs=ref_dirs,
)
# pop_size=455)
algorithm = CTAEA(ref_dirs=ref_dirs)


termination = get_termination("n_gen", number_generations)

res = pymoo.optimize.minimize(
    problem, algorithm, termination, seed=seed, save_history=True, verbose=True
)


f = pd.DataFrame(res.F)
x = pd.DataFrame(res.X)
p = np.array(pareto_values * x.shape[0]).reshape(x.shape[0], 4)

r = f - p
d = np.diag(1 / np.array(pareto_values))
percentage = r @ d * 100

dict_out = {"f": f, "x": x, "p": p}
path = "/home/manuel/Documents/VaccinationDistribution/code/objects/moo_dict.pkl"

with open(
    path,
    "wb",
) as output:
    out = dict_out
    pickle.dump(out, output, pickle.HIGHEST_PROTOCOL)


with open(
    path,
    "rb",
) as input:
    di = pickle.load(input)
f = np.array(di["f"]) / 10 ** 6
plot = Scatter(tight_layout=True, labels=["Belgium", "France", "Germany", "UK"])
plot.add(f, s=10)
# plot.add(F[10], s=30, color="red")
plot.show()

h = dict(zip(model.getParameterNames(), model.getParameters()))
