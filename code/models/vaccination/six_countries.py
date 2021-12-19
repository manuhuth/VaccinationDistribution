import numpy as np
import pickle
import pandas as pd
import multiprocessing as mp
import tqdm

from functions.tools import set_parameter_by_name
from functions.tools import get_start_splines_pop_based
from functions.tools import run_interventions_model
from functions.tools import get_infected_true
from functions.tools import create_R_rule_points
from functions.tools import get_spline_xx
from models.vaccination.spline_example import get_spline
from functions.plot_tools import plot_four_country_overview
from functions.tools import run_pareto_front
from functions.tools import create_constraints
from scipy.optimize import Bounds
from scipy.optimize import LinearConstraint
from scipy.optimize import minimize
from models.vaccination.create_model_vaccination import create_model_splines
from functions.create_vaccine_doses_inflow import create_inflow_from_data


create_model = False
model_name = "vaccination_multi_six_recovered"
path_sbml = f"stored_models/{model_name}/" + model_name
model_directory = "stored_models/" + model_name + "/vaccination_dir"

max_T = 70  # 112 works
number_yy = int(max_T / 14 + 1)
number_xx_R = int(112 / 14 + 1)
created_model = create_model_splines(
    create_model=create_model,
    number_areas=4,
    number_vaccines=2,
    vaccinated_compartments=["susceptible", "infectious", "recovered"],
    number_viruses=2,
    length_decision_period=max_T,
    number_yy=number_yy,
    number_xx_R=number_xx_R,
    model_name=model_name,
    path_sbml=path_sbml,
    model_directory=model_directory,
    # R0_periods=4,
)

model = created_model["model"]
solver = created_model["solver"]
solver.setMaxSteps(10000)
solver.setAbsoluteTolerance(1e-01)
solver.setRelativeTolerance(1e-01)
# solver.getSensitivityOrder()

areas = created_model["information"]["areas"]

# Countries:
# - Belgium, France, Germany, United Kingdom
# https://ec.europa.eu/eurostat/documents/2995521/11081093/3-10072020-AP-EN.pdf/d2f799bf-4412-05cc-a357-7b49b93615f1
infectious_t0 = {  # *1.2 due to undetected cases: https://www.cdc.gov/coronavirus/2019-ncov/cases-updates/burden.html
    # https://covid19.who.int/region/euro/country/be
    "infectious_countryA_vac0_virus1_t0": 67000,
    # https://covid19.who.int/region/euro/country/fr
    "infectious_countryB_vac0_virus1_t0": 257000,
    # https://covid19.who.int/region/euro/country/de
    "infectious_countryC_vac0_virus1_t0": 377000,
    # https://covid19.who.int/region/euro/country/gb
    "infectious_countryD_vac0_virus1_t0": 1233000,
    # mid april; https://www.nature.com/articles/d41586-021-01696-3 ; https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/
}

recovered_t0 = {
    "recovered_countryA_vac0_virus1_t0": 651225 * 1.2
    - 19834,  # https://covid19.who.int/region/euro/country/be
    # https://covid19.who.int/region/euro/country/fr
    "recovered_countryB_vac0_virus1_t0": 2608686 * 1.2 - 64543,
    # https://covid19.who.int/region/euro/country/de
    "recovered_countryC_vac0_virus1_t0": 1765666 * 1.2 - 34272,
    "recovered_countryD_vac0_virus1_t0": 2688021 * 1.2 - 74570,
}  # https://covid19.who.int/region/euro/country/gb


start_population = {
    "susceptible_countryA_vac0_t0": 11455.5 * 1000
    - infectious_t0["infectious_countryA_vac0_virus1_t0"]
    - recovered_t0["recovered_countryA_vac0_virus1_t0"],
    "susceptible_countryB_vac0_t0": 67012.9 * 1000
    - infectious_t0["infectious_countryB_vac0_virus1_t0"]
    - recovered_t0["recovered_countryB_vac0_virus1_t0"],
    "susceptible_countryC_vac0_t0": 83019.2 * 1000
    - infectious_t0["infectious_countryC_vac0_virus1_t0"]
    - recovered_t0["recovered_countryC_vac0_virus1_t0"],
    "susceptible_countryD_vac0_t0": 66647.1 * 1000
    - infectious_t0["infectious_countryD_vac0_virus1_t0"]
    - recovered_t0["recovered_countryD_vac0_virus1_t0"],
}

pop_sizes = np.array(list(infectious_t0.values())) + np.array(
    list(start_population.values())
)
relative_population = dict(zip(areas, pop_sizes / np.sum(pop_sizes)))

length = 420
grid_length = 1000
total_grid = int(length / max_T * grid_length)
time_grid = np.linspace(0, length, total_grid)

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
    zip(vaccine_supply_parameter_strings, vaccine_supply_parameter_values)
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
    "delta_vac1_virus2": 0.78,
    "delta_vac2_virus1": 0.62,
    "delta_vac2_virus2": 0.53,
}

omega = {
    "omega_vac1_virus1": 0.94,
    "omega_vac1_virus2": 0.89,
    "omega_vac2_virus1": 0.72,
    "omega_vac2_virus2": 0.50,
}


# https://www.cdc.gov/coronavirus/2019-ncov/variants/delta-variant.html
eta = {"eta_virus2": 2.3}


# R-values country

d = 1 / 1000
r = 1 / 1500
distances = {  # should be redefined
    "distance_countryA_countryB": d,
    "distance_countryB_countryA": d,
    "distance_countryA_countryC": d,
    "distance_countryC_countryA": d,
    "distance_countryA_countryD": r,
    "distance_countryD_countryA": d,
    "distance_countryB_countryC": d,
    "distance_countryC_countryB": d,
    "distance_countryB_countryD": d,
    "distance_countryD_countryB": d,
    "distance_countryC_countryD": r,
    "distance_countryD_countryC": r,
}


df_inf_true = get_infected_true()
countries = list(df_inf_true.columns)
end_data = list(df_inf_true.index)[-1]

spline_xx = get_spline_xx(
    leave_out=2, length=length, end_data=end_data, number_yy=number_yy, key_str="xx"
)
spline_xx_R = get_spline_xx(
    leave_out=2,
    length=length,
    end_data=end_data,
    number_yy=number_xx_R,
    key_str="xx_R0_",
)

maxi = 0.75
df_stringency = pd.read_csv("/home/manuel/Documents/VaccinationDistribution/stringency.csv")
df_stringency["Time"] = np.linspace(0, df_stringency.shape[0]-1, df_stringency.shape[0])
par_R = {}
for i in range(len(areas)):
    for j in range(number_xx_R):
        key = f"yR0_{areas[i]}_{j}"
        time_key = f"xx_R0_{j}"
        time_value = spline_xx_R[time_key]
        closest_val = df_stringency.loc[
                np.argmin(np.abs(df_stringency["Time"] - time_value)), countries[i].capitalize()
        ]
        par_R[key] = (1 - closest_val / 100 * maxi)  
        
#par_R = create_R_rule_points(areas, number_xx_R, spline_xx_R)

parameters = {
    **spline_xx,
    **spline_xx_R,
    **infectious_t0,
    **start_population,
    **general_parameters,
    **delta,
    **omega,
    **eta,
    **distances,
    **vaccine_inflow,
    **par_R,
}


interventionA = {
    "t": 46,  # first detected at day 53
    "parameter": {
        "infectious_countryD_vac0_virus2_t0": 20,
        "infectious_countryB_vac0_virus2_t0": 1,
        "infectious_countryC_vac0_virus2_t0": 20,
    },
}
interventionB = {
    "t": 76,  # first detected at day 53
    "parameter": {
        "infectious_countryA_vac0_virus2_t0": 20,
    },
}
interventions = [interventionA, interventionB]


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
rim["trajectories"]["t"] = time_grid[0 : (len(time_grid) - 1)] / 7
df_infected = pd.DataFrame(time_grid[0 : (len(time_grid) - 1)], columns=["Time"])
for i in areas:
    df_infected[i] = (
        rim["trajectories"][
            [x for x in rim["trajectories"].columns if ("infectious" in x) & (i in x)]
        ]
        .sum(axis=1)
        .reset_index(drop=True)
    )


vaccine_available = plot_four_country_overview(
    spline_xx,
    vaccine_inflow,
    number_yy,
    length,
    interventionA,
    interventionB,
    end_data,
    start_population,
    infectious_t0,
    recovered_t0,
    delta,
    omega,
    countries,
    areas,
    par_R,
    number_xx_R,
    total_grid,
    df_inf_true,
    df_infected,
    grid_data=np.linspace(0, end_data, int(total_grid * end_data / length)),
    grid_sim=np.linspace(
        end_data, length, total_grid - int(total_grid * end_data / length) - 1
    ),
    scale=10 ** 6,
    text_x=end_data / 7 / 3 - 2,
    text_y=0.07,
    ylim=[0, 0.8],
    text_str="Vaccination \nperiod",
    text_lockdown_x=((length - end_data) / 2 + end_data) / 7 - 7,
    text_lockdown_y=0.03,
    text_lockdown_str="Constant \nNPIs",
    color_prop=["C4", "C5", "C6", "C7"],
    label_vac1="mRNA",
    label_vac2="vector",
    color_vac1="C7",
    color_vac2="C8",
    title_vac="Available vaccines",
    title_setup="Set-up",
    position_start_vac=[0, 0.5],
    height_start_vac=0.3,
    letter_size=30,
    letter_y=1.06,
    size=(40, 20),
)

par_optimize = [x for x in model.getParameterNames() if int(x.rsplit("_", 1)[-1]) <= 3]


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


# dict_out = run_pareto_front(
#    par_optimize,
#    areas,
#    model,
#    number_yy,
#    interventions,
#    total_grid,
#    length,
#    parameters,
#   solver,
#    relative_population,
#    rim,
#    pop_size=100,
#    number_generations=500,
#    seed=1234,
#    crossover_eta=15,
#    crossover_prob=0.9,
#    mutation_eta=20,
#    lb=10 ** -15,
# )
# path = "/home/manuel/Documents/VaccinationDistribution/code/objects/manyoo_dict.pkl"

# with open(
#    path,
#    "wb",
# ) as output:
#    out = dict_out
#    pickle.dump(out, output, pickle.HIGHEST_PROTOCOL)

# start = np.array(list(get_start_splines_pop_based(model, areas, number_yy).values()))

# scipy
constraints = create_constraints(
    vaccines=created_model["information"]["vaccine_states"],
    par_optimize=par_optimize,
    areas=areas,
)

lb = 10 ** -9
ub = 1 - lb



theta_p = get_start_splines_pop_based(model, areas, number_yy)
keys = [x for x in theta_p.keys() if int(x.rsplit("_", 1)[-1]) <= 3]
dic = {}
for index in keys:
    dic[index] = theta_p[index]

z = np.array(list(dic.values()))
z_trans = 1 / (1 + np.exp(-z))


def get_optim_function_scipy(theta):
    return np.sum(get_optim_function_pareto(theta))


def draw_random_start(i_sim):
    save_sim = None
    np.random.seed(i_sim)
    B = constraints["lin_constraints"].copy()
    for i in B.index:
        for j in B.columns:
            col_sum = B.copy().replace(1, 0).sum(axis=1)[i]
            if B.loc[i][j] == 1:
                for c in range(len(areas)):
                    if areas[c] in j:
                        maxi = 1.2 * relative_population[areas[c]]
                B.loc[i, j] = np.random.uniform(lb, np.min([maxi, 1 - col_sum]))

    theta_par = B.sum()
    it_res = np.sum(get_optim_function_pareto(theta_par))
    #if it_res < 8 * 10 ** 5:
    save_sim = np.concatenate([theta_par, [it_res]])

    return save_sim


n_sim = 65000
with mp.Pool() as p:
    starts = list(
        tqdm.tqdm(p.imap_unordered(draw_random_start, range(0, n_sim)), total=n_sim)
    )

starts.append(np.concatenate([z_trans, [np.sum(get_optim_function_pareto(z_trans))]]))

start = pd.DataFrame(
    [np.concatenate([z_trans, [np.sum(get_optim_function_pareto(z_trans))]])],
    columns=par_optimize + ["f"],
)
for index in range(len(starts)):
    if all(np.repeat(starts[index] == starts[index], 2)):
        start = start.append(
            pd.DataFrame([starts[index]], columns=par_optimize + ["f"]),
            ignore_index=True,
        )

start = start.sort_values("f").reset_index(drop=True)

path = "/home/manuel/Documents/VaccinationDistribution/code/objects/starts.pkl"

with open(
    path,
    "rb",
) as input:
    start = pickle.load(input)

with open(
    path,
    "wb",
) as output:
    out = start
    pickle.dump(out, output, pickle.HIGHEST_PROTOCOL)


start_use = start.loc[0:49, :]
bounds = Bounds(np.repeat(lb, len(par_optimize)), np.repeat(ub, len(par_optimize)))
constraint = LinearConstraint(constraints["lin_constraints"], lb, ub)

def run_scipy_min_parallel(index):
    return minimize(
        fun=get_optim_function_scipy,
        x0=start_use[par_optimize].loc[index, :],
        method="trust-constr",
        constraints=constraint,
        bounds=bounds,
        options={"verbose": 2, "maxiter": 200},
    )


with mp.Pool() as p:
    results_unconstrained = list(
        tqdm.tqdm(
            p.imap_unordered(run_scipy_min_parallel, range(0, start_use.shape[0])),
            total=start_use.shape[0],
        )
    )

par_areas_vac1 = []
for i in range(len(areas)):
    par_areas_vac1.append(
        [x for x in par_optimize if (areas[i] in x) and ("vac1" in x)]
    )

par_areas_vac2 = []
for i in range(len(areas)):
    par_areas_vac2.append(
        [x for x in par_optimize if (areas[i] in x) and ("vac2" in x)]
    )

par_areas = {"vac1": par_areas_vac1, "vac2": par_areas_vac2}


r = results_unconstrained
par = par_optimize
df_results = pd.DataFrame(index=range(len(r)), columns=par + areas + ["f"])
ind = 0
for index in range(len(r)):
    paramet = r[index]["x"]
    par_d = dict(zip(par, paramet))
    for vac in ["vac1", "vac2"]:
        vac_ind = 0
        val = {}
        greater_zero = []
        for area in range(len(areas) - 1):
            pars = par_areas[vac][area]
            para = [par_d[x] for x in pars if 2 > 1]
            if any(np.array(para) <= 0) or any(np.array(para) >= 1):
                spline = np.repeat(1000, 5000)
            else:
                xx = np.array(list(spline_xx.values()))
                rel_pop_trans = relative_population[areas[area]]
                spline = get_spline(
                    array=np.concatenate(
                        [para, np.repeat(rel_pop_trans, len(xx) - len(pars))]
                    ),
                    periods=None,
                    length=None,
                    total_length=length,
                    grid_points=None,
                    transform=True,
                    x=xx,
                    grid=vaccine_available["t"],
                )
            val[areas[area]] = spline
        sum_df = 1 - pd.DataFrame(val).sum(axis=1)
        greater_zero.append(all(sum_df > 0))

    if all(greater_zero) is True:
        areas_results = get_optim_function_pareto(paramet)
        func_value = np.sum(areas_results)
        row = np.concatenate([paramet, areas_results, [func_value]])
        df_results.iloc[ind] = row
        ind += 1
df_results = df_results.dropna()



#pareto constraints
pareto = np.array(get_optim_function_pareto(z_trans))


def get_optim_function_scipy_constrained(theta):
    it = get_optim_function_pareto(theta)

    if all(it < pareto):
        return np.sum(get_optim_function_pareto(theta))
    else:
        bool_vec = it > pareto
        diff = it - pareto
        return np.sum(it) + 1000 * diff @ bool_vec


def draw_random_start_constr(i_sim):
    save_sim = None
    np.random.seed(i_sim)
    B = constraints["lin_constraints"].copy()
    for i in B.index:
        for j in B.columns:
            col_sum = B.copy().replace(1, 0).sum(axis=1)[i]
            if B.loc[i][j] == 1:
                for c in range(len(areas)):
                    if areas[c] in j:
                        maxi = 1.3 * relative_population[areas[c]]
                B.loc[i, j] = np.random.uniform(lb, np.min([maxi, 1 - col_sum]))

    theta_par = B.sum()
    it_res = get_optim_function_scipy_constrained(theta_par)
    
    save_sim = np.concatenate([theta_par, [it_res]])

    return save_sim


n_sim = 650000
with mp.Pool() as p:
    starts_constr = list(
        tqdm.tqdm(
            p.imap_unordered(draw_random_start_constr, range(0, n_sim)), total=n_sim
        )
    )

start_constr = pd.DataFrame(
    [np.concatenate([z_trans, [get_optim_function_scipy_constrained(z_trans)]])],
    columns=par_optimize + ["f"],
)
for index in range(len(starts_constr)):
    if not (starts_constr[index] is None):
        start_constr = start_constr.append(
            pd.DataFrame([starts_constr[index]], columns=par_optimize + ["f"]),
            ignore_index=True,
        )

start_constr = start_constr.sort_values("f").reset_index(drop=True)

path = "/home/manuel/Documents/VaccinationDistribution/code/objects/starts_constr.pkl"

with open(
    path,
    "rb",
) as input:
    start_constr = pickle.load(input)

with open(
    path,
    "wb",
) as output:
    out = start_constr
    pickle.dump(out, output, pickle.HIGHEST_PROTOCOL)

start_constr_use = start_constr.loc[0:49, :]

constraint = LinearConstraint(constraints["lin_constraints"], lb, ub)
def run_scipy_min_parallel_constrained(index):
    return minimize(
        fun=get_optim_function_scipy_constrained,
        x0=start_constr_use[par_optimize].loc[index, :],
        method="trust-constr",
        constraints=constraint,
        options={"verbose": 2, "maxiter": 100},
    )


with mp.Pool() as p:
    results_constrained = list(
        tqdm.tqdm(
            p.imap_unordered(
                run_scipy_min_parallel_constrained, range(0, start_constr_use.shape[0])
            ),
            total=start_constr_use.shape[0],
        )
    )

r = results_constrained
par = par_optimize
df_results_constr = pd.DataFrame(index=range(len(r)), columns=par + areas + ["f"])
ind = 0
for index in range(len(r)):
    paramet = r[index]["x"]
    par_d = dict(zip(par, paramet))
    for vac in ["vac1", "vac2"]:
        vac_ind = 0
        val = {}
        greater_zero = []
        for area in range(len(areas) - 1):
            pars = par_areas[vac][area]
            para = [par_d[x] for x in pars if 2 > 1]
            if any(np.array(para) <= 0) or any(np.array(para) >= 1):
                spline = np.repeat(1000, 5000)
            else:
                xx = np.array(list(spline_xx.values()))
                rel_pop_trans = relative_population[areas[area]]
                spline = get_spline(
                    array=np.concatenate(
                        [para, np.repeat(rel_pop_trans, len(xx) - len(pars))]
                    ),
                    periods=None,
                    length=None,
                    total_length=length,
                    grid_points=None,
                    transform=True,
                    x=xx,
                    grid=vaccine_available["t"],
                )
            val[areas[area]] = spline
        sum_df = 1 - pd.DataFrame(val).sum(axis=1)
        greater_zero.append(all(sum_df > 0))

    if all(greater_zero) is True:
        areas_results = get_optim_function_pareto(paramet)
        func_value = np.sum(areas_results)
        row = np.concatenate([paramet, areas_results, [func_value]])
        df_results_constr.iloc[ind] = row
        ind += 1
df_results_constr = df_results_constr.dropna()

df_results_constr_pareto = df_results_constr[
    (df_results_constr[areas[0]] <= pareto[0])
    & (df_results_constr[areas[1]] <= pareto[1])
    & (df_results_constr[areas[2]] <= pareto[2])
    & (df_results_constr[areas[3]] <= pareto[3])
].reset_index(drop=True)

df_pareto = pd.DataFrame(
    [list(df_results_constr_pareto.loc[np.argmin(df_results_constr_pareto["f"]), :])],
    columns=par_optimize + areas + ["f"],
)
df_pop = pd.DataFrame(
    [list(np.concatenate([z_trans, pareto, [np.sum(pareto)]]))],
    columns=par_optimize + areas + ["f"],
)
df_best = pd.DataFrame(
    [list(df_results.loc[np.argmin(df_results["f"]), :])],
    columns=par_optimize + areas + ["f"],
)

dfs = {"pop": df_pop, "constrained": df_pareto, "unconstrained": df_best}

vaccine_allocated = pd.DataFrame(vaccine_available["t"], columns=["t"])
for vac in ["vac1", "vac2"]:
    for type_allocation in ["pop", "constrained", "unconstrained"]:
        for area in range(len(areas[0 : (len(areas) - 1)])):
            pars = par_areas[vac][area]
            rel_pop_trans = relative_population[areas[area]]
            df = dfs[type_allocation]
            xx = np.array(list(spline_xx.values()))

            spline = get_spline(
                array=np.concatenate(
                    [df[pars].loc[0, :], np.repeat(rel_pop_trans, len(xx) - len(pars))]
                ),
                periods=None,
                length=None,
                total_length=length,
                grid_points=None,
                transform=True,
                x=xx,
                grid=vaccine_available["t"],
            )
            col_name = f"{type_allocation}_{areas[area]}_{vac}"
            vaccine_allocated[col_name] = spline

for vac in ["vac1", "vac2"]:
    for type_allocation in ["pop", "constrained", "unconstrained"]:
        cols = [x for x in vaccine_allocated if (vac in x) and (type_allocation in x)]
        if type_allocation == "constrained":
            cols = [x for x in cols if not ("unconstrained" in x)]
        c_D = 1 - vaccine_allocated[cols].sum(axis=1)
        col_name = f"{type_allocation}_countryD_{vac}"
        vaccine_allocated[col_name] = c_D

# trajectories
dfs = [df_pop, df_pareto, df_best]
names = ["pop_trajectories", "pareto_trajectories", "best_trajectories"]
dict_trajectories = {}
for i in range(len(dfs)):
    df = dfs[i]
    theta = df[par_optimize].loc[0, :]
    dict_theta = dict(zip(par_optimize, np.log(theta / (1 - theta))))
    set_parameter_by_name(model, dict_theta, fixed=False)
    trajectories = run_interventions_model(
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
    dict_trajectories[names[i]] = trajectories

# vaccine trajectories
sort_df_results = df_results.sort_values("f").reset_index(drop=True)
vaccine_all_results = pd.DataFrame(vaccine_available["t"], columns=["t"])
for vac in ["vac1", "vac2"]:
    for position_best in sort_df_results.index:
        for area in range(len(areas[0 : (len(areas) - 1)])):
            pars = par_areas[vac][area]
            rel_pop_trans = relative_population[areas[area]]
            df = sort_df_results
            xx = np.array(list(spline_xx.values()))

            spline = get_spline(
                array=np.array(
                    list(
                        np.concatenate(
                            [
                                df[pars].loc[position_best, :],
                                np.repeat(rel_pop_trans, len(xx) - len(pars)),
                            ]
                        )
                    )
                ),
                periods=None,
                length=None,
                total_length=length,
                grid_points=None,
                transform=True,
                x=xx,
                grid=vaccine_available["t"],
            )
            col_name = f"a{position_best}_{areas[area]}_{vac}"
            vaccine_all_results[col_name] = spline

for vac in ["vac1", "vac2"]:
    for position_best in sort_df_results.index:
        cols = [
            x
            for x in vaccine_all_results.columns
            if (vac in x) and (f"a{position_best}_" in x)
        ]
        c_D = 1 - vaccine_all_results[cols].sum(axis=1)
        col_name = f"a{position_best}_countryD_{vac}"
        vaccine_all_results[col_name] = c_D

dict_out = {
    "unconstr_best": df_results,
    "constr_best": df_pareto,
    "par_optimize": par_optimize,
    "pop_best": df_pop,
    "vaccine": vaccine_available,
    "allocated_best": vaccine_allocated,
    "constrained": df_results_constr_pareto,
    "results_pareto_optim_not_pareto_included": df_results_constr,
    "pareto_outcome": pareto,
    "vaccine_all_paths_best": vaccine_all_results,
    **dict_trajectories,
}

path = "/home/manuel/Documents/VaccinationDistribution/code/objects/results_trust_constr_4countries.pkl"
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
    dict_ = pickle.load(input)







# test-bang-bang:
names = []
for index_country in range(len(areas)):
    names.append(f"yy_{areas[index_country]}_vac1_0")
new_df1 = pd.DataFrame(np.diag(np.repeat(1, len(areas))), columns=names)

for index_vac in ["vac1"]:
    for index_period in range(1, 4):
        diags = []
        names = []
        for index_country in range(len(areas)):
            names.append(f"yy_{areas[index_country]}_{index_vac}_{index_period}")
        merge = pd.DataFrame(np.diag(np.repeat(1, len(areas))), columns=names)
        new_df1 = new_df1.merge(merge, how="cross")


names = []
for index_country in range(len(areas)):
    names.append(f"yy_{areas[index_country]}_vac2_0")
new_df2 = pd.DataFrame(np.diag(np.repeat(1, len(areas))), columns=names)
for index_vac in ["vac2"]:
    for index_period in range(1, 4):
        diags = []
        names = []
        for index_country in range(len(areas)):
            names.append(f"yy_{areas[index_country]}_{index_vac}_{index_period}")
        merge = pd.DataFrame(np.diag(np.repeat(1, len(areas))), columns=names)
        new_df2 = new_df2.merge(merge, how="cross")

new_df = new_df1.merge(new_df2, how="cross")
lb = 10 ** -6
start_bb = new_df[par_optimize]
start_bb = start_bb.replace(0, 0)
start_bb = start_bb.replace(1, 1 - lb)

results_bb = pd.DataFrame(
    columns=par_optimize + areas + ["f"], index=range(start_bb.shape[0])
)
for index in range(start_bb.shape[0]):
    it_res = get_optim_function_pareto(start_bb.loc[index, :])
    results_bb.loc[index, :] = np.concatenate(
        [start_bb.loc[index, :], it_res, [np.sum(it_res)]]
    )
    print(index)


def get_parallel_bb(index):
    it_res = get_optim_function_pareto(start_bb.loc[index, :])
    results = np.concatenate([start_bb.loc[index, :], it_res, [np.sum(it_res)]])

    return results


with mp.Pool() as p:
    results_bb = list(
        tqdm.tqdm(
            p.imap_unordered(get_parallel_bb, range(start_bb.shape[0])),
            total=start_bb.shape[0],
        )
    )

f_bb = []
for i in range(len(results_bb)):
    f_bb.append(results_bb[i][-1])
