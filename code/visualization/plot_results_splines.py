import pickle

from visualization.model_results import plot_states
from visualization.model_results import plot_observables
from visualization.model_results import get_substates
from visualization.model_results import get_observables_by_name

from functions.run_sbml import get_model_and_solver_from_sbml

# ------Set plot path----------------------------------------------------------
plot_path = "/home/manuel/Documents/VaccinationDistribution/paper/images/splines_"

# ----------Load optimization data---------------------------------------------
with open(
    "/home/manuel/Documents/VaccinationDistribution/code/objects/output_splines.pkl",
    "rb",
) as input:
    dict_out = pickle.load(input)
df_optimal_results = dict_out["df_optimal_results"]
model_optimal = dict_out["model_optimal"]
model_seperated = dict_out["model_seperated"]
model_current = dict_out["model_current"]
model_directory = dict_out["model_directory"]
path_sbml = dict_out["path_sbml"]
model_name = dict_out["model_name"]
compare_total = dict_out["compare_total"]
compare_A = dict_out["compare_A"]
compare_B = dict_out["compare_B"]
print(compare_total)
print(compare_A)
print(compare_B)


model_solver = model_and_solver = get_model_and_solver_from_sbml(
    path_sbml=path_sbml,
    model_name=model_name,
    model_directory=model_directory,
    only_import=True,
)
model = model_solver["model"]
observables = model.getObservableNames()

# ------------plots_model------------------------------------------------------
ylim_infectious = [0, 2.5 * 10 ** 7]
colors = ["C2", "C3"]
trajectory_states_optimal = model_optimal["states"]
trajectory_observables_optimal = model_optimal["observables"]

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
    observables, substrings=["nu_", "countryA"], include_all=True
)
observables_names_proportion_optimal = get_observables_by_name(
    observables, substrings=["proportion", "countryA"], include_all=True
)

plot_observables(
    results=model_optimal,
    model=model,
    # set_off_scientific_notation=True,
    observable_ids=observables_names_nu_optimal,
    time_name="t",
    title="Vaccination rate country A (Pareto optimal)",
    ylabel="Value",
    colors=colors,
    custom_label=[
        "Vaccine wild type A",
        "Vaccine wild type B",
        "Vaccine mutant A",
        "Vaccine mutant B",
    ],
)

plot_observables(
    results=model_optimal,
    model=model,
    observable_ids=["pline_countryA_vac1", "pline_countryA_vac2"],
    # set_off_scientific_notation=True,
    decimal_floats=4,
    title="Splines country A (Pareto optimal)",
    time_name="t",
    ylabel="Value",
    colors=colors,
    custom_label=["Vaccine wild type", "Vaccine mutant"],
)

# TODO create plot for proportion*vacc available (already implemented as
# observables but need to create model again)

# Proportion_observables
fig_obs, ax_obs = plot_observables(
    results=model_optimal,
    model=model,
    observable_ids=observables_names_proportion_optimal,
    # set_off_scientific_notation=True,
    decimal_floats=1,
    time_name="t",
    title="Vaccines assigned country A (Pareto optimal)",
    ylabel="Fraction",
    xlabel="Days",
    colors=colors,
    custom_label=["Vaccine wild type", "Vaccine mutant"],
    legend_next_to_plot=False,
)

# ax_obs.get_legend().remove()
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
    results=model_optimal,
    model=model,
    state_ids=substates_optimal_A,
    time_name="t",
    title="Unvaccinated infectious individuals (A) (Pareto optimal)",
    colors=["C2", "C3"],
    xlabel="Days",
    ylabel="Number of individuals",
    custom_label=["Virus wild type", "Virus mutant"],
    ylim=ylim_infectious,
    legend_next_to_plot=False,
)
# ax.get_legend().remove()
fig.savefig(plot_path + "infectious_A", bbox_inches="tight")
# -------

# ------- Infectious B
substates_optimal_B = get_substates(
    model=model, substrings=["vac0", "infectious", "countryB"], include_all=True
)

fig, ax = plot_states(
    results=model_optimal,
    model=model,
    state_ids=substates_optimal_B,
    time_name="t",
    title="Unvaccinated infectious individuals (B) (Pareto optimal)",
    colors=["C2", "C3"],
    xlabel="Days",
    ylabel="Number of individuals",
    custom_label=["Virus wild type", "Virus mutant"],
    ylim=ylim_infectious,
    legend_next_to_plot=False,
)
# ax.get_legend().remove()
fig.savefig(plot_path + "infectious_B", bbox_inches="tight")
# --------
# ------- Infectious A
substates_optimal_A = get_substates(
    model=model, substrings=["vac0", "countryA"], include_all=True
)

fig, ax = plot_states(
    results=model_optimal,
    model=model,
    state_ids=substates_optimal_A,
    time_name="t",
    title="Unvaccinated infectious individuals (A) (Pareto optimal)",
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
    # ylim = ylim_infectious,
    legend_next_to_plot=True,
)
# ax.get_legend().remove()
fig.savefig(plot_path + "unvaccinated_A", bbox_inches="tight")


# ------- Infectious B
substates_optimal_B = get_substates(
    model=model, substrings=["vac0", "countryB"], include_all=True
)

fig, ax = plot_states(
    results=model_optimal,
    model=model,
    state_ids=substates_optimal_B,
    time_name="t",
    title="Unvaccinated infectious individuals (B) (Pareto optimal)",
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
    # ylim = ylim_infectious,
    legend_next_to_plot=True,
)
# ax.get_legend().remove()
fig.savefig(plot_path + "unvaccinated_B", bbox_inches="tight")


# plot seperated
fig_obs, ax_obs = plot_observables(
    results=model_current,
    model=model,
    observable_ids=observables_names_proportion_optimal,
    set_off_scientific_notation=True,
    decimal_floats=4,
    time_name="t",
    title="Equiallocated fraction of vaccines assigned to country A (Current)",
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
    results=model_current,
    model=model,
    state_ids=substates_optimal_A,
    time_name="t",
    title="Unvaccinated infectious individuals (A) (Current)",
    colors=["C2", "C3"],
    xlabel="Days",
    ylabel="Number of individuals",
    custom_label=["Virus wild type", "Virus mutant"],
    ylim=ylim_infectious,
    legend_next_to_plot=False,
)
# ax.get_legend().remove()
fig.savefig(plot_path + "infectious_A_current", bbox_inches="tight")
# -------
# ------- Infectious B
substates_optimal_B = get_substates(
    model=model, substrings=["vac0", "infectious", "countryB"], include_all=True
)

fig, ax = plot_states(
    results=model_current,
    model=model,
    state_ids=substates_optimal_B,
    time_name="t",
    title="Unvaccinated infectious individuals (B) (Current)",
    colors=["C2", "C3"],
    xlabel="Days",
    ylabel="Number of individuals",
    custom_label=["Virus wild type", "Virus mutant"],
    ylim=ylim_infectious,
    legend_next_to_plot=False,
)
# ax.get_legend().remove()
fig.savefig(plot_path + "infectious_B_current", bbox_inches="tight")
