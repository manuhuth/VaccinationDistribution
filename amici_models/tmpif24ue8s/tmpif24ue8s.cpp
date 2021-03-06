#include "tmpif24ue8s.h"
#include <array>

namespace amici {

namespace model_tmpif24ue8s {

std::array<const char*, 22> parameterNames = {
    "yy_countryA_vac1_0",
"yy_countryA_vac1_1",
"yy_countryA_vac1_2",
"yy_countryA_vac1_3",
"yy_countryA_vac1_4",
"yy_countryA_vac1_5",
"yy_countryA_vac1_6",
"yy_countryA_vac1_7",
"yy_countryA_vac1_8",
"yy_countryA_vac1_9",
"yy_countryA_vac1_10",
"yy_countryA_vac2_0",
"yy_countryA_vac2_1",
"yy_countryA_vac2_2",
"yy_countryA_vac2_3",
"yy_countryA_vac2_4",
"yy_countryA_vac2_5",
"yy_countryA_vac2_6",
"yy_countryA_vac2_7",
"yy_countryA_vac2_8",
"yy_countryA_vac2_9",
"yy_countryA_vac2_10",
};

std::array<const char*, 93> fixedParameterNames = {
    "lambda1",
"prob_deceasing",
"gamma",
"beta",
"omega_vac1_virW",
"delta_vac1_virW",
"omega_vac2_virW",
"delta_vac2_virW",
"omega_vac1_virM",
"delta_vac1_virM",
"omega_vac2_virM",
"delta_vac2_virM",
"eta_virW",
"eta_virM",
"susceptible_countryA_vac0_t0",
"susceptible_countryA_vac1_t0",
"susceptible_countryA_vac2_t0",
"susceptible_countryB_vac0_t0",
"susceptible_countryB_vac1_t0",
"susceptible_countryB_vac2_t0",
"infectious_countryA_vac0_virW_t0",
"infectious_countryA_vac0_virM_t0",
"infectious_countryA_vac1_virW_t0",
"infectious_countryA_vac1_virM_t0",
"infectious_countryA_vac2_virW_t0",
"infectious_countryA_vac2_virM_t0",
"infectious_countryB_vac0_virW_t0",
"infectious_countryB_vac0_virM_t0",
"infectious_countryB_vac1_virW_t0",
"infectious_countryB_vac1_virM_t0",
"infectious_countryB_vac2_virW_t0",
"infectious_countryB_vac2_virM_t0",
"recovered_countryA_vac0_virW_t0",
"recovered_countryA_vac0_virM_t0",
"recovered_countryA_vac1_virW_t0",
"recovered_countryA_vac1_virM_t0",
"recovered_countryA_vac2_virW_t0",
"recovered_countryA_vac2_virM_t0",
"recovered_countryB_vac0_virW_t0",
"recovered_countryB_vac0_virM_t0",
"recovered_countryB_vac1_virW_t0",
"recovered_countryB_vac1_virM_t0",
"recovered_countryB_vac2_virW_t0",
"recovered_countryB_vac2_virM_t0",
"dead_countryA_vac0_virW_t0",
"dead_countryA_vac0_virM_t0",
"dead_countryA_vac1_virW_t0",
"dead_countryA_vac1_virM_t0",
"dead_countryA_vac2_virW_t0",
"dead_countryA_vac2_virM_t0",
"dead_countryB_vac0_virW_t0",
"dead_countryB_vac0_virM_t0",
"dead_countryB_vac1_virW_t0",
"dead_countryB_vac1_virM_t0",
"dead_countryB_vac2_virW_t0",
"dead_countryB_vac2_virM_t0",
"vaccine_supply_par_vac1_xx0",
"vaccine_supply_par_vac1_xx1",
"vaccine_supply_par_vac1_xx2",
"vaccine_supply_par_vac1_xx3",
"vaccine_supply_par_vac1_xx4",
"vaccine_supply_par_vac1_xx5",
"vaccine_supply_par_vac1_xx6",
"vaccine_supply_par_vac1_xx7",
"vaccine_supply_par_vac1_xx8",
"vaccine_supply_par_vac1_xx9",
"vaccine_supply_par_vac1_xx10",
"vaccine_supply_par_vac2_xx0",
"vaccine_supply_par_vac2_xx1",
"vaccine_supply_par_vac2_xx2",
"vaccine_supply_par_vac2_xx3",
"vaccine_supply_par_vac2_xx4",
"vaccine_supply_par_vac2_xx5",
"vaccine_supply_par_vac2_xx6",
"vaccine_supply_par_vac2_xx7",
"vaccine_supply_par_vac2_xx8",
"vaccine_supply_par_vac2_xx9",
"vaccine_supply_par_vac2_xx10",
"distance_countryA_countryA",
"distance_countryA_countryB",
"distance_countryB_countryA",
"distance_countryB_countryB",
"xx0",
"xx1",
"xx2",
"xx3",
"xx4",
"xx5",
"xx6",
"xx7",
"xx8",
"xx9",
"xx10",
};

std::array<const char*, 42> stateNames = {
    "susceptible_countryA_vac0",
"susceptible_countryA_vac1",
"susceptible_countryA_vac2",
"susceptible_countryB_vac0",
"susceptible_countryB_vac1",
"susceptible_countryB_vac2",
"infectious_countryA_vac0_virW",
"infectious_countryA_vac0_virM",
"infectious_countryA_vac1_virW",
"infectious_countryA_vac1_virM",
"infectious_countryA_vac2_virW",
"infectious_countryA_vac2_virM",
"infectious_countryB_vac0_virW",
"infectious_countryB_vac0_virM",
"infectious_countryB_vac1_virW",
"infectious_countryB_vac1_virM",
"infectious_countryB_vac2_virW",
"infectious_countryB_vac2_virM",
"recovered_countryA_vac0_virW",
"recovered_countryA_vac0_virM",
"recovered_countryA_vac1_virW",
"recovered_countryA_vac1_virM",
"recovered_countryA_vac2_virW",
"recovered_countryA_vac2_virM",
"recovered_countryB_vac0_virW",
"recovered_countryB_vac0_virM",
"recovered_countryB_vac1_virW",
"recovered_countryB_vac1_virM",
"recovered_countryB_vac2_virW",
"recovered_countryB_vac2_virM",
"dead_countryA_vac0_virW",
"dead_countryA_vac0_virM",
"dead_countryA_vac1_virW",
"dead_countryA_vac1_virM",
"dead_countryA_vac2_virW",
"dead_countryA_vac2_virM",
"dead_countryB_vac0_virW",
"dead_countryB_vac0_virM",
"dead_countryB_vac1_virW",
"dead_countryB_vac1_virM",
"dead_countryB_vac2_virW",
"dead_countryB_vac2_virM",
};

std::array<const char*, 1> observableNames = {
    "",
};

std::array<const char*, 22> parameterIds = {
    "yy_countryA_vac1_0",
"yy_countryA_vac1_1",
"yy_countryA_vac1_2",
"yy_countryA_vac1_3",
"yy_countryA_vac1_4",
"yy_countryA_vac1_5",
"yy_countryA_vac1_6",
"yy_countryA_vac1_7",
"yy_countryA_vac1_8",
"yy_countryA_vac1_9",
"yy_countryA_vac1_10",
"yy_countryA_vac2_0",
"yy_countryA_vac2_1",
"yy_countryA_vac2_2",
"yy_countryA_vac2_3",
"yy_countryA_vac2_4",
"yy_countryA_vac2_5",
"yy_countryA_vac2_6",
"yy_countryA_vac2_7",
"yy_countryA_vac2_8",
"yy_countryA_vac2_9",
"yy_countryA_vac2_10",
};

std::array<const char*, 93> fixedParameterIds = {
    "lambda1",
"prob_deceasing",
"gamma",
"beta",
"omega_vac1_virW",
"delta_vac1_virW",
"omega_vac2_virW",
"delta_vac2_virW",
"omega_vac1_virM",
"delta_vac1_virM",
"omega_vac2_virM",
"delta_vac2_virM",
"eta_virW",
"eta_virM",
"susceptible_countryA_vac0_t0",
"susceptible_countryA_vac1_t0",
"susceptible_countryA_vac2_t0",
"susceptible_countryB_vac0_t0",
"susceptible_countryB_vac1_t0",
"susceptible_countryB_vac2_t0",
"infectious_countryA_vac0_virW_t0",
"infectious_countryA_vac0_virM_t0",
"infectious_countryA_vac1_virW_t0",
"infectious_countryA_vac1_virM_t0",
"infectious_countryA_vac2_virW_t0",
"infectious_countryA_vac2_virM_t0",
"infectious_countryB_vac0_virW_t0",
"infectious_countryB_vac0_virM_t0",
"infectious_countryB_vac1_virW_t0",
"infectious_countryB_vac1_virM_t0",
"infectious_countryB_vac2_virW_t0",
"infectious_countryB_vac2_virM_t0",
"recovered_countryA_vac0_virW_t0",
"recovered_countryA_vac0_virM_t0",
"recovered_countryA_vac1_virW_t0",
"recovered_countryA_vac1_virM_t0",
"recovered_countryA_vac2_virW_t0",
"recovered_countryA_vac2_virM_t0",
"recovered_countryB_vac0_virW_t0",
"recovered_countryB_vac0_virM_t0",
"recovered_countryB_vac1_virW_t0",
"recovered_countryB_vac1_virM_t0",
"recovered_countryB_vac2_virW_t0",
"recovered_countryB_vac2_virM_t0",
"dead_countryA_vac0_virW_t0",
"dead_countryA_vac0_virM_t0",
"dead_countryA_vac1_virW_t0",
"dead_countryA_vac1_virM_t0",
"dead_countryA_vac2_virW_t0",
"dead_countryA_vac2_virM_t0",
"dead_countryB_vac0_virW_t0",
"dead_countryB_vac0_virM_t0",
"dead_countryB_vac1_virW_t0",
"dead_countryB_vac1_virM_t0",
"dead_countryB_vac2_virW_t0",
"dead_countryB_vac2_virM_t0",
"vaccine_supply_par_vac1_xx0",
"vaccine_supply_par_vac1_xx1",
"vaccine_supply_par_vac1_xx2",
"vaccine_supply_par_vac1_xx3",
"vaccine_supply_par_vac1_xx4",
"vaccine_supply_par_vac1_xx5",
"vaccine_supply_par_vac1_xx6",
"vaccine_supply_par_vac1_xx7",
"vaccine_supply_par_vac1_xx8",
"vaccine_supply_par_vac1_xx9",
"vaccine_supply_par_vac1_xx10",
"vaccine_supply_par_vac2_xx0",
"vaccine_supply_par_vac2_xx1",
"vaccine_supply_par_vac2_xx2",
"vaccine_supply_par_vac2_xx3",
"vaccine_supply_par_vac2_xx4",
"vaccine_supply_par_vac2_xx5",
"vaccine_supply_par_vac2_xx6",
"vaccine_supply_par_vac2_xx7",
"vaccine_supply_par_vac2_xx8",
"vaccine_supply_par_vac2_xx9",
"vaccine_supply_par_vac2_xx10",
"distance_countryA_countryA",
"distance_countryA_countryB",
"distance_countryB_countryA",
"distance_countryB_countryB",
"xx0",
"xx1",
"xx2",
"xx3",
"xx4",
"xx5",
"xx6",
"xx7",
"xx8",
"xx9",
"xx10",
};

std::array<const char*, 42> stateIds = {
    "susceptible_countryA_vac0",
"susceptible_countryA_vac1",
"susceptible_countryA_vac2",
"susceptible_countryB_vac0",
"susceptible_countryB_vac1",
"susceptible_countryB_vac2",
"infectious_countryA_vac0_virW",
"infectious_countryA_vac0_virM",
"infectious_countryA_vac1_virW",
"infectious_countryA_vac1_virM",
"infectious_countryA_vac2_virW",
"infectious_countryA_vac2_virM",
"infectious_countryB_vac0_virW",
"infectious_countryB_vac0_virM",
"infectious_countryB_vac1_virW",
"infectious_countryB_vac1_virM",
"infectious_countryB_vac2_virW",
"infectious_countryB_vac2_virM",
"recovered_countryA_vac0_virW",
"recovered_countryA_vac0_virM",
"recovered_countryA_vac1_virW",
"recovered_countryA_vac1_virM",
"recovered_countryA_vac2_virW",
"recovered_countryA_vac2_virM",
"recovered_countryB_vac0_virW",
"recovered_countryB_vac0_virM",
"recovered_countryB_vac1_virW",
"recovered_countryB_vac1_virM",
"recovered_countryB_vac2_virW",
"recovered_countryB_vac2_virM",
"dead_countryA_vac0_virW",
"dead_countryA_vac0_virM",
"dead_countryA_vac1_virW",
"dead_countryA_vac1_virM",
"dead_countryA_vac2_virW",
"dead_countryA_vac2_virM",
"dead_countryB_vac0_virW",
"dead_countryB_vac0_virM",
"dead_countryB_vac1_virW",
"dead_countryB_vac1_virM",
"dead_countryB_vac2_virW",
"dead_countryB_vac2_virM",
};

std::array<const char*, 1> observableIds = {
    "observable_deaths",
};


} // namespace model_tmpif24ue8s

} // namespace amici
