#include "test.h"
#include <array>

namespace amici {

namespace model_test {

std::array<const char*, 2> parameterNames = {
    "beta", // p[0]
"gamma", // p[1]
};

std::array<const char*, 0> fixedParameterNames = {
    
};

std::array<const char*, 3> stateNames = {
    "susceptible", // x_rdata[0]
"infectious", // x_rdata[1]
"recovered", // x_rdata[2]
};

std::array<const char*, 4> observableNames = {
    "susceptible", // y[0]
"infectious", // y[1]
"recovered", // y[2]
"A", // y[3]
};

std::array<const char*, 2> expressionNames = {
    "flux_r0", // w[0]
"flux_r1", // w[1]
};

std::array<const char*, 2> parameterIds = {
    "beta", // p[0]
"gamma", // p[1]
};

std::array<const char*, 0> fixedParameterIds = {
    
};

std::array<const char*, 3> stateIds = {
    "susceptible", // x_rdata[0]
"infectious", // x_rdata[1]
"recovered", // x_rdata[2]
};

std::array<const char*, 4> observableIds = {
    "ysusceptible", // y[0]
"yinfectious", // y[1]
"yrecovered", // y[2]
"yA", // y[3]
};

std::array<const char*, 2> expressionIds = {
    "flux_r0", // w[0]
"flux_r1", // w[1]
};

} // namespace model_test

} // namespace amici
