#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "basic_sir_x_rdata.h"

namespace amici {
namespace model_basic_sir {

void x_solver_basic_sir(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = susceptible;
    x_solver[1] = infectious;
    x_solver[2] = recovered;
}

} // namespace model_basic_sir
} // namespace amici