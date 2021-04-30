#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "basic_sir_p.h"

namespace amici {
namespace model_basic_sir {

void x0_basic_sir(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 1000.0;
    x0[1] = 9.9999999999999995e-7;
}

} // namespace model_basic_sir
} // namespace amici
