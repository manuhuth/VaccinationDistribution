#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "basic_sir_p.h"
#include "basic_sir_sigmay.h"

namespace amici {
namespace model_basic_sir {

void sigmay_basic_sir(realtype *sigmay, const realtype t, const realtype *p, const realtype *k){
    sigma_ysusceptible = 1.0;  // sigmay[0]
    sigma_yinfectious = 1.0;  // sigmay[1]
    sigma_yrecovered = 1.0;  // sigmay[2]
    sigma_yA = 1.0;  // sigmay[3]
}

} // namespace model_basic_sir
} // namespace amici
