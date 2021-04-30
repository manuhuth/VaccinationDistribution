#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "basic_sir_x.h"
#include "basic_sir_p.h"
#include "basic_sir_w.h"

namespace amici {
namespace model_basic_sir {

void y_basic_sir(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = susceptible;
    y[1] = infectious;
    y[2] = recovered;
    y[3] = 1.0;
}

} // namespace model_basic_sir
} // namespace amici
