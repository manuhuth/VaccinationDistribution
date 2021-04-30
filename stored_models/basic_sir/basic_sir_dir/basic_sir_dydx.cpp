#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "basic_sir_x.h"
#include "basic_sir_p.h"
#include "basic_sir_w.h"
#include "basic_sir_dwdx.h"

namespace amici {
namespace model_basic_sir {

void dydx_basic_sir(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[5] = 1;
    dydx[10] = 1;
}

} // namespace model_basic_sir
} // namespace amici
