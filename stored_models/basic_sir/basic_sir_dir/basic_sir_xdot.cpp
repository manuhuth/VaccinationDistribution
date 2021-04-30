#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "basic_sir_x.h"
#include "basic_sir_p.h"
#include "basic_sir_w.h"
#include "basic_sir_xdot.h"

namespace amici {
namespace model_basic_sir {

void xdot_basic_sir(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot0 = -1.0*flux_r0;  // xdot[0]
    xdot1 = 1.0*flux_r0 - 1.0*flux_r1;  // xdot[1]
    xdot2 = 1.0*flux_r1;  // xdot[2]
}

} // namespace model_basic_sir
} // namespace amici
