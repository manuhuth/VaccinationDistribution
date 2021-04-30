#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "basic_sir_x.h"
#include "basic_sir_p.h"
#include "basic_sir_w.h"
#include "basic_sir_dxdotdw.h"

namespace amici {
namespace model_basic_sir {

void dxdotdw_basic_sir(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dxdot0_dflux_r0 = -1.0;  // dxdotdw[0]
    dxdot1_dflux_r0 = 1.0;  // dxdotdw[1]
    dxdot1_dflux_r1 = -1.0;  // dxdotdw[2]
    dxdot2_dflux_r1 = 1.0;  // dxdotdw[3]
}

} // namespace model_basic_sir
} // namespace amici
