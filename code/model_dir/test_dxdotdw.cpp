#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "test_x.h"
#include "test_p.h"
#include "test_w.h"
#include "test_dxdotdw.h"

namespace amici {
namespace model_test {

void dxdotdw_test(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dxdot0_dflux_r0 = -1.0;  // dxdotdw[0]
    dxdot1_dflux_r0 = 1.0;  // dxdotdw[1]
    dxdot1_dflux_r1 = -1.0;  // dxdotdw[2]
    dxdot2_dflux_r1 = 1.0;  // dxdotdw[3]
}

} // namespace model_test
} // namespace amici
