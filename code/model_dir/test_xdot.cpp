#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "test_x.h"
#include "test_p.h"
#include "test_w.h"
#include "test_xdot.h"

namespace amici {
namespace model_test {

void xdot_test(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot0 = -1.0*flux_r0;  // xdot[0]
    xdot1 = 1.0*flux_r0 - 1.0*flux_r1;  // xdot[1]
    xdot2 = 1.0*flux_r1;  // xdot[2]
}

} // namespace model_test
} // namespace amici
