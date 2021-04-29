#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "test_x.h"
#include "test_p.h"
#include "test_w.h"
#include "test_dwdx.h"

namespace amici {
namespace model_test {

void dydx_test(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[5] = 1;
    dydx[10] = 1;
}

} // namespace model_test
} // namespace amici
