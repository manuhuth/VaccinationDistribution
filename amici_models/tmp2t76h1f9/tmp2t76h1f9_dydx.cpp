#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "tmp2t76h1f9_x.h"
#include "tmp2t76h1f9_p.h"
#include "tmp2t76h1f9_k.h"
#include "tmp2t76h1f9_h.h"
#include "tmp2t76h1f9_w.h"
#include "tmp2t76h1f9_dwdx.h"

namespace amici {
namespace model_tmp2t76h1f9 {

void dydx_tmp2t76h1f9(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[30] = 1;
    dydx[31] = 1;
    dydx[32] = 1;
    dydx[33] = 1;
    dydx[34] = 1;
    dydx[35] = 1;
    dydx[36] = 1;
    dydx[37] = 1;
    dydx[38] = 1;
    dydx[39] = 1;
    dydx[40] = 1;
    dydx[41] = 1;
}

} // namespace model_tmp2t76h1f9
} // namespace amici
