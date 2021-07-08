#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "p.h"
#include "k.h"
#include "w.h"
#include "x.h"
#include "dwdx.h"

namespace amici {
namespace model_tmpif24ue8s {

void dydx_tmpif24ue8s(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
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

} // namespace amici
} // namespace model_tmpif24ue8s