#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "vaccination_x.h"
#include "vaccination_p.h"
#include "vaccination_w.h"
#include "vaccination_dwdx.h"

namespace amici {
namespace model_vaccination {

void dydx_vaccination(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = 1;
    dydx[23] = 1;
    dydx[46] = 1;
    dydx[69] = 1;
    dydx[92] = 1;
    dydx[115] = 1;
    dydx[138] = 1;
    dydx[161] = 1;
    dydx[184] = 1;
    dydx[207] = 1;
    dydx[230] = 1;
    dydx[253] = 1;
    dydx[276] = 1;
    dydx[299] = 1;
    dydx[322] = 1;
    dydx[345] = 1;
    dydx[368] = 1;
    dydx[391] = 1;
    dydx[414] = 1;
    dydx[437] = 1;
    dydx[460] = 1;
}

} // namespace model_vaccination
} // namespace amici