#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "p.h"
#include "k.h"

namespace amici {
namespace model_vaccination {

void dspline_valuesdp_vaccination(realtype *dspline_valuesdp, const realtype *p, const realtype *k){
    dspline_valuesdp[94] = 1;
    dspline_valuesdp[211] = 1;
    dspline_valuesdp[328] = 1;
    dspline_valuesdp[445] = 1;
    dspline_valuesdp[562] = 1;
    dspline_valuesdp[679] = 1;
    dspline_valuesdp[796] = 1;
    dspline_valuesdp[913] = 1;
    dspline_valuesdp[1030] = 1;
    dspline_valuesdp[1147] = 1;
    dspline_valuesdp[1264] = 1;
    dspline_valuesdp[1381] = 1;
    dspline_valuesdp[1498] = 1;
    dspline_valuesdp[1615] = 1;
    dspline_valuesdp[1732] = 1;
    dspline_valuesdp[1849] = 1;
    dspline_valuesdp[1966] = 1;
    dspline_valuesdp[2083] = 1;
    dspline_valuesdp[2200] = 1;
    dspline_valuesdp[2317] = 1;
    dspline_valuesdp[2434] = 1;
    dspline_valuesdp[2551] = 1;
}

} // namespace amici
} // namespace model_vaccination