#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "p.h"
#include "k.h"

namespace amici {
namespace model_tmp_akih_pn {

void dspline_valuesdp_tmp_akih_pn(realtype *dspline_valuesdp, const realtype *p, const realtype *k){
    dspline_valuesdp[0] = 1;
    dspline_valuesdp[23] = 1;
    dspline_valuesdp[46] = 1;
    dspline_valuesdp[69] = 1;
    dspline_valuesdp[92] = 1;
    dspline_valuesdp[115] = 1;
    dspline_valuesdp[138] = 1;
    dspline_valuesdp[161] = 1;
    dspline_valuesdp[184] = 1;
    dspline_valuesdp[207] = 1;
    dspline_valuesdp[230] = 1;
    dspline_valuesdp[253] = 1;
    dspline_valuesdp[276] = 1;
    dspline_valuesdp[299] = 1;
    dspline_valuesdp[322] = 1;
    dspline_valuesdp[345] = 1;
    dspline_valuesdp[368] = 1;
    dspline_valuesdp[391] = 1;
    dspline_valuesdp[414] = 1;
    dspline_valuesdp[437] = 1;
    dspline_valuesdp[460] = 1;
    dspline_valuesdp[483] = 1;
}

} // namespace amici
} // namespace model_tmp_akih_pn