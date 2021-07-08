#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "p.h"
#include "k.h"

namespace amici {
namespace model_vaccination {

void dspline_valuesdp_vaccination(realtype *dspline_valuesdp, const realtype *p, const realtype *k){
    dspline_valuesdp[91] = 1;
    dspline_valuesdp[205] = 1;
    dspline_valuesdp[319] = 1;
    dspline_valuesdp[433] = 1;
    dspline_valuesdp[547] = 1;
    dspline_valuesdp[661] = 1;
    dspline_valuesdp[775] = 1;
    dspline_valuesdp[889] = 1;
    dspline_valuesdp[1003] = 1;
    dspline_valuesdp[1117] = 1;
    dspline_valuesdp[1231] = 1;
    dspline_valuesdp[1345] = 1;
    dspline_valuesdp[1459] = 1;
    dspline_valuesdp[1573] = 1;
    dspline_valuesdp[1687] = 1;
    dspline_valuesdp[1801] = 1;
    dspline_valuesdp[1915] = 1;
    dspline_valuesdp[2029] = 1;
    dspline_valuesdp[2143] = 1;
    dspline_valuesdp[2257] = 1;
    dspline_valuesdp[2371] = 1;
    dspline_valuesdp[2485] = 1;
}

} // namespace amici
} // namespace model_vaccination