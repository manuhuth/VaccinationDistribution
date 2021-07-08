#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "p.h"
#include "k.h"

namespace amici {
namespace model_vaccination_piecewise {

void sigmay_vaccination_piecewise(realtype *sigmay, const realtype t, const realtype *p, const realtype *k){
    sigmay[0] = 1.0;
    sigmay[1] = 1.0;
    sigmay[2] = 1.0;
    sigmay[3] = 1.0;
    sigmay[4] = 1.0;
    sigmay[5] = 1.0;
    sigmay[6] = 1.0;
    sigmay[7] = 1.0;
    sigmay[8] = 1.0;
    sigmay[9] = 1.0;
    sigmay[10] = 1.0;
    sigmay[11] = 1.0;
}

} // namespace amici
} // namespace model_vaccination_piecewise