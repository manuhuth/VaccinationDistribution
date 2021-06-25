#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "vaccination_partly_x.h"
#include "vaccination_partly_p.h"
#include "vaccination_partly_h.h"

namespace amici {
namespace model_vaccination_partly {

void root_vaccination_partly(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h){
    root[0] = amici_t - xx1;
    root[1] = -amici_t + xx1;
    root[2] = amici_t - xx2;
    root[3] = -amici_t + xx2;
    root[4] = amici_t - xx3;
    root[5] = -amici_t + xx3;
    root[6] = amici_t - xx4;
    root[7] = -amici_t + xx4;
    root[8] = amici_t - xx5;
    root[9] = -amici_t + xx5;
    root[10] = amici_t - xx6;
    root[11] = -amici_t + xx6;
    root[12] = amici_t - xx7;
    root[13] = -amici_t + xx7;
    root[14] = amici_t - xx8;
    root[15] = -amici_t + xx8;
    root[16] = amici_t - xx9;
    root[17] = -amici_t + xx9;
}

} // namespace model_vaccination_partly
} // namespace amici
