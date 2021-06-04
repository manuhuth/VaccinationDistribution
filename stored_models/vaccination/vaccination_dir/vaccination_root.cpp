#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "vaccination_x.h"
#include "vaccination_p.h"
#include "vaccination_h.h"

namespace amici {
namespace model_vaccination {

void root_vaccination(realtype *root, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h){
    root[0] = amici_t - 80.0/19.0;
    root[1] = 80.0/19.0 - amici_t;
    root[2] = amici_t - 160.0/19.0;
    root[3] = 160.0/19.0 - amici_t;
    root[4] = amici_t - 240.0/19.0;
    root[5] = 240.0/19.0 - amici_t;
    root[6] = amici_t - 320.0/19.0;
    root[7] = 320.0/19.0 - amici_t;
    root[8] = amici_t - 400.0/19.0;
    root[9] = 400.0/19.0 - amici_t;
    root[10] = amici_t - 480.0/19.0;
    root[11] = 480.0/19.0 - amici_t;
    root[12] = amici_t - 560.0/19.0;
    root[13] = 560.0/19.0 - amici_t;
    root[14] = amici_t - 640.0/19.0;
    root[15] = 640.0/19.0 - amici_t;
    root[16] = amici_t - 720.0/19.0;
    root[17] = 720.0/19.0 - amici_t;
    root[18] = amici_t - 800.0/19.0;
    root[19] = 800.0/19.0 - amici_t;
    root[20] = amici_t - 880.0/19.0;
    root[21] = 880.0/19.0 - amici_t;
    root[22] = amici_t - 960.0/19.0;
    root[23] = 960.0/19.0 - amici_t;
    root[24] = amici_t - 1040.0/19.0;
    root[25] = 1040.0/19.0 - amici_t;
    root[26] = amici_t - 1120.0/19.0;
    root[27] = 1120.0/19.0 - amici_t;
    root[28] = amici_t - 1200.0/19.0;
    root[29] = 1200.0/19.0 - amici_t;
    root[30] = amici_t - 1280.0/19.0;
    root[31] = 1280.0/19.0 - amici_t;
    root[32] = amici_t - 1360.0/19.0;
    root[33] = 1360.0/19.0 - amici_t;
    root[34] = amici_t - 1440.0/19.0;
    root[35] = 1440.0/19.0 - amici_t;
}

} // namespace model_vaccination
} // namespace amici
