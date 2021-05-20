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
    root[0] = amici_t;
    root[1] = -amici_t;
    root[2] = amici_t - 3;
    root[3] = 3 - amici_t;
    root[4] = amici_t - 9;
    root[5] = 9 - amici_t;
    root[6] = amici_t - 6;
    root[7] = 6 - amici_t;
    root[8] = amici_t - 12;
    root[9] = 12 - amici_t;
}

} // namespace model_vaccination
} // namespace amici
