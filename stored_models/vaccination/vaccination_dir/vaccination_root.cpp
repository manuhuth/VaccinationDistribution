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
    root[0] = amici_t - 0.20000000000000001;
    root[1] = 0.20000000000000001 - amici_t;
}

} // namespace model_vaccination
} // namespace amici
