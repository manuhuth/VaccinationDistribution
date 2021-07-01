#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "vaccination_x.h"
#include "vaccination_p.h"
#include "vaccination_h.h"
#include "vaccination_w.h"

namespace amici {
namespace model_vaccination {

void y_vaccination(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = amici_t;
    y[1] = number_vac1*proportion_countryA_vac1;
    y[2] = number_vac2*proportion_countryA_vac2;
    y[3] = number_vac1*proportion_countryB_vac1;
    y[4] = number_vac2*proportion_countryB_vac2;
}

} // namespace model_vaccination
} // namespace amici
