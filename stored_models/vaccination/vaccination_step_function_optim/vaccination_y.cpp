#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "vaccination_x.h"
#include "vaccination_p.h"
#include "vaccination_w.h"

namespace amici {
namespace model_vaccination {

void y_vaccination(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = nu_countryA_vac1;
    y[1] = nu_countryB_vac1;
    y[2] = proportion_countryA_vac1;
    y[3] = proportion_countryB_vac1;
}

} // namespace model_vaccination
} // namespace amici
