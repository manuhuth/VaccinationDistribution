#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "one_country_p.h"

namespace amici {
namespace model_one_country {

void x0_one_country(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 10000.0;
    x0[3] = 1.0;
    x0[4] = 1.0;
}

} // namespace model_one_country
} // namespace amici
