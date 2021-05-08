#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "vaccination_p.h"
#include "vaccination_sigmay.h"

namespace amici {
namespace model_vaccination {

void sigmay_vaccination(realtype *sigmay, const realtype t, const realtype *p, const realtype *k){
    sigma_observable_nu_countryA_vac1 = 1.0;  // sigmay[0]
    sigma_observable_nu_countryB_vac1 = 1.0;  // sigmay[1]
    sigma_observable_nu_countryA_vac2 = 1.0;  // sigmay[2]
    sigma_observable_nu_countryB_vac2 = 1.0;  // sigmay[3]
}

} // namespace model_vaccination
} // namespace amici
