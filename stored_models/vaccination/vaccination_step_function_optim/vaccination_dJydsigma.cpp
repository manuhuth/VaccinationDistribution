#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "vaccination_p.h"
#include "vaccination_y.h"
#include "vaccination_sigmay.h"
#include "vaccination_my.h"

namespace amici {
namespace model_vaccination {

void dJydsigma_vaccination(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_observable_nu_countryA_vac1 - 1.0*std::pow(-mobservable_nu_countryA_vac1 + observable_nu_countryA_vac1, 2)/std::pow(sigma_observable_nu_countryA_vac1, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_observable_nu_countryB_vac1 - 1.0*std::pow(-mobservable_nu_countryB_vac1 + observable_nu_countryB_vac1, 2)/std::pow(sigma_observable_nu_countryB_vac1, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_observable_proportion_countryA_vac1 - 1.0*std::pow(-mobservable_proportion_countryA_vac1 + observable_proportion_countryA_vac1, 2)/std::pow(sigma_observable_proportion_countryA_vac1, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_observable_proportion_countryB_vac1 - 1.0*std::pow(-mobservable_proportion_countryB_vac1 + observable_proportion_countryB_vac1, 2)/std::pow(sigma_observable_proportion_countryB_vac1, 3);
            break;
    }
}

} // namespace model_vaccination
} // namespace amici
