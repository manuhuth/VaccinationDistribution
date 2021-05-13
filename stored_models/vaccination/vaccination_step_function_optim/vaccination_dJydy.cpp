#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "vaccination_p.h"
#include "vaccination_y.h"
#include "vaccination_sigmay.h"
#include "vaccination_my.h"
#include "vaccination_dJydy.h"

namespace amici {
namespace model_vaccination {

void dJydy_vaccination(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mobservable_nu_countryA_vac1 + 1.0*observable_nu_countryA_vac1)/std::pow(sigma_observable_nu_countryA_vac1, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*mobservable_nu_countryB_vac1 + 1.0*observable_nu_countryB_vac1)/std::pow(sigma_observable_nu_countryB_vac1, 2);
            break;
        case 2:
            dJydy[0] = (-1.0*mobservable_proportion_countryA_vac1 + 1.0*observable_proportion_countryA_vac1)/std::pow(sigma_observable_proportion_countryA_vac1, 2);
            break;
        case 3:
            dJydy[0] = (-1.0*mobservable_proportion_countryB_vac1 + 1.0*observable_proportion_countryB_vac1)/std::pow(sigma_observable_proportion_countryB_vac1, 2);
            break;
    }
}

} // namespace model_vaccination
} // namespace amici
