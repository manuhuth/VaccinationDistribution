#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "one_country_p.h"
#include "one_country_sigmay.h"

namespace amici {
namespace model_one_country {

void sigmay_one_country(realtype *sigmay, const realtype t, const realtype *p, const realtype *k){
    sigma_ysusceptible_countryA_vac0 = 1.0;  // sigmay[0]
    sigma_ysusceptible_countryA_vac1 = 1.0;  // sigmay[1]
    sigma_ysusceptible_countryA_vac2 = 1.0;  // sigmay[2]
    sigma_yinfectious_countryA_vac0_virW = 1.0;  // sigmay[3]
    sigma_yinfectious_countryA_vac0_virM = 1.0;  // sigmay[4]
    sigma_yinfectious_countryA_vac1_virW = 1.0;  // sigmay[5]
    sigma_yinfectious_countryA_vac1_virM = 1.0;  // sigmay[6]
    sigma_yinfectious_countryA_vac2_virW = 1.0;  // sigmay[7]
    sigma_yinfectious_countryA_vac2_virM = 1.0;  // sigmay[8]
    sigma_yrecovered_countryA_vac0_virW = 1.0;  // sigmay[9]
    sigma_yrecovered_countryA_vac0_virM = 1.0;  // sigmay[10]
    sigma_yrecovered_countryA_vac1_virW = 1.0;  // sigmay[11]
    sigma_yrecovered_countryA_vac1_virM = 1.0;  // sigmay[12]
    sigma_yrecovered_countryA_vac2_virW = 1.0;  // sigmay[13]
    sigma_yrecovered_countryA_vac2_virM = 1.0;  // sigmay[14]
    sigma_ydead_countryA_vac0_virW = 1.0;  // sigmay[15]
    sigma_ydead_countryA_vac0_virM = 1.0;  // sigmay[16]
    sigma_ydead_countryA_vac1_virW = 1.0;  // sigmay[17]
    sigma_ydead_countryA_vac1_virM = 1.0;  // sigmay[18]
    sigma_ydead_countryA_vac2_virW = 1.0;  // sigmay[19]
    sigma_ydead_countryA_vac2_virM = 1.0;  // sigmay[20]
    sigma_ycountryA = 1.0;  // sigmay[21]
}

} // namespace model_one_country
} // namespace amici
