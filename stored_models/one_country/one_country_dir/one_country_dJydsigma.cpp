#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "one_country_p.h"
#include "one_country_y.h"
#include "one_country_sigmay.h"
#include "one_country_my.h"

namespace amici {
namespace model_one_country {

void dJydsigma_one_country(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_ysusceptible_countryA_vac0 - 1.0*std::pow(-mysusceptible_countryA_vac0 + ysusceptible_countryA_vac0, 2)/std::pow(sigma_ysusceptible_countryA_vac0, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_ysusceptible_countryA_vac1 - 1.0*std::pow(-mysusceptible_countryA_vac1 + ysusceptible_countryA_vac1, 2)/std::pow(sigma_ysusceptible_countryA_vac1, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_ysusceptible_countryA_vac2 - 1.0*std::pow(-mysusceptible_countryA_vac2 + ysusceptible_countryA_vac2, 2)/std::pow(sigma_ysusceptible_countryA_vac2, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_yinfectious_countryA_vac0_virW - 1.0*std::pow(-myinfectious_countryA_vac0_virW + yinfectious_countryA_vac0_virW, 2)/std::pow(sigma_yinfectious_countryA_vac0_virW, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_yinfectious_countryA_vac0_virM - 1.0*std::pow(-myinfectious_countryA_vac0_virM + yinfectious_countryA_vac0_virM, 2)/std::pow(sigma_yinfectious_countryA_vac0_virM, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_yinfectious_countryA_vac1_virW - 1.0*std::pow(-myinfectious_countryA_vac1_virW + yinfectious_countryA_vac1_virW, 2)/std::pow(sigma_yinfectious_countryA_vac1_virW, 3);
            break;
        case 6:
            dJydsigma[6] = 1.0/sigma_yinfectious_countryA_vac1_virM - 1.0*std::pow(-myinfectious_countryA_vac1_virM + yinfectious_countryA_vac1_virM, 2)/std::pow(sigma_yinfectious_countryA_vac1_virM, 3);
            break;
        case 7:
            dJydsigma[7] = 1.0/sigma_yinfectious_countryA_vac2_virW - 1.0*std::pow(-myinfectious_countryA_vac2_virW + yinfectious_countryA_vac2_virW, 2)/std::pow(sigma_yinfectious_countryA_vac2_virW, 3);
            break;
        case 8:
            dJydsigma[8] = 1.0/sigma_yinfectious_countryA_vac2_virM - 1.0*std::pow(-myinfectious_countryA_vac2_virM + yinfectious_countryA_vac2_virM, 2)/std::pow(sigma_yinfectious_countryA_vac2_virM, 3);
            break;
        case 9:
            dJydsigma[9] = 1.0/sigma_yrecovered_countryA_vac0_virW - 1.0*std::pow(-myrecovered_countryA_vac0_virW + yrecovered_countryA_vac0_virW, 2)/std::pow(sigma_yrecovered_countryA_vac0_virW, 3);
            break;
        case 10:
            dJydsigma[10] = 1.0/sigma_yrecovered_countryA_vac0_virM - 1.0*std::pow(-myrecovered_countryA_vac0_virM + yrecovered_countryA_vac0_virM, 2)/std::pow(sigma_yrecovered_countryA_vac0_virM, 3);
            break;
        case 11:
            dJydsigma[11] = 1.0/sigma_yrecovered_countryA_vac1_virW - 1.0*std::pow(-myrecovered_countryA_vac1_virW + yrecovered_countryA_vac1_virW, 2)/std::pow(sigma_yrecovered_countryA_vac1_virW, 3);
            break;
        case 12:
            dJydsigma[12] = 1.0/sigma_yrecovered_countryA_vac1_virM - 1.0*std::pow(-myrecovered_countryA_vac1_virM + yrecovered_countryA_vac1_virM, 2)/std::pow(sigma_yrecovered_countryA_vac1_virM, 3);
            break;
        case 13:
            dJydsigma[13] = 1.0/sigma_yrecovered_countryA_vac2_virW - 1.0*std::pow(-myrecovered_countryA_vac2_virW + yrecovered_countryA_vac2_virW, 2)/std::pow(sigma_yrecovered_countryA_vac2_virW, 3);
            break;
        case 14:
            dJydsigma[14] = 1.0/sigma_yrecovered_countryA_vac2_virM - 1.0*std::pow(-myrecovered_countryA_vac2_virM + yrecovered_countryA_vac2_virM, 2)/std::pow(sigma_yrecovered_countryA_vac2_virM, 3);
            break;
        case 15:
            dJydsigma[15] = 1.0/sigma_ydead_countryA_vac0_virW - 1.0*std::pow(-mydead_countryA_vac0_virW + ydead_countryA_vac0_virW, 2)/std::pow(sigma_ydead_countryA_vac0_virW, 3);
            break;
        case 16:
            dJydsigma[16] = 1.0/sigma_ydead_countryA_vac0_virM - 1.0*std::pow(-mydead_countryA_vac0_virM + ydead_countryA_vac0_virM, 2)/std::pow(sigma_ydead_countryA_vac0_virM, 3);
            break;
        case 17:
            dJydsigma[17] = 1.0/sigma_ydead_countryA_vac1_virW - 1.0*std::pow(-mydead_countryA_vac1_virW + ydead_countryA_vac1_virW, 2)/std::pow(sigma_ydead_countryA_vac1_virW, 3);
            break;
        case 18:
            dJydsigma[18] = 1.0/sigma_ydead_countryA_vac1_virM - 1.0*std::pow(-mydead_countryA_vac1_virM + ydead_countryA_vac1_virM, 2)/std::pow(sigma_ydead_countryA_vac1_virM, 3);
            break;
        case 19:
            dJydsigma[19] = 1.0/sigma_ydead_countryA_vac2_virW - 1.0*std::pow(-mydead_countryA_vac2_virW + ydead_countryA_vac2_virW, 2)/std::pow(sigma_ydead_countryA_vac2_virW, 3);
            break;
        case 20:
            dJydsigma[20] = 1.0/sigma_ydead_countryA_vac2_virM - 1.0*std::pow(-mydead_countryA_vac2_virM + ydead_countryA_vac2_virM, 2)/std::pow(sigma_ydead_countryA_vac2_virM, 3);
            break;
        case 21:
            dJydsigma[21] = 1.0/sigma_ycountryA - 1.0*std::pow(-mycountryA + ycountryA, 2)/std::pow(sigma_ycountryA, 3);
            break;
    }
}

} // namespace model_one_country
} // namespace amici
