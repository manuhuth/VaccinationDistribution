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

void Jy_vaccination(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ysusceptible_countryA_vac0, 2)) + 0.5*std::pow(-mysusceptible_countryA_vac0 + ysusceptible_countryA_vac0, 2)/std::pow(sigma_ysusceptible_countryA_vac0, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ysusceptible_countryA_vac1, 2)) + 0.5*std::pow(-mysusceptible_countryA_vac1 + ysusceptible_countryA_vac1, 2)/std::pow(sigma_ysusceptible_countryA_vac1, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ysusceptible_countryA_vac2, 2)) + 0.5*std::pow(-mysusceptible_countryA_vac2 + ysusceptible_countryA_vac2, 2)/std::pow(sigma_ysusceptible_countryA_vac2, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yinfectious_countryA_vac0_virW, 2)) + 0.5*std::pow(-myinfectious_countryA_vac0_virW + yinfectious_countryA_vac0_virW, 2)/std::pow(sigma_yinfectious_countryA_vac0_virW, 2);
            break;
        case 4:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yinfectious_countryA_vac0_virM, 2)) + 0.5*std::pow(-myinfectious_countryA_vac0_virM + yinfectious_countryA_vac0_virM, 2)/std::pow(sigma_yinfectious_countryA_vac0_virM, 2);
            break;
        case 5:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yinfectious_countryA_vac1_virW, 2)) + 0.5*std::pow(-myinfectious_countryA_vac1_virW + yinfectious_countryA_vac1_virW, 2)/std::pow(sigma_yinfectious_countryA_vac1_virW, 2);
            break;
        case 6:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yinfectious_countryA_vac1_virM, 2)) + 0.5*std::pow(-myinfectious_countryA_vac1_virM + yinfectious_countryA_vac1_virM, 2)/std::pow(sigma_yinfectious_countryA_vac1_virM, 2);
            break;
        case 7:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yinfectious_countryA_vac2_virW, 2)) + 0.5*std::pow(-myinfectious_countryA_vac2_virW + yinfectious_countryA_vac2_virW, 2)/std::pow(sigma_yinfectious_countryA_vac2_virW, 2);
            break;
        case 8:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yinfectious_countryA_vac2_virM, 2)) + 0.5*std::pow(-myinfectious_countryA_vac2_virM + yinfectious_countryA_vac2_virM, 2)/std::pow(sigma_yinfectious_countryA_vac2_virM, 2);
            break;
        case 9:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yrecovered_countryA_vac0_virW, 2)) + 0.5*std::pow(-myrecovered_countryA_vac0_virW + yrecovered_countryA_vac0_virW, 2)/std::pow(sigma_yrecovered_countryA_vac0_virW, 2);
            break;
        case 10:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yrecovered_countryA_vac0_virM, 2)) + 0.5*std::pow(-myrecovered_countryA_vac0_virM + yrecovered_countryA_vac0_virM, 2)/std::pow(sigma_yrecovered_countryA_vac0_virM, 2);
            break;
        case 11:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yrecovered_countryA_vac1_virW, 2)) + 0.5*std::pow(-myrecovered_countryA_vac1_virW + yrecovered_countryA_vac1_virW, 2)/std::pow(sigma_yrecovered_countryA_vac1_virW, 2);
            break;
        case 12:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yrecovered_countryA_vac1_virM, 2)) + 0.5*std::pow(-myrecovered_countryA_vac1_virM + yrecovered_countryA_vac1_virM, 2)/std::pow(sigma_yrecovered_countryA_vac1_virM, 2);
            break;
        case 13:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yrecovered_countryA_vac2_virW, 2)) + 0.5*std::pow(-myrecovered_countryA_vac2_virW + yrecovered_countryA_vac2_virW, 2)/std::pow(sigma_yrecovered_countryA_vac2_virW, 2);
            break;
        case 14:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yrecovered_countryA_vac2_virM, 2)) + 0.5*std::pow(-myrecovered_countryA_vac2_virM + yrecovered_countryA_vac2_virM, 2)/std::pow(sigma_yrecovered_countryA_vac2_virM, 2);
            break;
        case 15:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ydead_countryA_vac0_virW, 2)) + 0.5*std::pow(-mydead_countryA_vac0_virW + ydead_countryA_vac0_virW, 2)/std::pow(sigma_ydead_countryA_vac0_virW, 2);
            break;
        case 16:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ydead_countryA_vac0_virM, 2)) + 0.5*std::pow(-mydead_countryA_vac0_virM + ydead_countryA_vac0_virM, 2)/std::pow(sigma_ydead_countryA_vac0_virM, 2);
            break;
        case 17:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ydead_countryA_vac1_virW, 2)) + 0.5*std::pow(-mydead_countryA_vac1_virW + ydead_countryA_vac1_virW, 2)/std::pow(sigma_ydead_countryA_vac1_virW, 2);
            break;
        case 18:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ydead_countryA_vac1_virM, 2)) + 0.5*std::pow(-mydead_countryA_vac1_virM + ydead_countryA_vac1_virM, 2)/std::pow(sigma_ydead_countryA_vac1_virM, 2);
            break;
        case 19:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ydead_countryA_vac2_virW, 2)) + 0.5*std::pow(-mydead_countryA_vac2_virW + ydead_countryA_vac2_virW, 2)/std::pow(sigma_ydead_countryA_vac2_virW, 2);
            break;
        case 20:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ydead_countryA_vac2_virM, 2)) + 0.5*std::pow(-mydead_countryA_vac2_virM + ydead_countryA_vac2_virM, 2)/std::pow(sigma_ydead_countryA_vac2_virM, 2);
            break;
        case 21:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ycountryA, 2)) + 0.5*std::pow(-mycountryA + ycountryA, 2)/std::pow(sigma_ycountryA, 2);
            break;
    }
}

} // namespace model_vaccination
} // namespace amici
