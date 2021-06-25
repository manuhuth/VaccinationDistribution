#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "vaccination_partly_p.h"
#include "vaccination_partly_y.h"
#include "vaccination_partly_sigmay.h"
#include "vaccination_partly_my.h"
#include "vaccination_partly_dJydy.h"

namespace amici {
namespace model_vaccination_partly {

void dJydy_vaccination_partly(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mysusceptible_countryA_vac0 + 1.0*ysusceptible_countryA_vac0)/std::pow(sigma_ysusceptible_countryA_vac0, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*mysusceptible_countryA_vac1 + 1.0*ysusceptible_countryA_vac1)/std::pow(sigma_ysusceptible_countryA_vac1, 2);
            break;
        case 2:
            dJydy[0] = (-1.0*mysusceptible_countryA_vac2 + 1.0*ysusceptible_countryA_vac2)/std::pow(sigma_ysusceptible_countryA_vac2, 2);
            break;
        case 3:
            dJydy[0] = (-1.0*mysusceptible_countryB_vac0 + 1.0*ysusceptible_countryB_vac0)/std::pow(sigma_ysusceptible_countryB_vac0, 2);
            break;
        case 4:
            dJydy[0] = (-1.0*mysusceptible_countryB_vac1 + 1.0*ysusceptible_countryB_vac1)/std::pow(sigma_ysusceptible_countryB_vac1, 2);
            break;
        case 5:
            dJydy[0] = (-1.0*mysusceptible_countryB_vac2 + 1.0*ysusceptible_countryB_vac2)/std::pow(sigma_ysusceptible_countryB_vac2, 2);
            break;
        case 6:
            dJydy[0] = (-1.0*myinfectious_countryA_vac0_virW + 1.0*yinfectious_countryA_vac0_virW)/std::pow(sigma_yinfectious_countryA_vac0_virW, 2);
            break;
        case 7:
            dJydy[0] = (-1.0*myinfectious_countryA_vac0_virM + 1.0*yinfectious_countryA_vac0_virM)/std::pow(sigma_yinfectious_countryA_vac0_virM, 2);
            break;
        case 8:
            dJydy[0] = (-1.0*myinfectious_countryA_vac1_virW + 1.0*yinfectious_countryA_vac1_virW)/std::pow(sigma_yinfectious_countryA_vac1_virW, 2);
            break;
        case 9:
            dJydy[0] = (-1.0*myinfectious_countryA_vac1_virM + 1.0*yinfectious_countryA_vac1_virM)/std::pow(sigma_yinfectious_countryA_vac1_virM, 2);
            break;
        case 10:
            dJydy[0] = (-1.0*myinfectious_countryA_vac2_virW + 1.0*yinfectious_countryA_vac2_virW)/std::pow(sigma_yinfectious_countryA_vac2_virW, 2);
            break;
        case 11:
            dJydy[0] = (-1.0*myinfectious_countryA_vac2_virM + 1.0*yinfectious_countryA_vac2_virM)/std::pow(sigma_yinfectious_countryA_vac2_virM, 2);
            break;
        case 12:
            dJydy[0] = (-1.0*myinfectious_countryB_vac0_virW + 1.0*yinfectious_countryB_vac0_virW)/std::pow(sigma_yinfectious_countryB_vac0_virW, 2);
            break;
        case 13:
            dJydy[0] = (-1.0*myinfectious_countryB_vac0_virM + 1.0*yinfectious_countryB_vac0_virM)/std::pow(sigma_yinfectious_countryB_vac0_virM, 2);
            break;
        case 14:
            dJydy[0] = (-1.0*myinfectious_countryB_vac1_virW + 1.0*yinfectious_countryB_vac1_virW)/std::pow(sigma_yinfectious_countryB_vac1_virW, 2);
            break;
        case 15:
            dJydy[0] = (-1.0*myinfectious_countryB_vac1_virM + 1.0*yinfectious_countryB_vac1_virM)/std::pow(sigma_yinfectious_countryB_vac1_virM, 2);
            break;
        case 16:
            dJydy[0] = (-1.0*myinfectious_countryB_vac2_virW + 1.0*yinfectious_countryB_vac2_virW)/std::pow(sigma_yinfectious_countryB_vac2_virW, 2);
            break;
        case 17:
            dJydy[0] = (-1.0*myinfectious_countryB_vac2_virM + 1.0*yinfectious_countryB_vac2_virM)/std::pow(sigma_yinfectious_countryB_vac2_virM, 2);
            break;
        case 18:
            dJydy[0] = (-1.0*myrecovered_countryA_vac0_virW + 1.0*yrecovered_countryA_vac0_virW)/std::pow(sigma_yrecovered_countryA_vac0_virW, 2);
            break;
        case 19:
            dJydy[0] = (-1.0*myrecovered_countryA_vac0_virM + 1.0*yrecovered_countryA_vac0_virM)/std::pow(sigma_yrecovered_countryA_vac0_virM, 2);
            break;
        case 20:
            dJydy[0] = (-1.0*myrecovered_countryA_vac1_virW + 1.0*yrecovered_countryA_vac1_virW)/std::pow(sigma_yrecovered_countryA_vac1_virW, 2);
            break;
        case 21:
            dJydy[0] = (-1.0*myrecovered_countryA_vac1_virM + 1.0*yrecovered_countryA_vac1_virM)/std::pow(sigma_yrecovered_countryA_vac1_virM, 2);
            break;
        case 22:
            dJydy[0] = (-1.0*myrecovered_countryA_vac2_virW + 1.0*yrecovered_countryA_vac2_virW)/std::pow(sigma_yrecovered_countryA_vac2_virW, 2);
            break;
        case 23:
            dJydy[0] = (-1.0*myrecovered_countryA_vac2_virM + 1.0*yrecovered_countryA_vac2_virM)/std::pow(sigma_yrecovered_countryA_vac2_virM, 2);
            break;
        case 24:
            dJydy[0] = (-1.0*myrecovered_countryB_vac0_virW + 1.0*yrecovered_countryB_vac0_virW)/std::pow(sigma_yrecovered_countryB_vac0_virW, 2);
            break;
        case 25:
            dJydy[0] = (-1.0*myrecovered_countryB_vac0_virM + 1.0*yrecovered_countryB_vac0_virM)/std::pow(sigma_yrecovered_countryB_vac0_virM, 2);
            break;
        case 26:
            dJydy[0] = (-1.0*myrecovered_countryB_vac1_virW + 1.0*yrecovered_countryB_vac1_virW)/std::pow(sigma_yrecovered_countryB_vac1_virW, 2);
            break;
        case 27:
            dJydy[0] = (-1.0*myrecovered_countryB_vac1_virM + 1.0*yrecovered_countryB_vac1_virM)/std::pow(sigma_yrecovered_countryB_vac1_virM, 2);
            break;
        case 28:
            dJydy[0] = (-1.0*myrecovered_countryB_vac2_virW + 1.0*yrecovered_countryB_vac2_virW)/std::pow(sigma_yrecovered_countryB_vac2_virW, 2);
            break;
        case 29:
            dJydy[0] = (-1.0*myrecovered_countryB_vac2_virM + 1.0*yrecovered_countryB_vac2_virM)/std::pow(sigma_yrecovered_countryB_vac2_virM, 2);
            break;
        case 30:
            dJydy[0] = (-1.0*mydead_countryA_vac0_virW + 1.0*ydead_countryA_vac0_virW)/std::pow(sigma_ydead_countryA_vac0_virW, 2);
            break;
        case 31:
            dJydy[0] = (-1.0*mydead_countryA_vac0_virM + 1.0*ydead_countryA_vac0_virM)/std::pow(sigma_ydead_countryA_vac0_virM, 2);
            break;
        case 32:
            dJydy[0] = (-1.0*mydead_countryA_vac1_virW + 1.0*ydead_countryA_vac1_virW)/std::pow(sigma_ydead_countryA_vac1_virW, 2);
            break;
        case 33:
            dJydy[0] = (-1.0*mydead_countryA_vac1_virM + 1.0*ydead_countryA_vac1_virM)/std::pow(sigma_ydead_countryA_vac1_virM, 2);
            break;
        case 34:
            dJydy[0] = (-1.0*mydead_countryA_vac2_virW + 1.0*ydead_countryA_vac2_virW)/std::pow(sigma_ydead_countryA_vac2_virW, 2);
            break;
        case 35:
            dJydy[0] = (-1.0*mydead_countryA_vac2_virM + 1.0*ydead_countryA_vac2_virM)/std::pow(sigma_ydead_countryA_vac2_virM, 2);
            break;
        case 36:
            dJydy[0] = (-1.0*mydead_countryB_vac0_virW + 1.0*ydead_countryB_vac0_virW)/std::pow(sigma_ydead_countryB_vac0_virW, 2);
            break;
        case 37:
            dJydy[0] = (-1.0*mydead_countryB_vac0_virM + 1.0*ydead_countryB_vac0_virM)/std::pow(sigma_ydead_countryB_vac0_virM, 2);
            break;
        case 38:
            dJydy[0] = (-1.0*mydead_countryB_vac1_virW + 1.0*ydead_countryB_vac1_virW)/std::pow(sigma_ydead_countryB_vac1_virW, 2);
            break;
        case 39:
            dJydy[0] = (-1.0*mydead_countryB_vac1_virM + 1.0*ydead_countryB_vac1_virM)/std::pow(sigma_ydead_countryB_vac1_virM, 2);
            break;
        case 40:
            dJydy[0] = (-1.0*mydead_countryB_vac2_virW + 1.0*ydead_countryB_vac2_virW)/std::pow(sigma_ydead_countryB_vac2_virW, 2);
            break;
        case 41:
            dJydy[0] = (-1.0*mydead_countryB_vac2_virM + 1.0*ydead_countryB_vac2_virM)/std::pow(sigma_ydead_countryB_vac2_virM, 2);
            break;
        case 42:
            dJydy[0] = (-1.0*myt + 1.0*yt)/std::pow(sigma_yt, 2);
            break;
        case 43:
            dJydy[0] = (-1.0*myproportion_countryA_vac1 + 1.0*yproportion_countryA_vac1)/std::pow(sigma_yproportion_countryA_vac1, 2);
            break;
        case 44:
            dJydy[0] = (-1.0*myproportion_countryA_vac2 + 1.0*yproportion_countryA_vac2)/std::pow(sigma_yproportion_countryA_vac2, 2);
            break;
        case 45:
            dJydy[0] = (-1.0*myproportion_countryB_vac1 + 1.0*yproportion_countryB_vac1)/std::pow(sigma_yproportion_countryB_vac1, 2);
            break;
        case 46:
            dJydy[0] = (-1.0*myproportion_countryB_vac2 + 1.0*yproportion_countryB_vac2)/std::pow(sigma_yproportion_countryB_vac2, 2);
            break;
        case 47:
            dJydy[0] = (-1.0*mynu_countryA_vac1 + 1.0*ynu_countryA_vac1)/std::pow(sigma_ynu_countryA_vac1, 2);
            break;
        case 48:
            dJydy[0] = (-1.0*mynu_countryB_vac1 + 1.0*ynu_countryB_vac1)/std::pow(sigma_ynu_countryB_vac1, 2);
            break;
        case 49:
            dJydy[0] = (-1.0*mynu_countryA_vac2 + 1.0*ynu_countryA_vac2)/std::pow(sigma_ynu_countryA_vac2, 2);
            break;
        case 50:
            dJydy[0] = (-1.0*mynu_countryB_vac2 + 1.0*ynu_countryB_vac2)/std::pow(sigma_ynu_countryB_vac2, 2);
            break;
        case 51:
            dJydy[0] = (-1.0*myspline_countryA_vac1 + 1.0*yspline_countryA_vac1)/std::pow(sigma_yspline_countryA_vac1, 2);
            break;
        case 52:
            dJydy[0] = (-1.0*myspline_countryA_vac2 + 1.0*yspline_countryA_vac2)/std::pow(sigma_yspline_countryA_vac2, 2);
            break;
        case 53:
            dJydy[0] = (-1.0*mycountryA + 1.0*ycountryA)/std::pow(sigma_ycountryA, 2);
            break;
        case 54:
            dJydy[0] = (-1.0*mycountryB + 1.0*ycountryB)/std::pow(sigma_ycountryB, 2);
            break;
    }
}

} // namespace model_vaccination_partly
} // namespace amici
