#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "vaccination_partly_p.h"
#include "vaccination_partly_y.h"
#include "vaccination_partly_sigmay.h"
#include "vaccination_partly_my.h"

namespace amici {
namespace model_vaccination_partly {

void dJydsigma_vaccination_partly(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
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
            dJydsigma[3] = 1.0/sigma_ysusceptible_countryB_vac0 - 1.0*std::pow(-mysusceptible_countryB_vac0 + ysusceptible_countryB_vac0, 2)/std::pow(sigma_ysusceptible_countryB_vac0, 3);
            break;
        case 4:
            dJydsigma[4] = 1.0/sigma_ysusceptible_countryB_vac1 - 1.0*std::pow(-mysusceptible_countryB_vac1 + ysusceptible_countryB_vac1, 2)/std::pow(sigma_ysusceptible_countryB_vac1, 3);
            break;
        case 5:
            dJydsigma[5] = 1.0/sigma_ysusceptible_countryB_vac2 - 1.0*std::pow(-mysusceptible_countryB_vac2 + ysusceptible_countryB_vac2, 2)/std::pow(sigma_ysusceptible_countryB_vac2, 3);
            break;
        case 6:
            dJydsigma[6] = 1.0/sigma_yinfectious_countryA_vac0_virW - 1.0*std::pow(-myinfectious_countryA_vac0_virW + yinfectious_countryA_vac0_virW, 2)/std::pow(sigma_yinfectious_countryA_vac0_virW, 3);
            break;
        case 7:
            dJydsigma[7] = 1.0/sigma_yinfectious_countryA_vac0_virM - 1.0*std::pow(-myinfectious_countryA_vac0_virM + yinfectious_countryA_vac0_virM, 2)/std::pow(sigma_yinfectious_countryA_vac0_virM, 3);
            break;
        case 8:
            dJydsigma[8] = 1.0/sigma_yinfectious_countryA_vac1_virW - 1.0*std::pow(-myinfectious_countryA_vac1_virW + yinfectious_countryA_vac1_virW, 2)/std::pow(sigma_yinfectious_countryA_vac1_virW, 3);
            break;
        case 9:
            dJydsigma[9] = 1.0/sigma_yinfectious_countryA_vac1_virM - 1.0*std::pow(-myinfectious_countryA_vac1_virM + yinfectious_countryA_vac1_virM, 2)/std::pow(sigma_yinfectious_countryA_vac1_virM, 3);
            break;
        case 10:
            dJydsigma[10] = 1.0/sigma_yinfectious_countryA_vac2_virW - 1.0*std::pow(-myinfectious_countryA_vac2_virW + yinfectious_countryA_vac2_virW, 2)/std::pow(sigma_yinfectious_countryA_vac2_virW, 3);
            break;
        case 11:
            dJydsigma[11] = 1.0/sigma_yinfectious_countryA_vac2_virM - 1.0*std::pow(-myinfectious_countryA_vac2_virM + yinfectious_countryA_vac2_virM, 2)/std::pow(sigma_yinfectious_countryA_vac2_virM, 3);
            break;
        case 12:
            dJydsigma[12] = 1.0/sigma_yinfectious_countryB_vac0_virW - 1.0*std::pow(-myinfectious_countryB_vac0_virW + yinfectious_countryB_vac0_virW, 2)/std::pow(sigma_yinfectious_countryB_vac0_virW, 3);
            break;
        case 13:
            dJydsigma[13] = 1.0/sigma_yinfectious_countryB_vac0_virM - 1.0*std::pow(-myinfectious_countryB_vac0_virM + yinfectious_countryB_vac0_virM, 2)/std::pow(sigma_yinfectious_countryB_vac0_virM, 3);
            break;
        case 14:
            dJydsigma[14] = 1.0/sigma_yinfectious_countryB_vac1_virW - 1.0*std::pow(-myinfectious_countryB_vac1_virW + yinfectious_countryB_vac1_virW, 2)/std::pow(sigma_yinfectious_countryB_vac1_virW, 3);
            break;
        case 15:
            dJydsigma[15] = 1.0/sigma_yinfectious_countryB_vac1_virM - 1.0*std::pow(-myinfectious_countryB_vac1_virM + yinfectious_countryB_vac1_virM, 2)/std::pow(sigma_yinfectious_countryB_vac1_virM, 3);
            break;
        case 16:
            dJydsigma[16] = 1.0/sigma_yinfectious_countryB_vac2_virW - 1.0*std::pow(-myinfectious_countryB_vac2_virW + yinfectious_countryB_vac2_virW, 2)/std::pow(sigma_yinfectious_countryB_vac2_virW, 3);
            break;
        case 17:
            dJydsigma[17] = 1.0/sigma_yinfectious_countryB_vac2_virM - 1.0*std::pow(-myinfectious_countryB_vac2_virM + yinfectious_countryB_vac2_virM, 2)/std::pow(sigma_yinfectious_countryB_vac2_virM, 3);
            break;
        case 18:
            dJydsigma[18] = 1.0/sigma_yrecovered_countryA_vac0_virW - 1.0*std::pow(-myrecovered_countryA_vac0_virW + yrecovered_countryA_vac0_virW, 2)/std::pow(sigma_yrecovered_countryA_vac0_virW, 3);
            break;
        case 19:
            dJydsigma[19] = 1.0/sigma_yrecovered_countryA_vac0_virM - 1.0*std::pow(-myrecovered_countryA_vac0_virM + yrecovered_countryA_vac0_virM, 2)/std::pow(sigma_yrecovered_countryA_vac0_virM, 3);
            break;
        case 20:
            dJydsigma[20] = 1.0/sigma_yrecovered_countryA_vac1_virW - 1.0*std::pow(-myrecovered_countryA_vac1_virW + yrecovered_countryA_vac1_virW, 2)/std::pow(sigma_yrecovered_countryA_vac1_virW, 3);
            break;
        case 21:
            dJydsigma[21] = 1.0/sigma_yrecovered_countryA_vac1_virM - 1.0*std::pow(-myrecovered_countryA_vac1_virM + yrecovered_countryA_vac1_virM, 2)/std::pow(sigma_yrecovered_countryA_vac1_virM, 3);
            break;
        case 22:
            dJydsigma[22] = 1.0/sigma_yrecovered_countryA_vac2_virW - 1.0*std::pow(-myrecovered_countryA_vac2_virW + yrecovered_countryA_vac2_virW, 2)/std::pow(sigma_yrecovered_countryA_vac2_virW, 3);
            break;
        case 23:
            dJydsigma[23] = 1.0/sigma_yrecovered_countryA_vac2_virM - 1.0*std::pow(-myrecovered_countryA_vac2_virM + yrecovered_countryA_vac2_virM, 2)/std::pow(sigma_yrecovered_countryA_vac2_virM, 3);
            break;
        case 24:
            dJydsigma[24] = 1.0/sigma_yrecovered_countryB_vac0_virW - 1.0*std::pow(-myrecovered_countryB_vac0_virW + yrecovered_countryB_vac0_virW, 2)/std::pow(sigma_yrecovered_countryB_vac0_virW, 3);
            break;
        case 25:
            dJydsigma[25] = 1.0/sigma_yrecovered_countryB_vac0_virM - 1.0*std::pow(-myrecovered_countryB_vac0_virM + yrecovered_countryB_vac0_virM, 2)/std::pow(sigma_yrecovered_countryB_vac0_virM, 3);
            break;
        case 26:
            dJydsigma[26] = 1.0/sigma_yrecovered_countryB_vac1_virW - 1.0*std::pow(-myrecovered_countryB_vac1_virW + yrecovered_countryB_vac1_virW, 2)/std::pow(sigma_yrecovered_countryB_vac1_virW, 3);
            break;
        case 27:
            dJydsigma[27] = 1.0/sigma_yrecovered_countryB_vac1_virM - 1.0*std::pow(-myrecovered_countryB_vac1_virM + yrecovered_countryB_vac1_virM, 2)/std::pow(sigma_yrecovered_countryB_vac1_virM, 3);
            break;
        case 28:
            dJydsigma[28] = 1.0/sigma_yrecovered_countryB_vac2_virW - 1.0*std::pow(-myrecovered_countryB_vac2_virW + yrecovered_countryB_vac2_virW, 2)/std::pow(sigma_yrecovered_countryB_vac2_virW, 3);
            break;
        case 29:
            dJydsigma[29] = 1.0/sigma_yrecovered_countryB_vac2_virM - 1.0*std::pow(-myrecovered_countryB_vac2_virM + yrecovered_countryB_vac2_virM, 2)/std::pow(sigma_yrecovered_countryB_vac2_virM, 3);
            break;
        case 30:
            dJydsigma[30] = 1.0/sigma_ydead_countryA_vac0_virW - 1.0*std::pow(-mydead_countryA_vac0_virW + ydead_countryA_vac0_virW, 2)/std::pow(sigma_ydead_countryA_vac0_virW, 3);
            break;
        case 31:
            dJydsigma[31] = 1.0/sigma_ydead_countryA_vac0_virM - 1.0*std::pow(-mydead_countryA_vac0_virM + ydead_countryA_vac0_virM, 2)/std::pow(sigma_ydead_countryA_vac0_virM, 3);
            break;
        case 32:
            dJydsigma[32] = 1.0/sigma_ydead_countryA_vac1_virW - 1.0*std::pow(-mydead_countryA_vac1_virW + ydead_countryA_vac1_virW, 2)/std::pow(sigma_ydead_countryA_vac1_virW, 3);
            break;
        case 33:
            dJydsigma[33] = 1.0/sigma_ydead_countryA_vac1_virM - 1.0*std::pow(-mydead_countryA_vac1_virM + ydead_countryA_vac1_virM, 2)/std::pow(sigma_ydead_countryA_vac1_virM, 3);
            break;
        case 34:
            dJydsigma[34] = 1.0/sigma_ydead_countryA_vac2_virW - 1.0*std::pow(-mydead_countryA_vac2_virW + ydead_countryA_vac2_virW, 2)/std::pow(sigma_ydead_countryA_vac2_virW, 3);
            break;
        case 35:
            dJydsigma[35] = 1.0/sigma_ydead_countryA_vac2_virM - 1.0*std::pow(-mydead_countryA_vac2_virM + ydead_countryA_vac2_virM, 2)/std::pow(sigma_ydead_countryA_vac2_virM, 3);
            break;
        case 36:
            dJydsigma[36] = 1.0/sigma_ydead_countryB_vac0_virW - 1.0*std::pow(-mydead_countryB_vac0_virW + ydead_countryB_vac0_virW, 2)/std::pow(sigma_ydead_countryB_vac0_virW, 3);
            break;
        case 37:
            dJydsigma[37] = 1.0/sigma_ydead_countryB_vac0_virM - 1.0*std::pow(-mydead_countryB_vac0_virM + ydead_countryB_vac0_virM, 2)/std::pow(sigma_ydead_countryB_vac0_virM, 3);
            break;
        case 38:
            dJydsigma[38] = 1.0/sigma_ydead_countryB_vac1_virW - 1.0*std::pow(-mydead_countryB_vac1_virW + ydead_countryB_vac1_virW, 2)/std::pow(sigma_ydead_countryB_vac1_virW, 3);
            break;
        case 39:
            dJydsigma[39] = 1.0/sigma_ydead_countryB_vac1_virM - 1.0*std::pow(-mydead_countryB_vac1_virM + ydead_countryB_vac1_virM, 2)/std::pow(sigma_ydead_countryB_vac1_virM, 3);
            break;
        case 40:
            dJydsigma[40] = 1.0/sigma_ydead_countryB_vac2_virW - 1.0*std::pow(-mydead_countryB_vac2_virW + ydead_countryB_vac2_virW, 2)/std::pow(sigma_ydead_countryB_vac2_virW, 3);
            break;
        case 41:
            dJydsigma[41] = 1.0/sigma_ydead_countryB_vac2_virM - 1.0*std::pow(-mydead_countryB_vac2_virM + ydead_countryB_vac2_virM, 2)/std::pow(sigma_ydead_countryB_vac2_virM, 3);
            break;
        case 42:
            dJydsigma[42] = 1.0/sigma_yt - 1.0*std::pow(-myt + yt, 2)/std::pow(sigma_yt, 3);
            break;
        case 43:
            dJydsigma[43] = 1.0/sigma_yproportion_countryA_vac1 - 1.0*std::pow(-myproportion_countryA_vac1 + yproportion_countryA_vac1, 2)/std::pow(sigma_yproportion_countryA_vac1, 3);
            break;
        case 44:
            dJydsigma[44] = 1.0/sigma_yproportion_countryA_vac2 - 1.0*std::pow(-myproportion_countryA_vac2 + yproportion_countryA_vac2, 2)/std::pow(sigma_yproportion_countryA_vac2, 3);
            break;
        case 45:
            dJydsigma[45] = 1.0/sigma_yproportion_countryB_vac1 - 1.0*std::pow(-myproportion_countryB_vac1 + yproportion_countryB_vac1, 2)/std::pow(sigma_yproportion_countryB_vac1, 3);
            break;
        case 46:
            dJydsigma[46] = 1.0/sigma_yproportion_countryB_vac2 - 1.0*std::pow(-myproportion_countryB_vac2 + yproportion_countryB_vac2, 2)/std::pow(sigma_yproportion_countryB_vac2, 3);
            break;
        case 47:
            dJydsigma[47] = 1.0/sigma_ynu_countryA_vac1 - 1.0*std::pow(-mynu_countryA_vac1 + ynu_countryA_vac1, 2)/std::pow(sigma_ynu_countryA_vac1, 3);
            break;
        case 48:
            dJydsigma[48] = 1.0/sigma_ynu_countryB_vac1 - 1.0*std::pow(-mynu_countryB_vac1 + ynu_countryB_vac1, 2)/std::pow(sigma_ynu_countryB_vac1, 3);
            break;
        case 49:
            dJydsigma[49] = 1.0/sigma_ynu_countryA_vac2 - 1.0*std::pow(-mynu_countryA_vac2 + ynu_countryA_vac2, 2)/std::pow(sigma_ynu_countryA_vac2, 3);
            break;
        case 50:
            dJydsigma[50] = 1.0/sigma_ynu_countryB_vac2 - 1.0*std::pow(-mynu_countryB_vac2 + ynu_countryB_vac2, 2)/std::pow(sigma_ynu_countryB_vac2, 3);
            break;
        case 51:
            dJydsigma[51] = 1.0/sigma_yspline_countryA_vac1 - 1.0*std::pow(-myspline_countryA_vac1 + yspline_countryA_vac1, 2)/std::pow(sigma_yspline_countryA_vac1, 3);
            break;
        case 52:
            dJydsigma[52] = 1.0/sigma_yspline_countryA_vac2 - 1.0*std::pow(-myspline_countryA_vac2 + yspline_countryA_vac2, 2)/std::pow(sigma_yspline_countryA_vac2, 3);
            break;
        case 53:
            dJydsigma[53] = 1.0/sigma_ycountryA - 1.0*std::pow(-mycountryA + ycountryA, 2)/std::pow(sigma_ycountryA, 3);
            break;
        case 54:
            dJydsigma[54] = 1.0/sigma_ycountryB - 1.0*std::pow(-mycountryB + ycountryB, 2)/std::pow(sigma_ycountryB, 3);
            break;
    }
}

} // namespace model_vaccination_partly
} // namespace amici
