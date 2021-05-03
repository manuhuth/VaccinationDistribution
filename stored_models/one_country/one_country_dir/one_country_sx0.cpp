#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "one_country_x.h"
#include "one_country_p.h"

namespace amici {
namespace model_one_country {

void sx0_one_country(realtype *sx0, const realtype t,const realtype *x, const realtype *p, const realtype *k, const int ip){
    switch(ip) {
        case 16:
            sx0[0] = 1;
            break;
        case 17:
            sx0[1] = 1;
            break;
        case 18:
            sx0[2] = 1;
            break;
        case 19:
            sx0[3] = 1;
            break;
        case 20:
            sx0[4] = 1;
            break;
        case 21:
            sx0[5] = 1;
            break;
        case 22:
            sx0[6] = 1;
            break;
        case 23:
            sx0[7] = 1;
            break;
        case 24:
            sx0[8] = 1;
            break;
        case 25:
            sx0[9] = 1;
            break;
        case 26:
            sx0[10] = 1;
            break;
        case 27:
            sx0[11] = 1;
            break;
        case 28:
            sx0[12] = 1;
            break;
        case 29:
            sx0[13] = 1;
            break;
        case 30:
            sx0[14] = 1;
            break;
        case 31:
            sx0[15] = 1;
            break;
        case 32:
            sx0[16] = 1;
            break;
        case 33:
            sx0[17] = 1;
            break;
        case 34:
            sx0[18] = 1;
            break;
        case 35:
            sx0[19] = 1;
            break;
        case 36:
            sx0[20] = 1;
            break;
    }
}

} // namespace model_one_country
} // namespace amici
