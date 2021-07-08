#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "p.h"
#include "k.h"
#include "x.h"

namespace amici {
namespace model_vaccination {

void sx0_vaccination(realtype *sx0, const realtype t,const realtype *x, const realtype *p, const realtype *k, const int ip){
    switch(ip) {
        case 14:
            sx0[0] = 1;
            break;
        case 15:
            sx0[1] = 1;
            break;
        case 16:
            sx0[2] = 1;
            break;
        case 17:
            sx0[3] = 1;
            break;
        case 18:
            sx0[4] = 1;
            break;
        case 19:
            sx0[5] = 1;
            break;
        case 20:
            sx0[6] = 1;
            break;
        case 21:
            sx0[7] = 1;
            break;
        case 22:
            sx0[8] = 1;
            break;
        case 23:
            sx0[9] = 1;
            break;
        case 24:
            sx0[10] = 1;
            break;
        case 25:
            sx0[11] = 1;
            break;
        case 26:
            sx0[12] = 1;
            break;
        case 27:
            sx0[13] = 1;
            break;
        case 28:
            sx0[14] = 1;
            break;
        case 29:
            sx0[15] = 1;
            break;
        case 30:
            sx0[16] = 1;
            break;
        case 31:
            sx0[17] = 1;
            break;
        case 32:
            sx0[18] = 1;
            break;
        case 33:
            sx0[19] = 1;
            break;
        case 34:
            sx0[20] = 1;
            break;
        case 35:
            sx0[21] = 1;
            break;
        case 36:
            sx0[22] = 1;
            break;
        case 37:
            sx0[23] = 1;
            break;
        case 38:
            sx0[24] = 1;
            break;
        case 39:
            sx0[25] = 1;
            break;
        case 40:
            sx0[26] = 1;
            break;
        case 41:
            sx0[27] = 1;
            break;
        case 42:
            sx0[28] = 1;
            break;
        case 43:
            sx0[29] = 1;
            break;
        case 44:
            sx0[30] = 1;
            break;
        case 45:
            sx0[31] = 1;
            break;
        case 46:
            sx0[32] = 1;
            break;
        case 47:
            sx0[33] = 1;
            break;
        case 48:
            sx0[34] = 1;
            break;
        case 49:
            sx0[35] = 1;
            break;
        case 50:
            sx0[36] = 1;
            break;
        case 51:
            sx0[37] = 1;
            break;
        case 52:
            sx0[38] = 1;
            break;
        case 53:
            sx0[39] = 1;
            break;
        case 54:
            sx0[40] = 1;
            break;
        case 55:
            sx0[41] = 1;
            break;
    }
}

} // namespace amici
} // namespace model_vaccination