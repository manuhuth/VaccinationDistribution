#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "vaccination_x.h"
#include "vaccination_p.h"
#include "vaccination_h.h"
#include "vaccination_sx.h"

namespace amici {
namespace model_vaccination {

void stau_vaccination(realtype *stau, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *sx, const int ip, const int ie){
    switch(ie) {
        case 0:
            switch(ip) {
                case 0:
                    stau[0] = -1.0*sx42;
                    break;
                case 1:
                    stau[0] = -1.0*sx42;
                    break;
                case 2:
                    stau[0] = -1.0*sx42;
                    break;
                case 3:
                    stau[0] = -1.0*sx42;
                    break;
                case 4:
                    stau[0] = -1.0*sx42;
                    break;
                case 5:
                    stau[0] = -1.0*sx42;
                    break;
                case 6:
                    stau[0] = -1.0*sx42;
                    break;
                case 7:
                    stau[0] = -1.0*sx42;
                    break;
                case 8:
                    stau[0] = -1.0*sx42;
                    break;
                case 9:
                    stau[0] = -1.0*sx42;
                    break;
                case 10:
                    stau[0] = -1.0*sx42;
                    break;
                case 11:
                    stau[0] = -1.0*sx42;
                    break;
                case 12:
                    stau[0] = -1.0*sx42;
                    break;
                case 13:
                    stau[0] = -1.0*sx42;
                    break;
                case 14:
                    stau[0] = -1.0*sx42;
                    break;
                case 15:
                    stau[0] = -1.0*sx42;
                    break;
                case 16:
                    stau[0] = -1.0*sx42;
                    break;
                case 17:
                    stau[0] = -1.0*sx42;
                    break;
                case 18:
                    stau[0] = -1.0*sx42;
                    break;
                case 19:
                    stau[0] = -1.0*sx42;
                    break;
                case 20:
                    stau[0] = -1.0*sx42;
                    break;
                case 21:
                    stau[0] = -1.0*sx42;
                    break;
                case 22:
                    stau[0] = -1.0*sx42;
                    break;
                case 23:
                    stau[0] = -1.0*sx42;
                    break;
                case 24:
                    stau[0] = -1.0*sx42;
                    break;
                case 25:
                    stau[0] = -1.0*sx42;
                    break;
                case 26:
                    stau[0] = -1.0*sx42;
                    break;
                case 27:
                    stau[0] = -1.0*sx42;
                    break;
                case 28:
                    stau[0] = -1.0*sx42;
                    break;
                case 29:
                    stau[0] = -1.0*sx42;
                    break;
                case 30:
                    stau[0] = -1.0*sx42;
                    break;
                case 31:
                    stau[0] = -1.0*sx42;
                    break;
                case 32:
                    stau[0] = -1.0*sx42;
                    break;
                case 33:
                    stau[0] = -1.0*sx42;
                    break;
                case 34:
                    stau[0] = -1.0*sx42;
                    break;
                case 35:
                    stau[0] = -1.0*sx42;
                    break;
                case 36:
                    stau[0] = -1.0*sx42;
                    break;
                case 37:
                    stau[0] = -1.0*sx42;
                    break;
                case 38:
                    stau[0] = -1.0*sx42;
                    break;
                case 39:
                    stau[0] = -1.0*sx42;
                    break;
                case 40:
                    stau[0] = -1.0*sx42;
                    break;
                case 41:
                    stau[0] = -1.0*sx42;
                    break;
                case 42:
                    stau[0] = -1.0*sx42;
                    break;
                case 43:
                    stau[0] = -1.0*sx42;
                    break;
                case 44:
                    stau[0] = -1.0*sx42;
                    break;
                case 45:
                    stau[0] = -1.0*sx42;
                    break;
                case 46:
                    stau[0] = -1.0*sx42;
                    break;
                case 47:
                    stau[0] = -1.0*sx42;
                    break;
                case 48:
                    stau[0] = -1.0*sx42;
                    break;
                case 49:
                    stau[0] = -1.0*sx42;
                    break;
                case 50:
                    stau[0] = -1.0*sx42;
                    break;
                case 51:
                    stau[0] = -1.0*sx42;
                    break;
                case 52:
                    stau[0] = -1.0*sx42;
                    break;
                case 53:
                    stau[0] = -1.0*sx42;
                    break;
                case 54:
                    stau[0] = -1.0*sx42;
                    break;
                case 55:
                    stau[0] = -1.0*sx42;
                    break;
                case 56:
                    stau[0] = -1.0*sx42;
                    break;
                case 57:
                    stau[0] = -1.0*sx42;
                    break;
                case 58:
                    stau[0] = -1.0*sx42;
                    break;
                case 59:
                    stau[0] = -1.0*sx42;
                    break;
                case 60:
                    stau[0] = -1.0*sx42;
                    break;
                case 61:
                    stau[0] = -1.0*sx42;
                    break;
                case 62:
                    stau[0] = -1.0*sx42;
                    break;
            }
            break;
        case 1:
            switch(ip) {
                case 0:
                    stau[0] = -1.0*sx42;
                    break;
                case 1:
                    stau[0] = -1.0*sx42;
                    break;
                case 2:
                    stau[0] = -1.0*sx42;
                    break;
                case 3:
                    stau[0] = -1.0*sx42;
                    break;
                case 4:
                    stau[0] = -1.0*sx42;
                    break;
                case 5:
                    stau[0] = -1.0*sx42;
                    break;
                case 6:
                    stau[0] = -1.0*sx42;
                    break;
                case 7:
                    stau[0] = -1.0*sx42;
                    break;
                case 8:
                    stau[0] = -1.0*sx42;
                    break;
                case 9:
                    stau[0] = -1.0*sx42;
                    break;
                case 10:
                    stau[0] = -1.0*sx42;
                    break;
                case 11:
                    stau[0] = -1.0*sx42;
                    break;
                case 12:
                    stau[0] = -1.0*sx42;
                    break;
                case 13:
                    stau[0] = -1.0*sx42;
                    break;
                case 14:
                    stau[0] = -1.0*sx42;
                    break;
                case 15:
                    stau[0] = -1.0*sx42;
                    break;
                case 16:
                    stau[0] = -1.0*sx42;
                    break;
                case 17:
                    stau[0] = -1.0*sx42;
                    break;
                case 18:
                    stau[0] = -1.0*sx42;
                    break;
                case 19:
                    stau[0] = -1.0*sx42;
                    break;
                case 20:
                    stau[0] = -1.0*sx42;
                    break;
                case 21:
                    stau[0] = -1.0*sx42;
                    break;
                case 22:
                    stau[0] = -1.0*sx42;
                    break;
                case 23:
                    stau[0] = -1.0*sx42;
                    break;
                case 24:
                    stau[0] = -1.0*sx42;
                    break;
                case 25:
                    stau[0] = -1.0*sx42;
                    break;
                case 26:
                    stau[0] = -1.0*sx42;
                    break;
                case 27:
                    stau[0] = -1.0*sx42;
                    break;
                case 28:
                    stau[0] = -1.0*sx42;
                    break;
                case 29:
                    stau[0] = -1.0*sx42;
                    break;
                case 30:
                    stau[0] = -1.0*sx42;
                    break;
                case 31:
                    stau[0] = -1.0*sx42;
                    break;
                case 32:
                    stau[0] = -1.0*sx42;
                    break;
                case 33:
                    stau[0] = -1.0*sx42;
                    break;
                case 34:
                    stau[0] = -1.0*sx42;
                    break;
                case 35:
                    stau[0] = -1.0*sx42;
                    break;
                case 36:
                    stau[0] = -1.0*sx42;
                    break;
                case 37:
                    stau[0] = -1.0*sx42;
                    break;
                case 38:
                    stau[0] = -1.0*sx42;
                    break;
                case 39:
                    stau[0] = -1.0*sx42;
                    break;
                case 40:
                    stau[0] = -1.0*sx42;
                    break;
                case 41:
                    stau[0] = -1.0*sx42;
                    break;
                case 42:
                    stau[0] = -1.0*sx42;
                    break;
                case 43:
                    stau[0] = -1.0*sx42;
                    break;
                case 44:
                    stau[0] = -1.0*sx42;
                    break;
                case 45:
                    stau[0] = -1.0*sx42;
                    break;
                case 46:
                    stau[0] = -1.0*sx42;
                    break;
                case 47:
                    stau[0] = -1.0*sx42;
                    break;
                case 48:
                    stau[0] = -1.0*sx42;
                    break;
                case 49:
                    stau[0] = -1.0*sx42;
                    break;
                case 50:
                    stau[0] = -1.0*sx42;
                    break;
                case 51:
                    stau[0] = -1.0*sx42;
                    break;
                case 52:
                    stau[0] = -1.0*sx42;
                    break;
                case 53:
                    stau[0] = -1.0*sx42;
                    break;
                case 54:
                    stau[0] = -1.0*sx42;
                    break;
                case 55:
                    stau[0] = -1.0*sx42;
                    break;
                case 56:
                    stau[0] = -1.0*sx42;
                    break;
                case 57:
                    stau[0] = -1.0*sx42;
                    break;
                case 58:
                    stau[0] = -1.0*sx42;
                    break;
                case 59:
                    stau[0] = -1.0*sx42;
                    break;
                case 60:
                    stau[0] = -1.0*sx42;
                    break;
                case 61:
                    stau[0] = -1.0*sx42;
                    break;
                case 62:
                    stau[0] = -1.0*sx42;
                    break;
            }
            break;
    }
}

} // namespace model_vaccination
} // namespace amici
