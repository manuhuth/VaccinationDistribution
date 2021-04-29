#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "test_p.h"
#include "test_y.h"
#include "test_sigmay.h"
#include "test_my.h"
#include "test_dJydy.h"

namespace amici {
namespace model_test {

void dJydy_test(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mysusceptible + 1.0*ysusceptible)/std::pow(sigma_ysusceptible, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*myinfectious + 1.0*yinfectious)/std::pow(sigma_yinfectious, 2);
            break;
        case 2:
            dJydy[0] = (-1.0*myrecovered + 1.0*yrecovered)/std::pow(sigma_yrecovered, 2);
            break;
        case 3:
            dJydy[0] = (-1.0*myA + 1.0*yA)/std::pow(sigma_yA, 2);
            break;
    }
}

} // namespace model_test
} // namespace amici
