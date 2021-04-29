#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "test_p.h"
#include "test_y.h"
#include "test_sigmay.h"
#include "test_my.h"

namespace amici {
namespace model_test {

void dJydsigma_test(realtype *dJydsigma, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigma[0] = 1.0/sigma_ysusceptible - 1.0*std::pow(-mysusceptible + ysusceptible, 2)/std::pow(sigma_ysusceptible, 3);
            break;
        case 1:
            dJydsigma[1] = 1.0/sigma_yinfectious - 1.0*std::pow(-myinfectious + yinfectious, 2)/std::pow(sigma_yinfectious, 3);
            break;
        case 2:
            dJydsigma[2] = 1.0/sigma_yrecovered - 1.0*std::pow(-myrecovered + yrecovered, 2)/std::pow(sigma_yrecovered, 3);
            break;
        case 3:
            dJydsigma[3] = 1.0/sigma_yA - 1.0*std::pow(-myA + yA, 2)/std::pow(sigma_yA, 3);
            break;
    }
}

} // namespace model_test
} // namespace amici
