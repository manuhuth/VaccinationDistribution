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

void Jy_test(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_ysusceptible, 2)) + 0.5*std::pow(-mysusceptible + ysusceptible, 2)/std::pow(sigma_ysusceptible, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yinfectious, 2)) + 0.5*std::pow(-myinfectious + yinfectious, 2)/std::pow(sigma_yinfectious, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yrecovered, 2)) + 0.5*std::pow(-myrecovered + yrecovered, 2)/std::pow(sigma_yrecovered, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(2*amici::pi*std::pow(sigma_yA, 2)) + 0.5*std::pow(-myA + yA, 2)/std::pow(sigma_yA, 2);
            break;
    }
}

} // namespace model_test
} // namespace amici
