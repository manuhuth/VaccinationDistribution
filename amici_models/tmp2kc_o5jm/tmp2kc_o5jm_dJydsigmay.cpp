#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

namespace amici {
namespace model_tmp2kc_o5jm {

void dJydsigmay_tmp2kc_o5jm(realtype *dJydsigmay, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigmay[0] = 1.0/sigmaobservable_deaths - 1.0*std::pow(-mobservable_deaths + observable_deaths, 2)/std::pow(sigmaobservable_deaths, 3);
            break;
    }
}

} // namespace amici
} // namespace model_tmp2kc_o5jm