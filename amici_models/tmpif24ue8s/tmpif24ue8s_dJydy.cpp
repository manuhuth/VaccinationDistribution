#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

namespace amici {
namespace model_tmpif24ue8s {

void dJydy_tmpif24ue8s(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mobservable_deaths + 1.0*observable_deaths)/std::pow(sigmaobservable_deaths, 2);
            break;
    }
}

} // namespace amici
} // namespace model_tmpif24ue8s