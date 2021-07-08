#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

namespace amici {
namespace model_vaccination {

void dJydy_vaccination(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mobservable_quantity_countryA_vac1 + 1.0*observable_quantity_countryA_vac1)/std::pow(sigmaobservable_quantity_countryA_vac1, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*mobservable_quantity_countryA_vac2 + 1.0*observable_quantity_countryA_vac2)/std::pow(sigmaobservable_quantity_countryA_vac2, 2);
            break;
        case 2:
            dJydy[0] = (-1.0*mobservable_quantity_countryB_vac1 + 1.0*observable_quantity_countryB_vac1)/std::pow(sigmaobservable_quantity_countryB_vac1, 2);
            break;
        case 3:
            dJydy[0] = (-1.0*mobservable_quantity_countryB_vac2 + 1.0*observable_quantity_countryB_vac2)/std::pow(sigmaobservable_quantity_countryB_vac2, 2);
            break;
    }
}

} // namespace amici
} // namespace model_vaccination