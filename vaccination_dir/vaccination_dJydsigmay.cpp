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

void dJydsigmay_vaccination(realtype *dJydsigmay, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigmay[0] = 1.0/sigmaobservable_quantity_countryA_vac1 - 1.0*std::pow(-mobservable_quantity_countryA_vac1 + observable_quantity_countryA_vac1, 2)/std::pow(sigmaobservable_quantity_countryA_vac1, 3);
            break;
        case 1:
            dJydsigmay[1] = 1.0/sigmaobservable_quantity_countryA_vac2 - 1.0*std::pow(-mobservable_quantity_countryA_vac2 + observable_quantity_countryA_vac2, 2)/std::pow(sigmaobservable_quantity_countryA_vac2, 3);
            break;
        case 2:
            dJydsigmay[2] = 1.0/sigmaobservable_quantity_countryB_vac1 - 1.0*std::pow(-mobservable_quantity_countryB_vac1 + observable_quantity_countryB_vac1, 2)/std::pow(sigmaobservable_quantity_countryB_vac1, 3);
            break;
        case 3:
            dJydsigmay[3] = 1.0/sigmaobservable_quantity_countryB_vac2 - 1.0*std::pow(-mobservable_quantity_countryB_vac2 + observable_quantity_countryB_vac2, 2)/std::pow(sigmaobservable_quantity_countryB_vac2, 3);
            break;
    }
}

} // namespace amici
} // namespace model_vaccination