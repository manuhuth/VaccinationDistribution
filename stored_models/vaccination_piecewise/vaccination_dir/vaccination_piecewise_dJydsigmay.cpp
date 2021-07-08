#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

namespace amici {
namespace model_vaccination_piecewise {

void dJydsigmay_vaccination_piecewise(realtype *dJydsigmay, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydsigmay[0] = 1.0/sigmaobservable_nu_countryA_vac1 - 1.0*std::pow(-mobservable_nu_countryA_vac1 + observable_nu_countryA_vac1, 2)/std::pow(sigmaobservable_nu_countryA_vac1, 3);
            break;
        case 1:
            dJydsigmay[1] = 1.0/sigmaobservable_nu_countryB_vac1 - 1.0*std::pow(-mobservable_nu_countryB_vac1 + observable_nu_countryB_vac1, 2)/std::pow(sigmaobservable_nu_countryB_vac1, 3);
            break;
        case 2:
            dJydsigmay[2] = 1.0/sigmaobservable_nu_countryA_vac2 - 1.0*std::pow(-mobservable_nu_countryA_vac2 + observable_nu_countryA_vac2, 2)/std::pow(sigmaobservable_nu_countryA_vac2, 3);
            break;
        case 3:
            dJydsigmay[3] = 1.0/sigmaobservable_nu_countryB_vac2 - 1.0*std::pow(-mobservable_nu_countryB_vac2 + observable_nu_countryB_vac2, 2)/std::pow(sigmaobservable_nu_countryB_vac2, 3);
            break;
        case 4:
            dJydsigmay[4] = 1.0/sigmaobservable_proportion_countryA_vac1 - 1.0*std::pow(-mobservable_proportion_countryA_vac1 + observable_proportion_countryA_vac1, 2)/std::pow(sigmaobservable_proportion_countryA_vac1, 3);
            break;
        case 5:
            dJydsigmay[5] = 1.0/sigmaobservable_proportion_countryB_vac1 - 1.0*std::pow(-mobservable_proportion_countryB_vac1 + observable_proportion_countryB_vac1, 2)/std::pow(sigmaobservable_proportion_countryB_vac1, 3);
            break;
        case 6:
            dJydsigmay[6] = 1.0/sigmaobservable_proportion_countryA_vac2 - 1.0*std::pow(-mobservable_proportion_countryA_vac2 + observable_proportion_countryA_vac2, 2)/std::pow(sigmaobservable_proportion_countryA_vac2, 3);
            break;
        case 7:
            dJydsigmay[7] = 1.0/sigmaobservable_proportion_countryB_vac2 - 1.0*std::pow(-mobservable_proportion_countryB_vac2 + observable_proportion_countryB_vac2, 2)/std::pow(sigmaobservable_proportion_countryB_vac2, 3);
            break;
        case 8:
            dJydsigmay[8] = 1.0/sigmaobservable_quantity_countryA_vac1 - 1.0*std::pow(-mobservable_quantity_countryA_vac1 + observable_quantity_countryA_vac1, 2)/std::pow(sigmaobservable_quantity_countryA_vac1, 3);
            break;
        case 9:
            dJydsigmay[9] = 1.0/sigmaobservable_quantity_countryA_vac2 - 1.0*std::pow(-mobservable_quantity_countryA_vac2 + observable_quantity_countryA_vac2, 2)/std::pow(sigmaobservable_quantity_countryA_vac2, 3);
            break;
        case 10:
            dJydsigmay[10] = 1.0/sigmaobservable_quantity_countryB_vac1 - 1.0*std::pow(-mobservable_quantity_countryB_vac1 + observable_quantity_countryB_vac1, 2)/std::pow(sigmaobservable_quantity_countryB_vac1, 3);
            break;
        case 11:
            dJydsigmay[11] = 1.0/sigmaobservable_quantity_countryB_vac2 - 1.0*std::pow(-mobservable_quantity_countryB_vac2 + observable_quantity_countryB_vac2, 2)/std::pow(sigmaobservable_quantity_countryB_vac2, 3);
            break;
    }
}

} // namespace amici
} // namespace model_vaccination_piecewise