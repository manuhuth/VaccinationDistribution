#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "y.h"
#include "my.h"
#include "p.h"
#include "k.h"
#include "sigmay.h"

namespace amici {
namespace model_vaccination_neural_network {

void dJydy_vaccination_neural_network(realtype *dJydy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            dJydy[0] = (-1.0*mobservable_time + 1.0*observable_time)/std::pow(sigmaobservable_time, 2);
            break;
        case 1:
            dJydy[0] = (-1.0*mobservable_nu_countryA_vac1 + 1.0*observable_nu_countryA_vac1)/std::pow(sigmaobservable_nu_countryA_vac1, 2);
            break;
        case 2:
            dJydy[0] = (-1.0*mobservable_nu_countryB_vac1 + 1.0*observable_nu_countryB_vac1)/std::pow(sigmaobservable_nu_countryB_vac1, 2);
            break;
        case 3:
            dJydy[0] = (-1.0*mobservable_nu_countryA_vac2 + 1.0*observable_nu_countryA_vac2)/std::pow(sigmaobservable_nu_countryA_vac2, 2);
            break;
        case 4:
            dJydy[0] = (-1.0*mobservable_nu_countryB_vac2 + 1.0*observable_nu_countryB_vac2)/std::pow(sigmaobservable_nu_countryB_vac2, 2);
            break;
        case 5:
            dJydy[0] = (-1.0*mobservable_proportion_countryA_vac1 + 1.0*observable_proportion_countryA_vac1)/std::pow(sigmaobservable_proportion_countryA_vac1, 2);
            break;
        case 6:
            dJydy[0] = (-1.0*mobservable_proportion_countryB_vac1 + 1.0*observable_proportion_countryB_vac1)/std::pow(sigmaobservable_proportion_countryB_vac1, 2);
            break;
        case 7:
            dJydy[0] = (-1.0*mobservable_proportion_countryA_vac2 + 1.0*observable_proportion_countryA_vac2)/std::pow(sigmaobservable_proportion_countryA_vac2, 2);
            break;
        case 8:
            dJydy[0] = (-1.0*mobservable_proportion_countryB_vac2 + 1.0*observable_proportion_countryB_vac2)/std::pow(sigmaobservable_proportion_countryB_vac2, 2);
            break;
    }
}

} // namespace amici
} // namespace model_vaccination_neural_network