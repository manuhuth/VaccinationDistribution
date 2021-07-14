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

void Jy_vaccination_neural_network(realtype *Jy, const int iy, const realtype *p, const realtype *k, const realtype *y, const realtype *sigmay, const realtype *my){
    switch(iy) {
        case 0:
            Jy[0] = 0.5*std::log(std::pow(sigmaobservable_time, 2)) + 0.5*M_LN2 + 0.5*std::log(amici::pi) + 0.5*std::pow(-mobservable_time + observable_time, 2)/std::pow(sigmaobservable_time, 2);
            break;
        case 1:
            Jy[0] = 0.5*std::log(std::pow(sigmaobservable_nu_countryA_vac1, 2)) + 0.5*M_LN2 + 0.5*std::log(amici::pi) + 0.5*std::pow(-mobservable_nu_countryA_vac1 + observable_nu_countryA_vac1, 2)/std::pow(sigmaobservable_nu_countryA_vac1, 2);
            break;
        case 2:
            Jy[0] = 0.5*std::log(std::pow(sigmaobservable_nu_countryB_vac1, 2)) + 0.5*M_LN2 + 0.5*std::log(amici::pi) + 0.5*std::pow(-mobservable_nu_countryB_vac1 + observable_nu_countryB_vac1, 2)/std::pow(sigmaobservable_nu_countryB_vac1, 2);
            break;
        case 3:
            Jy[0] = 0.5*std::log(std::pow(sigmaobservable_nu_countryA_vac2, 2)) + 0.5*M_LN2 + 0.5*std::log(amici::pi) + 0.5*std::pow(-mobservable_nu_countryA_vac2 + observable_nu_countryA_vac2, 2)/std::pow(sigmaobservable_nu_countryA_vac2, 2);
            break;
        case 4:
            Jy[0] = 0.5*std::log(std::pow(sigmaobservable_nu_countryB_vac2, 2)) + 0.5*M_LN2 + 0.5*std::log(amici::pi) + 0.5*std::pow(-mobservable_nu_countryB_vac2 + observable_nu_countryB_vac2, 2)/std::pow(sigmaobservable_nu_countryB_vac2, 2);
            break;
        case 5:
            Jy[0] = 0.5*std::log(std::pow(sigmaobservable_proportion_countryA_vac1, 2)) + 0.5*M_LN2 + 0.5*std::log(amici::pi) + 0.5*std::pow(-mobservable_proportion_countryA_vac1 + observable_proportion_countryA_vac1, 2)/std::pow(sigmaobservable_proportion_countryA_vac1, 2);
            break;
        case 6:
            Jy[0] = 0.5*std::log(std::pow(sigmaobservable_proportion_countryB_vac1, 2)) + 0.5*M_LN2 + 0.5*std::log(amici::pi) + 0.5*std::pow(-mobservable_proportion_countryB_vac1 + observable_proportion_countryB_vac1, 2)/std::pow(sigmaobservable_proportion_countryB_vac1, 2);
            break;
        case 7:
            Jy[0] = 0.5*std::log(std::pow(sigmaobservable_proportion_countryA_vac2, 2)) + 0.5*M_LN2 + 0.5*std::log(amici::pi) + 0.5*std::pow(-mobservable_proportion_countryA_vac2 + observable_proportion_countryA_vac2, 2)/std::pow(sigmaobservable_proportion_countryA_vac2, 2);
            break;
        case 8:
            Jy[0] = 0.5*std::log(std::pow(sigmaobservable_proportion_countryB_vac2, 2)) + 0.5*M_LN2 + 0.5*std::log(amici::pi) + 0.5*std::pow(-mobservable_proportion_countryB_vac2 + observable_proportion_countryB_vac2, 2)/std::pow(sigmaobservable_proportion_countryB_vac2, 2);
            break;
    }
}

} // namespace amici
} // namespace model_vaccination_neural_network