#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "vaccination_x.h"
#include "vaccination_p.h"
#include "vaccination_h.h"
#include "vaccination_w.h"

namespace amici {
namespace model_vaccination {

void dydp_vaccination(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dtcldp){
    switch(ip) {
        case 56:
            dydp[0] = proportion_countryA_vac1/(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0);
            dydp[1] = proportion_countryB_vac1/(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0);
            break;
        case 60:
            dydp[0] = number_vac1*(1.0 - 1.0*amici::heaviside(amici_t - 200))*std::pow(amici::heaviside(amici_t), 2)/(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0);
            dydp[1] = number_vac1*(1.0*amici::heaviside(amici_t - 200) - 1.0)*std::pow(amici::heaviside(amici_t), 2)/(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0);
            dydp[2] = (1.0 - 1.0*amici::heaviside(amici_t - 200))*std::pow(amici::heaviside(amici_t), 2);
            dydp[3] = (1.0*amici::heaviside(amici_t - 200) - 1.0)*std::pow(amici::heaviside(amici_t), 2);
            break;
        case 61:
            dydp[0] = number_vac1*(-1.0*(amici::heaviside(amici_t - 200) - 1)*amici::heaviside(amici_t) - 1.0)*(amici::heaviside(amici_t - 400) - 1)*amici::heaviside(amici_t)*amici::heaviside(amici_t - 200)/(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0);
            dydp[1] = number_vac1*(1.0*(amici::heaviside(amici_t - 200) - 1)*amici::heaviside(amici_t) + 1.0)*(amici::heaviside(amici_t - 400) - 1)*amici::heaviside(amici_t)*amici::heaviside(amici_t - 200)/(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0);
            dydp[2] = (-1.0*(amici::heaviside(amici_t - 200) - 1)*amici::heaviside(amici_t) - 1.0)*(amici::heaviside(amici_t - 400) - 1)*amici::heaviside(amici_t)*amici::heaviside(amici_t - 200);
            dydp[3] = (1.0*(amici::heaviside(amici_t - 200) - 1)*amici::heaviside(amici_t) + 1.0)*(amici::heaviside(amici_t - 400) - 1)*amici::heaviside(amici_t)*amici::heaviside(amici_t - 200);
            break;
    }
}

} // namespace model_vaccination
} // namespace amici
