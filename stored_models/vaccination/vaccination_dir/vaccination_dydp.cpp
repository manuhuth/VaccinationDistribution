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
            dydp[45] = proportion_countryA_vac1/(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0);
            dydp[46] = proportion_countryB_vac1/(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0);
            break;
        case 57:
            dydp[47] = proportion_countryA_vac2/(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0);
            dydp[48] = proportion_countryB_vac2/(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0);
            break;
        case 58:
            dydp[47] = number_vac2/(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0);
            break;
        case 59:
            dydp[48] = number_vac2/(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0);
            break;
        case 60:
            dydp[43] = ((1.0 - 1.0*amici::heaviside(amici_t - 0.20000000000000001))*(5.0*amici_t*(5.0*amici_t*(2.5*amici_t - 0.5) - 1) + 1) + 1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*((99.5 - 2.5*amici_t)*(0.050505050505050504*amici_t - 0.010101010101010102) - 49.5)*amici::heaviside(amici_t - 0.20000000000000001))*std::exp(-spline_countryA_vac1)/std::pow(1 + std::exp(-spline_countryA_vac1), 2);
            dydp[44] = (-(1.0 - 1.0*amici::heaviside(amici_t - 0.20000000000000001))*(5.0*amici_t*(5.0*amici_t*(2.5*amici_t - 0.5) - 1) + 1) - 1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*((99.5 - 2.5*amici_t)*(0.050505050505050504*amici_t - 0.010101010101010102) - 49.5)*amici::heaviside(amici_t - 0.20000000000000001))*std::exp(-spline_countryA_vac1)/std::pow(1 + std::exp(-spline_countryA_vac1), 2);
            dydp[45] = number_vac1*((1.0 - 1.0*amici::heaviside(amici_t - 0.20000000000000001))*(5.0*amici_t*(5.0*amici_t*(2.5*amici_t - 0.5) - 1) + 1) + 1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*((99.5 - 2.5*amici_t)*(0.050505050505050504*amici_t - 0.010101010101010102) - 49.5)*amici::heaviside(amici_t - 0.20000000000000001))*std::exp(-spline_countryA_vac1)/(std::pow(1 + std::exp(-spline_countryA_vac1), 2)*(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0));
            dydp[46] = -number_vac1*((1.0 - 1.0*amici::heaviside(amici_t - 0.20000000000000001))*(5.0*amici_t*(5.0*amici_t*(2.5*amici_t - 0.5) - 1) + 1) + 1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*((99.5 - 2.5*amici_t)*(0.050505050505050504*amici_t - 0.010101010101010102) - 49.5)*amici::heaviside(amici_t - 0.20000000000000001))*std::exp(-spline_countryA_vac1)/(std::pow(1 + std::exp(-spline_countryA_vac1), 2)*(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0));
            dydp[49] = (1.0 - 1.0*amici::heaviside(amici_t - 0.20000000000000001))*(5.0*amici_t*(5.0*amici_t*(2.5*amici_t - 0.5) - 1) + 1) + 1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*((99.5 - 2.5*amici_t)*(0.050505050505050504*amici_t - 0.010101010101010102) - 49.5)*amici::heaviside(amici_t - 0.20000000000000001);
            break;
        case 61:
            dydp[43] = (5.0*amici_t*(1.0 - 1.0*amici::heaviside(amici_t - 0.20000000000000001))*(5.0*amici_t*(0.50505050505050497 - 2.5252525252525251*amici_t) + 1) + (1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*((0.050505050505050504*amici_t - 0.010101010101010102)*(2.5252525252525251*amici_t - 100.50505050505051) + 49) + 1.0)*amici::heaviside(amici_t - 0.20000000000000001))*std::exp(-spline_countryA_vac1)/std::pow(1 + std::exp(-spline_countryA_vac1), 2);
            dydp[44] = (-5.0*amici_t*(1.0 - 1.0*amici::heaviside(amici_t - 0.20000000000000001))*(5.0*amici_t*(0.50505050505050497 - 2.5252525252525251*amici_t) + 1) - (1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*((0.050505050505050504*amici_t - 0.010101010101010102)*(2.5252525252525251*amici_t - 100.50505050505051) + 49) + 1.0)*amici::heaviside(amici_t - 0.20000000000000001))*std::exp(-spline_countryA_vac1)/std::pow(1 + std::exp(-spline_countryA_vac1), 2);
            dydp[45] = number_vac1*(5.0*amici_t*(1.0 - 1.0*amici::heaviside(amici_t - 0.20000000000000001))*(5.0*amici_t*(0.50505050505050497 - 2.5252525252525251*amici_t) + 1) + (1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*((0.050505050505050504*amici_t - 0.010101010101010102)*(2.5252525252525251*amici_t - 100.50505050505051) + 49) + 1.0)*amici::heaviside(amici_t - 0.20000000000000001))*std::exp(-spline_countryA_vac1)/(std::pow(1 + std::exp(-spline_countryA_vac1), 2)*(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0));
            dydp[46] = -number_vac1*(5.0*amici_t*(1.0 - 1.0*amici::heaviside(amici_t - 0.20000000000000001))*(5.0*amici_t*(0.50505050505050497 - 2.5252525252525251*amici_t) + 1) + (1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*((0.050505050505050504*amici_t - 0.010101010101010102)*(2.5252525252525251*amici_t - 100.50505050505051) + 49) + 1.0)*amici::heaviside(amici_t - 0.20000000000000001))*std::exp(-spline_countryA_vac1)/(std::pow(1 + std::exp(-spline_countryA_vac1), 2)*(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0));
            dydp[49] = 5.0*amici_t*(1.0 - 1.0*amici::heaviside(amici_t - 0.20000000000000001))*(5.0*amici_t*(0.50505050505050497 - 2.5252525252525251*amici_t) + 1) + (1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*((0.050505050505050504*amici_t - 0.010101010101010102)*(2.5252525252525251*amici_t - 100.50505050505051) + 49) + 1.0)*amici::heaviside(amici_t - 0.20000000000000001);
            break;
        case 62:
            dydp[43] = (25.0*std::pow(amici_t, 2)*(1.0 - 1.0*amici::heaviside(amici_t - 0.20000000000000001))*(0.025252525252525249*amici_t - 0.0050505050505050501) + 1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*((1.005050505050505 - 0.025252525252525252*amici_t)*(0.050505050505050504*amici_t - 0.010101010101010102) + 0.5)*amici::heaviside(amici_t - 0.20000000000000001))*std::exp(-spline_countryA_vac1)/std::pow(1 + std::exp(-spline_countryA_vac1), 2);
            dydp[44] = (-25.0*std::pow(amici_t, 2)*(1.0 - 1.0*amici::heaviside(amici_t - 0.20000000000000001))*(0.025252525252525249*amici_t - 0.0050505050505050501) - 1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*((1.005050505050505 - 0.025252525252525252*amici_t)*(0.050505050505050504*amici_t - 0.010101010101010102) + 0.5)*amici::heaviside(amici_t - 0.20000000000000001))*std::exp(-spline_countryA_vac1)/std::pow(1 + std::exp(-spline_countryA_vac1), 2);
            dydp[45] = number_vac1*(25.0*std::pow(amici_t, 2)*(1.0 - 1.0*amici::heaviside(amici_t - 0.20000000000000001))*(0.025252525252525249*amici_t - 0.0050505050505050501) + 1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*((1.005050505050505 - 0.025252525252525252*amici_t)*(0.050505050505050504*amici_t - 0.010101010101010102) + 0.5)*amici::heaviside(amici_t - 0.20000000000000001))*std::exp(-spline_countryA_vac1)/(std::pow(1 + std::exp(-spline_countryA_vac1), 2)*(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0));
            dydp[46] = -number_vac1*(25.0*std::pow(amici_t, 2)*(1.0 - 1.0*amici::heaviside(amici_t - 0.20000000000000001))*(0.025252525252525249*amici_t - 0.0050505050505050501) + 1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*((1.005050505050505 - 0.025252525252525252*amici_t)*(0.050505050505050504*amici_t - 0.010101010101010102) + 0.5)*amici::heaviside(amici_t - 0.20000000000000001))*std::exp(-spline_countryA_vac1)/(std::pow(1 + std::exp(-spline_countryA_vac1), 2)*(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0));
            dydp[49] = 25.0*std::pow(amici_t, 2)*(1.0 - 1.0*amici::heaviside(amici_t - 0.20000000000000001))*(0.025252525252525249*amici_t - 0.0050505050505050501) + 1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*((1.005050505050505 - 0.025252525252525252*amici_t)*(0.050505050505050504*amici_t - 0.010101010101010102) + 0.5)*amici::heaviside(amici_t - 0.20000000000000001);
            break;
    }
}

} // namespace model_vaccination
} // namespace amici
