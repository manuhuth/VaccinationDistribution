#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "vaccination_x.h"
#include "vaccination_p.h"
#include "vaccination_h.h"
#include "vaccination_w.h"
#include "vaccination_dwdx.h"

namespace amici {
namespace model_vaccination {

void dydx_vaccination(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = -2*number_vac1*proportion_countryA_vac1/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[16] = -2*number_vac1*proportion_countryB_vac1/std::pow(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0, 2);
    dydx[90] = -number_vac1*proportion_countryA_vac1/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[95] = -number_vac1*proportion_countryA_vac1/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[121] = -number_vac1*proportion_countryB_vac1/std::pow(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0, 2);
    dydx[126] = -number_vac1*proportion_countryB_vac1/std::pow(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0, 2);
    dydx[210] = number_vac1*(2.0*proportion_par_countryA_vac1_0*(1 - amici::heaviside(amici_t - 200))*amici::dirac(amici_t)*amici::heaviside(amici_t) - 1.0*proportion_par_countryA_vac1_0*amici::dirac(amici_t - 200)*std::pow(amici::heaviside(amici_t), 2) - 1.0*proportion_par_countryA_vac1_200*((amici::heaviside(amici_t - 200) - 1)*amici::dirac(amici_t) + amici::dirac(amici_t - 200)*amici::heaviside(amici_t))*(amici::heaviside(amici_t - 400) - 1)*amici::heaviside(amici_t)*amici::heaviside(amici_t - 200) - 1.0*proportion_par_countryA_vac1_200*((amici::heaviside(amici_t - 200) - 1)*amici::heaviside(amici_t) + 1)*(amici::heaviside(amici_t - 400) - 1)*amici::dirac(amici_t)*amici::heaviside(amici_t - 200) - 1.0*proportion_par_countryA_vac1_200*((amici::heaviside(amici_t - 200) - 1)*amici::heaviside(amici_t) + 1)*(amici::heaviside(amici_t - 400) - 1)*amici::dirac(amici_t - 200)*amici::heaviside(amici_t) - 1.0*proportion_par_countryA_vac1_200*((amici::heaviside(amici_t - 200) - 1)*amici::heaviside(amici_t) + 1)*amici::dirac(amici_t - 400)*amici::heaviside(amici_t)*amici::heaviside(amici_t - 200))/(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0);
    dydx[211] = number_vac1*(2.0*(1 - proportion_par_countryA_vac1_0)*(1 - amici::heaviside(amici_t - 200))*amici::dirac(amici_t)*amici::heaviside(amici_t) - 1.0*(1 - proportion_par_countryA_vac1_0)*amici::dirac(amici_t - 200)*std::pow(amici::heaviside(amici_t), 2) + (1.0*proportion_par_countryA_vac1_200 - 1.0)*((amici::heaviside(amici_t - 200) - 1)*amici::dirac(amici_t) + amici::dirac(amici_t - 200)*amici::heaviside(amici_t))*(amici::heaviside(amici_t - 400) - 1)*amici::heaviside(amici_t)*amici::heaviside(amici_t - 200) + (1.0*proportion_par_countryA_vac1_200 - 1.0)*((amici::heaviside(amici_t - 200) - 1)*amici::heaviside(amici_t) + 1)*(amici::heaviside(amici_t - 400) - 1)*amici::dirac(amici_t)*amici::heaviside(amici_t - 200) + (1.0*proportion_par_countryA_vac1_200 - 1.0)*((amici::heaviside(amici_t - 200) - 1)*amici::heaviside(amici_t) + 1)*(amici::heaviside(amici_t - 400) - 1)*amici::dirac(amici_t - 200)*amici::heaviside(amici_t) + (1.0*proportion_par_countryA_vac1_200 - 1.0)*((amici::heaviside(amici_t - 200) - 1)*amici::heaviside(amici_t) + 1)*amici::dirac(amici_t - 400)*amici::heaviside(amici_t)*amici::heaviside(amici_t - 200))/(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0);
    dydx[212] = 2.0*proportion_par_countryA_vac1_0*(1 - amici::heaviside(amici_t - 200))*amici::dirac(amici_t)*amici::heaviside(amici_t) - 1.0*proportion_par_countryA_vac1_0*amici::dirac(amici_t - 200)*std::pow(amici::heaviside(amici_t), 2) - 1.0*proportion_par_countryA_vac1_200*((amici::heaviside(amici_t - 200) - 1)*amici::dirac(amici_t) + amici::dirac(amici_t - 200)*amici::heaviside(amici_t))*(amici::heaviside(amici_t - 400) - 1)*amici::heaviside(amici_t)*amici::heaviside(amici_t - 200) - 1.0*proportion_par_countryA_vac1_200*((amici::heaviside(amici_t - 200) - 1)*amici::heaviside(amici_t) + 1)*(amici::heaviside(amici_t - 400) - 1)*amici::dirac(amici_t)*amici::heaviside(amici_t - 200) - 1.0*proportion_par_countryA_vac1_200*((amici::heaviside(amici_t - 200) - 1)*amici::heaviside(amici_t) + 1)*(amici::heaviside(amici_t - 400) - 1)*amici::dirac(amici_t - 200)*amici::heaviside(amici_t) - 1.0*proportion_par_countryA_vac1_200*((amici::heaviside(amici_t - 200) - 1)*amici::heaviside(amici_t) + 1)*amici::dirac(amici_t - 400)*amici::heaviside(amici_t)*amici::heaviside(amici_t - 200);
    dydx[213] = 2.0*(1 - proportion_par_countryA_vac1_0)*(1 - amici::heaviside(amici_t - 200))*amici::dirac(amici_t)*amici::heaviside(amici_t) - 1.0*(1 - proportion_par_countryA_vac1_0)*amici::dirac(amici_t - 200)*std::pow(amici::heaviside(amici_t), 2) + (1.0*proportion_par_countryA_vac1_200 - 1.0)*((amici::heaviside(amici_t - 200) - 1)*amici::dirac(amici_t) + amici::dirac(amici_t - 200)*amici::heaviside(amici_t))*(amici::heaviside(amici_t - 400) - 1)*amici::heaviside(amici_t)*amici::heaviside(amici_t - 200) + (1.0*proportion_par_countryA_vac1_200 - 1.0)*((amici::heaviside(amici_t - 200) - 1)*amici::heaviside(amici_t) + 1)*(amici::heaviside(amici_t - 400) - 1)*amici::dirac(amici_t)*amici::heaviside(amici_t - 200) + (1.0*proportion_par_countryA_vac1_200 - 1.0)*((amici::heaviside(amici_t - 200) - 1)*amici::heaviside(amici_t) + 1)*(amici::heaviside(amici_t - 400) - 1)*amici::dirac(amici_t - 200)*amici::heaviside(amici_t) + (1.0*proportion_par_countryA_vac1_200 - 1.0)*((amici::heaviside(amici_t - 200) - 1)*amici::heaviside(amici_t) + 1)*amici::dirac(amici_t - 400)*amici::heaviside(amici_t)*amici::heaviside(amici_t - 200);
    dydx[214] = 1;
}

} // namespace model_vaccination
} // namespace amici
