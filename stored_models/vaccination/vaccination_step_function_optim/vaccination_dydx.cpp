#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "vaccination_x.h"
#include "vaccination_p.h"
#include "vaccination_w.h"
#include "vaccination_dwdx.h"

namespace amici {
namespace model_vaccination {

void dydx_vaccination(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = -2*number_vac1*proportion_countryA_vac1/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[2] = -2*number_vac2*proportion_countryA_vac2/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[25] = -2*number_vac1*proportion_countryB_vac1/std::pow(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0, 2);
    dydx[27] = -2*number_vac2*proportion_countryB_vac2/std::pow(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0, 2);
    dydx[144] = -number_vac1*proportion_countryA_vac1/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[146] = -number_vac2*proportion_countryA_vac2/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[152] = -number_vac1*proportion_countryA_vac1/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[154] = -number_vac2*proportion_countryA_vac2/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[193] = -number_vac1*proportion_countryB_vac1/std::pow(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0, 2);
    dydx[195] = -number_vac2*proportion_countryB_vac2/std::pow(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0, 2);
    dydx[201] = -number_vac1*proportion_countryB_vac1/std::pow(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0, 2);
    dydx[203] = -number_vac2*proportion_countryB_vac2/std::pow(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0, 2);
}

} // namespace model_vaccination
} // namespace amici
