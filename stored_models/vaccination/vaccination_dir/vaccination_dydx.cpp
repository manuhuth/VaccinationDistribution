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
    dydx[13] = -2*number_vac1*proportion_countryB_vac1/std::pow(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0, 2);
    dydx[15] = -2*number_vac2*proportion_countryB_vac2/std::pow(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0, 2);
    dydx[72] = -number_vac1*proportion_countryA_vac1/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[74] = -number_vac2*proportion_countryA_vac2/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[76] = -number_vac1*proportion_countryA_vac1/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[78] = -number_vac2*proportion_countryA_vac2/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[97] = -number_vac1*proportion_countryB_vac1/std::pow(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0, 2);
    dydx[99] = -number_vac2*proportion_countryB_vac2/std::pow(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0, 2);
    dydx[101] = -number_vac1*proportion_countryB_vac1/std::pow(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0, 2);
    dydx[103] = -number_vac2*proportion_countryB_vac2/std::pow(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0, 2);
}

} // namespace model_vaccination
} // namespace amici
