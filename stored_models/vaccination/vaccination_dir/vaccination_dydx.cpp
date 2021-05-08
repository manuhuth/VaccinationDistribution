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
    dydx[0] = 1;
    dydx[21] = -2*number_vac1*proportion_countryA_vac1/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[22] = -2*number_vac2*proportion_countryA_vac2/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[25] = 1;
    dydx[50] = 1;
    dydx[75] = 1;
    dydx[100] = 1;
    dydx[125] = 1;
    dydx[150] = 1;
    dydx[175] = 1;
    dydx[200] = 1;
    dydx[225] = 1;
    dydx[237] = -number_vac1*proportion_countryA_vac1/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[238] = -number_vac2*proportion_countryA_vac2/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[250] = 1;
    dydx[261] = -number_vac1*proportion_countryA_vac1/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[262] = -number_vac2*proportion_countryA_vac2/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[275] = 1;
    dydx[300] = 1;
    dydx[325] = 1;
    dydx[350] = 1;
    dydx[375] = 1;
    dydx[400] = 1;
    dydx[425] = 1;
    dydx[450] = 1;
    dydx[475] = 1;
    dydx[500] = 1;
}

} // namespace model_vaccination
} // namespace amici
