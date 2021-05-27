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
    dydx[0] = 1;
    dydx[45] = -2*number_vac1*proportion_countryA_vac1/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[47] = -2*number_vac2*proportion_countryA_vac2/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[53] = 1;
    dydx[106] = 1;
    dydx[159] = 1;
    dydx[202] = -2*number_vac1*proportion_countryB_vac1/std::pow(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0, 2);
    dydx[204] = -2*number_vac2*proportion_countryB_vac2/std::pow(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0, 2);
    dydx[212] = 1;
    dydx[265] = 1;
    dydx[318] = 1;
    dydx[371] = 1;
    dydx[424] = 1;
    dydx[477] = 1;
    dydx[530] = 1;
    dydx[583] = 1;
    dydx[636] = 1;
    dydx[689] = 1;
    dydx[742] = 1;
    dydx[795] = 1;
    dydx[848] = 1;
    dydx[901] = 1;
    dydx[954] = 1;
    dydx[981] = -number_vac1*proportion_countryA_vac1/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[983] = -number_vac2*proportion_countryA_vac2/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[1007] = 1;
    dydx[1033] = -number_vac1*proportion_countryA_vac1/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[1035] = -number_vac2*proportion_countryA_vac2/std::pow(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0, 2);
    dydx[1060] = 1;
    dydx[1113] = 1;
    dydx[1166] = 1;
    dydx[1219] = 1;
    dydx[1272] = 1;
    dydx[1294] = -number_vac1*proportion_countryB_vac1/std::pow(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0, 2);
    dydx[1296] = -number_vac2*proportion_countryB_vac2/std::pow(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0, 2);
    dydx[1325] = 1;
    dydx[1346] = -number_vac1*proportion_countryB_vac1/std::pow(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0, 2);
    dydx[1348] = -number_vac2*proportion_countryB_vac2/std::pow(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0, 2);
    dydx[1378] = 1;
    dydx[1431] = 1;
    dydx[1484] = 1;
    dydx[1537] = 1;
    dydx[1590] = 1;
    dydx[1643] = 1;
    dydx[1696] = 1;
    dydx[1749] = 1;
    dydx[1802] = 1;
    dydx[1855] = 1;
    dydx[1908] = 1;
    dydx[1961] = 1;
    dydx[2014] = 1;
    dydx[2067] = 1;
    dydx[2120] = 1;
    dydx[2173] = 1;
    dydx[2226] = 1;
    dydx[2227] = ((1.0 - 1.0*amici::heaviside(amici_t - 0.20000000000000001))*(25.0*amici_t*(5.0*amici_t*(0.5*y1 - 0.50505050505050497*y2 + 0.0050505050505050501*y3) - 0.5*y1 + 0.50505050505050497*y2 - 0.0050505050505050501*y3) + 5.0*amici_t*(25.0*amici_t*(0.5*y1 - 0.50505050505050497*y2 + 0.0050505050505050501*y3) + 5.0*amici_t*(2.5*y1 - 2.5252525252525251*y2 + 0.025252525252525249*y3) - 2.5*y1 + 2.5252525252525251*y2 - 0.025252525252525249*y3) - 5.0*y1 + 5.0*y2) + (1.0*y2 + 1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*(-49.5*y1 + 49*y2 + 0.5*y3 + (0.050505050505050504*amici_t - 0.010101010101010102)*(99*y1 - 100*y2 + y3 + (0.050505050505050504*amici_t - 0.010101010101010102)*(-49.5*y1 + 50*y2 - 0.5*y3))))*amici::dirac(amici_t - 0.20000000000000001) - 1.0*(5.0*amici_t*(5.0*amici_t*(5.0*amici_t*(0.5*y1 - 0.50505050505050497*y2 + 0.0050505050505050501*y3) - 0.5*y1 + 0.50505050505050497*y2 - 0.0050505050505050501*y3) - y1 + y2) + y1)*amici::dirac(amici_t - 0.20000000000000001) + (-2.5*y1 + 2.4747474747474749*y2 + 0.025252525252525252*y3 + 0.050505050505050504*(0.050505050505050504*amici_t - 0.010101010101010102)*(99*y1 - 100*y2 + y3 + (0.050505050505050504*amici_t - 0.010101010101010102)*(-49.5*y1 + 50*y2 - 0.5*y3)) + 1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*(5.0*y1 - 5.0505050505050502*y2 + 0.050505050505050504*y3 + 0.050505050505050504*(0.050505050505050504*amici_t - 0.010101010101010102)*(-49.5*y1 + 50*y2 - 0.5*y3) + (0.050505050505050504*amici_t - 0.010101010101010102)*(-2.5*y1 + 2.5252525252525251*y2 - 0.025252525252525252*y3)))*amici::heaviside(amici_t - 0.20000000000000001))*std::exp(-spline_countryA_vac1)/std::pow(1 + std::exp(-spline_countryA_vac1), 2);
    dydx[2228] = (-(1.0 - 1.0*amici::heaviside(amici_t - 0.20000000000000001))*(25.0*amici_t*(5.0*amici_t*(0.5*y1 - 0.50505050505050497*y2 + 0.0050505050505050501*y3) - 0.5*y1 + 0.50505050505050497*y2 - 0.0050505050505050501*y3) + 5.0*amici_t*(25.0*amici_t*(0.5*y1 - 0.50505050505050497*y2 + 0.0050505050505050501*y3) + 5.0*amici_t*(2.5*y1 - 2.5252525252525251*y2 + 0.025252525252525249*y3) - 2.5*y1 + 2.5252525252525251*y2 - 0.025252525252525249*y3) - 5.0*y1 + 5.0*y2) - (1.0*y2 + 1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*(-49.5*y1 + 49*y2 + 0.5*y3 + (0.050505050505050504*amici_t - 0.010101010101010102)*(99*y1 - 100*y2 + y3 + (0.050505050505050504*amici_t - 0.010101010101010102)*(-49.5*y1 + 50*y2 - 0.5*y3))))*amici::dirac(amici_t - 0.20000000000000001) + 1.0*(5.0*amici_t*(5.0*amici_t*(5.0*amici_t*(0.5*y1 - 0.50505050505050497*y2 + 0.0050505050505050501*y3) - 0.5*y1 + 0.50505050505050497*y2 - 0.0050505050505050501*y3) - y1 + y2) + y1)*amici::dirac(amici_t - 0.20000000000000001) - (-2.5*y1 + 2.4747474747474749*y2 + 0.025252525252525252*y3 + 0.050505050505050504*(0.050505050505050504*amici_t - 0.010101010101010102)*(99*y1 - 100*y2 + y3 + (0.050505050505050504*amici_t - 0.010101010101010102)*(-49.5*y1 + 50*y2 - 0.5*y3)) + 1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*(5.0*y1 - 5.0505050505050502*y2 + 0.050505050505050504*y3 + 0.050505050505050504*(0.050505050505050504*amici_t - 0.010101010101010102)*(-49.5*y1 + 50*y2 - 0.5*y3) + (0.050505050505050504*amici_t - 0.010101010101010102)*(-2.5*y1 + 2.5252525252525251*y2 - 0.025252525252525252*y3)))*amici::heaviside(amici_t - 0.20000000000000001))*std::exp(-spline_countryA_vac1)/std::pow(1 + std::exp(-spline_countryA_vac1), 2);
    dydx[2229] = number_vac1*((1.0 - 1.0*amici::heaviside(amici_t - 0.20000000000000001))*(25.0*amici_t*(5.0*amici_t*(0.5*y1 - 0.50505050505050497*y2 + 0.0050505050505050501*y3) - 0.5*y1 + 0.50505050505050497*y2 - 0.0050505050505050501*y3) + 5.0*amici_t*(25.0*amici_t*(0.5*y1 - 0.50505050505050497*y2 + 0.0050505050505050501*y3) + 5.0*amici_t*(2.5*y1 - 2.5252525252525251*y2 + 0.025252525252525249*y3) - 2.5*y1 + 2.5252525252525251*y2 - 0.025252525252525249*y3) - 5.0*y1 + 5.0*y2) + (1.0*y2 + 1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*(-49.5*y1 + 49*y2 + 0.5*y3 + (0.050505050505050504*amici_t - 0.010101010101010102)*(99*y1 - 100*y2 + y3 + (0.050505050505050504*amici_t - 0.010101010101010102)*(-49.5*y1 + 50*y2 - 0.5*y3))))*amici::dirac(amici_t - 0.20000000000000001) - 1.0*(5.0*amici_t*(5.0*amici_t*(5.0*amici_t*(0.5*y1 - 0.50505050505050497*y2 + 0.0050505050505050501*y3) - 0.5*y1 + 0.50505050505050497*y2 - 0.0050505050505050501*y3) - y1 + y2) + y1)*amici::dirac(amici_t - 0.20000000000000001) + (-2.5*y1 + 2.4747474747474749*y2 + 0.025252525252525252*y3 + 0.050505050505050504*(0.050505050505050504*amici_t - 0.010101010101010102)*(99*y1 - 100*y2 + y3 + (0.050505050505050504*amici_t - 0.010101010101010102)*(-49.5*y1 + 50*y2 - 0.5*y3)) + 1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*(5.0*y1 - 5.0505050505050502*y2 + 0.050505050505050504*y3 + 0.050505050505050504*(0.050505050505050504*amici_t - 0.010101010101010102)*(-49.5*y1 + 50*y2 - 0.5*y3) + (0.050505050505050504*amici_t - 0.010101010101010102)*(-2.5*y1 + 2.5252525252525251*y2 - 0.025252525252525252*y3)))*amici::heaviside(amici_t - 0.20000000000000001))*std::exp(-spline_countryA_vac1)/(std::pow(1 + std::exp(-spline_countryA_vac1), 2)*(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0));
    dydx[2230] = -number_vac1*((1.0 - 1.0*amici::heaviside(amici_t - 0.20000000000000001))*(25.0*amici_t*(5.0*amici_t*(0.5*y1 - 0.50505050505050497*y2 + 0.0050505050505050501*y3) - 0.5*y1 + 0.50505050505050497*y2 - 0.0050505050505050501*y3) + 5.0*amici_t*(25.0*amici_t*(0.5*y1 - 0.50505050505050497*y2 + 0.0050505050505050501*y3) + 5.0*amici_t*(2.5*y1 - 2.5252525252525251*y2 + 0.025252525252525249*y3) - 2.5*y1 + 2.5252525252525251*y2 - 0.025252525252525249*y3) - 5.0*y1 + 5.0*y2) + (1.0*y2 + 1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*(-49.5*y1 + 49*y2 + 0.5*y3 + (0.050505050505050504*amici_t - 0.010101010101010102)*(99*y1 - 100*y2 + y3 + (0.050505050505050504*amici_t - 0.010101010101010102)*(-49.5*y1 + 50*y2 - 0.5*y3))))*amici::dirac(amici_t - 0.20000000000000001) - 1.0*(5.0*amici_t*(5.0*amici_t*(5.0*amici_t*(0.5*y1 - 0.50505050505050497*y2 + 0.0050505050505050501*y3) - 0.5*y1 + 0.50505050505050497*y2 - 0.0050505050505050501*y3) - y1 + y2) + y1)*amici::dirac(amici_t - 0.20000000000000001) + (-2.5*y1 + 2.4747474747474749*y2 + 0.025252525252525252*y3 + 0.050505050505050504*(0.050505050505050504*amici_t - 0.010101010101010102)*(99*y1 - 100*y2 + y3 + (0.050505050505050504*amici_t - 0.010101010101010102)*(-49.5*y1 + 50*y2 - 0.5*y3)) + 1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*(5.0*y1 - 5.0505050505050502*y2 + 0.050505050505050504*y3 + 0.050505050505050504*(0.050505050505050504*amici_t - 0.010101010101010102)*(-49.5*y1 + 50*y2 - 0.5*y3) + (0.050505050505050504*amici_t - 0.010101010101010102)*(-2.5*y1 + 2.5252525252525251*y2 - 0.025252525252525252*y3)))*amici::heaviside(amici_t - 0.20000000000000001))*std::exp(-spline_countryA_vac1)/(std::pow(1 + std::exp(-spline_countryA_vac1), 2)*(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0));
    dydx[2233] = (1.0 - 1.0*amici::heaviside(amici_t - 0.20000000000000001))*(25.0*amici_t*(5.0*amici_t*(0.5*y1 - 0.50505050505050497*y2 + 0.0050505050505050501*y3) - 0.5*y1 + 0.50505050505050497*y2 - 0.0050505050505050501*y3) + 5.0*amici_t*(25.0*amici_t*(0.5*y1 - 0.50505050505050497*y2 + 0.0050505050505050501*y3) + 5.0*amici_t*(2.5*y1 - 2.5252525252525251*y2 + 0.025252525252525249*y3) - 2.5*y1 + 2.5252525252525251*y2 - 0.025252525252525249*y3) - 5.0*y1 + 5.0*y2) + (1.0*y2 + 1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*(-49.5*y1 + 49*y2 + 0.5*y3 + (0.050505050505050504*amici_t - 0.010101010101010102)*(99*y1 - 100*y2 + y3 + (0.050505050505050504*amici_t - 0.010101010101010102)*(-49.5*y1 + 50*y2 - 0.5*y3))))*amici::dirac(amici_t - 0.20000000000000001) - 1.0*(5.0*amici_t*(5.0*amici_t*(5.0*amici_t*(0.5*y1 - 0.50505050505050497*y2 + 0.0050505050505050501*y3) - 0.5*y1 + 0.50505050505050497*y2 - 0.0050505050505050501*y3) - y1 + y2) + y1)*amici::dirac(amici_t - 0.20000000000000001) + (-2.5*y1 + 2.4747474747474749*y2 + 0.025252525252525252*y3 + 0.050505050505050504*(0.050505050505050504*amici_t - 0.010101010101010102)*(99*y1 - 100*y2 + y3 + (0.050505050505050504*amici_t - 0.010101010101010102)*(-49.5*y1 + 50*y2 - 0.5*y3)) + 1.0*(0.050505050505050504*amici_t - 0.010101010101010102)*(5.0*y1 - 5.0505050505050502*y2 + 0.050505050505050504*y3 + 0.050505050505050504*(0.050505050505050504*amici_t - 0.010101010101010102)*(-49.5*y1 + 50*y2 - 0.5*y3) + (0.050505050505050504*amici_t - 0.010101010101010102)*(-2.5*y1 + 2.5252525252525251*y2 - 0.025252525252525252*y3)))*amici::heaviside(amici_t - 0.20000000000000001);
}

} // namespace model_vaccination
} // namespace amici
