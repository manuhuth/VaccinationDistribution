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

void y_vaccination(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = susceptible_countryA_vac0;
    y[1] = susceptible_countryA_vac1;
    y[2] = susceptible_countryA_vac2;
    y[3] = susceptible_countryB_vac0;
    y[4] = susceptible_countryB_vac1;
    y[5] = susceptible_countryB_vac2;
    y[6] = infectious_countryA_vac0_virW;
    y[7] = infectious_countryA_vac0_virM;
    y[8] = infectious_countryA_vac1_virW;
    y[9] = infectious_countryA_vac1_virM;
    y[10] = infectious_countryA_vac2_virW;
    y[11] = infectious_countryA_vac2_virM;
    y[12] = infectious_countryB_vac0_virW;
    y[13] = infectious_countryB_vac0_virM;
    y[14] = infectious_countryB_vac1_virW;
    y[15] = infectious_countryB_vac1_virM;
    y[16] = infectious_countryB_vac2_virW;
    y[17] = infectious_countryB_vac2_virM;
    y[18] = recovered_countryA_vac0_virW;
    y[19] = recovered_countryA_vac0_virM;
    y[20] = recovered_countryA_vac1_virW;
    y[21] = recovered_countryA_vac1_virM;
    y[22] = recovered_countryA_vac2_virW;
    y[23] = recovered_countryA_vac2_virM;
    y[24] = recovered_countryB_vac0_virW;
    y[25] = recovered_countryB_vac0_virM;
    y[26] = recovered_countryB_vac1_virW;
    y[27] = recovered_countryB_vac1_virM;
    y[28] = recovered_countryB_vac2_virW;
    y[29] = recovered_countryB_vac2_virM;
    y[30] = dead_countryA_vac0_virW;
    y[31] = dead_countryA_vac0_virM;
    y[32] = dead_countryA_vac1_virW;
    y[33] = dead_countryA_vac1_virM;
    y[34] = dead_countryA_vac2_virW;
    y[35] = dead_countryA_vac2_virM;
    y[36] = dead_countryB_vac0_virW;
    y[37] = dead_countryB_vac0_virM;
    y[38] = dead_countryB_vac1_virW;
    y[39] = dead_countryB_vac1_virM;
    y[40] = dead_countryB_vac2_virW;
    y[41] = dead_countryB_vac2_virM;
    y[42] = amici_t;
    y[43] = 1.0/(1 + std::exp(-spline_countryA_vac1));
    y[44] = 1 - proportion_countryA_vac1;
    y[45] = number_vac1*proportion_countryA_vac1/(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0);
    y[46] = number_vac1*proportion_countryB_vac1/(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0);
    y[47] = number_vac2*proportion_countryA_vac2/(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0);
    y[48] = number_vac2*proportion_countryB_vac2/(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0);
    y[49] = (1.0 - 1.0*amici::heaviside(amici_t - 0.20000000000000001))*(5.0*amici_t*(5.0*amici_t*(5.0*amici_t*(0.5*y1 - 0.50505050505050497*y2 + 0.0050505050505050501*y3) - 0.5*y1 + 0.50505050505050497*y2 - 0.0050505050505050501*y3) - y1 + y2) + y1) + 1.0*(y2 + (0.050505050505050504*amici_t - 0.010101010101010102)*(-49.5*y1 + 49*y2 + 0.5*y3 + (0.050505050505050504*amici_t - 0.010101010101010102)*(99*y1 - 100*y2 + y3 + (0.050505050505050504*amici_t - 0.010101010101010102)*(-49.5*y1 + 50*y2 - 0.5*y3))))*amici::heaviside(amici_t - 0.20000000000000001);
    y[50] = 1.0;
    y[51] = 1.0;
}

} // namespace model_vaccination
} // namespace amici
