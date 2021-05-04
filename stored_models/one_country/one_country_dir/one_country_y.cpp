#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "one_country_x.h"
#include "one_country_p.h"
#include "one_country_w.h"

namespace amici {
namespace model_one_country {

void y_one_country(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
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
    y[42] = 1.0;
    y[43] = 1.0;
}

} // namespace model_one_country
} // namespace amici
