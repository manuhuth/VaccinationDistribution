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
    y[3] = infectious_countryA_vac0_virW;
    y[4] = infectious_countryA_vac0_virM;
    y[5] = infectious_countryA_vac1_virW;
    y[6] = infectious_countryA_vac1_virM;
    y[7] = infectious_countryA_vac2_virW;
    y[8] = infectious_countryA_vac2_virM;
    y[9] = recovered_countryA_vac0_virW;
    y[10] = recovered_countryA_vac0_virM;
    y[11] = recovered_countryA_vac1_virW;
    y[12] = recovered_countryA_vac1_virM;
    y[13] = recovered_countryA_vac2_virW;
    y[14] = recovered_countryA_vac2_virM;
    y[15] = dead_countryA_vac0_virW;
    y[16] = dead_countryA_vac0_virM;
    y[17] = dead_countryA_vac1_virW;
    y[18] = dead_countryA_vac1_virM;
    y[19] = dead_countryA_vac2_virW;
    y[20] = dead_countryA_vac2_virM;
    y[21] = 1.0;
}

} // namespace model_one_country
} // namespace amici