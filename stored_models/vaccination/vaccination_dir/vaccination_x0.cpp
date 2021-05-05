#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "vaccination_p.h"

namespace amici {
namespace model_vaccination {

void x0_vaccination(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = susceptible_countryA_vac0_t0;
    x0[1] = susceptible_countryA_vac1_t0;
    x0[2] = susceptible_countryA_vac2_t0;
    x0[3] = infectious_countryA_vac0_virW_t0;
    x0[4] = infectious_countryA_vac0_virM_t0;
    x0[5] = infectious_countryA_vac1_virW_t0;
    x0[6] = infectious_countryA_vac1_virM_t0;
    x0[7] = infectious_countryA_vac2_virW_t0;
    x0[8] = infectious_countryA_vac2_virM_t0;
    x0[9] = recovered_countryA_vac0_virW_t0;
    x0[10] = recovered_countryA_vac0_virM_t0;
    x0[11] = recovered_countryA_vac1_virW_t0;
    x0[12] = recovered_countryA_vac1_virM_t0;
    x0[13] = recovered_countryA_vac2_virW_t0;
    x0[14] = recovered_countryA_vac2_virM_t0;
    x0[15] = dead_countryA_vac0_virW_t0;
    x0[16] = dead_countryA_vac0_virM_t0;
    x0[17] = dead_countryA_vac1_virW_t0;
    x0[18] = dead_countryA_vac1_virM_t0;
    x0[19] = dead_countryA_vac2_virW_t0;
    x0[20] = dead_countryA_vac2_virM_t0;
}

} // namespace model_vaccination
} // namespace amici
