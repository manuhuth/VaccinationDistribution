#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "tmp2t76h1f9_p.h"
#include "tmp2t76h1f9_k.h"

namespace amici {
namespace model_tmp2t76h1f9 {

void x0_tmp2t76h1f9(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = susceptible_countryA_vac0_t0;
    x0[1] = susceptible_countryA_vac1_t0;
    x0[2] = susceptible_countryA_vac2_t0;
    x0[3] = susceptible_countryB_vac0_t0;
    x0[4] = susceptible_countryB_vac1_t0;
    x0[5] = susceptible_countryB_vac2_t0;
    x0[6] = infectious_countryA_vac0_virW_t0;
    x0[7] = infectious_countryA_vac0_virM_t0;
    x0[8] = infectious_countryA_vac1_virW_t0;
    x0[9] = infectious_countryA_vac1_virM_t0;
    x0[10] = infectious_countryA_vac2_virW_t0;
    x0[11] = infectious_countryA_vac2_virM_t0;
    x0[12] = infectious_countryB_vac0_virW_t0;
    x0[13] = infectious_countryB_vac0_virM_t0;
    x0[14] = infectious_countryB_vac1_virW_t0;
    x0[15] = infectious_countryB_vac1_virM_t0;
    x0[16] = infectious_countryB_vac2_virW_t0;
    x0[17] = infectious_countryB_vac2_virM_t0;
    x0[18] = recovered_countryA_vac0_virW_t0;
    x0[19] = recovered_countryA_vac0_virM_t0;
    x0[20] = recovered_countryA_vac1_virW_t0;
    x0[21] = recovered_countryA_vac1_virM_t0;
    x0[22] = recovered_countryA_vac2_virW_t0;
    x0[23] = recovered_countryA_vac2_virM_t0;
    x0[24] = recovered_countryB_vac0_virW_t0;
    x0[25] = recovered_countryB_vac0_virM_t0;
    x0[26] = recovered_countryB_vac1_virW_t0;
    x0[27] = recovered_countryB_vac1_virM_t0;
    x0[28] = recovered_countryB_vac2_virW_t0;
    x0[29] = recovered_countryB_vac2_virM_t0;
    x0[30] = dead_countryA_vac0_virW_t0;
    x0[31] = dead_countryA_vac0_virM_t0;
    x0[32] = dead_countryA_vac1_virW_t0;
    x0[33] = dead_countryA_vac1_virM_t0;
    x0[34] = dead_countryA_vac2_virW_t0;
    x0[35] = dead_countryA_vac2_virM_t0;
    x0[36] = dead_countryB_vac0_virW_t0;
    x0[37] = dead_countryB_vac0_virM_t0;
    x0[38] = dead_countryB_vac1_virW_t0;
    x0[39] = dead_countryB_vac1_virM_t0;
    x0[40] = dead_countryB_vac2_virW_t0;
    x0[41] = dead_countryB_vac2_virM_t0;
}

} // namespace model_tmp2t76h1f9
} // namespace amici
