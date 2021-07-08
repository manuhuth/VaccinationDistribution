#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "p.h"
#include "k.h"

namespace amici {
namespace model_tmpif24ue8s {

void x0_fixedParameters_tmpif24ue8s(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k){
x0_fixedParameters[0] = susceptible_countryA_vac0_t0;
x0_fixedParameters[1] = susceptible_countryA_vac1_t0;
x0_fixedParameters[2] = susceptible_countryA_vac2_t0;
x0_fixedParameters[3] = susceptible_countryB_vac0_t0;
x0_fixedParameters[4] = susceptible_countryB_vac1_t0;
x0_fixedParameters[5] = susceptible_countryB_vac2_t0;
x0_fixedParameters[6] = infectious_countryA_vac0_virW_t0;
x0_fixedParameters[7] = infectious_countryA_vac0_virM_t0;
x0_fixedParameters[8] = infectious_countryA_vac1_virW_t0;
x0_fixedParameters[9] = infectious_countryA_vac1_virM_t0;
x0_fixedParameters[10] = infectious_countryA_vac2_virW_t0;
x0_fixedParameters[11] = infectious_countryA_vac2_virM_t0;
x0_fixedParameters[12] = infectious_countryB_vac0_virW_t0;
x0_fixedParameters[13] = infectious_countryB_vac0_virM_t0;
x0_fixedParameters[14] = infectious_countryB_vac1_virW_t0;
x0_fixedParameters[15] = infectious_countryB_vac1_virM_t0;
x0_fixedParameters[16] = infectious_countryB_vac2_virW_t0;
x0_fixedParameters[17] = infectious_countryB_vac2_virM_t0;
x0_fixedParameters[18] = recovered_countryA_vac0_virW_t0;
x0_fixedParameters[19] = recovered_countryA_vac0_virM_t0;
x0_fixedParameters[20] = recovered_countryA_vac1_virW_t0;
x0_fixedParameters[21] = recovered_countryA_vac1_virM_t0;
x0_fixedParameters[22] = recovered_countryA_vac2_virW_t0;
x0_fixedParameters[23] = recovered_countryA_vac2_virM_t0;
x0_fixedParameters[24] = recovered_countryB_vac0_virW_t0;
x0_fixedParameters[25] = recovered_countryB_vac0_virM_t0;
x0_fixedParameters[26] = recovered_countryB_vac1_virW_t0;
x0_fixedParameters[27] = recovered_countryB_vac1_virM_t0;
x0_fixedParameters[28] = recovered_countryB_vac2_virW_t0;
x0_fixedParameters[29] = recovered_countryB_vac2_virM_t0;
x0_fixedParameters[30] = dead_countryA_vac0_virW_t0;
x0_fixedParameters[31] = dead_countryA_vac0_virM_t0;
x0_fixedParameters[32] = dead_countryA_vac1_virW_t0;
x0_fixedParameters[33] = dead_countryA_vac1_virM_t0;
x0_fixedParameters[34] = dead_countryA_vac2_virW_t0;
x0_fixedParameters[35] = dead_countryA_vac2_virM_t0;
x0_fixedParameters[36] = dead_countryB_vac0_virW_t0;
x0_fixedParameters[37] = dead_countryB_vac0_virM_t0;
x0_fixedParameters[38] = dead_countryB_vac1_virW_t0;
x0_fixedParameters[39] = dead_countryB_vac1_virM_t0;
x0_fixedParameters[40] = dead_countryB_vac2_virW_t0;
x0_fixedParameters[41] = dead_countryB_vac2_virM_t0;
}

} // namespace amici
} // namespace model_tmpif24ue8s