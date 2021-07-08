#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "tmp0n6wklye_p.h"
#include "tmp0n6wklye_k.h"

namespace amici {
namespace model_tmp0n6wklye {

void x0_fixedParameters_tmp0n6wklye(realtype *x0_fixedParameters, const realtype t, const realtype *p, const realtype *k, gsl::span<const int> reinitialization_state_idxs){
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 0) != reinitialization_state_idxs.cend())
        x0_fixedParameters[0] = susceptible_countryA_vac0_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 1) != reinitialization_state_idxs.cend())
        x0_fixedParameters[1] = susceptible_countryA_vac1_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 2) != reinitialization_state_idxs.cend())
        x0_fixedParameters[2] = susceptible_countryA_vac2_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 3) != reinitialization_state_idxs.cend())
        x0_fixedParameters[3] = susceptible_countryB_vac0_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 4) != reinitialization_state_idxs.cend())
        x0_fixedParameters[4] = susceptible_countryB_vac1_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 5) != reinitialization_state_idxs.cend())
        x0_fixedParameters[5] = susceptible_countryB_vac2_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 6) != reinitialization_state_idxs.cend())
        x0_fixedParameters[6] = infectious_countryA_vac0_virW_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 7) != reinitialization_state_idxs.cend())
        x0_fixedParameters[7] = infectious_countryA_vac0_virM_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 8) != reinitialization_state_idxs.cend())
        x0_fixedParameters[8] = infectious_countryA_vac1_virW_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 9) != reinitialization_state_idxs.cend())
        x0_fixedParameters[9] = infectious_countryA_vac1_virM_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 10) != reinitialization_state_idxs.cend())
        x0_fixedParameters[10] = infectious_countryA_vac2_virW_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 11) != reinitialization_state_idxs.cend())
        x0_fixedParameters[11] = infectious_countryA_vac2_virM_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 12) != reinitialization_state_idxs.cend())
        x0_fixedParameters[12] = infectious_countryB_vac0_virW_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 13) != reinitialization_state_idxs.cend())
        x0_fixedParameters[13] = infectious_countryB_vac0_virM_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 14) != reinitialization_state_idxs.cend())
        x0_fixedParameters[14] = infectious_countryB_vac1_virW_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 15) != reinitialization_state_idxs.cend())
        x0_fixedParameters[15] = infectious_countryB_vac1_virM_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 16) != reinitialization_state_idxs.cend())
        x0_fixedParameters[16] = infectious_countryB_vac2_virW_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 17) != reinitialization_state_idxs.cend())
        x0_fixedParameters[17] = infectious_countryB_vac2_virM_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 18) != reinitialization_state_idxs.cend())
        x0_fixedParameters[18] = recovered_countryA_vac0_virW_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 19) != reinitialization_state_idxs.cend())
        x0_fixedParameters[19] = recovered_countryA_vac0_virM_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 20) != reinitialization_state_idxs.cend())
        x0_fixedParameters[20] = recovered_countryA_vac1_virW_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 21) != reinitialization_state_idxs.cend())
        x0_fixedParameters[21] = recovered_countryA_vac1_virM_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 22) != reinitialization_state_idxs.cend())
        x0_fixedParameters[22] = recovered_countryA_vac2_virW_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 23) != reinitialization_state_idxs.cend())
        x0_fixedParameters[23] = recovered_countryA_vac2_virM_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 24) != reinitialization_state_idxs.cend())
        x0_fixedParameters[24] = recovered_countryB_vac0_virW_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 25) != reinitialization_state_idxs.cend())
        x0_fixedParameters[25] = recovered_countryB_vac0_virM_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 26) != reinitialization_state_idxs.cend())
        x0_fixedParameters[26] = recovered_countryB_vac1_virW_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 27) != reinitialization_state_idxs.cend())
        x0_fixedParameters[27] = recovered_countryB_vac1_virM_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 28) != reinitialization_state_idxs.cend())
        x0_fixedParameters[28] = recovered_countryB_vac2_virW_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 29) != reinitialization_state_idxs.cend())
        x0_fixedParameters[29] = recovered_countryB_vac2_virM_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 30) != reinitialization_state_idxs.cend())
        x0_fixedParameters[30] = dead_countryA_vac0_virW_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 31) != reinitialization_state_idxs.cend())
        x0_fixedParameters[31] = dead_countryA_vac0_virM_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 32) != reinitialization_state_idxs.cend())
        x0_fixedParameters[32] = dead_countryA_vac1_virW_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 33) != reinitialization_state_idxs.cend())
        x0_fixedParameters[33] = dead_countryA_vac1_virM_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 34) != reinitialization_state_idxs.cend())
        x0_fixedParameters[34] = dead_countryA_vac2_virW_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 35) != reinitialization_state_idxs.cend())
        x0_fixedParameters[35] = dead_countryA_vac2_virM_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 36) != reinitialization_state_idxs.cend())
        x0_fixedParameters[36] = dead_countryB_vac0_virW_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 37) != reinitialization_state_idxs.cend())
        x0_fixedParameters[37] = dead_countryB_vac0_virM_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 38) != reinitialization_state_idxs.cend())
        x0_fixedParameters[38] = dead_countryB_vac1_virW_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 39) != reinitialization_state_idxs.cend())
        x0_fixedParameters[39] = dead_countryB_vac1_virM_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 40) != reinitialization_state_idxs.cend())
        x0_fixedParameters[40] = dead_countryB_vac2_virW_t0;
    if(std::find(reinitialization_state_idxs.cbegin(), reinitialization_state_idxs.cend(), 41) != reinitialization_state_idxs.cend())
        x0_fixedParameters[41] = dead_countryB_vac2_virM_t0;
}

} // namespace model_tmp0n6wklye
} // namespace amici
