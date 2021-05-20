#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "vaccination_x_rdata.h"

namespace amici {
namespace model_vaccination {

void x_solver_vaccination(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = susceptible_countryA_vac0;
    x_solver[1] = susceptible_countryA_vac1;
    x_solver[2] = susceptible_countryA_vac2;
    x_solver[3] = susceptible_countryB_vac0;
    x_solver[4] = susceptible_countryB_vac1;
    x_solver[5] = susceptible_countryB_vac2;
    x_solver[6] = infectious_countryA_vac0_virW;
    x_solver[7] = infectious_countryA_vac0_virM;
    x_solver[8] = infectious_countryA_vac1_virW;
    x_solver[9] = infectious_countryA_vac1_virM;
    x_solver[10] = infectious_countryA_vac2_virW;
    x_solver[11] = infectious_countryA_vac2_virM;
    x_solver[12] = infectious_countryB_vac0_virW;
    x_solver[13] = infectious_countryB_vac0_virM;
    x_solver[14] = infectious_countryB_vac1_virW;
    x_solver[15] = infectious_countryB_vac1_virM;
    x_solver[16] = infectious_countryB_vac2_virW;
    x_solver[17] = infectious_countryB_vac2_virM;
    x_solver[18] = recovered_countryA_vac0_virW;
    x_solver[19] = recovered_countryA_vac0_virM;
    x_solver[20] = recovered_countryA_vac1_virW;
    x_solver[21] = recovered_countryA_vac1_virM;
    x_solver[22] = recovered_countryA_vac2_virW;
    x_solver[23] = recovered_countryA_vac2_virM;
    x_solver[24] = recovered_countryB_vac0_virW;
    x_solver[25] = recovered_countryB_vac0_virM;
    x_solver[26] = recovered_countryB_vac1_virW;
    x_solver[27] = recovered_countryB_vac1_virM;
    x_solver[28] = recovered_countryB_vac2_virW;
    x_solver[29] = recovered_countryB_vac2_virM;
    x_solver[30] = dead_countryA_vac0_virW;
    x_solver[31] = dead_countryA_vac0_virM;
    x_solver[32] = dead_countryA_vac1_virW;
    x_solver[33] = dead_countryA_vac1_virM;
    x_solver[34] = dead_countryA_vac2_virW;
    x_solver[35] = dead_countryA_vac2_virM;
    x_solver[36] = dead_countryB_vac0_virW;
    x_solver[37] = dead_countryB_vac0_virM;
    x_solver[38] = dead_countryB_vac1_virW;
    x_solver[39] = dead_countryB_vac1_virM;
    x_solver[40] = dead_countryB_vac2_virW;
    x_solver[41] = dead_countryB_vac2_virM;
    x_solver[42] = amici_t;
}

} // namespace model_vaccination
} // namespace amici
