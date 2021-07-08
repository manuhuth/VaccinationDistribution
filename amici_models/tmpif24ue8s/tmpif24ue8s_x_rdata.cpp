#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "tcl.h"
#include "x.h"

namespace amici {
namespace model_tmpif24ue8s {

void x_rdata_tmpif24ue8s(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = susceptible_countryA_vac0;
    x_rdata[1] = susceptible_countryA_vac1;
    x_rdata[2] = susceptible_countryA_vac2;
    x_rdata[3] = susceptible_countryB_vac0;
    x_rdata[4] = susceptible_countryB_vac1;
    x_rdata[5] = susceptible_countryB_vac2;
    x_rdata[6] = infectious_countryA_vac0_virW;
    x_rdata[7] = infectious_countryA_vac0_virM;
    x_rdata[8] = infectious_countryA_vac1_virW;
    x_rdata[9] = infectious_countryA_vac1_virM;
    x_rdata[10] = infectious_countryA_vac2_virW;
    x_rdata[11] = infectious_countryA_vac2_virM;
    x_rdata[12] = infectious_countryB_vac0_virW;
    x_rdata[13] = infectious_countryB_vac0_virM;
    x_rdata[14] = infectious_countryB_vac1_virW;
    x_rdata[15] = infectious_countryB_vac1_virM;
    x_rdata[16] = infectious_countryB_vac2_virW;
    x_rdata[17] = infectious_countryB_vac2_virM;
    x_rdata[18] = recovered_countryA_vac0_virW;
    x_rdata[19] = recovered_countryA_vac0_virM;
    x_rdata[20] = recovered_countryA_vac1_virW;
    x_rdata[21] = recovered_countryA_vac1_virM;
    x_rdata[22] = recovered_countryA_vac2_virW;
    x_rdata[23] = recovered_countryA_vac2_virM;
    x_rdata[24] = recovered_countryB_vac0_virW;
    x_rdata[25] = recovered_countryB_vac0_virM;
    x_rdata[26] = recovered_countryB_vac1_virW;
    x_rdata[27] = recovered_countryB_vac1_virM;
    x_rdata[28] = recovered_countryB_vac2_virW;
    x_rdata[29] = recovered_countryB_vac2_virM;
    x_rdata[30] = dead_countryA_vac0_virW;
    x_rdata[31] = dead_countryA_vac0_virM;
    x_rdata[32] = dead_countryA_vac1_virW;
    x_rdata[33] = dead_countryA_vac1_virM;
    x_rdata[34] = dead_countryA_vac2_virW;
    x_rdata[35] = dead_countryA_vac2_virM;
    x_rdata[36] = dead_countryB_vac0_virW;
    x_rdata[37] = dead_countryB_vac0_virM;
    x_rdata[38] = dead_countryB_vac1_virW;
    x_rdata[39] = dead_countryB_vac1_virM;
    x_rdata[40] = dead_countryB_vac2_virW;
    x_rdata[41] = dead_countryB_vac2_virM;
}

} // namespace amici
} // namespace model_tmpif24ue8s