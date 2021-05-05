#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "vaccination_x.h"

namespace amici {
namespace model_vaccination {

void x_rdata_vaccination(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = susceptible_countryA_vac0;
    x_rdata[1] = susceptible_countryA_vac1;
    x_rdata[2] = susceptible_countryA_vac2;
    x_rdata[3] = infectious_countryA_vac0_virW;
    x_rdata[4] = infectious_countryA_vac0_virM;
    x_rdata[5] = infectious_countryA_vac1_virW;
    x_rdata[6] = infectious_countryA_vac1_virM;
    x_rdata[7] = infectious_countryA_vac2_virW;
    x_rdata[8] = infectious_countryA_vac2_virM;
    x_rdata[9] = recovered_countryA_vac0_virW;
    x_rdata[10] = recovered_countryA_vac0_virM;
    x_rdata[11] = recovered_countryA_vac1_virW;
    x_rdata[12] = recovered_countryA_vac1_virM;
    x_rdata[13] = recovered_countryA_vac2_virW;
    x_rdata[14] = recovered_countryA_vac2_virM;
    x_rdata[15] = dead_countryA_vac0_virW;
    x_rdata[16] = dead_countryA_vac0_virM;
    x_rdata[17] = dead_countryA_vac1_virW;
    x_rdata[18] = dead_countryA_vac1_virM;
    x_rdata[19] = dead_countryA_vac2_virW;
    x_rdata[20] = dead_countryA_vac2_virM;
}

} // namespace model_vaccination
} // namespace amici
