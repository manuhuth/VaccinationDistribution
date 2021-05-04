#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "one_country_x_rdata.h"

namespace amici {
namespace model_one_country {

void x_solver_one_country(realtype *x_solver, const realtype *x_rdata){
    x_solver[0] = susceptible_countryA_vac0;
    x_solver[1] = susceptible_countryA_vac1;
    x_solver[2] = susceptible_countryA_vac2;
    x_solver[3] = infectious_countryA_vac0_virW;
    x_solver[4] = infectious_countryA_vac0_virM;
    x_solver[5] = infectious_countryA_vac1_virW;
    x_solver[6] = infectious_countryA_vac1_virM;
    x_solver[7] = infectious_countryA_vac2_virW;
    x_solver[8] = infectious_countryA_vac2_virM;
    x_solver[9] = recovered_countryA_vac0_virW;
    x_solver[10] = recovered_countryA_vac0_virM;
    x_solver[11] = recovered_countryA_vac1_virW;
    x_solver[12] = recovered_countryA_vac1_virM;
    x_solver[13] = recovered_countryA_vac2_virW;
    x_solver[14] = recovered_countryA_vac2_virM;
    x_solver[15] = dead_countryA_vac0_virW;
    x_solver[16] = dead_countryA_vac0_virM;
    x_solver[17] = dead_countryA_vac1_virW;
    x_solver[18] = dead_countryA_vac1_virM;
    x_solver[19] = dead_countryA_vac2_virW;
    x_solver[20] = dead_countryA_vac2_virM;
}

} // namespace model_one_country
} // namespace amici
