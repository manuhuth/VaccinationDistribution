#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "vaccination_x.h"
#include "vaccination_p.h"
#include "vaccination_w.h"
#include "vaccination_dwdw.h"

namespace amici {
namespace model_vaccination {

void dwdw_vaccination(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dflux_r30_dnu_countryA_vac1 = susceptible_countryA_vac0;  // dwdw[0]
    dflux_r32_dnu_countryA_vac1 = recovered_countryA_vac0_virW;  // dwdw[1]
    dflux_r34_dnu_countryA_vac1 = recovered_countryA_vac0_virM;  // dwdw[2]
    dflux_r31_dnu_countryA_vac2 = susceptible_countryA_vac0;  // dwdw[3]
    dflux_r33_dnu_countryA_vac2 = recovered_countryA_vac0_virW;  // dwdw[4]
    dflux_r35_dnu_countryA_vac2 = recovered_countryA_vac0_virM;  // dwdw[5]
}

} // namespace model_vaccination
} // namespace amici
