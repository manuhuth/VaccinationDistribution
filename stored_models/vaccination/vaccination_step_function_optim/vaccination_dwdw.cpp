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
    dflux_r99_dnu_countryB_vac2 = susceptible_countryB_vac0;  // dwdw[0]
    dflux_r105_dnu_countryB_vac2 = recovered_countryB_vac0_virW;  // dwdw[1]
    dflux_r107_dnu_countryB_vac2 = recovered_countryB_vac0_virM;  // dwdw[2]
    dflux_r96_dnu_countryA_vac1 = susceptible_countryA_vac0;  // dwdw[3]
    dflux_r100_dnu_countryA_vac1 = recovered_countryA_vac0_virW;  // dwdw[4]
    dflux_r102_dnu_countryA_vac1 = recovered_countryA_vac0_virM;  // dwdw[5]
    dflux_r98_dnu_countryB_vac1 = susceptible_countryB_vac0;  // dwdw[6]
    dflux_r104_dnu_countryB_vac1 = recovered_countryB_vac0_virW;  // dwdw[7]
    dflux_r106_dnu_countryB_vac1 = recovered_countryB_vac0_virM;  // dwdw[8]
    dflux_r97_dnu_countryA_vac2 = susceptible_countryA_vac0;  // dwdw[9]
    dflux_r101_dnu_countryA_vac2 = recovered_countryA_vac0_virW;  // dwdw[10]
    dflux_r103_dnu_countryA_vac2 = recovered_countryA_vac0_virM;  // dwdw[11]
}

} // namespace model_vaccination
} // namespace amici
