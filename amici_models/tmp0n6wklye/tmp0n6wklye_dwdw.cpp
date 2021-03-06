#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "tmp0n6wklye_x.h"
#include "tmp0n6wklye_p.h"
#include "tmp0n6wklye_k.h"
#include "tmp0n6wklye_h.h"
#include "tmp0n6wklye_w.h"
#include "tmp0n6wklye_dwdw.h"

namespace amici {
namespace model_tmp0n6wklye {

void dwdw_tmp0n6wklye(realtype *dwdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dproportion_countryA_vac1_dspline_countryA_vac1 = std::exp(-spline_countryA_vac1)/std::pow(1 + std::exp(-spline_countryA_vac1), 2);  // dwdw[0]
    dnu_countryA_vac2_dnumber_vac2 = proportion_countryA_vac2/(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0);  // dwdw[1]
    dnu_countryB_vac2_dnumber_vac2 = proportion_countryB_vac2/(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0);  // dwdw[2]
    dnu_countryA_vac1_dnumber_vac1 = proportion_countryA_vac1/(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0);  // dwdw[3]
    dnu_countryB_vac1_dnumber_vac1 = proportion_countryB_vac1/(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0);  // dwdw[4]
    dproportion_countryA_vac2_dspline_countryA_vac2 = std::exp(-spline_countryA_vac2)/std::pow(1 + std::exp(-spline_countryA_vac2), 2);  // dwdw[5]
    dproportion_countryB_vac1_dproportion_countryA_vac1 = -1;  // dwdw[6]
    dnu_countryA_vac1_dproportion_countryA_vac1 = number_vac1/(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0);  // dwdw[7]
    dproportion_countryB_vac2_dproportion_countryA_vac2 = -1;  // dwdw[8]
    dnu_countryA_vac2_dproportion_countryA_vac2 = number_vac2/(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0);  // dwdw[9]
    dnu_countryB_vac2_dproportion_countryB_vac2 = number_vac2/(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0);  // dwdw[10]
    dnu_countryB_vac1_dproportion_countryB_vac1 = number_vac1/(recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + 2*susceptible_countryB_vac0);  // dwdw[11]
    dflux_r96_dnu_countryA_vac1 = susceptible_countryA_vac0;  // dwdw[12]
    dflux_r100_dnu_countryA_vac1 = recovered_countryA_vac0_virW;  // dwdw[13]
    dflux_r102_dnu_countryA_vac1 = recovered_countryA_vac0_virM;  // dwdw[14]
    dflux_r97_dnu_countryA_vac2 = susceptible_countryA_vac0;  // dwdw[15]
    dflux_r101_dnu_countryA_vac2 = recovered_countryA_vac0_virW;  // dwdw[16]
    dflux_r103_dnu_countryA_vac2 = recovered_countryA_vac0_virM;  // dwdw[17]
    dflux_r99_dnu_countryB_vac2 = susceptible_countryB_vac0;  // dwdw[18]
    dflux_r105_dnu_countryB_vac2 = recovered_countryB_vac0_virW;  // dwdw[19]
    dflux_r107_dnu_countryB_vac2 = recovered_countryB_vac0_virM;  // dwdw[20]
    dflux_r98_dnu_countryB_vac1 = susceptible_countryB_vac0;  // dwdw[21]
    dflux_r104_dnu_countryB_vac1 = recovered_countryB_vac0_virW;  // dwdw[22]
    dflux_r106_dnu_countryB_vac1 = recovered_countryB_vac0_virM;  // dwdw[23]
}

} // namespace model_tmp0n6wklye
} // namespace amici
