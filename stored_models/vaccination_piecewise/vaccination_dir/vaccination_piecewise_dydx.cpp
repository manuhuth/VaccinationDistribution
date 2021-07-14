#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "p.h"
#include "k.h"
#include "w.h"
#include "x.h"
#include "dwdx.h"

namespace amici {
namespace model_vaccination_piecewise {

void dydx_vaccination_piecewise(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    dydx[0] = -((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac1_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac1_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac1_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac1_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac1_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac1_xx5
)
: ((t >= xx6 && t < xx7) ? (
   vaccine_supply_par_vac1_xx6
)
: ((t >= xx7 && t < xx8) ? (
   vaccine_supply_par_vac1_xx7
)
: ((t >= xx8 && t < xx9) ? (
   vaccine_supply_par_vac1_xx8
)
: (
   vaccine_supply_par_vac1_xx9
))))))))))*((t < 0) ? (
   0
)
: ((t < 14) ? (
   proportion_par_countryA_vac1_0
)
: ((t < 28) ? (
   proportion_par_countryA_vac1_14
)
: ((t < 42) ? (
   proportion_par_countryA_vac1_28
)
: ((t < 56) ? (
   proportion_par_countryA_vac1_42
)
: ((t < 70) ? (
   proportion_par_countryA_vac1_56
)
: ((t < 84) ? (
   proportion_par_countryA_vac1_70
)
: ((t < 98) ? (
   proportion_par_countryA_vac1_84
)
: ((t < 112) ? (
   proportion_par_countryA_vac1_98
)
: ((t < 126) ? (
   proportion_par_countryA_vac1_112
)
: ((t < 140) ? (
   proportion_par_countryA_vac1_126
)
: (
   0
))))))))))))/std::pow(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0, 2);
    dydx[2] = -((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac2_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac2_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac2_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac2_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac2_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac2_xx5
)
: ((t >= xx6 && t < xx7) ? (
   vaccine_supply_par_vac2_xx6
)
: ((t >= xx7 && t < xx8) ? (
   vaccine_supply_par_vac2_xx7
)
: ((t >= xx8 && t < xx9) ? (
   vaccine_supply_par_vac2_xx8
)
: (
   vaccine_supply_par_vac2_xx9
))))))))))*((t < 0) ? (
   0
)
: ((t < 14) ? (
   proportion_par_countryA_vac2_0
)
: ((t < 28) ? (
   proportion_par_countryA_vac2_14
)
: ((t < 42) ? (
   proportion_par_countryA_vac2_28
)
: ((t < 56) ? (
   proportion_par_countryA_vac2_42
)
: ((t < 70) ? (
   proportion_par_countryA_vac2_56
)
: ((t < 84) ? (
   proportion_par_countryA_vac2_70
)
: ((t < 98) ? (
   proportion_par_countryA_vac2_84
)
: ((t < 112) ? (
   proportion_par_countryA_vac2_98
)
: ((t < 126) ? (
   proportion_par_countryA_vac2_112
)
: ((t < 140) ? (
   proportion_par_countryA_vac2_126
)
: (
   0
))))))))))))/std::pow(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0, 2);
    dydx[37] = -((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac1_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac1_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac1_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac1_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac1_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac1_xx5
)
: ((t >= xx6 && t < xx7) ? (
   vaccine_supply_par_vac1_xx6
)
: ((t >= xx7 && t < xx8) ? (
   vaccine_supply_par_vac1_xx7
)
: ((t >= xx8 && t < xx9) ? (
   vaccine_supply_par_vac1_xx8
)
: (
   vaccine_supply_par_vac1_xx9
))))))))))*((t < 0) ? (
   0
)
: ((t < 14) ? (
   1 - proportion_par_countryA_vac1_0
)
: ((t < 28) ? (
   1 - proportion_par_countryA_vac1_14
)
: ((t < 42) ? (
   1 - proportion_par_countryA_vac1_28
)
: ((t < 56) ? (
   1 - proportion_par_countryA_vac1_42
)
: ((t < 70) ? (
   1 - proportion_par_countryA_vac1_56
)
: ((t < 84) ? (
   1 - proportion_par_countryA_vac1_70
)
: ((t < 98) ? (
   1 - proportion_par_countryA_vac1_84
)
: ((t < 112) ? (
   1 - proportion_par_countryA_vac1_98
)
: ((t < 126) ? (
   1 - proportion_par_countryA_vac1_112
)
: ((t < 140) ? (
   1 - proportion_par_countryA_vac1_126
)
: (
   0
))))))))))))/std::pow(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0, 2);
    dydx[39] = -((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac2_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac2_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac2_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac2_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac2_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac2_xx5
)
: ((t >= xx6 && t < xx7) ? (
   vaccine_supply_par_vac2_xx6
)
: ((t >= xx7 && t < xx8) ? (
   vaccine_supply_par_vac2_xx7
)
: ((t >= xx8 && t < xx9) ? (
   vaccine_supply_par_vac2_xx8
)
: (
   vaccine_supply_par_vac2_xx9
))))))))))*((t < 0) ? (
   0
)
: ((t < 14) ? (
   1 - proportion_par_countryA_vac2_0
)
: ((t < 28) ? (
   1 - proportion_par_countryA_vac2_14
)
: ((t < 42) ? (
   1 - proportion_par_countryA_vac2_28
)
: ((t < 56) ? (
   1 - proportion_par_countryA_vac2_42
)
: ((t < 70) ? (
   1 - proportion_par_countryA_vac2_56
)
: ((t < 84) ? (
   1 - proportion_par_countryA_vac2_70
)
: ((t < 98) ? (
   1 - proportion_par_countryA_vac2_84
)
: ((t < 112) ? (
   1 - proportion_par_countryA_vac2_98
)
: ((t < 126) ? (
   1 - proportion_par_countryA_vac2_112
)
: ((t < 140) ? (
   1 - proportion_par_countryA_vac2_126
)
: (
   0
))))))))))))/std::pow(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0, 2);
    dydx[72] = -((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac1_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac1_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac1_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac1_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac1_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac1_xx5
)
: ((t >= xx6 && t < xx7) ? (
   vaccine_supply_par_vac1_xx6
)
: ((t >= xx7 && t < xx8) ? (
   vaccine_supply_par_vac1_xx7
)
: ((t >= xx8 && t < xx9) ? (
   vaccine_supply_par_vac1_xx8
)
: (
   vaccine_supply_par_vac1_xx9
))))))))))*((t < 0) ? (
   0
)
: ((t < 14) ? (
   proportion_par_countryA_vac1_0
)
: ((t < 28) ? (
   proportion_par_countryA_vac1_14
)
: ((t < 42) ? (
   proportion_par_countryA_vac1_28
)
: ((t < 56) ? (
   proportion_par_countryA_vac1_42
)
: ((t < 70) ? (
   proportion_par_countryA_vac1_56
)
: ((t < 84) ? (
   proportion_par_countryA_vac1_70
)
: ((t < 98) ? (
   proportion_par_countryA_vac1_84
)
: ((t < 112) ? (
   proportion_par_countryA_vac1_98
)
: ((t < 126) ? (
   proportion_par_countryA_vac1_112
)
: ((t < 140) ? (
   proportion_par_countryA_vac1_126
)
: (
   0
))))))))))))/std::pow(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0, 2);
    dydx[74] = -((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac2_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac2_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac2_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac2_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac2_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac2_xx5
)
: ((t >= xx6 && t < xx7) ? (
   vaccine_supply_par_vac2_xx6
)
: ((t >= xx7 && t < xx8) ? (
   vaccine_supply_par_vac2_xx7
)
: ((t >= xx8 && t < xx9) ? (
   vaccine_supply_par_vac2_xx8
)
: (
   vaccine_supply_par_vac2_xx9
))))))))))*((t < 0) ? (
   0
)
: ((t < 14) ? (
   proportion_par_countryA_vac2_0
)
: ((t < 28) ? (
   proportion_par_countryA_vac2_14
)
: ((t < 42) ? (
   proportion_par_countryA_vac2_28
)
: ((t < 56) ? (
   proportion_par_countryA_vac2_42
)
: ((t < 70) ? (
   proportion_par_countryA_vac2_56
)
: ((t < 84) ? (
   proportion_par_countryA_vac2_70
)
: ((t < 98) ? (
   proportion_par_countryA_vac2_84
)
: ((t < 112) ? (
   proportion_par_countryA_vac2_98
)
: ((t < 126) ? (
   proportion_par_countryA_vac2_112
)
: ((t < 140) ? (
   proportion_par_countryA_vac2_126
)
: (
   0
))))))))))))/std::pow(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0, 2);
    dydx[84] = -((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac1_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac1_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac1_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac1_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac1_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac1_xx5
)
: ((t >= xx6 && t < xx7) ? (
   vaccine_supply_par_vac1_xx6
)
: ((t >= xx7 && t < xx8) ? (
   vaccine_supply_par_vac1_xx7
)
: ((t >= xx8 && t < xx9) ? (
   vaccine_supply_par_vac1_xx8
)
: (
   vaccine_supply_par_vac1_xx9
))))))))))*((t < 0) ? (
   0
)
: ((t < 14) ? (
   proportion_par_countryA_vac1_0
)
: ((t < 28) ? (
   proportion_par_countryA_vac1_14
)
: ((t < 42) ? (
   proportion_par_countryA_vac1_28
)
: ((t < 56) ? (
   proportion_par_countryA_vac1_42
)
: ((t < 70) ? (
   proportion_par_countryA_vac1_56
)
: ((t < 84) ? (
   proportion_par_countryA_vac1_70
)
: ((t < 98) ? (
   proportion_par_countryA_vac1_84
)
: ((t < 112) ? (
   proportion_par_countryA_vac1_98
)
: ((t < 126) ? (
   proportion_par_countryA_vac1_112
)
: ((t < 140) ? (
   proportion_par_countryA_vac1_126
)
: (
   0
))))))))))))/std::pow(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0, 2);
    dydx[86] = -((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac2_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac2_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac2_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac2_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac2_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac2_xx5
)
: ((t >= xx6 && t < xx7) ? (
   vaccine_supply_par_vac2_xx6
)
: ((t >= xx7 && t < xx8) ? (
   vaccine_supply_par_vac2_xx7
)
: ((t >= xx8 && t < xx9) ? (
   vaccine_supply_par_vac2_xx8
)
: (
   vaccine_supply_par_vac2_xx9
))))))))))*((t < 0) ? (
   0
)
: ((t < 14) ? (
   proportion_par_countryA_vac2_0
)
: ((t < 28) ? (
   proportion_par_countryA_vac2_14
)
: ((t < 42) ? (
   proportion_par_countryA_vac2_28
)
: ((t < 56) ? (
   proportion_par_countryA_vac2_42
)
: ((t < 70) ? (
   proportion_par_countryA_vac2_56
)
: ((t < 84) ? (
   proportion_par_countryA_vac2_70
)
: ((t < 98) ? (
   proportion_par_countryA_vac2_84
)
: ((t < 112) ? (
   proportion_par_countryA_vac2_98
)
: ((t < 126) ? (
   proportion_par_countryA_vac2_112
)
: ((t < 140) ? (
   proportion_par_countryA_vac2_126
)
: (
   0
))))))))))))/std::pow(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0, 2);
    dydx[145] = -((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac1_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac1_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac1_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac1_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac1_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac1_xx5
)
: ((t >= xx6 && t < xx7) ? (
   vaccine_supply_par_vac1_xx6
)
: ((t >= xx7 && t < xx8) ? (
   vaccine_supply_par_vac1_xx7
)
: ((t >= xx8 && t < xx9) ? (
   vaccine_supply_par_vac1_xx8
)
: (
   vaccine_supply_par_vac1_xx9
))))))))))*((t < 0) ? (
   0
)
: ((t < 14) ? (
   1 - proportion_par_countryA_vac1_0
)
: ((t < 28) ? (
   1 - proportion_par_countryA_vac1_14
)
: ((t < 42) ? (
   1 - proportion_par_countryA_vac1_28
)
: ((t < 56) ? (
   1 - proportion_par_countryA_vac1_42
)
: ((t < 70) ? (
   1 - proportion_par_countryA_vac1_56
)
: ((t < 84) ? (
   1 - proportion_par_countryA_vac1_70
)
: ((t < 98) ? (
   1 - proportion_par_countryA_vac1_84
)
: ((t < 112) ? (
   1 - proportion_par_countryA_vac1_98
)
: ((t < 126) ? (
   1 - proportion_par_countryA_vac1_112
)
: ((t < 140) ? (
   1 - proportion_par_countryA_vac1_126
)
: (
   0
))))))))))))/std::pow(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0, 2);
    dydx[147] = -((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac2_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac2_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac2_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac2_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac2_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac2_xx5
)
: ((t >= xx6 && t < xx7) ? (
   vaccine_supply_par_vac2_xx6
)
: ((t >= xx7 && t < xx8) ? (
   vaccine_supply_par_vac2_xx7
)
: ((t >= xx8 && t < xx9) ? (
   vaccine_supply_par_vac2_xx8
)
: (
   vaccine_supply_par_vac2_xx9
))))))))))*((t < 0) ? (
   0
)
: ((t < 14) ? (
   1 - proportion_par_countryA_vac2_0
)
: ((t < 28) ? (
   1 - proportion_par_countryA_vac2_14
)
: ((t < 42) ? (
   1 - proportion_par_countryA_vac2_28
)
: ((t < 56) ? (
   1 - proportion_par_countryA_vac2_42
)
: ((t < 70) ? (
   1 - proportion_par_countryA_vac2_56
)
: ((t < 84) ? (
   1 - proportion_par_countryA_vac2_70
)
: ((t < 98) ? (
   1 - proportion_par_countryA_vac2_84
)
: ((t < 112) ? (
   1 - proportion_par_countryA_vac2_98
)
: ((t < 126) ? (
   1 - proportion_par_countryA_vac2_112
)
: ((t < 140) ? (
   1 - proportion_par_countryA_vac2_126
)
: (
   0
))))))))))))/std::pow(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0, 2);
    dydx[157] = -((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac1_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac1_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac1_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac1_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac1_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac1_xx5
)
: ((t >= xx6 && t < xx7) ? (
   vaccine_supply_par_vac1_xx6
)
: ((t >= xx7 && t < xx8) ? (
   vaccine_supply_par_vac1_xx7
)
: ((t >= xx8 && t < xx9) ? (
   vaccine_supply_par_vac1_xx8
)
: (
   vaccine_supply_par_vac1_xx9
))))))))))*((t < 0) ? (
   0
)
: ((t < 14) ? (
   1 - proportion_par_countryA_vac1_0
)
: ((t < 28) ? (
   1 - proportion_par_countryA_vac1_14
)
: ((t < 42) ? (
   1 - proportion_par_countryA_vac1_28
)
: ((t < 56) ? (
   1 - proportion_par_countryA_vac1_42
)
: ((t < 70) ? (
   1 - proportion_par_countryA_vac1_56
)
: ((t < 84) ? (
   1 - proportion_par_countryA_vac1_70
)
: ((t < 98) ? (
   1 - proportion_par_countryA_vac1_84
)
: ((t < 112) ? (
   1 - proportion_par_countryA_vac1_98
)
: ((t < 126) ? (
   1 - proportion_par_countryA_vac1_112
)
: ((t < 140) ? (
   1 - proportion_par_countryA_vac1_126
)
: (
   0
))))))))))))/std::pow(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0, 2);
    dydx[159] = -((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac2_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac2_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac2_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac2_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac2_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac2_xx5
)
: ((t >= xx6 && t < xx7) ? (
   vaccine_supply_par_vac2_xx6
)
: ((t >= xx7 && t < xx8) ? (
   vaccine_supply_par_vac2_xx7
)
: ((t >= xx8 && t < xx9) ? (
   vaccine_supply_par_vac2_xx8
)
: (
   vaccine_supply_par_vac2_xx9
))))))))))*((t < 0) ? (
   0
)
: ((t < 14) ? (
   1 - proportion_par_countryA_vac2_0
)
: ((t < 28) ? (
   1 - proportion_par_countryA_vac2_14
)
: ((t < 42) ? (
   1 - proportion_par_countryA_vac2_28
)
: ((t < 56) ? (
   1 - proportion_par_countryA_vac2_42
)
: ((t < 70) ? (
   1 - proportion_par_countryA_vac2_56
)
: ((t < 84) ? (
   1 - proportion_par_countryA_vac2_70
)
: ((t < 98) ? (
   1 - proportion_par_countryA_vac2_84
)
: ((t < 112) ? (
   1 - proportion_par_countryA_vac2_98
)
: ((t < 126) ? (
   1 - proportion_par_countryA_vac2_112
)
: ((t < 140) ? (
   1 - proportion_par_countryA_vac2_126
)
: (
   0
))))))))))))/std::pow(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0, 2);
    dydx[216] = -((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac1_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac1_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac1_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac1_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac1_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac1_xx5
)
: ((t >= xx6 && t < xx7) ? (
   vaccine_supply_par_vac1_xx6
)
: ((t >= xx7 && t < xx8) ? (
   vaccine_supply_par_vac1_xx7
)
: ((t >= xx8 && t < xx9) ? (
   vaccine_supply_par_vac1_xx8
)
: (
   vaccine_supply_par_vac1_xx9
))))))))))*((t < 0) ? (
   0
)
: ((t < 14) ? (
   proportion_par_countryA_vac1_0
)
: ((t < 28) ? (
   proportion_par_countryA_vac1_14
)
: ((t < 42) ? (
   proportion_par_countryA_vac1_28
)
: ((t < 56) ? (
   proportion_par_countryA_vac1_42
)
: ((t < 70) ? (
   proportion_par_countryA_vac1_56
)
: ((t < 84) ? (
   proportion_par_countryA_vac1_70
)
: ((t < 98) ? (
   proportion_par_countryA_vac1_84
)
: ((t < 112) ? (
   proportion_par_countryA_vac1_98
)
: ((t < 126) ? (
   proportion_par_countryA_vac1_112
)
: ((t < 140) ? (
   proportion_par_countryA_vac1_126
)
: (
   0
))))))))))))/std::pow(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0, 2);
    dydx[218] = -((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac2_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac2_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac2_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac2_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac2_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac2_xx5
)
: ((t >= xx6 && t < xx7) ? (
   vaccine_supply_par_vac2_xx6
)
: ((t >= xx7 && t < xx8) ? (
   vaccine_supply_par_vac2_xx7
)
: ((t >= xx8 && t < xx9) ? (
   vaccine_supply_par_vac2_xx8
)
: (
   vaccine_supply_par_vac2_xx9
))))))))))*((t < 0) ? (
   0
)
: ((t < 14) ? (
   proportion_par_countryA_vac2_0
)
: ((t < 28) ? (
   proportion_par_countryA_vac2_14
)
: ((t < 42) ? (
   proportion_par_countryA_vac2_28
)
: ((t < 56) ? (
   proportion_par_countryA_vac2_42
)
: ((t < 70) ? (
   proportion_par_countryA_vac2_56
)
: ((t < 84) ? (
   proportion_par_countryA_vac2_70
)
: ((t < 98) ? (
   proportion_par_countryA_vac2_84
)
: ((t < 112) ? (
   proportion_par_countryA_vac2_98
)
: ((t < 126) ? (
   proportion_par_countryA_vac2_112
)
: ((t < 140) ? (
   proportion_par_countryA_vac2_126
)
: (
   0
))))))))))))/std::pow(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0, 2);
    dydx[228] = -((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac1_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac1_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac1_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac1_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac1_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac1_xx5
)
: ((t >= xx6 && t < xx7) ? (
   vaccine_supply_par_vac1_xx6
)
: ((t >= xx7 && t < xx8) ? (
   vaccine_supply_par_vac1_xx7
)
: ((t >= xx8 && t < xx9) ? (
   vaccine_supply_par_vac1_xx8
)
: (
   vaccine_supply_par_vac1_xx9
))))))))))*((t < 0) ? (
   0
)
: ((t < 14) ? (
   proportion_par_countryA_vac1_0
)
: ((t < 28) ? (
   proportion_par_countryA_vac1_14
)
: ((t < 42) ? (
   proportion_par_countryA_vac1_28
)
: ((t < 56) ? (
   proportion_par_countryA_vac1_42
)
: ((t < 70) ? (
   proportion_par_countryA_vac1_56
)
: ((t < 84) ? (
   proportion_par_countryA_vac1_70
)
: ((t < 98) ? (
   proportion_par_countryA_vac1_84
)
: ((t < 112) ? (
   proportion_par_countryA_vac1_98
)
: ((t < 126) ? (
   proportion_par_countryA_vac1_112
)
: ((t < 140) ? (
   proportion_par_countryA_vac1_126
)
: (
   0
))))))))))))/std::pow(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0, 2);
    dydx[230] = -((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac2_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac2_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac2_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac2_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac2_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac2_xx5
)
: ((t >= xx6 && t < xx7) ? (
   vaccine_supply_par_vac2_xx6
)
: ((t >= xx7 && t < xx8) ? (
   vaccine_supply_par_vac2_xx7
)
: ((t >= xx8 && t < xx9) ? (
   vaccine_supply_par_vac2_xx8
)
: (
   vaccine_supply_par_vac2_xx9
))))))))))*((t < 0) ? (
   0
)
: ((t < 14) ? (
   proportion_par_countryA_vac2_0
)
: ((t < 28) ? (
   proportion_par_countryA_vac2_14
)
: ((t < 42) ? (
   proportion_par_countryA_vac2_28
)
: ((t < 56) ? (
   proportion_par_countryA_vac2_42
)
: ((t < 70) ? (
   proportion_par_countryA_vac2_56
)
: ((t < 84) ? (
   proportion_par_countryA_vac2_70
)
: ((t < 98) ? (
   proportion_par_countryA_vac2_84
)
: ((t < 112) ? (
   proportion_par_countryA_vac2_98
)
: ((t < 126) ? (
   proportion_par_countryA_vac2_112
)
: ((t < 140) ? (
   proportion_par_countryA_vac2_126
)
: (
   0
))))))))))))/std::pow(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0, 2);
    dydx[289] = -((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac1_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac1_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac1_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac1_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac1_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac1_xx5
)
: ((t >= xx6 && t < xx7) ? (
   vaccine_supply_par_vac1_xx6
)
: ((t >= xx7 && t < xx8) ? (
   vaccine_supply_par_vac1_xx7
)
: ((t >= xx8 && t < xx9) ? (
   vaccine_supply_par_vac1_xx8
)
: (
   vaccine_supply_par_vac1_xx9
))))))))))*((t < 0) ? (
   0
)
: ((t < 14) ? (
   1 - proportion_par_countryA_vac1_0
)
: ((t < 28) ? (
   1 - proportion_par_countryA_vac1_14
)
: ((t < 42) ? (
   1 - proportion_par_countryA_vac1_28
)
: ((t < 56) ? (
   1 - proportion_par_countryA_vac1_42
)
: ((t < 70) ? (
   1 - proportion_par_countryA_vac1_56
)
: ((t < 84) ? (
   1 - proportion_par_countryA_vac1_70
)
: ((t < 98) ? (
   1 - proportion_par_countryA_vac1_84
)
: ((t < 112) ? (
   1 - proportion_par_countryA_vac1_98
)
: ((t < 126) ? (
   1 - proportion_par_countryA_vac1_112
)
: ((t < 140) ? (
   1 - proportion_par_countryA_vac1_126
)
: (
   0
))))))))))))/std::pow(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0, 2);
    dydx[291] = -((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac2_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac2_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac2_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac2_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac2_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac2_xx5
)
: ((t >= xx6 && t < xx7) ? (
   vaccine_supply_par_vac2_xx6
)
: ((t >= xx7 && t < xx8) ? (
   vaccine_supply_par_vac2_xx7
)
: ((t >= xx8 && t < xx9) ? (
   vaccine_supply_par_vac2_xx8
)
: (
   vaccine_supply_par_vac2_xx9
))))))))))*((t < 0) ? (
   0
)
: ((t < 14) ? (
   1 - proportion_par_countryA_vac2_0
)
: ((t < 28) ? (
   1 - proportion_par_countryA_vac2_14
)
: ((t < 42) ? (
   1 - proportion_par_countryA_vac2_28
)
: ((t < 56) ? (
   1 - proportion_par_countryA_vac2_42
)
: ((t < 70) ? (
   1 - proportion_par_countryA_vac2_56
)
: ((t < 84) ? (
   1 - proportion_par_countryA_vac2_70
)
: ((t < 98) ? (
   1 - proportion_par_countryA_vac2_84
)
: ((t < 112) ? (
   1 - proportion_par_countryA_vac2_98
)
: ((t < 126) ? (
   1 - proportion_par_countryA_vac2_112
)
: ((t < 140) ? (
   1 - proportion_par_countryA_vac2_126
)
: (
   0
))))))))))))/std::pow(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0, 2);
    dydx[301] = -((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac1_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac1_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac1_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac1_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac1_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac1_xx5
)
: ((t >= xx6 && t < xx7) ? (
   vaccine_supply_par_vac1_xx6
)
: ((t >= xx7 && t < xx8) ? (
   vaccine_supply_par_vac1_xx7
)
: ((t >= xx8 && t < xx9) ? (
   vaccine_supply_par_vac1_xx8
)
: (
   vaccine_supply_par_vac1_xx9
))))))))))*((t < 0) ? (
   0
)
: ((t < 14) ? (
   1 - proportion_par_countryA_vac1_0
)
: ((t < 28) ? (
   1 - proportion_par_countryA_vac1_14
)
: ((t < 42) ? (
   1 - proportion_par_countryA_vac1_28
)
: ((t < 56) ? (
   1 - proportion_par_countryA_vac1_42
)
: ((t < 70) ? (
   1 - proportion_par_countryA_vac1_56
)
: ((t < 84) ? (
   1 - proportion_par_countryA_vac1_70
)
: ((t < 98) ? (
   1 - proportion_par_countryA_vac1_84
)
: ((t < 112) ? (
   1 - proportion_par_countryA_vac1_98
)
: ((t < 126) ? (
   1 - proportion_par_countryA_vac1_112
)
: ((t < 140) ? (
   1 - proportion_par_countryA_vac1_126
)
: (
   0
))))))))))))/std::pow(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0, 2);
    dydx[303] = -((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac2_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac2_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac2_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac2_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac2_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac2_xx5
)
: ((t >= xx6 && t < xx7) ? (
   vaccine_supply_par_vac2_xx6
)
: ((t >= xx7 && t < xx8) ? (
   vaccine_supply_par_vac2_xx7
)
: ((t >= xx8 && t < xx9) ? (
   vaccine_supply_par_vac2_xx8
)
: (
   vaccine_supply_par_vac2_xx9
))))))))))*((t < 0) ? (
   0
)
: ((t < 14) ? (
   1 - proportion_par_countryA_vac2_0
)
: ((t < 28) ? (
   1 - proportion_par_countryA_vac2_14
)
: ((t < 42) ? (
   1 - proportion_par_countryA_vac2_28
)
: ((t < 56) ? (
   1 - proportion_par_countryA_vac2_42
)
: ((t < 70) ? (
   1 - proportion_par_countryA_vac2_56
)
: ((t < 84) ? (
   1 - proportion_par_countryA_vac2_70
)
: ((t < 98) ? (
   1 - proportion_par_countryA_vac2_84
)
: ((t < 112) ? (
   1 - proportion_par_countryA_vac2_98
)
: ((t < 126) ? (
   1 - proportion_par_countryA_vac2_112
)
: ((t < 140) ? (
   1 - proportion_par_countryA_vac2_126
)
: (
   0
))))))))))))/std::pow(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0, 2);
}

} // namespace amici
} // namespace model_vaccination_piecewise