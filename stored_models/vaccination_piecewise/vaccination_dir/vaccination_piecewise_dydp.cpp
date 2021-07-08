#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "p.h"
#include "k.h"
#include "w.h"
#include "x.h"
#include "dwdp.h"

namespace amici {
namespace model_vaccination_piecewise {

void dydp_vaccination_piecewise(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dwdp){
    switch(ip) {
        case 56:
            dydp[0] = ((t < 0) ? (
   0
)
: ((t < 14) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[1] = ((t < 0) ? (
   0
)
: ((t < 14) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[4] = ((t < 0) ? (
   0
)
: ((t < 14) ? (
   1
)
: (
   0
)));
            dydp[5] = ((t < 0) ? (
   0
)
: ((t < 14) ? (
   -1
)
: (
   0
)));
            dydp[8] = ((t < 0) ? (
   0
)
: ((t < 14) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            dydp[10] = ((t < 0) ? (
   0
)
: ((t < 14) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            break;
        case 57:
            dydp[0] = ((t < 14) ? (
   0
)
: ((t < 28) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[1] = ((t < 14) ? (
   0
)
: ((t < 28) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[4] = ((t < 14) ? (
   0
)
: ((t < 28) ? (
   1
)
: (
   0
)));
            dydp[5] = ((t < 14) ? (
   0
)
: ((t < 28) ? (
   -1
)
: (
   0
)));
            dydp[8] = ((t < 14) ? (
   0
)
: ((t < 28) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            dydp[10] = ((t < 14) ? (
   0
)
: ((t < 28) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            break;
        case 58:
            dydp[0] = ((t < 28) ? (
   0
)
: ((t < 42) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[1] = ((t < 28) ? (
   0
)
: ((t < 42) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[4] = ((t < 28) ? (
   0
)
: ((t < 42) ? (
   1
)
: (
   0
)));
            dydp[5] = ((t < 28) ? (
   0
)
: ((t < 42) ? (
   -1
)
: (
   0
)));
            dydp[8] = ((t < 28) ? (
   0
)
: ((t < 42) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            dydp[10] = ((t < 28) ? (
   0
)
: ((t < 42) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            break;
        case 59:
            dydp[0] = ((t < 42) ? (
   0
)
: ((t < 56) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[1] = ((t < 42) ? (
   0
)
: ((t < 56) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[4] = ((t < 42) ? (
   0
)
: ((t < 56) ? (
   1
)
: (
   0
)));
            dydp[5] = ((t < 42) ? (
   0
)
: ((t < 56) ? (
   -1
)
: (
   0
)));
            dydp[8] = ((t < 42) ? (
   0
)
: ((t < 56) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            dydp[10] = ((t < 42) ? (
   0
)
: ((t < 56) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            break;
        case 60:
            dydp[0] = ((t < 56) ? (
   0
)
: ((t < 70) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[1] = ((t < 56) ? (
   0
)
: ((t < 70) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[4] = ((t < 56) ? (
   0
)
: ((t < 70) ? (
   1
)
: (
   0
)));
            dydp[5] = ((t < 56) ? (
   0
)
: ((t < 70) ? (
   -1
)
: (
   0
)));
            dydp[8] = ((t < 56) ? (
   0
)
: ((t < 70) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            dydp[10] = ((t < 56) ? (
   0
)
: ((t < 70) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            break;
        case 61:
            dydp[0] = ((t < 70) ? (
   0
)
: ((t < 84) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[1] = ((t < 70) ? (
   0
)
: ((t < 84) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[4] = ((t < 70) ? (
   0
)
: ((t < 84) ? (
   1
)
: (
   0
)));
            dydp[5] = ((t < 70) ? (
   0
)
: ((t < 84) ? (
   -1
)
: (
   0
)));
            dydp[8] = ((t < 70) ? (
   0
)
: ((t < 84) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            dydp[10] = ((t < 70) ? (
   0
)
: ((t < 84) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            break;
        case 62:
            dydp[0] = ((t < 84) ? (
   0
)
: ((t < 98) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[1] = ((t < 84) ? (
   0
)
: ((t < 98) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[4] = ((t < 84) ? (
   0
)
: ((t < 98) ? (
   1
)
: (
   0
)));
            dydp[5] = ((t < 84) ? (
   0
)
: ((t < 98) ? (
   -1
)
: (
   0
)));
            dydp[8] = ((t < 84) ? (
   0
)
: ((t < 98) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            dydp[10] = ((t < 84) ? (
   0
)
: ((t < 98) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            break;
        case 63:
            dydp[0] = ((t < 98) ? (
   0
)
: ((t < 112) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[1] = ((t < 98) ? (
   0
)
: ((t < 112) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[4] = ((t < 98) ? (
   0
)
: ((t < 112) ? (
   1
)
: (
   0
)));
            dydp[5] = ((t < 98) ? (
   0
)
: ((t < 112) ? (
   -1
)
: (
   0
)));
            dydp[8] = ((t < 98) ? (
   0
)
: ((t < 112) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            dydp[10] = ((t < 98) ? (
   0
)
: ((t < 112) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            break;
        case 64:
            dydp[0] = ((t < 112) ? (
   0
)
: ((t < 126) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[1] = ((t < 112) ? (
   0
)
: ((t < 126) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[4] = ((t < 112) ? (
   0
)
: ((t < 126) ? (
   1
)
: (
   0
)));
            dydp[5] = ((t < 112) ? (
   0
)
: ((t < 126) ? (
   -1
)
: (
   0
)));
            dydp[8] = ((t < 112) ? (
   0
)
: ((t < 126) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            dydp[10] = ((t < 112) ? (
   0
)
: ((t < 126) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            break;
        case 65:
            dydp[0] = ((t < 126) ? (
   0
)
: ((t < 140) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[1] = ((t < 126) ? (
   0
)
: ((t < 140) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[4] = ((t < 126) ? (
   0
)
: ((t < 140) ? (
   1
)
: (
   0
)));
            dydp[5] = ((t < 126) ? (
   0
)
: ((t < 140) ? (
   -1
)
: (
   0
)));
            dydp[8] = ((t < 126) ? (
   0
)
: ((t < 140) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            dydp[10] = ((t < 126) ? (
   0
)
: ((t < 140) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            break;
        case 66:
            dydp[0] = ((t < 140) ? (
   0
)
: ((t < 154) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[1] = ((t < 140) ? (
   0
)
: ((t < 154) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[4] = ((t < 140) ? (
   0
)
: ((t < 154) ? (
   1
)
: (
   0
)));
            dydp[5] = ((t < 140) ? (
   0
)
: ((t < 154) ? (
   -1
)
: (
   0
)));
            dydp[8] = ((t < 140) ? (
   0
)
: ((t < 154) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            dydp[10] = ((t < 140) ? (
   0
)
: ((t < 154) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            break;
        case 67:
            dydp[2] = ((t < 0) ? (
   0
)
: ((t < 14) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[3] = ((t < 0) ? (
   0
)
: ((t < 14) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[6] = ((t < 0) ? (
   0
)
: ((t < 14) ? (
   1
)
: (
   0
)));
            dydp[7] = ((t < 0) ? (
   0
)
: ((t < 14) ? (
   -1
)
: (
   0
)));
            dydp[9] = ((t < 0) ? (
   0
)
: ((t < 14) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            dydp[11] = ((t < 0) ? (
   0
)
: ((t < 14) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            break;
        case 68:
            dydp[2] = ((t < 14) ? (
   0
)
: ((t < 28) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[3] = ((t < 14) ? (
   0
)
: ((t < 28) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[6] = ((t < 14) ? (
   0
)
: ((t < 28) ? (
   1
)
: (
   0
)));
            dydp[7] = ((t < 14) ? (
   0
)
: ((t < 28) ? (
   -1
)
: (
   0
)));
            dydp[9] = ((t < 14) ? (
   0
)
: ((t < 28) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            dydp[11] = ((t < 14) ? (
   0
)
: ((t < 28) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            break;
        case 69:
            dydp[2] = ((t < 28) ? (
   0
)
: ((t < 42) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[3] = ((t < 28) ? (
   0
)
: ((t < 42) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[6] = ((t < 28) ? (
   0
)
: ((t < 42) ? (
   1
)
: (
   0
)));
            dydp[7] = ((t < 28) ? (
   0
)
: ((t < 42) ? (
   -1
)
: (
   0
)));
            dydp[9] = ((t < 28) ? (
   0
)
: ((t < 42) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            dydp[11] = ((t < 28) ? (
   0
)
: ((t < 42) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            break;
        case 70:
            dydp[2] = ((t < 42) ? (
   0
)
: ((t < 56) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[3] = ((t < 42) ? (
   0
)
: ((t < 56) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[6] = ((t < 42) ? (
   0
)
: ((t < 56) ? (
   1
)
: (
   0
)));
            dydp[7] = ((t < 42) ? (
   0
)
: ((t < 56) ? (
   -1
)
: (
   0
)));
            dydp[9] = ((t < 42) ? (
   0
)
: ((t < 56) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            dydp[11] = ((t < 42) ? (
   0
)
: ((t < 56) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            break;
        case 71:
            dydp[2] = ((t < 56) ? (
   0
)
: ((t < 70) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[3] = ((t < 56) ? (
   0
)
: ((t < 70) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[6] = ((t < 56) ? (
   0
)
: ((t < 70) ? (
   1
)
: (
   0
)));
            dydp[7] = ((t < 56) ? (
   0
)
: ((t < 70) ? (
   -1
)
: (
   0
)));
            dydp[9] = ((t < 56) ? (
   0
)
: ((t < 70) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            dydp[11] = ((t < 56) ? (
   0
)
: ((t < 70) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            break;
        case 72:
            dydp[2] = ((t < 70) ? (
   0
)
: ((t < 84) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[3] = ((t < 70) ? (
   0
)
: ((t < 84) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[6] = ((t < 70) ? (
   0
)
: ((t < 84) ? (
   1
)
: (
   0
)));
            dydp[7] = ((t < 70) ? (
   0
)
: ((t < 84) ? (
   -1
)
: (
   0
)));
            dydp[9] = ((t < 70) ? (
   0
)
: ((t < 84) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            dydp[11] = ((t < 70) ? (
   0
)
: ((t < 84) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            break;
        case 73:
            dydp[2] = ((t < 84) ? (
   0
)
: ((t < 98) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[3] = ((t < 84) ? (
   0
)
: ((t < 98) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[6] = ((t < 84) ? (
   0
)
: ((t < 98) ? (
   1
)
: (
   0
)));
            dydp[7] = ((t < 84) ? (
   0
)
: ((t < 98) ? (
   -1
)
: (
   0
)));
            dydp[9] = ((t < 84) ? (
   0
)
: ((t < 98) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            dydp[11] = ((t < 84) ? (
   0
)
: ((t < 98) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            break;
        case 74:
            dydp[2] = ((t < 98) ? (
   0
)
: ((t < 112) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[3] = ((t < 98) ? (
   0
)
: ((t < 112) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[6] = ((t < 98) ? (
   0
)
: ((t < 112) ? (
   1
)
: (
   0
)));
            dydp[7] = ((t < 98) ? (
   0
)
: ((t < 112) ? (
   -1
)
: (
   0
)));
            dydp[9] = ((t < 98) ? (
   0
)
: ((t < 112) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            dydp[11] = ((t < 98) ? (
   0
)
: ((t < 112) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            break;
        case 75:
            dydp[2] = ((t < 112) ? (
   0
)
: ((t < 126) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[3] = ((t < 112) ? (
   0
)
: ((t < 126) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[6] = ((t < 112) ? (
   0
)
: ((t < 126) ? (
   1
)
: (
   0
)));
            dydp[7] = ((t < 112) ? (
   0
)
: ((t < 126) ? (
   -1
)
: (
   0
)));
            dydp[9] = ((t < 112) ? (
   0
)
: ((t < 126) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            dydp[11] = ((t < 112) ? (
   0
)
: ((t < 126) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            break;
        case 76:
            dydp[2] = ((t < 126) ? (
   0
)
: ((t < 140) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[3] = ((t < 126) ? (
   0
)
: ((t < 140) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[6] = ((t < 126) ? (
   0
)
: ((t < 140) ? (
   1
)
: (
   0
)));
            dydp[7] = ((t < 126) ? (
   0
)
: ((t < 140) ? (
   -1
)
: (
   0
)));
            dydp[9] = ((t < 126) ? (
   0
)
: ((t < 140) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            dydp[11] = ((t < 126) ? (
   0
)
: ((t < 140) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            break;
        case 77:
            dydp[2] = ((t < 140) ? (
   0
)
: ((t < 154) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[3] = ((t < 140) ? (
   0
)
: ((t < 154) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[6] = ((t < 140) ? (
   0
)
: ((t < 154) ? (
   1
)
: (
   0
)));
            dydp[7] = ((t < 140) ? (
   0
)
: ((t < 154) ? (
   -1
)
: (
   0
)));
            dydp[9] = ((t < 140) ? (
   0
)
: ((t < 154) ? (
   1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            dydp[11] = ((t < 140) ? (
   0
)
: ((t < 154) ? (
   -1
)
: (
   0
)))*((t >= xx0 && t < xx1) ? (
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
   0
))))))))));
            break;
        case 88:
            dydp[0] = ((t >= xx0 && t < xx1) ? (
   1
)
: (
   0
))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[1] = ((t >= xx0 && t < xx1) ? (
   1
)
: (
   0
))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[8] = ((t >= xx0 && t < xx1) ? (
   1
)
: (
   0
))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))));
            dydp[10] = ((t >= xx0 && t < xx1) ? (
   1
)
: (
   0
))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))));
            break;
        case 89:
            dydp[0] = ((t >= xx0 && t < xx1) ? (
   0
)
: ((t >= xx1 && t < xx2) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[1] = ((t >= xx0 && t < xx1) ? (
   0
)
: ((t >= xx1 && t < xx2) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[8] = ((t >= xx0 && t < xx1) ? (
   0
)
: ((t >= xx1 && t < xx2) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))));
            dydp[10] = ((t >= xx0 && t < xx1) ? (
   0
)
: ((t >= xx1 && t < xx2) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))));
            break;
        case 90:
            dydp[0] = (((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
   0
)
: ((t >= xx2 && t < xx3) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[1] = (((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
   0
)
: ((t >= xx2 && t < xx3) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[8] = (((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
   0
)
: ((t >= xx2 && t < xx3) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))));
            dydp[10] = (((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
   0
)
: ((t >= xx2 && t < xx3) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))));
            break;
        case 91:
            dydp[0] = (((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
   0
)
: ((t >= xx3 && t < xx4) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[1] = (((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
   0
)
: ((t >= xx3 && t < xx4) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[8] = (((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
   0
)
: ((t >= xx3 && t < xx4) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))));
            dydp[10] = (((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
   0
)
: ((t >= xx3 && t < xx4) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))));
            break;
        case 92:
            dydp[0] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
   0
)
: ((t >= xx4 && t < xx5) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[1] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
   0
)
: ((t >= xx4 && t < xx5) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[8] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
   0
)
: ((t >= xx4 && t < xx5) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))));
            dydp[10] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
   0
)
: ((t >= xx4 && t < xx5) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))));
            break;
        case 93:
            dydp[0] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
   0
)
: ((t >= xx5 && t < xx6) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[1] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
   0
)
: ((t >= xx5 && t < xx6) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[8] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
   0
)
: ((t >= xx5 && t < xx6) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))));
            dydp[10] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
   0
)
: ((t >= xx5 && t < xx6) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))));
            break;
        case 94:
            dydp[0] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
   0
)
: ((t >= xx6 && t < xx7) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[1] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
   0
)
: ((t >= xx6 && t < xx7) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[8] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
   0
)
: ((t >= xx6 && t < xx7) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))));
            dydp[10] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
   0
)
: ((t >= xx6 && t < xx7) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))));
            break;
        case 95:
            dydp[0] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
   0
)
: ((t >= xx7 && t < xx8) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[1] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
   0
)
: ((t >= xx7 && t < xx8) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[8] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
   0
)
: ((t >= xx7 && t < xx8) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))));
            dydp[10] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
   0
)
: ((t >= xx7 && t < xx8) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))));
            break;
        case 96:
            dydp[0] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
   0
)
: ((t >= xx8 && t < xx9) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[1] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
   0
)
: ((t >= xx8 && t < xx9) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[8] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
   0
)
: ((t >= xx8 && t < xx9) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))));
            dydp[10] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
   0
)
: ((t >= xx8 && t < xx9) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac1_140
)
: (
   0
)))))))))))));
            break;
        case 98:
            dydp[2] = ((t >= xx0 && t < xx1) ? (
   1
)
: (
   0
))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[3] = ((t >= xx0 && t < xx1) ? (
   1
)
: (
   0
))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[9] = ((t >= xx0 && t < xx1) ? (
   1
)
: (
   0
))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))));
            dydp[11] = ((t >= xx0 && t < xx1) ? (
   1
)
: (
   0
))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))));
            break;
        case 99:
            dydp[2] = ((t >= xx0 && t < xx1) ? (
   0
)
: ((t >= xx1 && t < xx2) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[3] = ((t >= xx0 && t < xx1) ? (
   0
)
: ((t >= xx1 && t < xx2) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[9] = ((t >= xx0 && t < xx1) ? (
   0
)
: ((t >= xx1 && t < xx2) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))));
            dydp[11] = ((t >= xx0 && t < xx1) ? (
   0
)
: ((t >= xx1 && t < xx2) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))));
            break;
        case 100:
            dydp[2] = (((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
   0
)
: ((t >= xx2 && t < xx3) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[3] = (((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
   0
)
: ((t >= xx2 && t < xx3) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[9] = (((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
   0
)
: ((t >= xx2 && t < xx3) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))));
            dydp[11] = (((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
   0
)
: ((t >= xx2 && t < xx3) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))));
            break;
        case 101:
            dydp[2] = (((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
   0
)
: ((t >= xx3 && t < xx4) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[3] = (((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
   0
)
: ((t >= xx3 && t < xx4) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[9] = (((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
   0
)
: ((t >= xx3 && t < xx4) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))));
            dydp[11] = (((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
   0
)
: ((t >= xx3 && t < xx4) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))));
            break;
        case 102:
            dydp[2] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
   0
)
: ((t >= xx4 && t < xx5) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[3] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
   0
)
: ((t >= xx4 && t < xx5) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[9] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
   0
)
: ((t >= xx4 && t < xx5) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))));
            dydp[11] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
   0
)
: ((t >= xx4 && t < xx5) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))));
            break;
        case 103:
            dydp[2] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
   0
)
: ((t >= xx5 && t < xx6) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[3] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
   0
)
: ((t >= xx5 && t < xx6) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[9] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
   0
)
: ((t >= xx5 && t < xx6) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))));
            dydp[11] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
   0
)
: ((t >= xx5 && t < xx6) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))));
            break;
        case 104:
            dydp[2] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
   0
)
: ((t >= xx6 && t < xx7) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[3] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
   0
)
: ((t >= xx6 && t < xx7) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[9] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
   0
)
: ((t >= xx6 && t < xx7) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))));
            dydp[11] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
   0
)
: ((t >= xx6 && t < xx7) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))));
            break;
        case 105:
            dydp[2] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
   0
)
: ((t >= xx7 && t < xx8) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[3] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
   0
)
: ((t >= xx7 && t < xx8) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[9] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
   0
)
: ((t >= xx7 && t < xx8) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))));
            dydp[11] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
   0
)
: ((t >= xx7 && t < xx8) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))));
            break;
        case 106:
            dydp[2] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
   0
)
: ((t >= xx8 && t < xx9) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
            dydp[3] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
   0
)
: ((t >= xx8 && t < xx9) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[9] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
   0
)
: ((t >= xx8 && t < xx9) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))));
            dydp[11] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
   0
)
: ((t >= xx8 && t < xx9) ? (
   1
)
: (
   0
)))*((t < 0) ? (
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
: ((t < 154) ? (
   1 - proportion_par_countryA_vac2_140
)
: (
   0
)))))))))))));
            break;
    }
}

} // namespace amici
} // namespace model_vaccination_piecewise