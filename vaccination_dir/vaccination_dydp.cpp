#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "p.h"
#include "k.h"
#include "w.h"
#include "x.h"
#include "dwdp.h"

namespace amici {
namespace model_vaccination {

void dydp_vaccination(realtype *dydp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const int ip, const realtype *w, const realtype *dwdp){
    switch(ip) {
        case 4:
            dydp[0] = dwdp168*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[1] = dwdp169*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            dydp[2] = -dwdp168*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[3] = -dwdp169*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            break;
        case 57:
            dydp[0] = ((t >= xx0 && t < xx1) ? (
   1
)
: (
   0
))/(1 + std::pow(exponentiale, -spline_countryA_vac1));
            dydp[2] = (1 - 1/(1 + std::pow(exponentiale, -spline_countryA_vac1)))*((t >= xx0 && t < xx1) ? (
   1
)
: (
   0
));
            break;
        case 58:
            dydp[0] = ((t >= xx0 && t < xx1) ? (
   0
)
: ((t >= xx1 && t < xx2) ? (
   1
)
: (
   0
)))/(1 + std::pow(exponentiale, -spline_countryA_vac1));
            dydp[2] = (1 - 1/(1 + std::pow(exponentiale, -spline_countryA_vac1)))*((t >= xx0 && t < xx1) ? (
   0
)
: ((t >= xx1 && t < xx2) ? (
   1
)
: (
   0
)));
            break;
        case 59:
            dydp[0] = (((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
   0
)
: ((t >= xx2 && t < xx3) ? (
   1
)
: (
   0
)))/(1 + std::pow(exponentiale, -spline_countryA_vac1));
            dydp[2] = (1 - 1/(1 + std::pow(exponentiale, -spline_countryA_vac1)))*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
   0
)
: ((t >= xx2 && t < xx3) ? (
   1
)
: (
   0
)));
            break;
        case 60:
            dydp[0] = (((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
   0
)
: ((t >= xx3 && t < xx4) ? (
   1
)
: (
   0
)))/(1 + std::pow(exponentiale, -spline_countryA_vac1));
            dydp[2] = (1 - 1/(1 + std::pow(exponentiale, -spline_countryA_vac1)))*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
   0
)
: ((t >= xx3 && t < xx4) ? (
   1
)
: (
   0
)));
            break;
        case 61:
            dydp[0] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
   0
)
: ((t >= xx4 && t < xx5) ? (
   1
)
: (
   0
)))/(1 + std::pow(exponentiale, -spline_countryA_vac1));
            dydp[2] = (1 - 1/(1 + std::pow(exponentiale, -spline_countryA_vac1)))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
   0
)
: ((t >= xx4 && t < xx5) ? (
   1
)
: (
   0
)));
            break;
        case 62:
            dydp[0] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
   0
)
: ((t >= xx5 && t < xx6) ? (
   1
)
: (
   0
)))/(1 + std::pow(exponentiale, -spline_countryA_vac1));
            dydp[2] = (1 - 1/(1 + std::pow(exponentiale, -spline_countryA_vac1)))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
   0
)
: ((t >= xx5 && t < xx6) ? (
   1
)
: (
   0
)));
            break;
        case 63:
            dydp[0] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
   0
)
: ((t >= xx6 && t < xx7) ? (
   1
)
: (
   0
)))/(1 + std::pow(exponentiale, -spline_countryA_vac1));
            dydp[2] = (1 - 1/(1 + std::pow(exponentiale, -spline_countryA_vac1)))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
   0
)
: ((t >= xx6 && t < xx7) ? (
   1
)
: (
   0
)));
            break;
        case 64:
            dydp[0] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
   0
)
: ((t >= xx7 && t < xx8) ? (
   1
)
: (
   0
)))/(1 + std::pow(exponentiale, -spline_countryA_vac1));
            dydp[2] = (1 - 1/(1 + std::pow(exponentiale, -spline_countryA_vac1)))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
   0
)
: ((t >= xx7 && t < xx8) ? (
   1
)
: (
   0
)));
            break;
        case 65:
            dydp[0] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
   0
)
: ((t >= xx8 && t < xx9) ? (
   1
)
: (
   0
)))/(1 + std::pow(exponentiale, -spline_countryA_vac1));
            dydp[2] = (1 - 1/(1 + std::pow(exponentiale, -spline_countryA_vac1)))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
   0
)
: ((t >= xx8 && t < xx9) ? (
   1
)
: (
   0
)));
            break;
        case 66:
            dydp[0] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: ((t >= xx9 && t < xx10) ? (
   1
)
: (
   0
)))/(1 + std::pow(exponentiale, -spline_countryA_vac1));
            dydp[2] = (1 - 1/(1 + std::pow(exponentiale, -spline_countryA_vac1)))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: ((t >= xx9 && t < xx10) ? (
   1
)
: (
   0
)));
            break;
        case 68:
            dydp[1] = ((t >= xx0 && t < xx1) ? (
   1
)
: (
   0
))/(1 + std::pow(exponentiale, -spline_countryA_vac2));
            dydp[3] = (1 - 1/(1 + std::pow(exponentiale, -spline_countryA_vac2)))*((t >= xx0 && t < xx1) ? (
   1
)
: (
   0
));
            break;
        case 69:
            dydp[1] = ((t >= xx0 && t < xx1) ? (
   0
)
: ((t >= xx1 && t < xx2) ? (
   1
)
: (
   0
)))/(1 + std::pow(exponentiale, -spline_countryA_vac2));
            dydp[3] = (1 - 1/(1 + std::pow(exponentiale, -spline_countryA_vac2)))*((t >= xx0 && t < xx1) ? (
   0
)
: ((t >= xx1 && t < xx2) ? (
   1
)
: (
   0
)));
            break;
        case 70:
            dydp[1] = (((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
   0
)
: ((t >= xx2 && t < xx3) ? (
   1
)
: (
   0
)))/(1 + std::pow(exponentiale, -spline_countryA_vac2));
            dydp[3] = (1 - 1/(1 + std::pow(exponentiale, -spline_countryA_vac2)))*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
   0
)
: ((t >= xx2 && t < xx3) ? (
   1
)
: (
   0
)));
            break;
        case 71:
            dydp[1] = (((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
   0
)
: ((t >= xx3 && t < xx4) ? (
   1
)
: (
   0
)))/(1 + std::pow(exponentiale, -spline_countryA_vac2));
            dydp[3] = (1 - 1/(1 + std::pow(exponentiale, -spline_countryA_vac2)))*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
   0
)
: ((t >= xx3 && t < xx4) ? (
   1
)
: (
   0
)));
            break;
        case 72:
            dydp[1] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
   0
)
: ((t >= xx4 && t < xx5) ? (
   1
)
: (
   0
)))/(1 + std::pow(exponentiale, -spline_countryA_vac2));
            dydp[3] = (1 - 1/(1 + std::pow(exponentiale, -spline_countryA_vac2)))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
   0
)
: ((t >= xx4 && t < xx5) ? (
   1
)
: (
   0
)));
            break;
        case 73:
            dydp[1] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
   0
)
: ((t >= xx5 && t < xx6) ? (
   1
)
: (
   0
)))/(1 + std::pow(exponentiale, -spline_countryA_vac2));
            dydp[3] = (1 - 1/(1 + std::pow(exponentiale, -spline_countryA_vac2)))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
   0
)
: ((t >= xx5 && t < xx6) ? (
   1
)
: (
   0
)));
            break;
        case 74:
            dydp[1] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
   0
)
: ((t >= xx6 && t < xx7) ? (
   1
)
: (
   0
)))/(1 + std::pow(exponentiale, -spline_countryA_vac2));
            dydp[3] = (1 - 1/(1 + std::pow(exponentiale, -spline_countryA_vac2)))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
   0
)
: ((t >= xx6 && t < xx7) ? (
   1
)
: (
   0
)));
            break;
        case 75:
            dydp[1] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
   0
)
: ((t >= xx7 && t < xx8) ? (
   1
)
: (
   0
)))/(1 + std::pow(exponentiale, -spline_countryA_vac2));
            dydp[3] = (1 - 1/(1 + std::pow(exponentiale, -spline_countryA_vac2)))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
   0
)
: ((t >= xx7 && t < xx8) ? (
   1
)
: (
   0
)));
            break;
        case 76:
            dydp[1] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
   0
)
: ((t >= xx8 && t < xx9) ? (
   1
)
: (
   0
)))/(1 + std::pow(exponentiale, -spline_countryA_vac2));
            dydp[3] = (1 - 1/(1 + std::pow(exponentiale, -spline_countryA_vac2)))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
   0
)
: ((t >= xx8 && t < xx9) ? (
   1
)
: (
   0
)));
            break;
        case 77:
            dydp[1] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: ((t >= xx9 && t < xx10) ? (
   1
)
: (
   0
)))/(1 + std::pow(exponentiale, -spline_countryA_vac2));
            dydp[3] = (1 - 1/(1 + std::pow(exponentiale, -spline_countryA_vac2)))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: ((t >= xx9 && t < xx10) ? (
   1
)
: (
   0
)));
            break;
        case 83:
            dydp[0] = dwdp598*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[1] = dwdp599*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            dydp[2] = -dwdp598*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[3] = -dwdp599*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            break;
        case 84:
            dydp[0] = dwdp620*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[1] = dwdp621*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            dydp[2] = -dwdp620*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[3] = -dwdp621*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            break;
        case 85:
            dydp[0] = dwdp642*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[1] = dwdp643*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            dydp[2] = -dwdp642*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[3] = -dwdp643*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            break;
        case 86:
            dydp[0] = dwdp664*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[1] = dwdp665*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            dydp[2] = -dwdp664*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[3] = -dwdp665*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            break;
        case 87:
            dydp[0] = dwdp686*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[1] = dwdp687*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            dydp[2] = -dwdp686*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[3] = -dwdp687*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            break;
        case 88:
            dydp[0] = dwdp708*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[1] = dwdp709*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            dydp[2] = -dwdp708*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[3] = -dwdp709*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            break;
        case 89:
            dydp[0] = dwdp730*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[1] = dwdp731*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            dydp[2] = -dwdp730*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[3] = -dwdp731*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            break;
        case 90:
            dydp[0] = dwdp752*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[1] = dwdp753*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            dydp[2] = -dwdp752*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[3] = -dwdp753*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            break;
        case 91:
            dydp[0] = dwdp774*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[1] = dwdp775*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            dydp[2] = -dwdp774*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[3] = -dwdp775*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            break;
        case 92:
            dydp[0] = dwdp796*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[1] = dwdp797*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            dydp[2] = -dwdp796*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[3] = -dwdp797*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            break;
        case 93:
            dydp[0] = dwdp818*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[1] = dwdp819*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            dydp[2] = -dwdp818*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[3] = -dwdp819*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            break;
        case 94:
            dydp[0] = dwdp840*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[2] = -dwdp840*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            break;
        case 95:
            dydp[0] = dwdp851*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[2] = -dwdp851*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            break;
        case 96:
            dydp[0] = dwdp862*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[2] = -dwdp862*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            break;
        case 97:
            dydp[0] = dwdp873*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[2] = -dwdp873*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            break;
        case 98:
            dydp[0] = dwdp884*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[2] = -dwdp884*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            break;
        case 99:
            dydp[0] = dwdp895*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[2] = -dwdp895*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            break;
        case 100:
            dydp[0] = dwdp906*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[2] = -dwdp906*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            break;
        case 101:
            dydp[0] = dwdp917*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[2] = -dwdp917*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            break;
        case 102:
            dydp[0] = dwdp928*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[2] = -dwdp928*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            break;
        case 103:
            dydp[0] = dwdp939*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[2] = -dwdp939*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            break;
        case 104:
            dydp[0] = dwdp950*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            dydp[2] = -dwdp950*std::pow(exponentiale, -spline_countryA_vac1)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac1_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac1), 2);
            break;
        case 105:
            dydp[1] = dwdp961*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            dydp[3] = -dwdp961*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            break;
        case 106:
            dydp[1] = dwdp972*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            dydp[3] = -dwdp972*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            break;
        case 107:
            dydp[1] = dwdp983*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            dydp[3] = -dwdp983*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            break;
        case 108:
            dydp[1] = dwdp994*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            dydp[3] = -dwdp994*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            break;
        case 109:
            dydp[1] = dwdp1005*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            dydp[3] = -dwdp1005*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            break;
        case 110:
            dydp[1] = dwdp1016*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            dydp[3] = -dwdp1016*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            break;
        case 111:
            dydp[1] = dwdp1027*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            dydp[3] = -dwdp1027*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            break;
        case 112:
            dydp[1] = dwdp1038*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            dydp[3] = -dwdp1038*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            break;
        case 113:
            dydp[1] = dwdp1049*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            dydp[3] = -dwdp1049*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            break;
        case 114:
            dydp[1] = dwdp1060*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            dydp[3] = -dwdp1060*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            break;
        case 115:
            dydp[1] = dwdp1071*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            dydp[3] = -dwdp1071*std::pow(exponentiale, -spline_countryA_vac2)*((t >= xx0 && t < xx1) ? (
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
: ((t >= xx9 && t < xx10) ? (
   vaccine_supply_par_vac2_xx9
)
: (
   0
)))))))))))*std::log(exponentiale)/std::pow(1 + std::pow(exponentiale, -spline_countryA_vac2), 2);
            break;
    }
}

} // namespace amici
} // namespace model_vaccination