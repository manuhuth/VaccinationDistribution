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
        case 56:
            dydp[0] = ((t >= xx0 && t < xx1) ? (
   1
)
: (
   0
))/((std::pow(0.36787944123356736, spline_countryA_vac1) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1))*((t >= xx0 && t < xx1) ? (
   1
)
: (
   0
))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[10] = ((t >= xx0 && t < xx1) ? (
   1
)
: (
   0
))/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1);
            dydp[12] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1))*((t >= xx0 && t < xx1) ? (
   1
)
: (
   0
));
            break;
        case 57:
            dydp[0] = ((t >= xx0 && t < xx1) ? (
   0
)
: ((t >= xx1 && t < xx2) ? (
   1
)
: (
   0
)))/((std::pow(0.36787944123356736, spline_countryA_vac1) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1))*((t >= xx0 && t < xx1) ? (
   0
)
: ((t >= xx1 && t < xx2) ? (
   1
)
: (
   0
)))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[10] = ((t >= xx0 && t < xx1) ? (
   0
)
: ((t >= xx1 && t < xx2) ? (
   1
)
: (
   0
)))/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1);
            dydp[12] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1))*((t >= xx0 && t < xx1) ? (
   0
)
: ((t >= xx1 && t < xx2) ? (
   1
)
: (
   0
)));
            break;
        case 58:
            dydp[0] = (((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
   0
)
: ((t >= xx2 && t < xx3) ? (
   1
)
: (
   0
)))/((std::pow(0.36787944123356736, spline_countryA_vac1) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1))*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
   0
)
: ((t >= xx2 && t < xx3) ? (
   1
)
: (
   0
)))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[10] = (((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
   0
)
: ((t >= xx2 && t < xx3) ? (
   1
)
: (
   0
)))/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1);
            dydp[12] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1))*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
   0
)
: ((t >= xx2 && t < xx3) ? (
   1
)
: (
   0
)));
            break;
        case 59:
            dydp[0] = (((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
   0
)
: ((t >= xx3 && t < xx4) ? (
   1
)
: (
   0
)))/((std::pow(0.36787944123356736, spline_countryA_vac1) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
   0
)
: ((t >= xx3 && t < xx4) ? (
   1
)
: (
   0
)))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[10] = (((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
   0
)
: ((t >= xx3 && t < xx4) ? (
   1
)
: (
   0
)))/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1);
            dydp[12] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
   0
)
: ((t >= xx3 && t < xx4) ? (
   1
)
: (
   0
)));
            break;
        case 60:
            dydp[0] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
   0
)
: ((t >= xx4 && t < xx5) ? (
   1
)
: (
   0
)))/((std::pow(0.36787944123356736, spline_countryA_vac1) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
   0
)
: ((t >= xx4 && t < xx5) ? (
   1
)
: (
   0
)))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[10] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
   0
)
: ((t >= xx4 && t < xx5) ? (
   1
)
: (
   0
)))/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1);
            dydp[12] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
   0
)
: ((t >= xx4 && t < xx5) ? (
   1
)
: (
   0
)));
            break;
        case 61:
            dydp[0] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
   0
)
: ((t >= xx5 && t < xx6) ? (
   1
)
: (
   0
)))/((std::pow(0.36787944123356736, spline_countryA_vac1) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
   0
)
: ((t >= xx5 && t < xx6) ? (
   1
)
: (
   0
)))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[10] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
   0
)
: ((t >= xx5 && t < xx6) ? (
   1
)
: (
   0
)))/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1);
            dydp[12] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
   0
)
: ((t >= xx5 && t < xx6) ? (
   1
)
: (
   0
)));
            break;
        case 62:
            dydp[0] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
   0
)
: ((t >= xx6 && t < xx7) ? (
   1
)
: (
   0
)))/((std::pow(0.36787944123356736, spline_countryA_vac1) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
   0
)
: ((t >= xx6 && t < xx7) ? (
   1
)
: (
   0
)))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[10] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
   0
)
: ((t >= xx6 && t < xx7) ? (
   1
)
: (
   0
)))/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1);
            dydp[12] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
   0
)
: ((t >= xx6 && t < xx7) ? (
   1
)
: (
   0
)));
            break;
        case 63:
            dydp[0] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
   0
)
: ((t >= xx7 && t < xx8) ? (
   1
)
: (
   0
)))/((std::pow(0.36787944123356736, spline_countryA_vac1) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
   0
)
: ((t >= xx7 && t < xx8) ? (
   1
)
: (
   0
)))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[10] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
   0
)
: ((t >= xx7 && t < xx8) ? (
   1
)
: (
   0
)))/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1);
            dydp[12] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
   0
)
: ((t >= xx7 && t < xx8) ? (
   1
)
: (
   0
)));
            break;
        case 64:
            dydp[0] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
   0
)
: ((t >= xx8 && t < xx9) ? (
   1
)
: (
   0
)))/((std::pow(0.36787944123356736, spline_countryA_vac1) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
   0
)
: ((t >= xx8 && t < xx9) ? (
   1
)
: (
   0
)))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[10] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
   0
)
: ((t >= xx8 && t < xx9) ? (
   1
)
: (
   0
)))/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1);
            dydp[12] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
   0
)
: ((t >= xx8 && t < xx9) ? (
   1
)
: (
   0
)));
            break;
        case 65:
            dydp[0] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
))/((std::pow(0.36787944123356736, spline_countryA_vac1) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[10] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
))/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1);
            dydp[12] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac1) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
));
            break;
        case 66:
            dydp[2] = ((t >= xx0 && t < xx1) ? (
   1
)
: (
   0
))/((std::pow(0.36787944123356736, spline_countryA_vac2) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1))*((t >= xx0 && t < xx1) ? (
   1
)
: (
   0
))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[11] = ((t >= xx0 && t < xx1) ? (
   1
)
: (
   0
))/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1);
            dydp[13] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1))*((t >= xx0 && t < xx1) ? (
   1
)
: (
   0
));
            break;
        case 67:
            dydp[2] = ((t >= xx0 && t < xx1) ? (
   0
)
: ((t >= xx1 && t < xx2) ? (
   1
)
: (
   0
)))/((std::pow(0.36787944123356736, spline_countryA_vac2) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1))*((t >= xx0 && t < xx1) ? (
   0
)
: ((t >= xx1 && t < xx2) ? (
   1
)
: (
   0
)))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[11] = ((t >= xx0 && t < xx1) ? (
   0
)
: ((t >= xx1 && t < xx2) ? (
   1
)
: (
   0
)))/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1);
            dydp[13] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1))*((t >= xx0 && t < xx1) ? (
   0
)
: ((t >= xx1 && t < xx2) ? (
   1
)
: (
   0
)));
            break;
        case 68:
            dydp[2] = (((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
   0
)
: ((t >= xx2 && t < xx3) ? (
   1
)
: (
   0
)))/((std::pow(0.36787944123356736, spline_countryA_vac2) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1))*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
   0
)
: ((t >= xx2 && t < xx3) ? (
   1
)
: (
   0
)))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[11] = (((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
   0
)
: ((t >= xx2 && t < xx3) ? (
   1
)
: (
   0
)))/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1);
            dydp[13] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1))*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
   0
)
: ((t >= xx2 && t < xx3) ? (
   1
)
: (
   0
)));
            break;
        case 69:
            dydp[2] = (((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
   0
)
: ((t >= xx3 && t < xx4) ? (
   1
)
: (
   0
)))/((std::pow(0.36787944123356736, spline_countryA_vac2) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
   0
)
: ((t >= xx3 && t < xx4) ? (
   1
)
: (
   0
)))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[11] = (((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
   0
)
: ((t >= xx3 && t < xx4) ? (
   1
)
: (
   0
)))/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1);
            dydp[13] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
   0
)
: ((t >= xx3 && t < xx4) ? (
   1
)
: (
   0
)));
            break;
        case 70:
            dydp[2] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
   0
)
: ((t >= xx4 && t < xx5) ? (
   1
)
: (
   0
)))/((std::pow(0.36787944123356736, spline_countryA_vac2) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
   0
)
: ((t >= xx4 && t < xx5) ? (
   1
)
: (
   0
)))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[11] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
   0
)
: ((t >= xx4 && t < xx5) ? (
   1
)
: (
   0
)))/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1);
            dydp[13] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
   0
)
: ((t >= xx4 && t < xx5) ? (
   1
)
: (
   0
)));
            break;
        case 71:
            dydp[2] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
   0
)
: ((t >= xx5 && t < xx6) ? (
   1
)
: (
   0
)))/((std::pow(0.36787944123356736, spline_countryA_vac2) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
   0
)
: ((t >= xx5 && t < xx6) ? (
   1
)
: (
   0
)))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[11] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
   0
)
: ((t >= xx5 && t < xx6) ? (
   1
)
: (
   0
)))/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1);
            dydp[13] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
   0
)
: ((t >= xx5 && t < xx6) ? (
   1
)
: (
   0
)));
            break;
        case 72:
            dydp[2] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
   0
)
: ((t >= xx6 && t < xx7) ? (
   1
)
: (
   0
)))/((std::pow(0.36787944123356736, spline_countryA_vac2) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
   0
)
: ((t >= xx6 && t < xx7) ? (
   1
)
: (
   0
)))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[11] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
   0
)
: ((t >= xx6 && t < xx7) ? (
   1
)
: (
   0
)))/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1);
            dydp[13] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
   0
)
: ((t >= xx6 && t < xx7) ? (
   1
)
: (
   0
)));
            break;
        case 73:
            dydp[2] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
   0
)
: ((t >= xx7 && t < xx8) ? (
   1
)
: (
   0
)))/((std::pow(0.36787944123356736, spline_countryA_vac2) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
   0
)
: ((t >= xx7 && t < xx8) ? (
   1
)
: (
   0
)))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[11] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
   0
)
: ((t >= xx7 && t < xx8) ? (
   1
)
: (
   0
)))/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1);
            dydp[13] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
   0
)
: ((t >= xx7 && t < xx8) ? (
   1
)
: (
   0
)));
            break;
        case 74:
            dydp[2] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
   0
)
: ((t >= xx8 && t < xx9) ? (
   1
)
: (
   0
)))/((std::pow(0.36787944123356736, spline_countryA_vac2) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
   0
)
: ((t >= xx8 && t < xx9) ? (
   1
)
: (
   0
)))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[11] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
   0
)
: ((t >= xx8 && t < xx9) ? (
   1
)
: (
   0
)))/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1);
            dydp[13] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
   0
)
: ((t >= xx8 && t < xx9) ? (
   1
)
: (
   0
)));
            break;
        case 75:
            dydp[2] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
))/((std::pow(0.36787944123356736, spline_countryA_vac2) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
            dydp[11] = (((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
))/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1);
            dydp[13] = (1 - 1/(std::pow(0.36787944123356736, spline_countryA_vac2) + 1))*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
));
            break;
        case 80:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp576*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp576*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp577*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp577*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp576/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp576/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp577/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp577/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[8] = dwdp576;
            dydp[9] = dwdp577;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp576*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp577*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp576*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp577*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 81:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp598*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp598*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp599*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp599*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp598/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp598/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp599/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp599/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[8] = dwdp598;
            dydp[9] = dwdp599;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp598*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp599*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp598*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp599*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 82:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp620*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp620*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp621*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp621*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp620/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp620/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp621/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp621/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[8] = dwdp620;
            dydp[9] = dwdp621;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp620*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp621*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp620*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp621*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 83:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp642*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp642*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp643*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp643*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp642/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp642/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp643/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp643/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[8] = dwdp642;
            dydp[9] = dwdp643;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp642*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp643*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp642*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp643*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 84:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp664*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp664*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp665*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp665*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp664/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp664/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp665/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp665/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[8] = dwdp664;
            dydp[9] = dwdp665;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp664*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp665*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp664*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp665*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 85:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp686*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp686*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp687*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp687*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp686/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp686/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp687/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp687/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[8] = dwdp686;
            dydp[9] = dwdp687;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp686*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp687*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp686*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp687*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 86:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp708*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp708*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp709*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp709*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp708/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp708/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp709/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp709/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[8] = dwdp708;
            dydp[9] = dwdp709;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp708*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp709*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp708*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp709*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 87:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp730*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp730*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp731*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp731*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp730/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp730/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp731/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp731/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[8] = dwdp730;
            dydp[9] = dwdp731;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp730*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp731*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp730*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp731*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 88:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp752*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp752*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp753*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp753*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp752/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp752/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp753/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp753/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[8] = dwdp752;
            dydp[9] = dwdp753;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp752*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp753*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp752*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp753*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 89:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp774*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp774*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp775*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp775*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp774/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp774/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp775/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp775/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[8] = dwdp774;
            dydp[9] = dwdp775;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp774*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp775*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp774*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp775*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 90:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp796*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp796*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp797*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp797*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp796/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp796/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp797/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp797/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[8] = dwdp796;
            dydp[9] = dwdp797;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp796*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp797*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp796*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp797*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 91:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp818*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp818*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp818/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp818/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[8] = dwdp818;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp818*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp818*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            break;
        case 92:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp829*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp829*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp829/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp829/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[8] = dwdp829;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp829*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp829*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            break;
        case 93:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp840*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp840*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp840/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp840/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[8] = dwdp840;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp840*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp840*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            break;
        case 94:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp851*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp851*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp851/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp851/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[8] = dwdp851;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp851*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp851*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            break;
        case 95:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp862*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp862*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp862/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp862/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[8] = dwdp862;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp862*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp862*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            break;
        case 96:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp873*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp873*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp873/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp873/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[8] = dwdp873;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp873*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp873*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            break;
        case 97:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp884*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp884*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp884/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp884/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[8] = dwdp884;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp884*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp884*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            break;
        case 98:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp895*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp895*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp895/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp895/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[8] = dwdp895;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp895*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp895*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            break;
        case 99:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp906*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp906*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp906/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp906/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[8] = dwdp906;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp906*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp906*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            break;
        case 100:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp917*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp917*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp917/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp917/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[8] = dwdp917;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp917*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp917*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            break;
        case 101:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp928*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp928*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp928/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp928/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[8] = dwdp928;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp928*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp928*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            break;
        case 102:
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp939*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp939*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp939/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp939/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[9] = dwdp939;
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp939*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp939*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 103:
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp950*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp950*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp950/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp950/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[9] = dwdp950;
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp950*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp950*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 104:
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp961*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp961*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp961/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp961/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[9] = dwdp961;
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp961*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp961*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 105:
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp972*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp972*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp972/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp972/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[9] = dwdp972;
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp972*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp972*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 106:
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp983*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp983*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp983/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp983/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[9] = dwdp983;
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp983*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp983*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 107:
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp994*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp994*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp994/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp994/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[9] = dwdp994;
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp994*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp994*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 108:
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1005*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1005*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1005/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1005/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[9] = dwdp1005;
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1005*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1005*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 109:
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1016*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1016*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1016/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1016/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[9] = dwdp1016;
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1016*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1016*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 110:
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1027*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1027*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1027/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1027/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[9] = dwdp1027;
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1027*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1027*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 111:
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1038*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1038*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1038/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1038/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[9] = dwdp1038;
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1038*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1038*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 112:
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1049*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1049*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1049/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1049/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[9] = dwdp1049;
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1049*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1049*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
    }
}

} // namespace amici
} // namespace model_vaccination