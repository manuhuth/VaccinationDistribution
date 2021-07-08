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
        case 80:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp556*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp556*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp557*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp557*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp556/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp556/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp557/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp557/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[8] = dwdp556;
            dydp[9] = dwdp557;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp556*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp557*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp556*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp557*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 81:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp578*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp578*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp579*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp579*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp578/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp578/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp579/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp579/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[8] = dwdp578;
            dydp[9] = dwdp579;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp578*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp579*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp578*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp579*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 82:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp600*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp600*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp601*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp601*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp600/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp600/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp601/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp601/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[8] = dwdp600;
            dydp[9] = dwdp601;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp600*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp601*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp600*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp601*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 83:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp622*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp622*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp623*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp623*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp622/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp622/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp623/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp623/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[8] = dwdp622;
            dydp[9] = dwdp623;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp622*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp623*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp622*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp623*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 84:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp644*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp644*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp645*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp645*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp644/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp644/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp645/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp645/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[8] = dwdp644;
            dydp[9] = dwdp645;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp644*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp645*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp644*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp645*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 85:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp666*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp666*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp667*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp667*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp666/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp666/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp667/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp667/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[8] = dwdp666;
            dydp[9] = dwdp667;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp666*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp667*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp666*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp667*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 86:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp688*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp688*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp689*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp689*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp688/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp688/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp689/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp689/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[8] = dwdp688;
            dydp[9] = dwdp689;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp688*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp689*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp688*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp689*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 87:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp710*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp710*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp711*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp711*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp710/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp710/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp711/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp711/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[8] = dwdp710;
            dydp[9] = dwdp711;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp710*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp711*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp710*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp711*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 88:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp732*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp732*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp733*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp733*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp732/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp732/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp733/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp733/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[8] = dwdp732;
            dydp[9] = dwdp733;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp732*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp733*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp732*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp733*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 89:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp754*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp754*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp755*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp755*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp754/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp754/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp755/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp755/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[8] = dwdp754;
            dydp[9] = dwdp755;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp754*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp755*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp754*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp755*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 90:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp776*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp776*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp777*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp777*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp776/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp776/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp777/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp777/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[8] = dwdp776;
            dydp[9] = dwdp777;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp776*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp777*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp776*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp777*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 91:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp798*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp798*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp798/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp798/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[8] = dwdp798;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp798*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp798*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            break;
        case 92:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp809*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp809*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp809/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp809/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[8] = dwdp809;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp809*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp809*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            break;
        case 93:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp820*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp820*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp820/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp820/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[8] = dwdp820;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp820*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp820*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            break;
        case 94:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp831*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp831*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp831/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp831/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[8] = dwdp831;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp831*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp831*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            break;
        case 95:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp842*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp842*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp842/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp842/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[8] = dwdp842;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp842*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp842*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            break;
        case 96:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp853*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp853*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp853/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp853/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[8] = dwdp853;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp853*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp853*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            break;
        case 97:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp864*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp864*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp864/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp864/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[8] = dwdp864;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp864*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp864*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            break;
        case 98:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp875*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp875*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp875/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp875/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[8] = dwdp875;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp875*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp875*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            break;
        case 99:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp886*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp886*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp886/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp886/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[8] = dwdp886;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp886*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp886*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            break;
        case 100:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp897*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp897*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp897/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp897/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[8] = dwdp897;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp897*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp897*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            break;
        case 101:
            dydp[0] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp908*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[1] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp908*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[4] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp908/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[5] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp908/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[8] = dwdp908;
            dydp[10] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp908*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            dydp[12] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac1)*dwdp908*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac1) + 1, 2);
            break;
        case 102:
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp919*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp919*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp919/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp919/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[9] = dwdp919;
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp919*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp919*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 103:
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp930*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp930*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp930/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp930/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[9] = dwdp930;
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp930*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp930*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 104:
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp941*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp941*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp941/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp941/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[9] = dwdp941;
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp941*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp941*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 105:
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp952*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp952*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp952/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp952/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[9] = dwdp952;
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp952*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp952*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 106:
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp963*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp963*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp963/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp963/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[9] = dwdp963;
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp963*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp963*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 107:
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp974*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp974*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp974/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp974/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[9] = dwdp974;
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp974*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp974*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 108:
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp985*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp985*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp985/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp985/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[9] = dwdp985;
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp985*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp985*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 109:
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp996*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp996*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp996/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp996/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[9] = dwdp996;
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp996*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp996*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 110:
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1007*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1007*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1007/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1007/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[9] = dwdp1007;
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1007*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1007*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 111:
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1018*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1018*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1018/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1018/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[9] = dwdp1018;
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1018*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1018*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
        case 112:
            dydp[2] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1029*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
            dydp[3] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1029*((t >= xx0 && t < xx1) ? (
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
))))))))))/(std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0));
            dydp[6] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1029/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[7] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1029/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[9] = dwdp1029;
            dydp[11] = 0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1029*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            dydp[13] = -0.99999999983112664*std::pow(0.36787944123356736, spline_countryA_vac2)*dwdp1029*((t >= xx0 && t < xx1) ? (
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
))))))))))/std::pow(std::pow(0.36787944123356736, spline_countryA_vac2) + 1, 2);
            break;
    }
}

} // namespace amici
} // namespace model_vaccination