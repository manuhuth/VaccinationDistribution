#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "tcl.h"
#include "p.h"
#include "k.h"
#include "w.h"
#include "spl.h"
#include "sspl.h"
#include "x.h"

namespace amici {
namespace model_vaccination_piecewise {

void dwdp_vaccination_piecewise(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp,const realtype *spl, const realtype *sspl){
    dwdp[0] = infectious_countryA_vac0_virW*(1 - prob_deceasing);
    dwdp[1] = infectious_countryA_vac0_virW*prob_deceasing;
    dwdp[2] = infectious_countryA_vac0_virM*(1 - prob_deceasing);
    dwdp[3] = infectious_countryA_vac0_virM*prob_deceasing;
    dwdp[4] = infectious_countryA_vac1_virW*(-prob_deceasing*(1 - omega_vac1_virW) + 1);
    dwdp[5] = infectious_countryA_vac1_virW*prob_deceasing*(1 - omega_vac1_virW);
    dwdp[6] = infectious_countryA_vac1_virM*(-prob_deceasing*(1 - omega_vac1_virM) + 1);
    dwdp[7] = infectious_countryA_vac1_virM*prob_deceasing*(1 - omega_vac1_virM);
    dwdp[8] = infectious_countryA_vac2_virW*(-prob_deceasing*(1 - omega_vac2_virW) + 1);
    dwdp[9] = infectious_countryA_vac2_virW*prob_deceasing*(1 - omega_vac2_virW);
    dwdp[10] = infectious_countryA_vac2_virM*(-prob_deceasing*(1 - omega_vac2_virM) + 1);
    dwdp[11] = infectious_countryA_vac2_virM*prob_deceasing*(1 - omega_vac2_virM);
    dwdp[12] = infectious_countryB_vac0_virW*(1 - prob_deceasing);
    dwdp[13] = infectious_countryB_vac0_virW*prob_deceasing;
    dwdp[14] = infectious_countryB_vac0_virM*(1 - prob_deceasing);
    dwdp[15] = infectious_countryB_vac0_virM*prob_deceasing;
    dwdp[16] = infectious_countryB_vac1_virW*(-prob_deceasing*(1 - omega_vac1_virW) + 1);
    dwdp[17] = infectious_countryB_vac1_virW*prob_deceasing*(1 - omega_vac1_virW);
    dwdp[18] = infectious_countryB_vac1_virM*(-prob_deceasing*(1 - omega_vac1_virM) + 1);
    dwdp[19] = infectious_countryB_vac1_virM*prob_deceasing*(1 - omega_vac1_virM);
    dwdp[20] = infectious_countryB_vac2_virW*(-prob_deceasing*(1 - omega_vac2_virW) + 1);
    dwdp[21] = infectious_countryB_vac2_virW*prob_deceasing*(1 - omega_vac2_virW);
    dwdp[22] = infectious_countryB_vac2_virM*(-prob_deceasing*(1 - omega_vac2_virM) + 1);
    dwdp[23] = infectious_countryB_vac2_virM*prob_deceasing*(1 - omega_vac2_virM);
    dwdp[24] = -infectious_countryA_vac0_virW*lambda1;
    dwdp[25] = infectious_countryA_vac0_virW*lambda1;
    dwdp[26] = -infectious_countryA_vac0_virM*lambda1;
    dwdp[27] = infectious_countryA_vac0_virM*lambda1;
    dwdp[28] = infectious_countryA_vac1_virW*lambda1*(omega_vac1_virW - 1);
    dwdp[29] = infectious_countryA_vac1_virW*lambda1*(1 - omega_vac1_virW);
    dwdp[30] = infectious_countryA_vac1_virM*lambda1*(omega_vac1_virM - 1);
    dwdp[31] = infectious_countryA_vac1_virM*lambda1*(1 - omega_vac1_virM);
    dwdp[32] = infectious_countryA_vac2_virW*lambda1*(omega_vac2_virW - 1);
    dwdp[33] = infectious_countryA_vac2_virW*lambda1*(1 - omega_vac2_virW);
    dwdp[34] = infectious_countryA_vac2_virM*lambda1*(omega_vac2_virM - 1);
    dwdp[35] = infectious_countryA_vac2_virM*lambda1*(1 - omega_vac2_virM);
    dwdp[36] = -infectious_countryB_vac0_virW*lambda1;
    dwdp[37] = infectious_countryB_vac0_virW*lambda1;
    dwdp[38] = -infectious_countryB_vac0_virM*lambda1;
    dwdp[39] = infectious_countryB_vac0_virM*lambda1;
    dwdp[40] = infectious_countryB_vac1_virW*lambda1*(omega_vac1_virW - 1);
    dwdp[41] = infectious_countryB_vac1_virW*lambda1*(1 - omega_vac1_virW);
    dwdp[42] = infectious_countryB_vac1_virM*lambda1*(omega_vac1_virM - 1);
    dwdp[43] = infectious_countryB_vac1_virM*lambda1*(1 - omega_vac1_virM);
    dwdp[44] = infectious_countryB_vac2_virW*lambda1*(omega_vac2_virW - 1);
    dwdp[45] = infectious_countryB_vac2_virW*lambda1*(1 - omega_vac2_virW);
    dwdp[46] = infectious_countryB_vac2_virM*lambda1*(omega_vac2_virM - 1);
    dwdp[47] = infectious_countryB_vac2_virM*lambda1*(1 - omega_vac2_virM);
    dwdp[48] = -beta*eta_virW*infectious_countryA_vac1_virW*susceptible_countryA_vac0*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[49] = -beta*eta_virW*infectious_countryA_vac1_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[50] = -beta*eta_virW*infectious_countryA_vac1_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[51] = -beta*distance_countryA_countryB*eta_virW*infectious_countryA_vac1_virW*susceptible_countryB_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[52] = -beta*distance_countryA_countryB*eta_virW*infectious_countryA_vac1_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[53] = -beta*distance_countryA_countryB*eta_virW*infectious_countryA_vac1_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[54] = -beta*eta_virM*infectious_countryA_vac1_virM*susceptible_countryA_vac0*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[55] = -beta*eta_virM*infectious_countryA_vac1_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[56] = -beta*eta_virM*infectious_countryA_vac1_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[57] = -beta*distance_countryA_countryB*eta_virM*infectious_countryA_vac1_virM*susceptible_countryB_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[58] = -beta*distance_countryA_countryB*eta_virM*infectious_countryA_vac1_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[59] = -beta*distance_countryA_countryB*eta_virM*infectious_countryA_vac1_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[60] = -beta*eta_virW*infectious_countryA_vac2_virW*susceptible_countryA_vac0*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[61] = -beta*eta_virW*infectious_countryA_vac2_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[62] = -beta*eta_virW*infectious_countryA_vac2_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[63] = -beta*distance_countryA_countryB*eta_virW*infectious_countryA_vac2_virW*susceptible_countryB_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[64] = -beta*distance_countryA_countryB*eta_virW*infectious_countryA_vac2_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[65] = -beta*distance_countryA_countryB*eta_virW*infectious_countryA_vac2_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[66] = -beta*eta_virM*infectious_countryA_vac2_virM*susceptible_countryA_vac0*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[67] = -beta*eta_virM*infectious_countryA_vac2_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[68] = -beta*eta_virM*infectious_countryA_vac2_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[69] = -beta*distance_countryA_countryB*eta_virM*infectious_countryA_vac2_virM*susceptible_countryB_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[70] = -beta*distance_countryA_countryB*eta_virM*infectious_countryA_vac2_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[71] = -beta*distance_countryA_countryB*eta_virM*infectious_countryA_vac2_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[72] = -beta*distance_countryB_countryA*eta_virW*infectious_countryB_vac1_virW*susceptible_countryA_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[73] = -beta*distance_countryB_countryA*eta_virW*infectious_countryB_vac1_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[74] = -beta*distance_countryB_countryA*eta_virW*infectious_countryB_vac1_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[75] = -beta*eta_virW*infectious_countryB_vac1_virW*susceptible_countryB_vac0*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[76] = -beta*eta_virW*infectious_countryB_vac1_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[77] = -beta*eta_virW*infectious_countryB_vac1_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[78] = -beta*distance_countryB_countryA*eta_virM*infectious_countryB_vac1_virM*susceptible_countryA_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[79] = -beta*distance_countryB_countryA*eta_virM*infectious_countryB_vac1_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[80] = -beta*distance_countryB_countryA*eta_virM*infectious_countryB_vac1_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[81] = -beta*eta_virM*infectious_countryB_vac1_virM*susceptible_countryB_vac0*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[82] = -beta*eta_virM*infectious_countryB_vac1_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[83] = -beta*eta_virM*infectious_countryB_vac1_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[84] = -beta*distance_countryB_countryA*eta_virW*infectious_countryB_vac2_virW*susceptible_countryA_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[85] = -beta*distance_countryB_countryA*eta_virW*infectious_countryB_vac2_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[86] = -beta*distance_countryB_countryA*eta_virW*infectious_countryB_vac2_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[87] = -beta*eta_virW*infectious_countryB_vac2_virW*susceptible_countryB_vac0*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[88] = -beta*eta_virW*infectious_countryB_vac2_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[89] = -beta*eta_virW*infectious_countryB_vac2_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[90] = -beta*distance_countryB_countryA*eta_virM*infectious_countryB_vac2_virM*susceptible_countryA_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[91] = -beta*distance_countryB_countryA*eta_virM*infectious_countryB_vac2_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[92] = -beta*distance_countryB_countryA*eta_virM*infectious_countryB_vac2_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[93] = -beta*eta_virM*infectious_countryB_vac2_virM*susceptible_countryB_vac0*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[94] = -beta*eta_virM*infectious_countryB_vac2_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[95] = -beta*eta_virM*infectious_countryB_vac2_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[96] = eta_virW*infectious_countryA_vac0_virW*susceptible_countryA_vac0*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[97] = eta_virW*infectious_countryA_vac0_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[98] = eta_virW*infectious_countryA_vac0_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[99] = distance_countryA_countryB*eta_virW*infectious_countryA_vac0_virW*susceptible_countryB_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[100] = distance_countryA_countryB*eta_virW*infectious_countryA_vac0_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[101] = distance_countryA_countryB*eta_virW*infectious_countryA_vac0_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[102] = eta_virM*infectious_countryA_vac0_virM*susceptible_countryA_vac0*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[103] = eta_virM*infectious_countryA_vac0_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[104] = eta_virM*infectious_countryA_vac0_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[105] = distance_countryA_countryB*eta_virM*infectious_countryA_vac0_virM*susceptible_countryB_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[106] = distance_countryA_countryB*eta_virM*infectious_countryA_vac0_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[107] = distance_countryA_countryB*eta_virM*infectious_countryA_vac0_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[108] = eta_virW*infectious_countryA_vac1_virW*susceptible_countryA_vac0*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[109] = eta_virW*infectious_countryA_vac1_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[110] = eta_virW*infectious_countryA_vac1_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[111] = distance_countryA_countryB*eta_virW*infectious_countryA_vac1_virW*susceptible_countryB_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[112] = distance_countryA_countryB*eta_virW*infectious_countryA_vac1_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[113] = distance_countryA_countryB*eta_virW*infectious_countryA_vac1_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[114] = eta_virM*infectious_countryA_vac1_virM*susceptible_countryA_vac0*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[115] = eta_virM*infectious_countryA_vac1_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[116] = eta_virM*infectious_countryA_vac1_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[117] = distance_countryA_countryB*eta_virM*infectious_countryA_vac1_virM*susceptible_countryB_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[118] = distance_countryA_countryB*eta_virM*infectious_countryA_vac1_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[119] = distance_countryA_countryB*eta_virM*infectious_countryA_vac1_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[120] = eta_virW*infectious_countryA_vac2_virW*susceptible_countryA_vac0*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[121] = eta_virW*infectious_countryA_vac2_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[122] = eta_virW*infectious_countryA_vac2_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[123] = distance_countryA_countryB*eta_virW*infectious_countryA_vac2_virW*susceptible_countryB_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[124] = distance_countryA_countryB*eta_virW*infectious_countryA_vac2_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[125] = distance_countryA_countryB*eta_virW*infectious_countryA_vac2_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[126] = eta_virM*infectious_countryA_vac2_virM*susceptible_countryA_vac0*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[127] = eta_virM*infectious_countryA_vac2_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[128] = eta_virM*infectious_countryA_vac2_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[129] = distance_countryA_countryB*eta_virM*infectious_countryA_vac2_virM*susceptible_countryB_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[130] = distance_countryA_countryB*eta_virM*infectious_countryA_vac2_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[131] = distance_countryA_countryB*eta_virM*infectious_countryA_vac2_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[132] = distance_countryB_countryA*eta_virW*infectious_countryB_vac0_virW*susceptible_countryA_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[133] = distance_countryB_countryA*eta_virW*infectious_countryB_vac0_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[134] = distance_countryB_countryA*eta_virW*infectious_countryB_vac0_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[135] = eta_virW*infectious_countryB_vac0_virW*susceptible_countryB_vac0*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[136] = eta_virW*infectious_countryB_vac0_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[137] = eta_virW*infectious_countryB_vac0_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[138] = distance_countryB_countryA*eta_virM*infectious_countryB_vac0_virM*susceptible_countryA_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[139] = distance_countryB_countryA*eta_virM*infectious_countryB_vac0_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[140] = distance_countryB_countryA*eta_virM*infectious_countryB_vac0_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[141] = eta_virM*infectious_countryB_vac0_virM*susceptible_countryB_vac0*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[142] = eta_virM*infectious_countryB_vac0_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[143] = eta_virM*infectious_countryB_vac0_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[144] = distance_countryB_countryA*eta_virW*infectious_countryB_vac1_virW*susceptible_countryA_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[145] = distance_countryB_countryA*eta_virW*infectious_countryB_vac1_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[146] = distance_countryB_countryA*eta_virW*infectious_countryB_vac1_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[147] = eta_virW*infectious_countryB_vac1_virW*susceptible_countryB_vac0*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[148] = eta_virW*infectious_countryB_vac1_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[149] = eta_virW*infectious_countryB_vac1_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[150] = distance_countryB_countryA*eta_virM*infectious_countryB_vac1_virM*susceptible_countryA_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[151] = distance_countryB_countryA*eta_virM*infectious_countryB_vac1_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[152] = distance_countryB_countryA*eta_virM*infectious_countryB_vac1_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[153] = eta_virM*infectious_countryB_vac1_virM*susceptible_countryB_vac0*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[154] = eta_virM*infectious_countryB_vac1_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[155] = eta_virM*infectious_countryB_vac1_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[156] = distance_countryB_countryA*eta_virW*infectious_countryB_vac2_virW*susceptible_countryA_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[157] = distance_countryB_countryA*eta_virW*infectious_countryB_vac2_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[158] = distance_countryB_countryA*eta_virW*infectious_countryB_vac2_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[159] = eta_virW*infectious_countryB_vac2_virW*susceptible_countryB_vac0*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[160] = eta_virW*infectious_countryB_vac2_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[161] = eta_virW*infectious_countryB_vac2_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[162] = distance_countryB_countryA*eta_virM*infectious_countryB_vac2_virM*susceptible_countryA_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[163] = distance_countryB_countryA*eta_virM*infectious_countryB_vac2_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[164] = distance_countryB_countryA*eta_virM*infectious_countryB_vac2_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[165] = eta_virM*infectious_countryB_vac2_virM*susceptible_countryB_vac0*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[166] = eta_virM*infectious_countryB_vac2_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[167] = eta_virM*infectious_countryB_vac2_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[168] = infectious_countryA_vac1_virW*lambda1*prob_deceasing;
    dwdp[169] = -infectious_countryA_vac1_virW*lambda1*prob_deceasing;
    dwdp[170] = infectious_countryB_vac1_virW*lambda1*prob_deceasing;
    dwdp[171] = -infectious_countryB_vac1_virW*lambda1*prob_deceasing;
    dwdp[172] = -beta*eta_virW*infectious_countryA_vac0_virW*susceptible_countryA_vac1*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[173] = -beta*distance_countryA_countryB*eta_virW*infectious_countryA_vac0_virW*susceptible_countryB_vac1/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[174] = -beta*eta_virW*infectious_countryA_vac1_virW*susceptible_countryA_vac1*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[175] = -beta*distance_countryA_countryB*eta_virW*infectious_countryA_vac1_virW*susceptible_countryB_vac1*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[176] = -beta*eta_virW*infectious_countryA_vac2_virW*susceptible_countryA_vac1*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[177] = -beta*distance_countryA_countryB*eta_virW*infectious_countryA_vac2_virW*susceptible_countryB_vac1*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[178] = -beta*distance_countryB_countryA*eta_virW*infectious_countryB_vac0_virW*susceptible_countryA_vac1/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[179] = -beta*eta_virW*infectious_countryB_vac0_virW*susceptible_countryB_vac1*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[180] = -beta*distance_countryB_countryA*eta_virW*infectious_countryB_vac1_virW*susceptible_countryA_vac1*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[181] = -beta*eta_virW*infectious_countryB_vac1_virW*susceptible_countryB_vac1*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[182] = -beta*distance_countryB_countryA*eta_virW*infectious_countryB_vac2_virW*susceptible_countryA_vac1*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[183] = -beta*eta_virW*infectious_countryB_vac2_virW*susceptible_countryB_vac1*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[184] = infectious_countryA_vac2_virW*lambda1*prob_deceasing;
    dwdp[185] = -infectious_countryA_vac2_virW*lambda1*prob_deceasing;
    dwdp[186] = infectious_countryB_vac2_virW*lambda1*prob_deceasing;
    dwdp[187] = -infectious_countryB_vac2_virW*lambda1*prob_deceasing;
    dwdp[188] = -beta*eta_virW*infectious_countryA_vac0_virW*susceptible_countryA_vac2*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[189] = -beta*distance_countryA_countryB*eta_virW*infectious_countryA_vac0_virW*susceptible_countryB_vac2/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[190] = -beta*eta_virW*infectious_countryA_vac1_virW*susceptible_countryA_vac2*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[191] = -beta*distance_countryA_countryB*eta_virW*infectious_countryA_vac1_virW*susceptible_countryB_vac2*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[192] = -beta*eta_virW*infectious_countryA_vac2_virW*susceptible_countryA_vac2*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[193] = -beta*distance_countryA_countryB*eta_virW*infectious_countryA_vac2_virW*susceptible_countryB_vac2*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[194] = -beta*distance_countryB_countryA*eta_virW*infectious_countryB_vac0_virW*susceptible_countryA_vac2/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[195] = -beta*eta_virW*infectious_countryB_vac0_virW*susceptible_countryB_vac2*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[196] = -beta*distance_countryB_countryA*eta_virW*infectious_countryB_vac1_virW*susceptible_countryA_vac2*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[197] = -beta*eta_virW*infectious_countryB_vac1_virW*susceptible_countryB_vac2*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[198] = -beta*distance_countryB_countryA*eta_virW*infectious_countryB_vac2_virW*susceptible_countryA_vac2*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[199] = -beta*eta_virW*infectious_countryB_vac2_virW*susceptible_countryB_vac2*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[200] = infectious_countryA_vac1_virM*lambda1*prob_deceasing;
    dwdp[201] = -infectious_countryA_vac1_virM*lambda1*prob_deceasing;
    dwdp[202] = infectious_countryB_vac1_virM*lambda1*prob_deceasing;
    dwdp[203] = -infectious_countryB_vac1_virM*lambda1*prob_deceasing;
    dwdp[204] = -beta*eta_virM*infectious_countryA_vac0_virM*susceptible_countryA_vac1*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[205] = -beta*distance_countryA_countryB*eta_virM*infectious_countryA_vac0_virM*susceptible_countryB_vac1/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[206] = -beta*eta_virM*infectious_countryA_vac1_virM*susceptible_countryA_vac1*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[207] = -beta*distance_countryA_countryB*eta_virM*infectious_countryA_vac1_virM*susceptible_countryB_vac1*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[208] = -beta*eta_virM*infectious_countryA_vac2_virM*susceptible_countryA_vac1*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[209] = -beta*distance_countryA_countryB*eta_virM*infectious_countryA_vac2_virM*susceptible_countryB_vac1*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[210] = -beta*distance_countryB_countryA*eta_virM*infectious_countryB_vac0_virM*susceptible_countryA_vac1/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[211] = -beta*eta_virM*infectious_countryB_vac0_virM*susceptible_countryB_vac1*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[212] = -beta*distance_countryB_countryA*eta_virM*infectious_countryB_vac1_virM*susceptible_countryA_vac1*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[213] = -beta*eta_virM*infectious_countryB_vac1_virM*susceptible_countryB_vac1*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[214] = -beta*distance_countryB_countryA*eta_virM*infectious_countryB_vac2_virM*susceptible_countryA_vac1*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[215] = -beta*eta_virM*infectious_countryB_vac2_virM*susceptible_countryB_vac1*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[216] = infectious_countryA_vac2_virM*lambda1*prob_deceasing;
    dwdp[217] = -infectious_countryA_vac2_virM*lambda1*prob_deceasing;
    dwdp[218] = infectious_countryB_vac2_virM*lambda1*prob_deceasing;
    dwdp[219] = -infectious_countryB_vac2_virM*lambda1*prob_deceasing;
    dwdp[220] = -beta*eta_virM*infectious_countryA_vac0_virM*susceptible_countryA_vac2*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[221] = -beta*distance_countryA_countryB*eta_virM*infectious_countryA_vac0_virM*susceptible_countryB_vac2/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[222] = -beta*eta_virM*infectious_countryA_vac1_virM*susceptible_countryA_vac2*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[223] = -beta*distance_countryA_countryB*eta_virM*infectious_countryA_vac1_virM*susceptible_countryB_vac2*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[224] = -beta*eta_virM*infectious_countryA_vac2_virM*susceptible_countryA_vac2*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[225] = -beta*distance_countryA_countryB*eta_virM*infectious_countryA_vac2_virM*susceptible_countryB_vac2*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[226] = -beta*distance_countryB_countryA*eta_virM*infectious_countryB_vac0_virM*susceptible_countryA_vac2/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[227] = -beta*eta_virM*infectious_countryB_vac0_virM*susceptible_countryB_vac2*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[228] = -beta*distance_countryB_countryA*eta_virM*infectious_countryB_vac1_virM*susceptible_countryA_vac2*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[229] = -beta*eta_virM*infectious_countryB_vac1_virM*susceptible_countryB_vac2*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[230] = -beta*distance_countryB_countryA*eta_virM*infectious_countryB_vac2_virM*susceptible_countryA_vac2*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[231] = -beta*eta_virM*infectious_countryB_vac2_virM*susceptible_countryB_vac2*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[232] = beta*infectious_countryA_vac0_virW*susceptible_countryA_vac0*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[233] = beta*infectious_countryA_vac0_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[234] = beta*infectious_countryA_vac0_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[235] = beta*distance_countryA_countryB*infectious_countryA_vac0_virW*susceptible_countryB_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[236] = beta*distance_countryA_countryB*infectious_countryA_vac0_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[237] = beta*distance_countryA_countryB*infectious_countryA_vac0_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[238] = beta*infectious_countryA_vac1_virW*susceptible_countryA_vac0*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[239] = beta*infectious_countryA_vac1_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[240] = beta*infectious_countryA_vac1_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[241] = beta*distance_countryA_countryB*infectious_countryA_vac1_virW*susceptible_countryB_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[242] = beta*distance_countryA_countryB*infectious_countryA_vac1_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[243] = beta*distance_countryA_countryB*infectious_countryA_vac1_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[244] = beta*infectious_countryA_vac2_virW*susceptible_countryA_vac0*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[245] = beta*infectious_countryA_vac2_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[246] = beta*infectious_countryA_vac2_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[247] = beta*distance_countryA_countryB*infectious_countryA_vac2_virW*susceptible_countryB_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[248] = beta*distance_countryA_countryB*infectious_countryA_vac2_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[249] = beta*distance_countryA_countryB*infectious_countryA_vac2_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[250] = beta*distance_countryB_countryA*infectious_countryB_vac0_virW*susceptible_countryA_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[251] = beta*distance_countryB_countryA*infectious_countryB_vac0_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[252] = beta*distance_countryB_countryA*infectious_countryB_vac0_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[253] = beta*infectious_countryB_vac0_virW*susceptible_countryB_vac0*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[254] = beta*infectious_countryB_vac0_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[255] = beta*infectious_countryB_vac0_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[256] = beta*distance_countryB_countryA*infectious_countryB_vac1_virW*susceptible_countryA_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[257] = beta*distance_countryB_countryA*infectious_countryB_vac1_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[258] = beta*distance_countryB_countryA*infectious_countryB_vac1_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[259] = beta*infectious_countryB_vac1_virW*susceptible_countryB_vac0*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[260] = beta*infectious_countryB_vac1_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[261] = beta*infectious_countryB_vac1_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[262] = beta*distance_countryB_countryA*infectious_countryB_vac2_virW*susceptible_countryA_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[263] = beta*distance_countryB_countryA*infectious_countryB_vac2_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[264] = beta*distance_countryB_countryA*infectious_countryB_vac2_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[265] = beta*infectious_countryB_vac2_virW*susceptible_countryB_vac0*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[266] = beta*infectious_countryB_vac2_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[267] = beta*infectious_countryB_vac2_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[268] = beta*infectious_countryA_vac0_virM*susceptible_countryA_vac0*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[269] = beta*infectious_countryA_vac0_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[270] = beta*infectious_countryA_vac0_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[271] = beta*distance_countryA_countryB*infectious_countryA_vac0_virM*susceptible_countryB_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[272] = beta*distance_countryA_countryB*infectious_countryA_vac0_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[273] = beta*distance_countryA_countryB*infectious_countryA_vac0_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[274] = beta*infectious_countryA_vac1_virM*susceptible_countryA_vac0*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[275] = beta*infectious_countryA_vac1_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[276] = beta*infectious_countryA_vac1_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[277] = beta*distance_countryA_countryB*infectious_countryA_vac1_virM*susceptible_countryB_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[278] = beta*distance_countryA_countryB*infectious_countryA_vac1_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[279] = beta*distance_countryA_countryB*infectious_countryA_vac1_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[280] = beta*infectious_countryA_vac2_virM*susceptible_countryA_vac0*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[281] = beta*infectious_countryA_vac2_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[282] = beta*infectious_countryA_vac2_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    dwdp[283] = beta*distance_countryA_countryB*infectious_countryA_vac2_virM*susceptible_countryB_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[284] = beta*distance_countryA_countryB*infectious_countryA_vac2_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[285] = beta*distance_countryA_countryB*infectious_countryA_vac2_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[286] = beta*distance_countryB_countryA*infectious_countryB_vac0_virM*susceptible_countryA_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[287] = beta*distance_countryB_countryA*infectious_countryB_vac0_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[288] = beta*distance_countryB_countryA*infectious_countryB_vac0_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[289] = beta*infectious_countryB_vac0_virM*susceptible_countryB_vac0*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[290] = beta*infectious_countryB_vac0_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[291] = beta*infectious_countryB_vac0_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[292] = beta*distance_countryB_countryA*infectious_countryB_vac1_virM*susceptible_countryA_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[293] = beta*distance_countryB_countryA*infectious_countryB_vac1_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[294] = beta*distance_countryB_countryA*infectious_countryB_vac1_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[295] = beta*infectious_countryB_vac1_virM*susceptible_countryB_vac0*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[296] = beta*infectious_countryB_vac1_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[297] = beta*infectious_countryB_vac1_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[298] = beta*distance_countryB_countryA*infectious_countryB_vac2_virM*susceptible_countryA_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[299] = beta*distance_countryB_countryA*infectious_countryB_vac2_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[300] = beta*distance_countryB_countryA*infectious_countryB_vac2_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[301] = beta*infectious_countryB_vac2_virM*susceptible_countryB_vac0*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[302] = beta*infectious_countryB_vac2_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[303] = beta*infectious_countryB_vac2_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[304] = susceptible_countryA_vac0*((t < 0) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[305] = susceptible_countryB_vac0*((t < 0) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[306] = infectious_countryA_vac0_virW*((t < 0) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[307] = infectious_countryA_vac0_virM*((t < 0) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[308] = infectious_countryB_vac0_virW*((t < 0) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[309] = infectious_countryB_vac0_virM*((t < 0) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[310] = recovered_countryA_vac0_virW*((t < 0) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[311] = recovered_countryA_vac0_virM*((t < 0) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[312] = recovered_countryB_vac0_virW*((t < 0) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[313] = recovered_countryB_vac0_virM*((t < 0) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[314] = susceptible_countryA_vac0*((t < 14) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[315] = susceptible_countryB_vac0*((t < 14) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[316] = infectious_countryA_vac0_virW*((t < 14) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[317] = infectious_countryA_vac0_virM*((t < 14) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[318] = infectious_countryB_vac0_virW*((t < 14) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[319] = infectious_countryB_vac0_virM*((t < 14) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[320] = recovered_countryA_vac0_virW*((t < 14) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[321] = recovered_countryA_vac0_virM*((t < 14) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[322] = recovered_countryB_vac0_virW*((t < 14) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[323] = recovered_countryB_vac0_virM*((t < 14) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[324] = susceptible_countryA_vac0*((t < 28) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[325] = susceptible_countryB_vac0*((t < 28) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[326] = infectious_countryA_vac0_virW*((t < 28) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[327] = infectious_countryA_vac0_virM*((t < 28) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[328] = infectious_countryB_vac0_virW*((t < 28) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[329] = infectious_countryB_vac0_virM*((t < 28) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[330] = recovered_countryA_vac0_virW*((t < 28) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[331] = recovered_countryA_vac0_virM*((t < 28) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[332] = recovered_countryB_vac0_virW*((t < 28) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[333] = recovered_countryB_vac0_virM*((t < 28) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[334] = susceptible_countryA_vac0*((t < 42) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[335] = susceptible_countryB_vac0*((t < 42) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[336] = infectious_countryA_vac0_virW*((t < 42) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[337] = infectious_countryA_vac0_virM*((t < 42) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[338] = infectious_countryB_vac0_virW*((t < 42) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[339] = infectious_countryB_vac0_virM*((t < 42) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[340] = recovered_countryA_vac0_virW*((t < 42) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[341] = recovered_countryA_vac0_virM*((t < 42) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[342] = recovered_countryB_vac0_virW*((t < 42) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[343] = recovered_countryB_vac0_virM*((t < 42) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[344] = susceptible_countryA_vac0*((t < 56) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[345] = susceptible_countryB_vac0*((t < 56) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[346] = infectious_countryA_vac0_virW*((t < 56) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[347] = infectious_countryA_vac0_virM*((t < 56) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[348] = infectious_countryB_vac0_virW*((t < 56) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[349] = infectious_countryB_vac0_virM*((t < 56) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[350] = recovered_countryA_vac0_virW*((t < 56) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[351] = recovered_countryA_vac0_virM*((t < 56) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[352] = recovered_countryB_vac0_virW*((t < 56) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[353] = recovered_countryB_vac0_virM*((t < 56) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[354] = susceptible_countryA_vac0*((t < 70) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[355] = susceptible_countryB_vac0*((t < 70) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[356] = infectious_countryA_vac0_virW*((t < 70) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[357] = infectious_countryA_vac0_virM*((t < 70) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[358] = infectious_countryB_vac0_virW*((t < 70) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[359] = infectious_countryB_vac0_virM*((t < 70) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[360] = recovered_countryA_vac0_virW*((t < 70) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[361] = recovered_countryA_vac0_virM*((t < 70) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[362] = recovered_countryB_vac0_virW*((t < 70) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[363] = recovered_countryB_vac0_virM*((t < 70) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[364] = susceptible_countryA_vac0*((t < 84) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[365] = susceptible_countryB_vac0*((t < 84) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[366] = infectious_countryA_vac0_virW*((t < 84) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[367] = infectious_countryA_vac0_virM*((t < 84) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[368] = infectious_countryB_vac0_virW*((t < 84) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[369] = infectious_countryB_vac0_virM*((t < 84) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[370] = recovered_countryA_vac0_virW*((t < 84) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[371] = recovered_countryA_vac0_virM*((t < 84) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[372] = recovered_countryB_vac0_virW*((t < 84) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[373] = recovered_countryB_vac0_virM*((t < 84) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[374] = susceptible_countryA_vac0*((t < 98) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[375] = susceptible_countryB_vac0*((t < 98) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[376] = infectious_countryA_vac0_virW*((t < 98) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[377] = infectious_countryA_vac0_virM*((t < 98) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[378] = infectious_countryB_vac0_virW*((t < 98) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[379] = infectious_countryB_vac0_virM*((t < 98) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[380] = recovered_countryA_vac0_virW*((t < 98) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[381] = recovered_countryA_vac0_virM*((t < 98) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[382] = recovered_countryB_vac0_virW*((t < 98) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[383] = recovered_countryB_vac0_virM*((t < 98) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[384] = susceptible_countryA_vac0*((t < 112) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[385] = susceptible_countryB_vac0*((t < 112) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[386] = infectious_countryA_vac0_virW*((t < 112) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[387] = infectious_countryA_vac0_virM*((t < 112) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[388] = infectious_countryB_vac0_virW*((t < 112) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[389] = infectious_countryB_vac0_virM*((t < 112) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[390] = recovered_countryA_vac0_virW*((t < 112) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[391] = recovered_countryA_vac0_virM*((t < 112) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[392] = recovered_countryB_vac0_virW*((t < 112) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[393] = recovered_countryB_vac0_virM*((t < 112) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[394] = susceptible_countryA_vac0*((t < 126) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[395] = susceptible_countryB_vac0*((t < 126) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[396] = infectious_countryA_vac0_virW*((t < 126) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[397] = infectious_countryA_vac0_virM*((t < 126) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[398] = infectious_countryB_vac0_virW*((t < 126) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[399] = infectious_countryB_vac0_virM*((t < 126) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[400] = recovered_countryA_vac0_virW*((t < 126) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[401] = recovered_countryA_vac0_virM*((t < 126) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[402] = recovered_countryB_vac0_virW*((t < 126) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[403] = recovered_countryB_vac0_virM*((t < 126) ? (
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
   vaccine_supply_par_vac1_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[404] = susceptible_countryA_vac0*((t < 0) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[405] = susceptible_countryB_vac0*((t < 0) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[406] = infectious_countryA_vac0_virW*((t < 0) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[407] = infectious_countryA_vac0_virM*((t < 0) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[408] = infectious_countryB_vac0_virW*((t < 0) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[409] = infectious_countryB_vac0_virM*((t < 0) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[410] = recovered_countryA_vac0_virW*((t < 0) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[411] = recovered_countryA_vac0_virM*((t < 0) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[412] = recovered_countryB_vac0_virW*((t < 0) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[413] = recovered_countryB_vac0_virM*((t < 0) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[414] = susceptible_countryA_vac0*((t < 14) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[415] = susceptible_countryB_vac0*((t < 14) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[416] = infectious_countryA_vac0_virW*((t < 14) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[417] = infectious_countryA_vac0_virM*((t < 14) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[418] = infectious_countryB_vac0_virW*((t < 14) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[419] = infectious_countryB_vac0_virM*((t < 14) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[420] = recovered_countryA_vac0_virW*((t < 14) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[421] = recovered_countryA_vac0_virM*((t < 14) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[422] = recovered_countryB_vac0_virW*((t < 14) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[423] = recovered_countryB_vac0_virM*((t < 14) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[424] = susceptible_countryA_vac0*((t < 28) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[425] = susceptible_countryB_vac0*((t < 28) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[426] = infectious_countryA_vac0_virW*((t < 28) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[427] = infectious_countryA_vac0_virM*((t < 28) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[428] = infectious_countryB_vac0_virW*((t < 28) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[429] = infectious_countryB_vac0_virM*((t < 28) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[430] = recovered_countryA_vac0_virW*((t < 28) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[431] = recovered_countryA_vac0_virM*((t < 28) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[432] = recovered_countryB_vac0_virW*((t < 28) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[433] = recovered_countryB_vac0_virM*((t < 28) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[434] = susceptible_countryA_vac0*((t < 42) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[435] = susceptible_countryB_vac0*((t < 42) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[436] = infectious_countryA_vac0_virW*((t < 42) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[437] = infectious_countryA_vac0_virM*((t < 42) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[438] = infectious_countryB_vac0_virW*((t < 42) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[439] = infectious_countryB_vac0_virM*((t < 42) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[440] = recovered_countryA_vac0_virW*((t < 42) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[441] = recovered_countryA_vac0_virM*((t < 42) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[442] = recovered_countryB_vac0_virW*((t < 42) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[443] = recovered_countryB_vac0_virM*((t < 42) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[444] = susceptible_countryA_vac0*((t < 56) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[445] = susceptible_countryB_vac0*((t < 56) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[446] = infectious_countryA_vac0_virW*((t < 56) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[447] = infectious_countryA_vac0_virM*((t < 56) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[448] = infectious_countryB_vac0_virW*((t < 56) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[449] = infectious_countryB_vac0_virM*((t < 56) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[450] = recovered_countryA_vac0_virW*((t < 56) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[451] = recovered_countryA_vac0_virM*((t < 56) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[452] = recovered_countryB_vac0_virW*((t < 56) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[453] = recovered_countryB_vac0_virM*((t < 56) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[454] = susceptible_countryA_vac0*((t < 70) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[455] = susceptible_countryB_vac0*((t < 70) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[456] = infectious_countryA_vac0_virW*((t < 70) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[457] = infectious_countryA_vac0_virM*((t < 70) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[458] = infectious_countryB_vac0_virW*((t < 70) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[459] = infectious_countryB_vac0_virM*((t < 70) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[460] = recovered_countryA_vac0_virW*((t < 70) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[461] = recovered_countryA_vac0_virM*((t < 70) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[462] = recovered_countryB_vac0_virW*((t < 70) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[463] = recovered_countryB_vac0_virM*((t < 70) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[464] = susceptible_countryA_vac0*((t < 84) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[465] = susceptible_countryB_vac0*((t < 84) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[466] = infectious_countryA_vac0_virW*((t < 84) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[467] = infectious_countryA_vac0_virM*((t < 84) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[468] = infectious_countryB_vac0_virW*((t < 84) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[469] = infectious_countryB_vac0_virM*((t < 84) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[470] = recovered_countryA_vac0_virW*((t < 84) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[471] = recovered_countryA_vac0_virM*((t < 84) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[472] = recovered_countryB_vac0_virW*((t < 84) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[473] = recovered_countryB_vac0_virM*((t < 84) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[474] = susceptible_countryA_vac0*((t < 98) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[475] = susceptible_countryB_vac0*((t < 98) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[476] = infectious_countryA_vac0_virW*((t < 98) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[477] = infectious_countryA_vac0_virM*((t < 98) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[478] = infectious_countryB_vac0_virW*((t < 98) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[479] = infectious_countryB_vac0_virM*((t < 98) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[480] = recovered_countryA_vac0_virW*((t < 98) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[481] = recovered_countryA_vac0_virM*((t < 98) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[482] = recovered_countryB_vac0_virW*((t < 98) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[483] = recovered_countryB_vac0_virM*((t < 98) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[484] = susceptible_countryA_vac0*((t < 112) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[485] = susceptible_countryB_vac0*((t < 112) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[486] = infectious_countryA_vac0_virW*((t < 112) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[487] = infectious_countryA_vac0_virM*((t < 112) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[488] = infectious_countryB_vac0_virW*((t < 112) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[489] = infectious_countryB_vac0_virM*((t < 112) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[490] = recovered_countryA_vac0_virW*((t < 112) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[491] = recovered_countryA_vac0_virM*((t < 112) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[492] = recovered_countryB_vac0_virW*((t < 112) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[493] = recovered_countryB_vac0_virM*((t < 112) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[494] = susceptible_countryA_vac0*((t < 126) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[495] = susceptible_countryB_vac0*((t < 126) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[496] = infectious_countryA_vac0_virW*((t < 126) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[497] = infectious_countryA_vac0_virM*((t < 126) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[498] = infectious_countryB_vac0_virW*((t < 126) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[499] = infectious_countryB_vac0_virM*((t < 126) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[500] = recovered_countryA_vac0_virW*((t < 126) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[501] = recovered_countryA_vac0_virM*((t < 126) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[502] = recovered_countryB_vac0_virW*((t < 126) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[503] = recovered_countryB_vac0_virM*((t < 126) ? (
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
   vaccine_supply_par_vac2_xx9
))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[504] = susceptible_countryA_vac0*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[505] = susceptible_countryB_vac0*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[506] = infectious_countryA_vac0_virW*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[507] = infectious_countryA_vac0_virM*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[508] = infectious_countryB_vac0_virW*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[509] = infectious_countryB_vac0_virM*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[510] = recovered_countryA_vac0_virW*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[511] = recovered_countryA_vac0_virM*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[512] = recovered_countryB_vac0_virW*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[513] = recovered_countryB_vac0_virM*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[514] = susceptible_countryA_vac0*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[515] = susceptible_countryB_vac0*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[516] = infectious_countryA_vac0_virW*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[517] = infectious_countryA_vac0_virM*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[518] = infectious_countryB_vac0_virW*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[519] = infectious_countryB_vac0_virM*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[520] = recovered_countryA_vac0_virW*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[521] = recovered_countryA_vac0_virM*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[522] = recovered_countryB_vac0_virW*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[523] = recovered_countryB_vac0_virM*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[524] = susceptible_countryA_vac0*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[525] = susceptible_countryB_vac0*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[526] = infectious_countryA_vac0_virW*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[527] = infectious_countryA_vac0_virM*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[528] = infectious_countryB_vac0_virW*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[529] = infectious_countryB_vac0_virM*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[530] = recovered_countryA_vac0_virW*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[531] = recovered_countryA_vac0_virM*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[532] = recovered_countryB_vac0_virW*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[533] = recovered_countryB_vac0_virM*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[534] = susceptible_countryA_vac0*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[535] = susceptible_countryB_vac0*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[536] = infectious_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[537] = infectious_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[538] = infectious_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[539] = infectious_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[540] = recovered_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[541] = recovered_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[542] = recovered_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[543] = recovered_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[544] = susceptible_countryA_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[545] = susceptible_countryB_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[546] = infectious_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[547] = infectious_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[548] = infectious_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[549] = infectious_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[550] = recovered_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[551] = recovered_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[552] = recovered_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[553] = recovered_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[554] = susceptible_countryA_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[555] = susceptible_countryB_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[556] = infectious_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[557] = infectious_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[558] = infectious_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[559] = infectious_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[560] = recovered_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[561] = recovered_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[562] = recovered_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[563] = recovered_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[564] = susceptible_countryA_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[565] = susceptible_countryB_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[566] = infectious_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[567] = infectious_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[568] = infectious_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[569] = infectious_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[570] = recovered_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[571] = recovered_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[572] = recovered_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[573] = recovered_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[574] = susceptible_countryA_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[575] = susceptible_countryB_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[576] = infectious_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[577] = infectious_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[578] = infectious_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[579] = infectious_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[580] = recovered_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[581] = recovered_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[582] = recovered_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[583] = recovered_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[584] = susceptible_countryA_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[585] = susceptible_countryB_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[586] = infectious_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[587] = infectious_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[588] = infectious_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[589] = infectious_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[590] = recovered_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[591] = recovered_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[592] = recovered_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[593] = recovered_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[594] = susceptible_countryA_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[595] = susceptible_countryB_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[596] = infectious_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[597] = infectious_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[598] = infectious_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[599] = infectious_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[600] = recovered_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[601] = recovered_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[602] = recovered_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[603] = recovered_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[604] = susceptible_countryA_vac0*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[605] = susceptible_countryB_vac0*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[606] = infectious_countryA_vac0_virW*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[607] = infectious_countryA_vac0_virM*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[608] = infectious_countryB_vac0_virW*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[609] = infectious_countryB_vac0_virM*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[610] = recovered_countryA_vac0_virW*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[611] = recovered_countryA_vac0_virM*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[612] = recovered_countryB_vac0_virW*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[613] = recovered_countryB_vac0_virM*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[614] = susceptible_countryA_vac0*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[615] = susceptible_countryB_vac0*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[616] = infectious_countryA_vac0_virW*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[617] = infectious_countryA_vac0_virM*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[618] = infectious_countryB_vac0_virW*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[619] = infectious_countryB_vac0_virM*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[620] = recovered_countryA_vac0_virW*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[621] = recovered_countryA_vac0_virM*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[622] = recovered_countryB_vac0_virW*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[623] = recovered_countryB_vac0_virM*((t >= xx0 && t < xx1) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[624] = susceptible_countryA_vac0*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[625] = susceptible_countryB_vac0*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[626] = infectious_countryA_vac0_virW*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[627] = infectious_countryA_vac0_virM*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[628] = infectious_countryB_vac0_virW*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[629] = infectious_countryB_vac0_virM*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[630] = recovered_countryA_vac0_virW*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[631] = recovered_countryA_vac0_virM*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[632] = recovered_countryB_vac0_virW*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[633] = recovered_countryB_vac0_virM*(((t >= xx0 || t >= xx1) && (t >= xx0 || t < xx2) && (t < xx1 || t < xx2)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[634] = susceptible_countryA_vac0*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[635] = susceptible_countryB_vac0*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[636] = infectious_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[637] = infectious_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[638] = infectious_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[639] = infectious_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[640] = recovered_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[641] = recovered_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[642] = recovered_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[643] = recovered_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2) && (t >= xx0 || t >= xx1 || t < xx3) && (t >= xx0 || t < xx2 || t < xx3) && (t < xx1 || t < xx2 || t < xx3)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[644] = susceptible_countryA_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[645] = susceptible_countryB_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[646] = infectious_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[647] = infectious_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[648] = infectious_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[649] = infectious_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[650] = recovered_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[651] = recovered_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[652] = recovered_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[653] = recovered_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4) && (t < xx1 || t < xx2 || t < xx3 || t < xx4)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[654] = susceptible_countryA_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[655] = susceptible_countryB_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[656] = infectious_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[657] = infectious_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[658] = infectious_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[659] = infectious_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[660] = recovered_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[661] = recovered_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[662] = recovered_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[663] = recovered_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[664] = susceptible_countryA_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[665] = susceptible_countryB_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[666] = infectious_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[667] = infectious_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[668] = infectious_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[669] = infectious_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[670] = recovered_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[671] = recovered_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[672] = recovered_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[673] = recovered_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[674] = susceptible_countryA_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[675] = susceptible_countryB_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[676] = infectious_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[677] = infectious_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[678] = infectious_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[679] = infectious_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[680] = recovered_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[681] = recovered_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[682] = recovered_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[683] = recovered_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[684] = susceptible_countryA_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[685] = susceptible_countryB_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[686] = infectious_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[687] = infectious_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[688] = infectious_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[689] = infectious_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[690] = recovered_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[691] = recovered_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[692] = recovered_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[693] = recovered_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8)) ? (
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[694] = susceptible_countryA_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[695] = susceptible_countryB_vac0*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[696] = infectious_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[697] = infectious_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[698] = infectious_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[699] = infectious_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[700] = recovered_countryA_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[701] = recovered_countryA_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
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
: (
   0
))))))))))))/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0);
    dwdp[702] = recovered_countryB_vac0_virW*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[703] = recovered_countryB_vac0_virM*(((t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t >= xx8) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t >= xx7 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t >= xx6 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t >= xx5 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t >= xx4 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t >= xx3 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t >= xx2 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t >= xx1 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t >= xx0 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9) && (t < xx1 || t < xx2 || t < xx3 || t < xx4 || t < xx5 || t < xx6 || t < xx7 || t < xx8 || t < xx9)) ? (
   0
)
: (
   1
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
: (
   0
))))))))))))/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + susceptible_countryB_vac0);
    dwdp[704] = -beta*eta_virW*infectious_countryA_vac0_virW*susceptible_countryA_vac0*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/((infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[705] = -beta*eta_virW*infectious_countryA_vac0_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/((infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[706] = -beta*eta_virW*infectious_countryA_vac0_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/((infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[707] = beta*eta_virW*infectious_countryA_vac0_virW*susceptible_countryB_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[708] = beta*eta_virW*infectious_countryA_vac0_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[709] = beta*eta_virW*infectious_countryA_vac0_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[710] = -beta*eta_virM*infectious_countryA_vac0_virM*susceptible_countryA_vac0*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/((infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[711] = -beta*eta_virM*infectious_countryA_vac0_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/((infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[712] = -beta*eta_virM*infectious_countryA_vac0_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/((infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[713] = beta*eta_virM*infectious_countryA_vac0_virM*susceptible_countryB_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[714] = beta*eta_virM*infectious_countryA_vac0_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[715] = beta*eta_virM*infectious_countryA_vac0_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[716] = -beta*eta_virW*infectious_countryA_vac1_virW*susceptible_countryA_vac0*(1 - gamma)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/((infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[717] = -beta*eta_virW*infectious_countryA_vac1_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(1 - gamma)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/((infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[718] = -beta*eta_virW*infectious_countryA_vac1_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(1 - gamma)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/((infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[719] = beta*eta_virW*infectious_countryA_vac1_virW*susceptible_countryB_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[720] = beta*eta_virW*infectious_countryA_vac1_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[721] = beta*eta_virW*infectious_countryA_vac1_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[722] = -beta*eta_virM*infectious_countryA_vac1_virM*susceptible_countryA_vac0*(1 - gamma)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/((infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[723] = -beta*eta_virM*infectious_countryA_vac1_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(1 - gamma)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/((infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[724] = -beta*eta_virM*infectious_countryA_vac1_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(1 - gamma)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/((infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[725] = beta*eta_virM*infectious_countryA_vac1_virM*susceptible_countryB_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[726] = beta*eta_virM*infectious_countryA_vac1_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[727] = beta*eta_virM*infectious_countryA_vac1_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[728] = -beta*eta_virW*infectious_countryA_vac2_virW*susceptible_countryA_vac0*(1 - gamma)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/((infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[729] = -beta*eta_virW*infectious_countryA_vac2_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(1 - gamma)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/((infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[730] = -beta*eta_virW*infectious_countryA_vac2_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(1 - gamma)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/((infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[731] = beta*eta_virW*infectious_countryA_vac2_virW*susceptible_countryB_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[732] = beta*eta_virW*infectious_countryA_vac2_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[733] = beta*eta_virW*infectious_countryA_vac2_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[734] = -beta*eta_virM*infectious_countryA_vac2_virM*susceptible_countryA_vac0*(1 - gamma)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/((infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[735] = -beta*eta_virM*infectious_countryA_vac2_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(1 - gamma)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/((infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[736] = -beta*eta_virM*infectious_countryA_vac2_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(1 - gamma)*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/((infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[737] = beta*eta_virM*infectious_countryA_vac2_virM*susceptible_countryB_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[738] = beta*eta_virM*infectious_countryA_vac2_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[739] = beta*eta_virM*infectious_countryA_vac2_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[740] = beta*eta_virW*infectious_countryB_vac0_virW*susceptible_countryA_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[741] = beta*eta_virW*infectious_countryB_vac0_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[742] = beta*eta_virW*infectious_countryB_vac0_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[743] = -beta*eta_virW*infectious_countryB_vac0_virW*susceptible_countryB_vac0*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/((infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[744] = -beta*eta_virW*infectious_countryB_vac0_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/((infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[745] = -beta*eta_virW*infectious_countryB_vac0_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/((infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[746] = beta*eta_virM*infectious_countryB_vac0_virM*susceptible_countryA_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[747] = beta*eta_virM*infectious_countryB_vac0_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[748] = beta*eta_virM*infectious_countryB_vac0_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[749] = -beta*eta_virM*infectious_countryB_vac0_virM*susceptible_countryB_vac0*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/((infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[750] = -beta*eta_virM*infectious_countryB_vac0_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/((infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[751] = -beta*eta_virM*infectious_countryB_vac0_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/((infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[752] = beta*eta_virW*infectious_countryB_vac1_virW*susceptible_countryA_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[753] = beta*eta_virW*infectious_countryB_vac1_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[754] = beta*eta_virW*infectious_countryB_vac1_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[755] = -beta*eta_virW*infectious_countryB_vac1_virW*susceptible_countryB_vac0*(1 - gamma)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/((infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[756] = -beta*eta_virW*infectious_countryB_vac1_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)*(1 - gamma)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/((infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[757] = -beta*eta_virW*infectious_countryB_vac1_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)*(1 - gamma)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/((infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[758] = beta*eta_virM*infectious_countryB_vac1_virM*susceptible_countryA_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[759] = beta*eta_virM*infectious_countryB_vac1_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[760] = beta*eta_virM*infectious_countryB_vac1_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[761] = -beta*eta_virM*infectious_countryB_vac1_virM*susceptible_countryB_vac0*(1 - gamma)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/((infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[762] = -beta*eta_virM*infectious_countryB_vac1_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)*(1 - gamma)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/((infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[763] = -beta*eta_virM*infectious_countryB_vac1_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)*(1 - gamma)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/((infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[764] = beta*eta_virW*infectious_countryB_vac2_virW*susceptible_countryA_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[765] = beta*eta_virW*infectious_countryB_vac2_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[766] = beta*eta_virW*infectious_countryB_vac2_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[767] = -beta*eta_virW*infectious_countryB_vac2_virW*susceptible_countryB_vac0*(1 - gamma)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/((infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[768] = -beta*eta_virW*infectious_countryB_vac2_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)*(1 - gamma)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/((infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[769] = -beta*eta_virW*infectious_countryB_vac2_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)*(1 - gamma)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/((infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[770] = beta*eta_virM*infectious_countryB_vac2_virM*susceptible_countryA_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[771] = beta*eta_virM*infectious_countryB_vac2_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[772] = beta*eta_virM*infectious_countryB_vac2_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    dwdp[773] = -beta*eta_virM*infectious_countryB_vac2_virM*susceptible_countryB_vac0*(1 - gamma)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/((infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[774] = -beta*eta_virM*infectious_countryB_vac2_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)*(1 - gamma)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/((infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
    dwdp[775] = -beta*eta_virM*infectious_countryB_vac2_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)*(1 - gamma)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/((infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2));
}

} // namespace amici
} // namespace model_vaccination_piecewise