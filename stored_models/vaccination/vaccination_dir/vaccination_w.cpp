#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "tcl.h"
#include "p.h"
#include "k.h"
#include "spl.h"
#include "x.h"

namespace amici {
namespace model_vaccination {

void w_vaccination(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl, const realtype *spl){
    w[0] = spl_0;
    w[1] = spl_1;
    w[2] = infectious_countryA_vac0_virW*lambda1*(1 - prob_deceasing);
    w[3] = infectious_countryA_vac0_virW*lambda1*prob_deceasing;
    w[4] = infectious_countryA_vac0_virM*lambda1*(1 - prob_deceasing);
    w[5] = infectious_countryA_vac0_virM*lambda1*prob_deceasing;
    w[6] = infectious_countryA_vac1_virW*lambda1*(-prob_deceasing*(1 - omega_vac1_virW) + 1);
    w[7] = infectious_countryA_vac1_virW*lambda1*prob_deceasing*(1 - omega_vac1_virW);
    w[8] = infectious_countryA_vac1_virM*lambda1*(-prob_deceasing*(1 - omega_vac1_virM) + 1);
    w[9] = infectious_countryA_vac1_virM*lambda1*prob_deceasing*(1 - omega_vac1_virM);
    w[10] = infectious_countryA_vac2_virW*lambda1*(-prob_deceasing*(1 - omega_vac2_virW) + 1);
    w[11] = infectious_countryA_vac2_virW*lambda1*prob_deceasing*(1 - omega_vac2_virW);
    w[12] = infectious_countryA_vac2_virM*lambda1*(-prob_deceasing*(1 - omega_vac2_virM) + 1);
    w[13] = infectious_countryA_vac2_virM*lambda1*prob_deceasing*(1 - omega_vac2_virM);
    w[14] = infectious_countryB_vac0_virW*lambda1*(1 - prob_deceasing);
    w[15] = infectious_countryB_vac0_virW*lambda1*prob_deceasing;
    w[16] = infectious_countryB_vac0_virM*lambda1*(1 - prob_deceasing);
    w[17] = infectious_countryB_vac0_virM*lambda1*prob_deceasing;
    w[18] = infectious_countryB_vac1_virW*lambda1*(-prob_deceasing*(1 - omega_vac1_virW) + 1);
    w[19] = infectious_countryB_vac1_virW*lambda1*prob_deceasing*(1 - omega_vac1_virW);
    w[20] = infectious_countryB_vac1_virM*lambda1*(-prob_deceasing*(1 - omega_vac1_virM) + 1);
    w[21] = infectious_countryB_vac1_virM*lambda1*prob_deceasing*(1 - omega_vac1_virM);
    w[22] = infectious_countryB_vac2_virW*lambda1*(-prob_deceasing*(1 - omega_vac2_virW) + 1);
    w[23] = infectious_countryB_vac2_virW*lambda1*prob_deceasing*(1 - omega_vac2_virW);
    w[24] = infectious_countryB_vac2_virM*lambda1*(-prob_deceasing*(1 - omega_vac2_virM) + 1);
    w[25] = infectious_countryB_vac2_virM*lambda1*prob_deceasing*(1 - omega_vac2_virM);
    w[26] = beta*eta_virW*infectious_countryA_vac0_virW*susceptible_countryA_vac0*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    w[27] = beta*eta_virW*infectious_countryA_vac0_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    w[28] = beta*eta_virW*infectious_countryA_vac0_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    w[29] = beta*distance_countryA_countryB*eta_virW*infectious_countryA_vac0_virW*susceptible_countryB_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[30] = beta*distance_countryA_countryB*eta_virW*infectious_countryA_vac0_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[31] = beta*distance_countryA_countryB*eta_virW*infectious_countryA_vac0_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[32] = beta*eta_virM*infectious_countryA_vac0_virM*susceptible_countryA_vac0*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    w[33] = beta*eta_virM*infectious_countryA_vac0_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    w[34] = beta*eta_virM*infectious_countryA_vac0_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    w[35] = beta*distance_countryA_countryB*eta_virM*infectious_countryA_vac0_virM*susceptible_countryB_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[36] = beta*distance_countryA_countryB*eta_virM*infectious_countryA_vac0_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[37] = beta*distance_countryA_countryB*eta_virM*infectious_countryA_vac0_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[38] = beta*eta_virW*infectious_countryA_vac1_virW*susceptible_countryA_vac0*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    w[39] = beta*eta_virW*infectious_countryA_vac1_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    w[40] = beta*eta_virW*infectious_countryA_vac1_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    w[41] = beta*distance_countryA_countryB*eta_virW*infectious_countryA_vac1_virW*susceptible_countryB_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[42] = beta*distance_countryA_countryB*eta_virW*infectious_countryA_vac1_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[43] = beta*distance_countryA_countryB*eta_virW*infectious_countryA_vac1_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[44] = beta*eta_virM*infectious_countryA_vac1_virM*susceptible_countryA_vac0*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    w[45] = beta*eta_virM*infectious_countryA_vac1_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    w[46] = beta*eta_virM*infectious_countryA_vac1_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    w[47] = beta*distance_countryA_countryB*eta_virM*infectious_countryA_vac1_virM*susceptible_countryB_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[48] = beta*distance_countryA_countryB*eta_virM*infectious_countryA_vac1_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[49] = beta*distance_countryA_countryB*eta_virM*infectious_countryA_vac1_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[50] = beta*eta_virW*infectious_countryA_vac2_virW*susceptible_countryA_vac0*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    w[51] = beta*eta_virW*infectious_countryA_vac2_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    w[52] = beta*eta_virW*infectious_countryA_vac2_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    w[53] = beta*distance_countryA_countryB*eta_virW*infectious_countryA_vac2_virW*susceptible_countryB_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[54] = beta*distance_countryA_countryB*eta_virW*infectious_countryA_vac2_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[55] = beta*distance_countryA_countryB*eta_virW*infectious_countryA_vac2_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[56] = beta*eta_virM*infectious_countryA_vac2_virM*susceptible_countryA_vac0*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    w[57] = beta*eta_virM*infectious_countryA_vac2_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    w[58] = beta*eta_virM*infectious_countryA_vac2_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(1 - gamma)*(-distance_countryA_countryB*(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);
    w[59] = beta*distance_countryA_countryB*eta_virM*infectious_countryA_vac2_virM*susceptible_countryB_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[60] = beta*distance_countryA_countryB*eta_virM*infectious_countryA_vac2_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[61] = beta*distance_countryA_countryB*eta_virM*infectious_countryA_vac2_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[62] = beta*distance_countryB_countryA*eta_virW*infectious_countryB_vac0_virW*susceptible_countryA_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[63] = beta*distance_countryB_countryA*eta_virW*infectious_countryB_vac0_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[64] = beta*distance_countryB_countryA*eta_virW*infectious_countryB_vac0_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[65] = beta*eta_virW*infectious_countryB_vac0_virW*susceptible_countryB_vac0*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[66] = beta*eta_virW*infectious_countryB_vac0_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[67] = beta*eta_virW*infectious_countryB_vac0_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[68] = beta*distance_countryB_countryA*eta_virM*infectious_countryB_vac0_virM*susceptible_countryA_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[69] = beta*distance_countryB_countryA*eta_virM*infectious_countryB_vac0_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[70] = beta*distance_countryB_countryA*eta_virM*infectious_countryB_vac0_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[71] = beta*eta_virM*infectious_countryB_vac0_virM*susceptible_countryB_vac0*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[72] = beta*eta_virM*infectious_countryB_vac0_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[73] = beta*eta_virM*infectious_countryB_vac0_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[74] = beta*distance_countryB_countryA*eta_virW*infectious_countryB_vac1_virW*susceptible_countryA_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[75] = beta*distance_countryB_countryA*eta_virW*infectious_countryB_vac1_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[76] = beta*distance_countryB_countryA*eta_virW*infectious_countryB_vac1_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[77] = beta*eta_virW*infectious_countryB_vac1_virW*susceptible_countryB_vac0*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[78] = beta*eta_virW*infectious_countryB_vac1_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[79] = beta*eta_virW*infectious_countryB_vac1_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[80] = beta*distance_countryB_countryA*eta_virM*infectious_countryB_vac1_virM*susceptible_countryA_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[81] = beta*distance_countryB_countryA*eta_virM*infectious_countryB_vac1_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[82] = beta*distance_countryB_countryA*eta_virM*infectious_countryB_vac1_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[83] = beta*eta_virM*infectious_countryB_vac1_virM*susceptible_countryB_vac0*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[84] = beta*eta_virM*infectious_countryB_vac1_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[85] = beta*eta_virM*infectious_countryB_vac1_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[86] = beta*distance_countryB_countryA*eta_virW*infectious_countryB_vac2_virW*susceptible_countryA_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[87] = beta*distance_countryB_countryA*eta_virW*infectious_countryB_vac2_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[88] = beta*distance_countryB_countryA*eta_virW*infectious_countryB_vac2_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[89] = beta*eta_virW*infectious_countryB_vac2_virW*susceptible_countryB_vac0*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[90] = beta*eta_virW*infectious_countryB_vac2_virW*susceptible_countryB_vac1*(1 - delta_vac1_virW)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[91] = beta*eta_virW*infectious_countryB_vac2_virW*susceptible_countryB_vac2*(1 - delta_vac2_virW)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[92] = beta*distance_countryB_countryA*eta_virM*infectious_countryB_vac2_virM*susceptible_countryA_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[93] = beta*distance_countryB_countryA*eta_virM*infectious_countryB_vac2_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[94] = beta*distance_countryB_countryA*eta_virM*infectious_countryB_vac2_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[95] = beta*eta_virM*infectious_countryB_vac2_virM*susceptible_countryB_vac0*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[96] = beta*eta_virM*infectious_countryB_vac2_virM*susceptible_countryB_vac1*(1 - delta_vac1_virM)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[97] = beta*eta_virM*infectious_countryB_vac2_virM*susceptible_countryB_vac2*(1 - delta_vac2_virM)*(1 - gamma)*(-distance_countryB_countryA*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2 + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2) + 1)/(infectious_countryB_vac0_virM + infectious_countryB_vac0_virW + infectious_countryB_vac1_virM + infectious_countryB_vac1_virW + infectious_countryB_vac2_virM + infectious_countryB_vac2_virW + recovered_countryB_vac0_virM + recovered_countryB_vac0_virW + recovered_countryB_vac1_virM + recovered_countryB_vac1_virW + recovered_countryB_vac2_virM + recovered_countryB_vac2_virW + susceptible_countryB_vac0 + susceptible_countryB_vac1 + susceptible_countryB_vac2);
    w[98] = susceptible_countryA_vac0*((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac1_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac1_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac1_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac1_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac1_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac1_xx5
)
: ((t >= xx6 && t < xx7) ? (
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
))))))))))/((std::pow(0.36787944123356736, spl_0) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
    w[99] = susceptible_countryA_vac0*((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac2_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac2_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac2_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac2_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac2_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac2_xx5
)
: ((t >= xx6 && t < xx7) ? (
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
))))))))))/((std::pow(0.36787944123356736, spl_1) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
    w[100] = susceptible_countryB_vac0*(1 - 1/(std::pow(0.36787944123356736, spl_0) + 1))*((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac1_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac1_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac1_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac1_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac1_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac1_xx5
)
: ((t >= xx6 && t < xx7) ? (
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
    w[101] = susceptible_countryB_vac0*(1 - 1/(std::pow(0.36787944123356736, spl_1) + 1))*((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac2_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac2_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac2_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac2_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac2_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac2_xx5
)
: ((t >= xx6 && t < xx7) ? (
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
    w[102] = infectious_countryA_vac0_virW*((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac1_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac1_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac1_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac1_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac1_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac1_xx5
)
: ((t >= xx6 && t < xx7) ? (
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
))))))))))/((std::pow(0.36787944123356736, spl_0) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
    w[103] = infectious_countryA_vac0_virW*((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac2_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac2_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac2_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac2_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac2_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac2_xx5
)
: ((t >= xx6 && t < xx7) ? (
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
))))))))))/((std::pow(0.36787944123356736, spl_1) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
    w[104] = infectious_countryA_vac0_virM*((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac1_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac1_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac1_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac1_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac1_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac1_xx5
)
: ((t >= xx6 && t < xx7) ? (
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
))))))))))/((std::pow(0.36787944123356736, spl_0) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
    w[105] = infectious_countryA_vac0_virM*((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac2_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac2_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac2_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac2_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac2_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac2_xx5
)
: ((t >= xx6 && t < xx7) ? (
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
))))))))))/((std::pow(0.36787944123356736, spl_1) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
    w[106] = infectious_countryB_vac0_virW*(1 - 1/(std::pow(0.36787944123356736, spl_0) + 1))*((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac1_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac1_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac1_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac1_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac1_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac1_xx5
)
: ((t >= xx6 && t < xx7) ? (
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
    w[107] = infectious_countryB_vac0_virW*(1 - 1/(std::pow(0.36787944123356736, spl_1) + 1))*((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac2_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac2_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac2_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac2_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac2_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac2_xx5
)
: ((t >= xx6 && t < xx7) ? (
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
    w[108] = infectious_countryB_vac0_virM*(1 - 1/(std::pow(0.36787944123356736, spl_0) + 1))*((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac1_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac1_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac1_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac1_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac1_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac1_xx5
)
: ((t >= xx6 && t < xx7) ? (
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
    w[109] = infectious_countryB_vac0_virM*(1 - 1/(std::pow(0.36787944123356736, spl_1) + 1))*((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac2_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac2_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac2_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac2_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac2_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac2_xx5
)
: ((t >= xx6 && t < xx7) ? (
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
    w[110] = recovered_countryA_vac0_virW*((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac1_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac1_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac1_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac1_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac1_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac1_xx5
)
: ((t >= xx6 && t < xx7) ? (
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
))))))))))/((std::pow(0.36787944123356736, spl_0) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
    w[111] = recovered_countryA_vac0_virW*((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac2_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac2_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac2_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac2_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac2_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac2_xx5
)
: ((t >= xx6 && t < xx7) ? (
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
))))))))))/((std::pow(0.36787944123356736, spl_1) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
    w[112] = recovered_countryA_vac0_virM*((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac1_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac1_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac1_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac1_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac1_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac1_xx5
)
: ((t >= xx6 && t < xx7) ? (
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
))))))))))/((std::pow(0.36787944123356736, spl_0) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
    w[113] = recovered_countryA_vac0_virM*((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac2_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac2_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac2_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac2_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac2_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac2_xx5
)
: ((t >= xx6 && t < xx7) ? (
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
))))))))))/((std::pow(0.36787944123356736, spl_1) + 1)*(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + susceptible_countryA_vac0));
    w[114] = recovered_countryB_vac0_virW*(1 - 1/(std::pow(0.36787944123356736, spl_0) + 1))*((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac1_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac1_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac1_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac1_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac1_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac1_xx5
)
: ((t >= xx6 && t < xx7) ? (
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
    w[115] = recovered_countryB_vac0_virW*(1 - 1/(std::pow(0.36787944123356736, spl_1) + 1))*((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac2_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac2_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac2_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac2_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac2_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac2_xx5
)
: ((t >= xx6 && t < xx7) ? (
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
    w[116] = recovered_countryB_vac0_virM*(1 - 1/(std::pow(0.36787944123356736, spl_0) + 1))*((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac1_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac1_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac1_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac1_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac1_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac1_xx5
)
: ((t >= xx6 && t < xx7) ? (
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
    w[117] = recovered_countryB_vac0_virM*(1 - 1/(std::pow(0.36787944123356736, spl_1) + 1))*((t >= xx0 && t < xx1) ? (
   vaccine_supply_par_vac2_xx0
)
: ((t >= xx1 && t < xx2) ? (
   vaccine_supply_par_vac2_xx1
)
: ((t >= xx2 && t < xx3) ? (
   vaccine_supply_par_vac2_xx2
)
: ((t >= xx3 && t < xx4) ? (
   vaccine_supply_par_vac2_xx3
)
: ((t >= xx4 && t < xx5) ? (
   vaccine_supply_par_vac2_xx4
)
: ((t >= xx5 && t < xx6) ? (
   vaccine_supply_par_vac2_xx5
)
: ((t >= xx6 && t < xx7) ? (
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
}

} // namespace amici
} // namespace model_vaccination