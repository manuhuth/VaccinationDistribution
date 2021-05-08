#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "vaccination_x.h"
#include "vaccination_p.h"
#include "vaccination_w.h"

namespace amici {
namespace model_vaccination {

void w_vaccination(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    nu_countryA_vac1 = number_vac1*proportion_countryA_vac1/(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0);  // w[0]
    nu_countryA_vac2 = number_vac2*proportion_countryA_vac2/(recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + 2*susceptible_countryA_vac0);  // w[1]
    flux_r0 = infectious_countryA_vac0_virW + lambda1*(1 - amici_p);  // w[2]
    flux_r1 = amici_p*infectious_countryA_vac0_virW*lambda1;  // w[3]
    flux_r2 = infectious_countryA_vac0_virM + lambda1*(1 - amici_p);  // w[4]
    flux_r3 = amici_p*infectious_countryA_vac0_virM*lambda1;  // w[5]
    flux_r4 = infectious_countryA_vac1_virW + lambda1*(-amici_p*(1 - omega_vac1_virW) + 1);  // w[6]
    flux_r5 = amici_p*infectious_countryA_vac1_virW*lambda1*(1 - omega_vac1_virW);  // w[7]
    flux_r6 = infectious_countryA_vac1_virM + lambda1*(-amici_p*(1 - omega_vac1_virM) + 1);  // w[8]
    flux_r7 = amici_p*infectious_countryA_vac1_virM*lambda1*(1 - omega_vac1_virM);  // w[9]
    flux_r8 = infectious_countryA_vac2_virW + lambda1*(-amici_p*(1 - omega_vac2_virW) + 1);  // w[10]
    flux_r9 = amici_p*infectious_countryA_vac2_virW*lambda1*(1 - omega_vac2_virW);  // w[11]
    flux_r10 = infectious_countryA_vac2_virM + lambda1*(-amici_p*(1 - omega_vac2_virM) + 1);  // w[12]
    flux_r11 = amici_p*infectious_countryA_vac2_virM*lambda1*(1 - omega_vac2_virM);  // w[13]
    flux_r12 = beta*eta_virW*infectious_countryA_vac0_virW*susceptible_countryA_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);  // w[14]
    flux_r13 = beta*eta_virW*infectious_countryA_vac0_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);  // w[15]
    flux_r14 = beta*eta_virW*infectious_countryA_vac0_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);  // w[16]
    flux_r15 = beta*eta_virM*infectious_countryA_vac0_virM*susceptible_countryA_vac0/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);  // w[17]
    flux_r16 = beta*eta_virM*infectious_countryA_vac0_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);  // w[18]
    flux_r17 = beta*eta_virM*infectious_countryA_vac0_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);  // w[19]
    flux_r18 = beta*eta_virW*infectious_countryA_vac1_virW*susceptible_countryA_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);  // w[20]
    flux_r19 = beta*eta_virW*infectious_countryA_vac1_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);  // w[21]
    flux_r20 = beta*eta_virW*infectious_countryA_vac1_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);  // w[22]
    flux_r21 = beta*eta_virM*infectious_countryA_vac1_virM*susceptible_countryA_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);  // w[23]
    flux_r22 = beta*eta_virM*infectious_countryA_vac1_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);  // w[24]
    flux_r23 = beta*eta_virM*infectious_countryA_vac1_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);  // w[25]
    flux_r24 = beta*eta_virW*infectious_countryA_vac2_virW*susceptible_countryA_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);  // w[26]
    flux_r25 = beta*eta_virW*infectious_countryA_vac2_virW*susceptible_countryA_vac1*(1 - delta_vac1_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);  // w[27]
    flux_r26 = beta*eta_virW*infectious_countryA_vac2_virW*susceptible_countryA_vac2*(1 - delta_vac2_virW)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);  // w[28]
    flux_r27 = beta*eta_virM*infectious_countryA_vac2_virM*susceptible_countryA_vac0*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);  // w[29]
    flux_r28 = beta*eta_virM*infectious_countryA_vac2_virM*susceptible_countryA_vac1*(1 - delta_vac1_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);  // w[30]
    flux_r29 = beta*eta_virM*infectious_countryA_vac2_virM*susceptible_countryA_vac2*(1 - delta_vac2_virM)*(1 - gamma)/(infectious_countryA_vac0_virM + infectious_countryA_vac0_virW + infectious_countryA_vac1_virM + infectious_countryA_vac1_virW + infectious_countryA_vac2_virM + infectious_countryA_vac2_virW + recovered_countryA_vac0_virM + recovered_countryA_vac0_virW + recovered_countryA_vac1_virM + recovered_countryA_vac1_virW + recovered_countryA_vac2_virM + recovered_countryA_vac2_virW + susceptible_countryA_vac0 + susceptible_countryA_vac1 + susceptible_countryA_vac2);  // w[31]
    flux_r30 = nu_countryA_vac1*susceptible_countryA_vac0;  // w[32]
    flux_r31 = nu_countryA_vac2*susceptible_countryA_vac0;  // w[33]
    flux_r32 = nu_countryA_vac1*recovered_countryA_vac0_virW;  // w[34]
    flux_r33 = nu_countryA_vac2*recovered_countryA_vac0_virW;  // w[35]
    flux_r34 = nu_countryA_vac1*recovered_countryA_vac0_virM;  // w[36]
    flux_r35 = nu_countryA_vac2*recovered_countryA_vac0_virM;  // w[37]
}

} // namespace model_vaccination
} // namespace amici
