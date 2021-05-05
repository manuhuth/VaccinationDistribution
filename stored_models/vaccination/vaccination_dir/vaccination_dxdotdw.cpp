#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "vaccination_x.h"
#include "vaccination_p.h"
#include "vaccination_w.h"
#include "vaccination_dxdotdw.h"

namespace amici {
namespace model_vaccination {

void dxdotdw_vaccination(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dxdot3_dflux_r0 = -1.0;  // dxdotdw[0]
    dxdot9_dflux_r0 = 1.0;  // dxdotdw[1]
    dxdot3_dflux_r1 = -1.0;  // dxdotdw[2]
    dxdot15_dflux_r1 = 1.0;  // dxdotdw[3]
    dxdot4_dflux_r2 = -1.0;  // dxdotdw[4]
    dxdot10_dflux_r2 = 1.0;  // dxdotdw[5]
    dxdot4_dflux_r3 = -1.0;  // dxdotdw[6]
    dxdot16_dflux_r3 = 1.0;  // dxdotdw[7]
    dxdot5_dflux_r4 = -1.0;  // dxdotdw[8]
    dxdot11_dflux_r4 = 1.0;  // dxdotdw[9]
    dxdot5_dflux_r5 = -1.0;  // dxdotdw[10]
    dxdot17_dflux_r5 = 1.0;  // dxdotdw[11]
    dxdot6_dflux_r6 = -1.0;  // dxdotdw[12]
    dxdot12_dflux_r6 = 1.0;  // dxdotdw[13]
    dxdot6_dflux_r7 = -1.0;  // dxdotdw[14]
    dxdot18_dflux_r7 = 1.0;  // dxdotdw[15]
    dxdot7_dflux_r8 = -1.0;  // dxdotdw[16]
    dxdot13_dflux_r8 = 1.0;  // dxdotdw[17]
    dxdot7_dflux_r9 = -1.0;  // dxdotdw[18]
    dxdot19_dflux_r9 = 1.0;  // dxdotdw[19]
    dxdot8_dflux_r10 = -1.0;  // dxdotdw[20]
    dxdot14_dflux_r10 = 1.0;  // dxdotdw[21]
    dxdot8_dflux_r11 = -1.0;  // dxdotdw[22]
    dxdot20_dflux_r11 = 1.0;  // dxdotdw[23]
    dxdot0_dflux_r12 = -1.0;  // dxdotdw[24]
    dxdot3_dflux_r12 = 1.0;  // dxdotdw[25]
    dxdot1_dflux_r13 = -1.0;  // dxdotdw[26]
    dxdot3_dflux_r13 = 1.0;  // dxdotdw[27]
    dxdot2_dflux_r14 = -1.0;  // dxdotdw[28]
    dxdot3_dflux_r14 = 1.0;  // dxdotdw[29]
    dxdot0_dflux_r15 = -1.0;  // dxdotdw[30]
    dxdot4_dflux_r15 = 1.0;  // dxdotdw[31]
    dxdot1_dflux_r16 = -1.0;  // dxdotdw[32]
    dxdot4_dflux_r16 = 1.0;  // dxdotdw[33]
    dxdot2_dflux_r17 = -1.0;  // dxdotdw[34]
    dxdot4_dflux_r17 = 1.0;  // dxdotdw[35]
    dxdot0_dflux_r18 = -1.0;  // dxdotdw[36]
    dxdot5_dflux_r18 = 1.0;  // dxdotdw[37]
    dxdot1_dflux_r19 = -1.0;  // dxdotdw[38]
    dxdot5_dflux_r19 = 1.0;  // dxdotdw[39]
    dxdot2_dflux_r20 = -1.0;  // dxdotdw[40]
    dxdot5_dflux_r20 = 1.0;  // dxdotdw[41]
    dxdot0_dflux_r21 = -1.0;  // dxdotdw[42]
    dxdot6_dflux_r21 = 1.0;  // dxdotdw[43]
    dxdot1_dflux_r22 = -1.0;  // dxdotdw[44]
    dxdot6_dflux_r22 = 1.0;  // dxdotdw[45]
    dxdot2_dflux_r23 = -1.0;  // dxdotdw[46]
    dxdot6_dflux_r23 = 1.0;  // dxdotdw[47]
    dxdot0_dflux_r24 = -1.0;  // dxdotdw[48]
    dxdot7_dflux_r24 = 1.0;  // dxdotdw[49]
    dxdot1_dflux_r25 = -1.0;  // dxdotdw[50]
    dxdot7_dflux_r25 = 1.0;  // dxdotdw[51]
    dxdot2_dflux_r26 = -1.0;  // dxdotdw[52]
    dxdot7_dflux_r26 = 1.0;  // dxdotdw[53]
    dxdot0_dflux_r27 = -1.0;  // dxdotdw[54]
    dxdot8_dflux_r27 = 1.0;  // dxdotdw[55]
    dxdot1_dflux_r28 = -1.0;  // dxdotdw[56]
    dxdot8_dflux_r28 = 1.0;  // dxdotdw[57]
    dxdot2_dflux_r29 = -1.0;  // dxdotdw[58]
    dxdot8_dflux_r29 = 1.0;  // dxdotdw[59]
    dxdot0_dflux_r30 = -1.0;  // dxdotdw[60]
    dxdot1_dflux_r30 = 1.0;  // dxdotdw[61]
    dxdot0_dflux_r31 = -1.0;  // dxdotdw[62]
    dxdot2_dflux_r31 = 1.0;  // dxdotdw[63]
    dxdot9_dflux_r32 = -1.0;  // dxdotdw[64]
    dxdot11_dflux_r32 = 1.0;  // dxdotdw[65]
    dxdot9_dflux_r33 = -1.0;  // dxdotdw[66]
    dxdot13_dflux_r33 = 1.0;  // dxdotdw[67]
    dxdot10_dflux_r34 = -1.0;  // dxdotdw[68]
    dxdot12_dflux_r34 = 1.0;  // dxdotdw[69]
    dxdot10_dflux_r35 = -1.0;  // dxdotdw[70]
    dxdot14_dflux_r35 = 1.0;  // dxdotdw[71]
}

} // namespace model_vaccination
} // namespace amici
