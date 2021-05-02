#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "one_country_x.h"
#include "one_country_p.h"
#include "one_country_w.h"
#include "one_country_xdot.h"

namespace amici {
namespace model_one_country {

void xdot_one_country(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot0 = -1.0*flux_r12 - 1.0*flux_r15 - 1.0*flux_r18 - 1.0*flux_r21 - 1.0*flux_r24 - 1.0*flux_r27 - 1.0*flux_r30 - 1.0*flux_r31;  // xdot[0]
    xdot1 = -1.0*flux_r13 - 1.0*flux_r16 - 1.0*flux_r19 - 1.0*flux_r22 - 1.0*flux_r25 - 1.0*flux_r28 + 1.0*flux_r30;  // xdot[1]
    xdot2 = -1.0*flux_r14 - 1.0*flux_r17 - 1.0*flux_r20 - 1.0*flux_r23 - 1.0*flux_r26 - 1.0*flux_r29 + 1.0*flux_r31;  // xdot[2]
    xdot3 = -1.0*flux_r0 - 1.0*flux_r1 + 1.0*flux_r12 + 1.0*flux_r13 + 1.0*flux_r14;  // xdot[3]
    xdot4 = 1.0*flux_r15 + 1.0*flux_r16 + 1.0*flux_r17 - 1.0*flux_r2 - 1.0*flux_r3;  // xdot[4]
    xdot5 = 1.0*flux_r18 + 1.0*flux_r19 + 1.0*flux_r20 - 1.0*flux_r4 - 1.0*flux_r5;  // xdot[5]
    xdot6 = 1.0*flux_r21 + 1.0*flux_r22 + 1.0*flux_r23 - 1.0*flux_r6 - 1.0*flux_r7;  // xdot[6]
    xdot7 = 1.0*flux_r24 + 1.0*flux_r25 + 1.0*flux_r26 - 1.0*flux_r8 - 1.0*flux_r9;  // xdot[7]
    xdot8 = -1.0*flux_r10 - 1.0*flux_r11 + 1.0*flux_r27 + 1.0*flux_r28 + 1.0*flux_r29;  // xdot[8]
    xdot9 = 1.0*flux_r0 - 1.0*flux_r32 - 1.0*flux_r33;  // xdot[9]
    xdot10 = 1.0*flux_r2 - 1.0*flux_r34 - 1.0*flux_r35;  // xdot[10]
    xdot11 = 1.0*flux_r32 + 1.0*flux_r4;  // xdot[11]
    xdot12 = 1.0*flux_r34 + 1.0*flux_r6;  // xdot[12]
    xdot13 = 1.0*flux_r33 + 1.0*flux_r8;  // xdot[13]
    xdot14 = 1.0*flux_r10 + 1.0*flux_r35;  // xdot[14]
    xdot15 = 1.0*flux_r1;  // xdot[15]
    xdot16 = 1.0*flux_r3;  // xdot[16]
    xdot17 = 1.0*flux_r5;  // xdot[17]
    xdot18 = 1.0*flux_r7;  // xdot[18]
    xdot19 = 1.0*flux_r9;  // xdot[19]
    xdot20 = 1.0*flux_r11;  // xdot[20]
}

} // namespace model_one_country
} // namespace amici
