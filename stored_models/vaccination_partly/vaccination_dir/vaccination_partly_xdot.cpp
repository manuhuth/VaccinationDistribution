#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "vaccination_partly_x.h"
#include "vaccination_partly_p.h"
#include "vaccination_partly_h.h"
#include "vaccination_partly_w.h"
#include "vaccination_partly_xdot.h"

namespace amici {
namespace model_vaccination_partly {

void xdot_vaccination_partly(realtype *xdot, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    xdot0 = -1.0*flux_r24 - 1.0*flux_r30 - 1.0*flux_r36 - 1.0*flux_r42 - 1.0*flux_r48 - 1.0*flux_r54 - 1.0*flux_r60 - 1.0*flux_r66 - 1.0*flux_r72 - 1.0*flux_r78 - 1.0*flux_r84 - 1.0*flux_r90 - 1.0*flux_r96 - 1.0*flux_r97;  // xdot[0]
    xdot1 = -1.0*flux_r25 - 1.0*flux_r31 - 1.0*flux_r37 - 1.0*flux_r43 - 1.0*flux_r49 - 1.0*flux_r55 - 1.0*flux_r61 - 1.0*flux_r67 - 1.0*flux_r73 - 1.0*flux_r79 - 1.0*flux_r85 - 1.0*flux_r91 + 1.0*flux_r96;  // xdot[1]
    xdot2 = -1.0*flux_r26 - 1.0*flux_r32 - 1.0*flux_r38 - 1.0*flux_r44 - 1.0*flux_r50 - 1.0*flux_r56 - 1.0*flux_r62 - 1.0*flux_r68 - 1.0*flux_r74 - 1.0*flux_r80 - 1.0*flux_r86 - 1.0*flux_r92 + 1.0*flux_r97;  // xdot[2]
    xdot3 = -1.0*flux_r27 - 1.0*flux_r33 - 1.0*flux_r39 - 1.0*flux_r45 - 1.0*flux_r51 - 1.0*flux_r57 - 1.0*flux_r63 - 1.0*flux_r69 - 1.0*flux_r75 - 1.0*flux_r81 - 1.0*flux_r87 - 1.0*flux_r93 - 1.0*flux_r98 - 1.0*flux_r99;  // xdot[3]
    xdot4 = -1.0*flux_r28 - 1.0*flux_r34 - 1.0*flux_r40 - 1.0*flux_r46 - 1.0*flux_r52 - 1.0*flux_r58 - 1.0*flux_r64 - 1.0*flux_r70 - 1.0*flux_r76 - 1.0*flux_r82 - 1.0*flux_r88 - 1.0*flux_r94 + 1.0*flux_r98;  // xdot[4]
    xdot5 = -1.0*flux_r29 - 1.0*flux_r35 - 1.0*flux_r41 - 1.0*flux_r47 - 1.0*flux_r53 - 1.0*flux_r59 - 1.0*flux_r65 - 1.0*flux_r71 - 1.0*flux_r77 - 1.0*flux_r83 - 1.0*flux_r89 - 1.0*flux_r95 + 1.0*flux_r99;  // xdot[5]
    xdot6 = -1.0*flux_r0 - 1.0*flux_r1 + 1.0*flux_r24 + 1.0*flux_r36 + 1.0*flux_r48 + 1.0*flux_r60 + 1.0*flux_r72 + 1.0*flux_r84;  // xdot[6]
    xdot7 = -1.0*flux_r2 - 1.0*flux_r3 + 1.0*flux_r30 + 1.0*flux_r42 + 1.0*flux_r54 + 1.0*flux_r66 + 1.0*flux_r78 + 1.0*flux_r90;  // xdot[7]
    xdot8 = 1.0*flux_r25 + 1.0*flux_r37 - 1.0*flux_r4 + 1.0*flux_r49 - 1.0*flux_r5 + 1.0*flux_r61 + 1.0*flux_r73 + 1.0*flux_r85;  // xdot[8]
    xdot9 = 1.0*flux_r31 + 1.0*flux_r43 + 1.0*flux_r55 - 1.0*flux_r6 + 1.0*flux_r67 - 1.0*flux_r7 + 1.0*flux_r79 + 1.0*flux_r91;  // xdot[9]
    xdot10 = 1.0*flux_r26 + 1.0*flux_r38 + 1.0*flux_r50 + 1.0*flux_r62 + 1.0*flux_r74 - 1.0*flux_r8 + 1.0*flux_r86 - 1.0*flux_r9;  // xdot[10]
    xdot11 = -1.0*flux_r10 - 1.0*flux_r11 + 1.0*flux_r32 + 1.0*flux_r44 + 1.0*flux_r56 + 1.0*flux_r68 + 1.0*flux_r80 + 1.0*flux_r92;  // xdot[11]
    xdot12 = -1.0*flux_r12 - 1.0*flux_r13 + 1.0*flux_r27 + 1.0*flux_r39 + 1.0*flux_r51 + 1.0*flux_r63 + 1.0*flux_r75 + 1.0*flux_r87;  // xdot[12]
    xdot13 = -1.0*flux_r14 - 1.0*flux_r15 + 1.0*flux_r33 + 1.0*flux_r45 + 1.0*flux_r57 + 1.0*flux_r69 + 1.0*flux_r81 + 1.0*flux_r93;  // xdot[13]
    xdot14 = -1.0*flux_r16 - 1.0*flux_r17 + 1.0*flux_r28 + 1.0*flux_r40 + 1.0*flux_r52 + 1.0*flux_r64 + 1.0*flux_r76 + 1.0*flux_r88;  // xdot[14]
    xdot15 = -1.0*flux_r18 - 1.0*flux_r19 + 1.0*flux_r34 + 1.0*flux_r46 + 1.0*flux_r58 + 1.0*flux_r70 + 1.0*flux_r82 + 1.0*flux_r94;  // xdot[15]
    xdot16 = -1.0*flux_r20 - 1.0*flux_r21 + 1.0*flux_r29 + 1.0*flux_r41 + 1.0*flux_r53 + 1.0*flux_r65 + 1.0*flux_r77 + 1.0*flux_r89;  // xdot[16]
    xdot17 = -1.0*flux_r22 - 1.0*flux_r23 + 1.0*flux_r35 + 1.0*flux_r47 + 1.0*flux_r59 + 1.0*flux_r71 + 1.0*flux_r83 + 1.0*flux_r95;  // xdot[17]
    xdot18 = 1.0*flux_r0 - 1.0*flux_r100 - 1.0*flux_r101;  // xdot[18]
    xdot19 = -1.0*flux_r102 - 1.0*flux_r103 + 1.0*flux_r2;  // xdot[19]
    xdot20 = 1.0*flux_r100 + 1.0*flux_r4;  // xdot[20]
    xdot21 = 1.0*flux_r102 + 1.0*flux_r6;  // xdot[21]
    xdot22 = 1.0*flux_r101 + 1.0*flux_r8;  // xdot[22]
    xdot23 = 1.0*flux_r10 + 1.0*flux_r103;  // xdot[23]
    xdot24 = -1.0*flux_r104 - 1.0*flux_r105 + 1.0*flux_r12;  // xdot[24]
    xdot25 = -1.0*flux_r106 - 1.0*flux_r107 + 1.0*flux_r14;  // xdot[25]
    xdot26 = 1.0*flux_r104 + 1.0*flux_r16;  // xdot[26]
    xdot27 = 1.0*flux_r106 + 1.0*flux_r18;  // xdot[27]
    xdot28 = 1.0*flux_r105 + 1.0*flux_r20;  // xdot[28]
    xdot29 = 1.0*flux_r107 + 1.0*flux_r22;  // xdot[29]
    xdot30 = 1.0*flux_r1;  // xdot[30]
    xdot31 = 1.0*flux_r3;  // xdot[31]
    xdot32 = 1.0*flux_r5;  // xdot[32]
    xdot33 = 1.0*flux_r7;  // xdot[33]
    xdot34 = 1.0*flux_r9;  // xdot[34]
    xdot35 = 1.0*flux_r11;  // xdot[35]
    xdot36 = 1.0*flux_r13;  // xdot[36]
    xdot37 = 1.0*flux_r15;  // xdot[37]
    xdot38 = 1.0*flux_r17;  // xdot[38]
    xdot39 = 1.0*flux_r19;  // xdot[39]
    xdot40 = 1.0*flux_r21;  // xdot[40]
    xdot41 = 1.0*flux_r23;  // xdot[41]
    xdot42 = 1.0;  // xdot[42]
}

} // namespace model_vaccination_partly
} // namespace amici
