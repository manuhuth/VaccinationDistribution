#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "vaccination_x.h"
#include "vaccination_p.h"
#include "vaccination_h.h"
#include "vaccination_w.h"
#include "vaccination_dxdotdw.h"

namespace amici {
namespace model_vaccination {

void dxdotdw_vaccination(realtype *dxdotdw, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    dxdot6_dflux_r0 = -1.0;  // dxdotdw[0]
    dxdot18_dflux_r0 = 1.0;  // dxdotdw[1]
    dxdot6_dflux_r1 = -1.0;  // dxdotdw[2]
    dxdot30_dflux_r1 = 1.0;  // dxdotdw[3]
    dxdot7_dflux_r2 = -1.0;  // dxdotdw[4]
    dxdot19_dflux_r2 = 1.0;  // dxdotdw[5]
    dxdot7_dflux_r3 = -1.0;  // dxdotdw[6]
    dxdot31_dflux_r3 = 1.0;  // dxdotdw[7]
    dxdot8_dflux_r4 = -1.0;  // dxdotdw[8]
    dxdot20_dflux_r4 = 1.0;  // dxdotdw[9]
    dxdot8_dflux_r5 = -1.0;  // dxdotdw[10]
    dxdot32_dflux_r5 = 1.0;  // dxdotdw[11]
    dxdot9_dflux_r6 = -1.0;  // dxdotdw[12]
    dxdot21_dflux_r6 = 1.0;  // dxdotdw[13]
    dxdot9_dflux_r7 = -1.0;  // dxdotdw[14]
    dxdot33_dflux_r7 = 1.0;  // dxdotdw[15]
    dxdot10_dflux_r8 = -1.0;  // dxdotdw[16]
    dxdot22_dflux_r8 = 1.0;  // dxdotdw[17]
    dxdot10_dflux_r9 = -1.0;  // dxdotdw[18]
    dxdot34_dflux_r9 = 1.0;  // dxdotdw[19]
    dxdot11_dflux_r10 = -1.0;  // dxdotdw[20]
    dxdot23_dflux_r10 = 1.0;  // dxdotdw[21]
    dxdot11_dflux_r11 = -1.0;  // dxdotdw[22]
    dxdot35_dflux_r11 = 1.0;  // dxdotdw[23]
    dxdot12_dflux_r12 = -1.0;  // dxdotdw[24]
    dxdot24_dflux_r12 = 1.0;  // dxdotdw[25]
    dxdot12_dflux_r13 = -1.0;  // dxdotdw[26]
    dxdot36_dflux_r13 = 1.0;  // dxdotdw[27]
    dxdot13_dflux_r14 = -1.0;  // dxdotdw[28]
    dxdot25_dflux_r14 = 1.0;  // dxdotdw[29]
    dxdot13_dflux_r15 = -1.0;  // dxdotdw[30]
    dxdot37_dflux_r15 = 1.0;  // dxdotdw[31]
    dxdot14_dflux_r16 = -1.0;  // dxdotdw[32]
    dxdot26_dflux_r16 = 1.0;  // dxdotdw[33]
    dxdot14_dflux_r17 = -1.0;  // dxdotdw[34]
    dxdot38_dflux_r17 = 1.0;  // dxdotdw[35]
    dxdot15_dflux_r18 = -1.0;  // dxdotdw[36]
    dxdot27_dflux_r18 = 1.0;  // dxdotdw[37]
    dxdot15_dflux_r19 = -1.0;  // dxdotdw[38]
    dxdot39_dflux_r19 = 1.0;  // dxdotdw[39]
    dxdot16_dflux_r20 = -1.0;  // dxdotdw[40]
    dxdot28_dflux_r20 = 1.0;  // dxdotdw[41]
    dxdot16_dflux_r21 = -1.0;  // dxdotdw[42]
    dxdot40_dflux_r21 = 1.0;  // dxdotdw[43]
    dxdot17_dflux_r22 = -1.0;  // dxdotdw[44]
    dxdot29_dflux_r22 = 1.0;  // dxdotdw[45]
    dxdot17_dflux_r23 = -1.0;  // dxdotdw[46]
    dxdot41_dflux_r23 = 1.0;  // dxdotdw[47]
    dxdot0_dflux_r24 = -1.0;  // dxdotdw[48]
    dxdot6_dflux_r24 = 1.0;  // dxdotdw[49]
    dxdot1_dflux_r25 = -1.0;  // dxdotdw[50]
    dxdot8_dflux_r25 = 1.0;  // dxdotdw[51]
    dxdot2_dflux_r26 = -1.0;  // dxdotdw[52]
    dxdot10_dflux_r26 = 1.0;  // dxdotdw[53]
    dxdot3_dflux_r27 = -1.0;  // dxdotdw[54]
    dxdot12_dflux_r27 = 1.0;  // dxdotdw[55]
    dxdot4_dflux_r28 = -1.0;  // dxdotdw[56]
    dxdot14_dflux_r28 = 1.0;  // dxdotdw[57]
    dxdot5_dflux_r29 = -1.0;  // dxdotdw[58]
    dxdot16_dflux_r29 = 1.0;  // dxdotdw[59]
    dxdot0_dflux_r30 = -1.0;  // dxdotdw[60]
    dxdot7_dflux_r30 = 1.0;  // dxdotdw[61]
    dxdot1_dflux_r31 = -1.0;  // dxdotdw[62]
    dxdot9_dflux_r31 = 1.0;  // dxdotdw[63]
    dxdot2_dflux_r32 = -1.0;  // dxdotdw[64]
    dxdot11_dflux_r32 = 1.0;  // dxdotdw[65]
    dxdot3_dflux_r33 = -1.0;  // dxdotdw[66]
    dxdot13_dflux_r33 = 1.0;  // dxdotdw[67]
    dxdot4_dflux_r34 = -1.0;  // dxdotdw[68]
    dxdot15_dflux_r34 = 1.0;  // dxdotdw[69]
    dxdot5_dflux_r35 = -1.0;  // dxdotdw[70]
    dxdot17_dflux_r35 = 1.0;  // dxdotdw[71]
    dxdot0_dflux_r36 = -1.0;  // dxdotdw[72]
    dxdot6_dflux_r36 = 1.0;  // dxdotdw[73]
    dxdot1_dflux_r37 = -1.0;  // dxdotdw[74]
    dxdot8_dflux_r37 = 1.0;  // dxdotdw[75]
    dxdot2_dflux_r38 = -1.0;  // dxdotdw[76]
    dxdot10_dflux_r38 = 1.0;  // dxdotdw[77]
    dxdot3_dflux_r39 = -1.0;  // dxdotdw[78]
    dxdot12_dflux_r39 = 1.0;  // dxdotdw[79]
    dxdot4_dflux_r40 = -1.0;  // dxdotdw[80]
    dxdot14_dflux_r40 = 1.0;  // dxdotdw[81]
    dxdot5_dflux_r41 = -1.0;  // dxdotdw[82]
    dxdot16_dflux_r41 = 1.0;  // dxdotdw[83]
    dxdot0_dflux_r42 = -1.0;  // dxdotdw[84]
    dxdot7_dflux_r42 = 1.0;  // dxdotdw[85]
    dxdot1_dflux_r43 = -1.0;  // dxdotdw[86]
    dxdot9_dflux_r43 = 1.0;  // dxdotdw[87]
    dxdot2_dflux_r44 = -1.0;  // dxdotdw[88]
    dxdot11_dflux_r44 = 1.0;  // dxdotdw[89]
    dxdot3_dflux_r45 = -1.0;  // dxdotdw[90]
    dxdot13_dflux_r45 = 1.0;  // dxdotdw[91]
    dxdot4_dflux_r46 = -1.0;  // dxdotdw[92]
    dxdot15_dflux_r46 = 1.0;  // dxdotdw[93]
    dxdot5_dflux_r47 = -1.0;  // dxdotdw[94]
    dxdot17_dflux_r47 = 1.0;  // dxdotdw[95]
    dxdot0_dflux_r48 = -1.0;  // dxdotdw[96]
    dxdot6_dflux_r48 = 1.0;  // dxdotdw[97]
    dxdot1_dflux_r49 = -1.0;  // dxdotdw[98]
    dxdot8_dflux_r49 = 1.0;  // dxdotdw[99]
    dxdot2_dflux_r50 = -1.0;  // dxdotdw[100]
    dxdot10_dflux_r50 = 1.0;  // dxdotdw[101]
    dxdot3_dflux_r51 = -1.0;  // dxdotdw[102]
    dxdot12_dflux_r51 = 1.0;  // dxdotdw[103]
    dxdot4_dflux_r52 = -1.0;  // dxdotdw[104]
    dxdot14_dflux_r52 = 1.0;  // dxdotdw[105]
    dxdot5_dflux_r53 = -1.0;  // dxdotdw[106]
    dxdot16_dflux_r53 = 1.0;  // dxdotdw[107]
    dxdot0_dflux_r54 = -1.0;  // dxdotdw[108]
    dxdot7_dflux_r54 = 1.0;  // dxdotdw[109]
    dxdot1_dflux_r55 = -1.0;  // dxdotdw[110]
    dxdot9_dflux_r55 = 1.0;  // dxdotdw[111]
    dxdot2_dflux_r56 = -1.0;  // dxdotdw[112]
    dxdot11_dflux_r56 = 1.0;  // dxdotdw[113]
    dxdot3_dflux_r57 = -1.0;  // dxdotdw[114]
    dxdot13_dflux_r57 = 1.0;  // dxdotdw[115]
    dxdot4_dflux_r58 = -1.0;  // dxdotdw[116]
    dxdot15_dflux_r58 = 1.0;  // dxdotdw[117]
    dxdot5_dflux_r59 = -1.0;  // dxdotdw[118]
    dxdot17_dflux_r59 = 1.0;  // dxdotdw[119]
    dxdot0_dflux_r60 = -1.0;  // dxdotdw[120]
    dxdot6_dflux_r60 = 1.0;  // dxdotdw[121]
    dxdot1_dflux_r61 = -1.0;  // dxdotdw[122]
    dxdot8_dflux_r61 = 1.0;  // dxdotdw[123]
    dxdot2_dflux_r62 = -1.0;  // dxdotdw[124]
    dxdot10_dflux_r62 = 1.0;  // dxdotdw[125]
    dxdot3_dflux_r63 = -1.0;  // dxdotdw[126]
    dxdot12_dflux_r63 = 1.0;  // dxdotdw[127]
    dxdot4_dflux_r64 = -1.0;  // dxdotdw[128]
    dxdot14_dflux_r64 = 1.0;  // dxdotdw[129]
    dxdot5_dflux_r65 = -1.0;  // dxdotdw[130]
    dxdot16_dflux_r65 = 1.0;  // dxdotdw[131]
    dxdot0_dflux_r66 = -1.0;  // dxdotdw[132]
    dxdot7_dflux_r66 = 1.0;  // dxdotdw[133]
    dxdot1_dflux_r67 = -1.0;  // dxdotdw[134]
    dxdot9_dflux_r67 = 1.0;  // dxdotdw[135]
    dxdot2_dflux_r68 = -1.0;  // dxdotdw[136]
    dxdot11_dflux_r68 = 1.0;  // dxdotdw[137]
    dxdot3_dflux_r69 = -1.0;  // dxdotdw[138]
    dxdot13_dflux_r69 = 1.0;  // dxdotdw[139]
    dxdot4_dflux_r70 = -1.0;  // dxdotdw[140]
    dxdot15_dflux_r70 = 1.0;  // dxdotdw[141]
    dxdot5_dflux_r71 = -1.0;  // dxdotdw[142]
    dxdot17_dflux_r71 = 1.0;  // dxdotdw[143]
    dxdot0_dflux_r72 = -1.0;  // dxdotdw[144]
    dxdot6_dflux_r72 = 1.0;  // dxdotdw[145]
    dxdot1_dflux_r73 = -1.0;  // dxdotdw[146]
    dxdot8_dflux_r73 = 1.0;  // dxdotdw[147]
    dxdot2_dflux_r74 = -1.0;  // dxdotdw[148]
    dxdot10_dflux_r74 = 1.0;  // dxdotdw[149]
    dxdot3_dflux_r75 = -1.0;  // dxdotdw[150]
    dxdot12_dflux_r75 = 1.0;  // dxdotdw[151]
    dxdot4_dflux_r76 = -1.0;  // dxdotdw[152]
    dxdot14_dflux_r76 = 1.0;  // dxdotdw[153]
    dxdot5_dflux_r77 = -1.0;  // dxdotdw[154]
    dxdot16_dflux_r77 = 1.0;  // dxdotdw[155]
    dxdot0_dflux_r78 = -1.0;  // dxdotdw[156]
    dxdot7_dflux_r78 = 1.0;  // dxdotdw[157]
    dxdot1_dflux_r79 = -1.0;  // dxdotdw[158]
    dxdot9_dflux_r79 = 1.0;  // dxdotdw[159]
    dxdot2_dflux_r80 = -1.0;  // dxdotdw[160]
    dxdot11_dflux_r80 = 1.0;  // dxdotdw[161]
    dxdot3_dflux_r81 = -1.0;  // dxdotdw[162]
    dxdot13_dflux_r81 = 1.0;  // dxdotdw[163]
    dxdot4_dflux_r82 = -1.0;  // dxdotdw[164]
    dxdot15_dflux_r82 = 1.0;  // dxdotdw[165]
    dxdot5_dflux_r83 = -1.0;  // dxdotdw[166]
    dxdot17_dflux_r83 = 1.0;  // dxdotdw[167]
    dxdot0_dflux_r84 = -1.0;  // dxdotdw[168]
    dxdot6_dflux_r84 = 1.0;  // dxdotdw[169]
    dxdot1_dflux_r85 = -1.0;  // dxdotdw[170]
    dxdot8_dflux_r85 = 1.0;  // dxdotdw[171]
    dxdot2_dflux_r86 = -1.0;  // dxdotdw[172]
    dxdot10_dflux_r86 = 1.0;  // dxdotdw[173]
    dxdot3_dflux_r87 = -1.0;  // dxdotdw[174]
    dxdot12_dflux_r87 = 1.0;  // dxdotdw[175]
    dxdot4_dflux_r88 = -1.0;  // dxdotdw[176]
    dxdot14_dflux_r88 = 1.0;  // dxdotdw[177]
    dxdot5_dflux_r89 = -1.0;  // dxdotdw[178]
    dxdot16_dflux_r89 = 1.0;  // dxdotdw[179]
    dxdot0_dflux_r90 = -1.0;  // dxdotdw[180]
    dxdot7_dflux_r90 = 1.0;  // dxdotdw[181]
    dxdot1_dflux_r91 = -1.0;  // dxdotdw[182]
    dxdot9_dflux_r91 = 1.0;  // dxdotdw[183]
    dxdot2_dflux_r92 = -1.0;  // dxdotdw[184]
    dxdot11_dflux_r92 = 1.0;  // dxdotdw[185]
    dxdot3_dflux_r93 = -1.0;  // dxdotdw[186]
    dxdot13_dflux_r93 = 1.0;  // dxdotdw[187]
    dxdot4_dflux_r94 = -1.0;  // dxdotdw[188]
    dxdot15_dflux_r94 = 1.0;  // dxdotdw[189]
    dxdot5_dflux_r95 = -1.0;  // dxdotdw[190]
    dxdot17_dflux_r95 = 1.0;  // dxdotdw[191]
    dxdot0_dflux_r96 = -1.0;  // dxdotdw[192]
    dxdot1_dflux_r96 = 1.0;  // dxdotdw[193]
    dxdot0_dflux_r97 = -1.0;  // dxdotdw[194]
    dxdot2_dflux_r97 = 1.0;  // dxdotdw[195]
    dxdot3_dflux_r98 = -1.0;  // dxdotdw[196]
    dxdot4_dflux_r98 = 1.0;  // dxdotdw[197]
    dxdot3_dflux_r99 = -1.0;  // dxdotdw[198]
    dxdot5_dflux_r99 = 1.0;  // dxdotdw[199]
    dxdot18_dflux_r100 = -1.0;  // dxdotdw[200]
    dxdot20_dflux_r100 = 1.0;  // dxdotdw[201]
    dxdot18_dflux_r101 = -1.0;  // dxdotdw[202]
    dxdot22_dflux_r101 = 1.0;  // dxdotdw[203]
    dxdot19_dflux_r102 = -1.0;  // dxdotdw[204]
    dxdot21_dflux_r102 = 1.0;  // dxdotdw[205]
    dxdot19_dflux_r103 = -1.0;  // dxdotdw[206]
    dxdot23_dflux_r103 = 1.0;  // dxdotdw[207]
    dxdot24_dflux_r104 = -1.0;  // dxdotdw[208]
    dxdot26_dflux_r104 = 1.0;  // dxdotdw[209]
    dxdot24_dflux_r105 = -1.0;  // dxdotdw[210]
    dxdot28_dflux_r105 = 1.0;  // dxdotdw[211]
    dxdot25_dflux_r106 = -1.0;  // dxdotdw[212]
    dxdot27_dflux_r106 = 1.0;  // dxdotdw[213]
    dxdot25_dflux_r107 = -1.0;  // dxdotdw[214]
    dxdot29_dflux_r107 = 1.0;  // dxdotdw[215]
}

} // namespace model_vaccination
} // namespace amici
