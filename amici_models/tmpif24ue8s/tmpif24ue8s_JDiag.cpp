#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "p.h"
#include "k.h"
#include "w.h"
#include "x.h"
#include "dwdx.h"

namespace amici {
namespace model_tmpif24ue8s {

void JDiag_tmpif24ue8s(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0 - 1.0*dwdx12 - 1.0*dwdx18 - 1.0*dwdx24 - 1.0*dwdx30 - 1.0*dwdx36 - 1.0*dwdx42 - 1.0*dwdx48 - 1.0*dwdx54 - 1.0*dwdx6 - 1.0*dwdx60 - 1.0*dwdx66 - 1.0*dwdx72 - 1.0*dwdx73;
    JDiag[1] = -1.0*dwdx101 - 1.0*dwdx107 - 1.0*dwdx113 - 1.0*dwdx119 - 1.0*dwdx125 - 1.0*dwdx131 - 1.0*dwdx137 - 1.0*dwdx143 - 1.0*dwdx149 - 1.0*dwdx83 - 1.0*dwdx89 - 1.0*dwdx95;
    JDiag[2] = -1.0*dwdx156 - 1.0*dwdx162 - 1.0*dwdx168 - 1.0*dwdx174 - 1.0*dwdx180 - 1.0*dwdx186 - 1.0*dwdx192 - 1.0*dwdx198 - 1.0*dwdx204 - 1.0*dwdx210 - 1.0*dwdx216 - 1.0*dwdx222;
    JDiag[3] = -1.0*dwdx229 - 1.0*dwdx235 - 1.0*dwdx241 - 1.0*dwdx247 - 1.0*dwdx253 - 1.0*dwdx259 - 1.0*dwdx265 - 1.0*dwdx271 - 1.0*dwdx277 - 1.0*dwdx283 - 1.0*dwdx289 - 1.0*dwdx295 - 1.0*dwdx298 - 1.0*dwdx299;
    JDiag[4] = -1.0*dwdx312 - 1.0*dwdx318 - 1.0*dwdx324 - 1.0*dwdx330 - 1.0*dwdx336 - 1.0*dwdx342 - 1.0*dwdx348 - 1.0*dwdx354 - 1.0*dwdx360 - 1.0*dwdx366 - 1.0*dwdx372 - 1.0*dwdx378;
    JDiag[5] = -1.0*dwdx385 - 1.0*dwdx391 - 1.0*dwdx397 - 1.0*dwdx403 - 1.0*dwdx409 - 1.0*dwdx415 - 1.0*dwdx421 - 1.0*dwdx427 - 1.0*dwdx433 - 1.0*dwdx439 - 1.0*dwdx445 - 1.0*dwdx451;
    JDiag[6] = -1.0*dwdx452 - 1.0*dwdx453 + 1.0*dwdx454 + 1.0*dwdx466 + 1.0*dwdx478 + 1.0*dwdx490 + 1.0*dwdx502 + 1.0*dwdx514 - 1.0*dwdx528 - 1.0*dwdx529;
    JDiag[7] = -1.0*dwdx536 - 1.0*dwdx537 + 1.0*dwdx544 + 1.0*dwdx556 + 1.0*dwdx568 + 1.0*dwdx580 + 1.0*dwdx592 + 1.0*dwdx604 - 1.0*dwdx614 - 1.0*dwdx615;
    JDiag[8] = -1.0*dwdx620 - 1.0*dwdx621 + 1.0*dwdx623 + 1.0*dwdx635 + 1.0*dwdx647 + 1.0*dwdx659 + 1.0*dwdx671 + 1.0*dwdx683;
    JDiag[9] = -1.0*dwdx694 - 1.0*dwdx695 + 1.0*dwdx703 + 1.0*dwdx715 + 1.0*dwdx727 + 1.0*dwdx739 + 1.0*dwdx751 + 1.0*dwdx763;
    JDiag[10] = -1.0*dwdx768 - 1.0*dwdx769 + 1.0*dwdx772 + 1.0*dwdx784 + 1.0*dwdx796 + 1.0*dwdx808 + 1.0*dwdx820 + 1.0*dwdx832;
    JDiag[11] = -1.0*dwdx842 - 1.0*dwdx843 + 1.0*dwdx852 + 1.0*dwdx864 + 1.0*dwdx876 + 1.0*dwdx888 + 1.0*dwdx900 + 1.0*dwdx912;
    JDiag[12] = -1.0*dwdx916 - 1.0*dwdx917 + 1.0*dwdx921 + 1.0*dwdx933 + 1.0*dwdx945 + 1.0*dwdx957 + 1.0*dwdx969 + 1.0*dwdx981 - 1.0*dwdx992 - 1.0*dwdx993;
    JDiag[13] = -1.0*dwdx1000 - 1.0*dwdx1001 + 1.0*dwdx1011 + 1.0*dwdx1023 + 1.0*dwdx1035 + 1.0*dwdx1047 + 1.0*dwdx1059 + 1.0*dwdx1071 - 1.0*dwdx1078 - 1.0*dwdx1079;
    JDiag[14] = -1.0*dwdx1084 - 1.0*dwdx1085 + 1.0*dwdx1090 + 1.0*dwdx1102 + 1.0*dwdx1114 + 1.0*dwdx1126 + 1.0*dwdx1138 + 1.0*dwdx1150;
    JDiag[15] = -1.0*dwdx1158 - 1.0*dwdx1159 + 1.0*dwdx1170 + 1.0*dwdx1182 + 1.0*dwdx1194 + 1.0*dwdx1206 + 1.0*dwdx1218 + 1.0*dwdx1230;
    JDiag[16] = -1.0*dwdx1232 - 1.0*dwdx1233 + 1.0*dwdx1239 + 1.0*dwdx1251 + 1.0*dwdx1263 + 1.0*dwdx1275 + 1.0*dwdx1287 + 1.0*dwdx1299;
    JDiag[17] = -1.0*dwdx1306 - 1.0*dwdx1307 + 1.0*dwdx1319 + 1.0*dwdx1331 + 1.0*dwdx1343 + 1.0*dwdx1355 + 1.0*dwdx1367 + 1.0*dwdx1379;
    JDiag[18] = -1.0*dwdx1458 - 1.0*dwdx1459;
    JDiag[19] = -1.0*dwdx1542 - 1.0*dwdx1543;
    JDiag[24] = -1.0*dwdx1910 - 1.0*dwdx1911;
    JDiag[25] = -1.0*dwdx1994 - 1.0*dwdx1995;
}

} // namespace amici
} // namespace model_tmpif24ue8s