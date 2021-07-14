#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "p.h"
#include "k.h"
#include "w.h"
#include "x.h"
#include "dwdx.h"

namespace amici {
namespace model_vaccination_neural_network {

void JDiag_vaccination_neural_network(realtype *JDiag, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
    JDiag[0] = -1.0*dwdx0 - 1.0*dwdx12 - 1.0*dwdx18 - 1.0*dwdx24 - 1.0*dwdx30 - 1.0*dwdx36 - 1.0*dwdx42 - 1.0*dwdx48 - 1.0*dwdx54 - 1.0*dwdx6 - 1.0*dwdx60 - 1.0*dwdx66;
    JDiag[1] = -1.0*dwdx103 - 1.0*dwdx109 - 1.0*dwdx115 - 1.0*dwdx121 - 1.0*dwdx127 - 1.0*dwdx133 - 1.0*dwdx139 - 1.0*dwdx73 - 1.0*dwdx79 - 1.0*dwdx85 - 1.0*dwdx91 - 1.0*dwdx97;
    JDiag[2] = -1.0*dwdx146 - 1.0*dwdx152 - 1.0*dwdx158 - 1.0*dwdx164 - 1.0*dwdx170 - 1.0*dwdx176 - 1.0*dwdx182 - 1.0*dwdx188 - 1.0*dwdx194 - 1.0*dwdx200 - 1.0*dwdx206 - 1.0*dwdx212;
    JDiag[3] = -1.0*dwdx219 - 1.0*dwdx225 - 1.0*dwdx231 - 1.0*dwdx237 - 1.0*dwdx243 - 1.0*dwdx249 - 1.0*dwdx255 - 1.0*dwdx261 - 1.0*dwdx267 - 1.0*dwdx273 - 1.0*dwdx279 - 1.0*dwdx285;
    JDiag[4] = -1.0*dwdx292 - 1.0*dwdx298 - 1.0*dwdx304 - 1.0*dwdx310 - 1.0*dwdx316 - 1.0*dwdx322 - 1.0*dwdx328 - 1.0*dwdx334 - 1.0*dwdx340 - 1.0*dwdx346 - 1.0*dwdx352 - 1.0*dwdx358;
    JDiag[5] = -1.0*dwdx365 - 1.0*dwdx371 - 1.0*dwdx377 - 1.0*dwdx383 - 1.0*dwdx389 - 1.0*dwdx395 - 1.0*dwdx401 - 1.0*dwdx407 - 1.0*dwdx413 - 1.0*dwdx419 - 1.0*dwdx425 - 1.0*dwdx431;
    JDiag[6] = -1.0*dwdx432 - 1.0*dwdx433 + 1.0*dwdx434 + 1.0*dwdx446 + 1.0*dwdx458 + 1.0*dwdx470 + 1.0*dwdx482 + 1.0*dwdx494;
    JDiag[7] = -1.0*dwdx506 - 1.0*dwdx507 + 1.0*dwdx514 + 1.0*dwdx526 + 1.0*dwdx538 + 1.0*dwdx550 + 1.0*dwdx562 + 1.0*dwdx574;
    JDiag[8] = -1.0*dwdx580 - 1.0*dwdx581 + 1.0*dwdx583 + 1.0*dwdx595 + 1.0*dwdx607 + 1.0*dwdx619 + 1.0*dwdx631 + 1.0*dwdx643;
    JDiag[9] = -1.0*dwdx654 - 1.0*dwdx655 + 1.0*dwdx663 + 1.0*dwdx675 + 1.0*dwdx687 + 1.0*dwdx699 + 1.0*dwdx711 + 1.0*dwdx723;
    JDiag[10] = -1.0*dwdx728 - 1.0*dwdx729 + 1.0*dwdx732 + 1.0*dwdx744 + 1.0*dwdx756 + 1.0*dwdx768 + 1.0*dwdx780 + 1.0*dwdx792;
    JDiag[11] = -1.0*dwdx802 - 1.0*dwdx803 + 1.0*dwdx812 + 1.0*dwdx824 + 1.0*dwdx836 + 1.0*dwdx848 + 1.0*dwdx860 + 1.0*dwdx872;
    JDiag[12] = -1.0*dwdx876 - 1.0*dwdx877 + 1.0*dwdx881 + 1.0*dwdx893 + 1.0*dwdx905 + 1.0*dwdx917 + 1.0*dwdx929 + 1.0*dwdx941;
    JDiag[13] = 1.0*dwdx1009 + 1.0*dwdx1021 - 1.0*dwdx950 - 1.0*dwdx951 + 1.0*dwdx961 + 1.0*dwdx973 + 1.0*dwdx985 + 1.0*dwdx997;
    JDiag[14] = -1.0*dwdx1024 - 1.0*dwdx1025 + 1.0*dwdx1030 + 1.0*dwdx1042 + 1.0*dwdx1054 + 1.0*dwdx1066 + 1.0*dwdx1078 + 1.0*dwdx1090;
    JDiag[15] = -1.0*dwdx1098 - 1.0*dwdx1099 + 1.0*dwdx1110 + 1.0*dwdx1122 + 1.0*dwdx1134 + 1.0*dwdx1146 + 1.0*dwdx1158 + 1.0*dwdx1170;
    JDiag[16] = -1.0*dwdx1172 - 1.0*dwdx1173 + 1.0*dwdx1179 + 1.0*dwdx1191 + 1.0*dwdx1203 + 1.0*dwdx1215 + 1.0*dwdx1227 + 1.0*dwdx1239;
    JDiag[17] = -1.0*dwdx1246 - 1.0*dwdx1247 + 1.0*dwdx1259 + 1.0*dwdx1271 + 1.0*dwdx1283 + 1.0*dwdx1295 + 1.0*dwdx1307 + 1.0*dwdx1319;
}

} // namespace amici
} // namespace model_vaccination_neural_network