#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "p.h"
#include "k.h"
#include "w.h"
#include "x.h"

namespace amici {
namespace model_tmp2kc_o5jm {

void y_tmp2kc_o5jm(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = dead_countryA_vac0_virM + dead_countryA_vac0_virW + dead_countryA_vac1_virM + dead_countryA_vac1_virW + dead_countryA_vac2_virM + dead_countryA_vac2_virW + dead_countryB_vac0_virM + dead_countryB_vac0_virW + dead_countryB_vac1_virM + dead_countryB_vac1_virW + dead_countryB_vac2_virM + dead_countryB_vac2_virW;
}

} // namespace amici
} // namespace model_tmp2kc_o5jm