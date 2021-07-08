#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "tmp0n6wklye_x.h"
#include "tmp0n6wklye_p.h"
#include "tmp0n6wklye_k.h"
#include "tmp0n6wklye_h.h"
#include "tmp0n6wklye_w.h"

namespace amici {
namespace model_tmp0n6wklye {

void y_tmp0n6wklye(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = dead_countryA_vac0_virM + dead_countryA_vac0_virW + dead_countryA_vac1_virM + dead_countryA_vac1_virW + dead_countryA_vac2_virM + dead_countryA_vac2_virW + dead_countryB_vac0_virM + dead_countryB_vac0_virW + dead_countryB_vac1_virM + dead_countryB_vac1_virW + dead_countryB_vac2_virM + dead_countryB_vac2_virW;
}

} // namespace model_tmp0n6wklye
} // namespace amici
