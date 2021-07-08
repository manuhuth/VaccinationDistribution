#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "tmp0n6wklye_p.h"
#include "tmp0n6wklye_k.h"
#include "tmp0n6wklye_sigmay.h"

namespace amici {
namespace model_tmp0n6wklye {

void sigmay_tmp0n6wklye(realtype *sigmay, const realtype t, const realtype *p, const realtype *k){
    sigma_observable_deaths = 1.0;  // sigmay[0]
}

} // namespace model_tmp0n6wklye
} // namespace amici
