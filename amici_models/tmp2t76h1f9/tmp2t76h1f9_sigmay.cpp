#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "tmp2t76h1f9_p.h"
#include "tmp2t76h1f9_k.h"
#include "tmp2t76h1f9_sigmay.h"

namespace amici {
namespace model_tmp2t76h1f9 {

void sigmay_tmp2t76h1f9(realtype *sigmay, const realtype t, const realtype *p, const realtype *k){
    sigma_observable_deaths = 1.0;  // sigmay[0]
}

} // namespace model_tmp2t76h1f9
} // namespace amici
