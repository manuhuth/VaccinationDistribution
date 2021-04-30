#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "basic_sir_x.h"
#include "basic_sir_p.h"
#include "basic_sir_w.h"

namespace amici {
namespace model_basic_sir {

void w_basic_sir(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    flux_r0 = beta*infectious*susceptible/(infectious + recovered + susceptible);  // w[0]
    flux_r1 = gamma*infectious;  // w[1]
}

} // namespace model_basic_sir
} // namespace amici
