#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "test_x.h"
#include "test_p.h"
#include "test_w.h"

namespace amici {
namespace model_test {

void w_test(realtype *w, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *tcl){
    flux_r0 = beta*infectious*susceptible/(infectious + recovered + susceptible);  // w[0]
    flux_r1 = gamma*infectious;  // w[1]
}

} // namespace model_test
} // namespace amici
