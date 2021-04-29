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

void y_test(realtype *y, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w){
    y[0] = susceptible;
    y[1] = infectious;
    y[2] = recovered;
    y[3] = 1.0;
}

} // namespace model_test
} // namespace amici
