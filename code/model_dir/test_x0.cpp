#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "test_p.h"

namespace amici {
namespace model_test {

void x0_test(realtype *x0, const realtype t, const realtype *p, const realtype *k){
    x0[0] = 100.0;
    x0[1] = 0.0001;
}

} // namespace model_test
} // namespace amici
