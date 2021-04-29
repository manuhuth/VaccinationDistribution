#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "test_x.h"

namespace amici {
namespace model_test {

void x_rdata_test(realtype *x_rdata, const realtype *x, const realtype *tcl){
    x_rdata[0] = susceptible;
    x_rdata[1] = infectious;
    x_rdata[2] = recovered;
}

} // namespace model_test
} // namespace amici
