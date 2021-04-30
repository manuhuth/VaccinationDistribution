#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "basic_sir_x.h"
#include "basic_sir_p.h"
#include "basic_sir_w.h"
#include "basic_sir_dwdp.h"

namespace amici {
namespace model_basic_sir {

void dwdp_basic_sir(realtype *dwdp, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl, const realtype *dtcldp){
    dflux_r0_dbeta = infectious*susceptible/(infectious + recovered + susceptible);  // dwdp[0]
    dflux_r1_dgamma = infectious;  // dwdp[1]
}

} // namespace model_basic_sir
} // namespace amici
