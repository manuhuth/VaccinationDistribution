#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include <gsl/gsl-lite.hpp>
#include <array>

#include "basic_sir_x.h"
#include "basic_sir_p.h"
#include "basic_sir_w.h"
#include "basic_sir_dwdx.h"

namespace amici {
namespace model_basic_sir {

void dwdx_basic_sir(realtype *dwdx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *tcl){
    dflux_r0_dsusceptible = -beta*infectious*susceptible/std::pow(infectious + recovered + susceptible, 2) + beta*infectious/(infectious + recovered + susceptible);  // dwdx[0]
    dflux_r0_dinfectious = -beta*infectious*susceptible/std::pow(infectious + recovered + susceptible, 2) + beta*susceptible/(infectious + recovered + susceptible);  // dwdx[1]
    dflux_r1_dinfectious = gamma;  // dwdx[2]
    dflux_r0_drecovered = -beta*infectious*susceptible/std::pow(infectious + recovered + susceptible, 2);  // dwdx[3]
}

} // namespace model_basic_sir
} // namespace amici
