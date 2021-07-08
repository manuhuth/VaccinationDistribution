#include "amici/symbolic_functions.h"
#include "amici/defines.h"
#include "sundials/sundials_types.h"

#include "p.h"
#include "k.h"
#include "w.h"
#include "x.h"
#include "dwdx.h"

namespace amici {
namespace model_vaccination {

void dydx_vaccination(realtype *dydx, const realtype t, const realtype *x, const realtype *p, const realtype *k, const realtype *h, const realtype *w, const realtype *dwdx){
}

} // namespace amici
} // namespace model_vaccination